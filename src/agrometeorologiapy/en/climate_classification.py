"""Climate classification: Köppen-Geiger, Thornthwaite (1948), Camargo (1991)
modified by Maluf (2000) and Holdridge life zones.

Each method has two versions:

- `<method>_classification(...)`: one location (lists/arrays of 12 months).
- `<method>_classification_grid(...)`: many locations at once (e.g. every
  pixel of a raster), with the 12 months on the first axis, shape (12, ...).

Locations without data (NaN in any month) get id 0 and class ''.
"""

import numpy as np

from .water_balance import climatic_water_balance_grid

__all__ = [
    "koppen_classification",
    "koppen_classification_grid",
    "thornthwaite_classification",
    "thornthwaite_classification_grid",
    "camargo_classification",
    "camargo_classification_grid",
    "holdridge_classification",
    "holdridge_classification_grid",
]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def _monthly(x, name):
    x = np.asarray(x, dtype=float)
    if x.shape[0] != 12:
        raise ValueError(f"{name} must have 12 months on the first axis (shape (12, ...)).")
    return x


def _spatial(x, shape):
    return np.broadcast_to(np.asarray(x, dtype=float), shape)


def _bands(x, limits, codes):
    """Band code [limits[i], limits[i+1]) of x; NaN -> 0."""
    idx = np.digitize(np.nan_to_num(x, nan=0.0), limits)
    out = np.asarray(codes, dtype=np.uint8)[idx]
    return np.where(np.isfinite(x), out, 0).astype(np.uint8)


def _labels(code, legend):
    """Array of integer codes -> array of strings (0 -> '')."""
    table = np.array([''] + [legend[k] for k in range(1, max(legend) + 1)], dtype=object)
    return np.asarray(table[code], dtype=object)


def _scalar(out):
    """Converts the output of a _grid function for a single location into scalars."""
    return {k: (v.item() if isinstance(v, np.ndarray) else v) for k, v in out.items()}


def _single_site(*monthly):
    for x in monthly:
        if x is not None and np.ndim(x) != 1:
            raise ValueError("For many locations use the _grid version.")


# Weighted astronomical seasons (southern hemisphere; index 0 = January):
# summer = 1/3 Dec + Jan + Feb + 2/3 Mar; autumn = 1/3 Mar + Apr + May + 2/3 Jun;
# winter = 1/3 Jun + Jul + Aug + 2/3 Sep; spring = 1/3 Sep + Oct + Nov + 2/3 Dec.
# In the northern hemisphere, summer <-> winter and autumn <-> spring.
_SOUTH_WEIGHTS = {
    'summer': {11: 1 / 3, 0: 1, 1: 1, 2: 2 / 3},
    'autumn': {2: 1 / 3, 3: 1, 4: 1, 5: 2 / 3},
    'winter': {5: 1 / 3, 6: 1, 7: 1, 8: 2 / 3},
    'spring': {8: 1 / 3, 9: 1, 10: 1, 11: 2 / 3},
}
_OPPOSITE = {'summer': 'winter', 'winter': 'summer', 'autumn': 'spring', 'spring': 'autumn'}


def _weights(season):
    w = np.zeros(12)
    for m, p in _SOUTH_WEIGHTS[season].items():
        w[m] = p
    return w


def _season_sum(x, season, north):
    """Weighted sum of x (12, ...) over the season, switching hemisphere where north."""
    extra = (1,) * (x.ndim - 1)
    south = (x * _weights(season).reshape((12,) + extra)).sum(0)
    nor = (x * _weights(_OPPOSITE[season]).reshape((12,) + extra)).sum(0)
    return np.where(north, nor, south)


# ---------------------------------------------------------------------------
# Köppen-Geiger (Alvares et al., 2013; seasonality of Kottek et al., 2006)
# ---------------------------------------------------------------------------
KOPPEN_CLASSES = {
    1: "Af", 2: "Am", 3: "As", 4: "Aw",
    5: "BSh", 6: "BSk", 7: "BWh", 8: "BWk",
    9: "Cfa", 10: "Cfb", 11: "Cfc",
    12: "Csa", 13: "Csb", 14: "Csc",
    15: "Cwa", 16: "Cwb", 17: "Cwc",
    18: "Dfa", 19: "Dfb", 20: "Dfc", 21: "Dfd",
    22: "Dsa", 23: "Dsb", 24: "Dsc", 25: "Dsd",
    26: "Dwa", 27: "Dwb", 28: "Dwc", 29: "Dwd",
    30: "ET", 31: "EF",
}

# Half-years Oct-Mar and Apr-Sep (index 0 = January).
_OCT_MAR = [9, 10, 11, 0, 1, 2]
_APR_SEP = [3, 4, 5, 6, 7, 8]


def koppen_classification_grid(T, P, lat):
    """
    Köppen-Geiger climate classification for many locations, following
    Alvares et al. (2013), with the f/s/w seasonality of groups C and D of
    Kottek et al. (2006).

    Summer is the Oct-Mar half-year in the southern hemisphere and Apr-Sep
    in the northern.

    Parameters
    ----------
    T : array_like, shape (12, ...)
        Monthly mean air temperature, in °C (January to December).
    P : array_like, shape (12, ...)
        Monthly precipitation, in mm/month.
    lat : float or array_like
        Latitude, in degrees (negative in the southern hemisphere); scalar
        or with the shape of the spatial dimensions.

    Returns
    -------
    dict of numpy.ndarray
        'id' (1-31, see KOPPEN_CLASSES; 0 = no data) and 'class' (e.g. 'Cfa').
    """
    T = _monthly(T, "T")
    P = _monthly(P, "P")
    if T.shape != P.shape:
        raise ValueError("T and P must have the same shape.")
    shape = T.shape[1:]
    north = _spatial(lat, shape) > 0
    valid = np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)

    tann = T.mean(0)
    tcold = T.min(0)
    thot = T.max(0)
    tmon10 = (T > 10).sum(0)
    rann = P.sum(0)
    rdry = P.min(0)

    om, aset = P[_OCT_MAR], P[_APR_SEP]
    p_summer = np.where(north, aset.sum(0), om.sum(0))
    p_winter = np.where(north, om.sum(0), aset.sum(0))
    psdry = np.where(north, aset.min(0), om.min(0))
    pwdry = np.where(north, om.min(0), aset.min(0))
    pswet = np.where(north, aset.max(0), om.max(0))
    pwwet = np.where(north, om.max(0), aset.max(0))

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_summer = np.where(rann > 0, p_summer / rann, 0.0)
        pct_winter = np.where(rann > 0, p_winter / rann, 0.0)

    # Aridity threshold (group B), in cm: 2T + 14, 2T (winter rain) or 2T + 28 (summer rain).
    p_th = 2 * tann + 14
    p_th = np.where(pct_winter >= 0.7, 2 * tann, p_th)
    p_th = np.where(pct_summer >= 0.7, 2 * tann + 28, p_th)

    group_a = tcold >= 18
    group_c = (thot > 10) & (tcold > -3) & (tcold < 18)
    group_d = (thot > 10) & (tcold <= -3)
    group_e = thot <= 10
    group_b = rann < 10 * p_th

    # A: f, m, s, w
    am_lim = 100 - rann / 25
    a_dry = group_a & (rdry < 60) & (rdry < am_lim)

    # C and D: s (dry summer), w (dry winter) and f (neither s nor w)
    is_s = (psdry < pwdry) & (pwwet > 3 * psdry) & (psdry < 40)
    is_w = (pwdry < psdry) & (pswet > 10 * pwdry)
    is_f = ~is_s & ~is_w
    hot = thot >= 22
    temperate = ~hot & (tmon10 >= 4)
    short = ~hot & ~temperate & (tmon10 >= 1)
    very_cold = tcold < -38

    code = np.zeros(shape, dtype=np.uint8)
    rules = [
        (group_a & (rdry >= 60), 1),
        (group_a & (rdry < 60) & (rdry >= am_lim), 2),
        (a_dry & (psdry < pwdry), 3),
        (a_dry & (psdry >= pwdry), 4),
    ]
    for group, base in ((group_c, 9), (group_d, 18)):
        step = 3 if base == 9 else 4
        for k, season_type in enumerate((is_f, is_s, is_w)):
            b = base + k * step
            rules += [(group & season_type & hot, b), (group & season_type & temperate, b + 1),
                       (group & season_type & short, b + 2)]
            if base == 18:
                rules.append((group & season_type & very_cold, b + 3))
    rules += [
        (group_e & (thot > 0), 30),
        (group_e & (thot <= 0), 31),
        # B last: overrides any other group.
        (group_b & (rann >= 5 * p_th) & (tann >= 18), 5),
        (group_b & (rann >= 5 * p_th) & (tann < 18), 6),
        (group_b & (rann < 5 * p_th) & (tann >= 18), 7),
        (group_b & (rann < 5 * p_th) & (tann < 18), 8),
    ]
    for cond, c in rules:
        code[cond] = c
    code[~valid] = 0
    return {'id': code, 'class': _labels(code, KOPPEN_CLASSES).astype(str)}


def koppen_classification(T, P, lat):
    """
    Köppen-Geiger climate classification of one location (Alvares et al.,
    2013; Kottek et al., 2006). See `koppen_classification_grid`.

    Parameters
    ----------
    T : list or array of 12 values
        Monthly mean air temperature (Jan-Dec), in °C.
    P : list or array of 12 values
        Monthly precipitation (Jan-Dec), in mm/month.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).

    Returns
    -------
    dict
        'id' (int) and 'class' (str, e.g. 'Cfa').
    """
    _single_site(T, P)
    return _scalar(koppen_classification_grid(T, P, lat))


# ---------------------------------------------------------------------------
# Thornthwaite (1948)
# ---------------------------------------------------------------------------
TH_MOISTURE = {1: "A", 2: "B4", 3: "B3", 4: "B2", 5: "B1", 6: "C2", 7: "C1", 8: "D", 9: "E"}
TH_MOISTURE_NAMES = {1: "Perhumid", 2: "Humid", 3: "Humid", 4: "Humid", 5: "Humid",
                     6: "Moist subhumid", 7: "Dry subhumid", 8: "Semiarid", 9: "Arid"}
TH_SUBTYPE = {1: "r", 2: "s", 3: "w", 4: "s2", 5: "w2",   # humid (Im >= 0), by the deficit
              6: "d", 7: "s", 8: "w", 9: "s2", 10: "w2"}  # dry (Im < 0), by the surplus
TH_THERMAL = {1: "A'", 2: "B'4", 3: "B'3", 4: "B'2", 5: "B'1", 6: "C'2", 7: "C'1", 8: "D'", 9: "E'"}
TH_THERMAL_NAMES = {1: "Megathermal", 2: "Mesothermal", 3: "Mesothermal", 4: "Mesothermal",
                    5: "Mesothermal", 6: "Microthermal", 7: "Microthermal", 8: "Tundra",
                    9: "Frost"}
TH_CONCENTRATION = {1: "a'", 2: "b'4", 3: "b'3", 4: "b'2", 5: "b'1", 6: "c'2", 7: "c'1", 8: "d'"}


def thornthwaite_classification_grid(P, PET, lat, AWC=100.0):
    """
    Thornthwaite (1948) climate classification for many locations, from the
    climatic water balance of Thornthwaite & Mather (1955).

    Tables of Aparecido et al. (2016), with three corrections from the
    original Thornthwaite (1948): B'3/B'2 limit at 855 mm (the paper has
    885); in dry climates, s2/w2 from Ih >= 20 (the paper has 33.3); and a
    larger surplus in summer = w, in winter = s (the paper has it reversed).

    Indices (%): Ih = 100 SUR/PET, Ia = 100 DEF/PET, Im = Ih - 0.6 Ia.
    Summer concentration: summer PET (1/3 Dec + Jan + Feb + 2/3 Mar, in the
    southern hemisphere) as % of the annual PET.

    Parameters
    ----------
    P, PET : array_like, shape (12, ...)
        Monthly precipitation and potential evapotranspiration, in mm/month.
    lat : float or array_like
        Latitude, in degrees (sets the hemisphere for summer and winter).
    AWC : float or array_like, optional
        Available water capacity of the soil, in mm. Default 100.

    Returns
    -------
    dict of numpy.ndarray
        'class' (e.g. "B1rA'a'"), 'moisture', 'subtype', 'thermal',
        'concentration' (integer codes, 0 = no data), 'Ih', 'Ia', 'Im',
        'PET_annual', 'DEF_annual', 'SUR_annual' (mm) and 'PET_summer_pct'.
    """
    P = _monthly(P, "P")
    PET = _monthly(PET, "ETP")
    shape = P.shape[1:]
    north = _spatial(lat, shape) > 0
    AWC = _spatial(AWC, shape)
    valid = (np.all(np.isfinite(P), axis=0) & np.all(np.isfinite(PET), axis=0)
              & np.isfinite(AWC) & (AWC > 0))

    wb = climatic_water_balance_grid(P, PET, np.where(valid, AWC, 1.0))
    pet_annual = PET.sum(0)
    pet_safe = np.maximum(pet_annual, 1e-9)
    def_annual, sur_annual = wb['DEF'].sum(0), wb['SUR'].sum(0)
    ih = 100 * sur_annual / pet_safe
    ia = 100 * def_annual / pet_safe
    im = ih - 0.6 * ia
    petr = 100 * _season_sum(PET, 'summer', north) / pet_safe

    moisture = _bands(im, [-40, -20, 0, 20, 40, 60, 80, 100], [9, 8, 7, 6, 5, 4, 3, 2, 1])

    def_s, def_w = _season_sum(wb['DEF'], 'summer', north), _season_sum(wb['DEF'], 'winter', north)
    sur_s, sur_w = _season_sum(wb['SUR'], 'summer', north), _season_sum(wb['SUR'], 'winter', north)
    humid, dry = im >= 0, im < 0
    def_summer, sur_summer = def_s > def_w, sur_s > sur_w
    subtype = np.zeros(shape, dtype=np.uint8)
    for cond, c in [
        (humid & (ia < 16.7), 1),
        (humid & (ia >= 16.7) & (ia < 33.3) & def_summer, 2),
        (humid & (ia >= 16.7) & (ia < 33.3) & ~def_summer, 3),
        (humid & (ia >= 33.3) & def_summer, 4),
        (humid & (ia >= 33.3) & ~def_summer, 5),
        (dry & (ih < 10), 6),
        (dry & (ih >= 10) & (ih < 20) & ~sur_summer, 7),
        (dry & (ih >= 10) & (ih < 20) & sur_summer, 8),
        (dry & (ih >= 20) & ~sur_summer, 9),
        (dry & (ih >= 20) & sur_summer, 10),
    ]:
        subtype[cond] = c

    thermal = _bands(pet_annual, [142, 285, 427, 570, 712, 855, 997, 1140], [9, 8, 7, 6, 5, 4, 3, 2, 1])
    concentration = _bands(petr, [48, 51.9, 56.3, 61.6, 68, 76.3, 88], [1, 2, 3, 4, 5, 6, 7, 8])

    for c in (moisture, subtype, thermal, concentration):
        c[~valid] = 0
    klass = (_labels(moisture, TH_MOISTURE) + _labels(subtype, TH_SUBTYPE)
              + _labels(thermal, TH_THERMAL) + _labels(concentration, TH_CONCENTRATION))
    klass = np.where(valid, klass, '').astype(str)

    nan = np.where(valid, 1.0, np.nan)
    return {
        'class': klass,
        'moisture': moisture, 'subtype': subtype, 'thermal': thermal, 'concentration': concentration,
        'Ih': ih * nan, 'Ia': ia * nan, 'Im': im * nan,
        'PET_annual': pet_annual * nan, 'DEF_annual': def_annual * nan, 'SUR_annual': sur_annual * nan,
        'PET_summer_pct': petr * nan,
    }


def thornthwaite_classification(P, PET, lat, AWC=100.0):
    """
    Thornthwaite (1948) climate classification of one location. See
    `thornthwaite_classification_grid`.

    Parameters
    ----------
    P, PET : list or array of 12 values
        Monthly precipitation and potential evapotranspiration (Jan-Dec), in mm/month.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    AWC : float, optional
        Available water capacity of the soil, in mm. Default 100.

    Returns
    -------
    dict
        'class' (e.g. "B1rA'a'"), 'description', codes and indices (Ih, Ia, Im...).
    """
    _single_site(P, PET)
    out = _scalar(thornthwaite_classification_grid(P, PET, lat, AWC))
    if out['class']:
        out['description'] = (f"{TH_MOISTURE_NAMES[out['moisture']]}, "
                              f"{TH_THERMAL_NAMES[out['thermal']].lower()}")
    else:
        out['description'] = ''
    return out


# ---------------------------------------------------------------------------
# Camargo (1991) modified by Maluf (2000)
# ---------------------------------------------------------------------------
CAMARGO_THERMAL = {1: 'GL', 2: 'FR', 3: 'CO', 4: 'TE', 5: 'STE', 6: 'ST', 7: 'TR', 8: 'EQ'}
CAMARGO_THERMAL_NAMES = {1: 'Glacial', 2: 'Cold', 3: 'Moderately cold', 4: 'Temperate',
                         5: 'Subtemperate', 6: 'Subtropical', 7: 'Tropical', 8: 'Equatorial'}
CAMARGO_MOISTURE = {1: 'DE', 2: 'AR', 3: 'SE', 4: 'MO', 5: 'SB', 6: 'UM', 7: 'PU', 8: 'SU'}
CAMARGO_MOISTURE_NAMES = {1: 'Desert', 2: 'Arid', 3: 'Dry', 4: 'Monsoonal', 5: 'Subhumid',
                          6: 'Humid', 7: 'Perhumid', 8: 'Extremely humid'}
CAMARGO_DRY_SEASON = {1: 'v', 2: 'o', 3: 'i', 4: 'p'}  # Portuguese initials of summer, autumn, winter, spring

_ZERO_MM = 0.05  # DEF and SUR below this count as zero (they round to 0.0 mm)


def camargo_classification_grid(T, P, PET, lat, AWC=100.0, thermal_table='consistent'):
    """
    Camargo (1991) climate classification modified by Maluf (2000), for
    many locations (Tables 6-8 of Aparecido et al., 2016).

    Thermal class from the annual mean temperature (Ty) and that of the
    coldest month (Tcold); moisture class from the annual deficit (DEF) and
    surplus (SUR) of the climatic water balance; letter of the season with
    the largest deficit (v, o, i, p: Portuguese initials of summer, autumn,
    winter and spring) in classes SE, MO, SB and UM.

    Table 6 of the paper is ambiguous for Tcold > 20 °C. With
    `thermal_table='consistent'` (default), Tcold only subdivides the range
    18 < Ty <= 22 (STE <= 13 < ST <= 20 < TR) and Ty > 25 is always EQ. With
    'literal', Tcold > 20 °C gives TR in any case (as in the GEE script of
    the climas_brasil repository), and class EQ practically never occurs.

    Table 7 does not cover 0 < DEF <= 150 mm with SUR = 0, nor DEF = 0 with
    SUR <= 200 mm; those cases go to SB (0 <= DEF <= 150 and 0 <= SUR <= 200).

    Parameters
    ----------
    T : array_like, shape (12, ...)
        Monthly mean air temperature, in °C.
    P, PET : array_like, shape (12, ...)
        Monthly precipitation and potential evapotranspiration, in mm/month.
    lat : float or array_like
        Latitude, in degrees (sets the hemisphere of the seasons).
    AWC : float or array_like, optional
        Available water capacity of the soil, in mm. Default 100.
    thermal_table : {'consistent', 'literal'}, optional
        Reading of Table 6 (see above). Default 'consistent'.

    Returns
    -------
    dict of numpy.ndarray
        'class' (e.g. 'ST-UMi'), 'thermal', 'moisture', 'dry_season'
        (integer codes, 0 = no data/no letter), 'T_annual', 'T_coldest',
        'DEF_annual' and 'SUR_annual'.
    """
    if thermal_table not in ('consistent', 'literal'):
        raise ValueError("thermal_table must be 'consistent' or 'literal'.")
    T = _monthly(T, "T")
    P = _monthly(P, "P")
    PET = _monthly(PET, "ETP")
    shape = T.shape[1:]
    north = _spatial(lat, shape) > 0
    AWC = _spatial(AWC, shape)
    valid = (np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)
              & np.all(np.isfinite(PET), axis=0) & np.isfinite(AWC) & (AWC > 0))

    ty = T.mean(0)
    tcold = T.min(0)
    wb = climatic_water_balance_grid(P, PET, np.where(valid, AWC, 1.0))
    def_annual, sur_annual = wb['DEF'].sum(0), wb['SUR'].sum(0)
    DEF = np.where(def_annual < _ZERO_MM, 0.0, def_annual)
    SUR = np.where(sur_annual < _ZERO_MM, 0.0, sur_annual)

    thermal = np.zeros(shape, dtype=np.uint8)
    for cond, c in [
        (ty <= 3, 1),
        ((ty > 3) & (ty <= 7), 2),
        ((ty > 7) & (ty <= 12), 3),
        ((ty > 12) & (ty <= 18), 4),
        ((ty > 18) & (ty <= 22) & (tcold <= 13), 5),
        ((ty > 18) & (ty <= 22) & (tcold > 13) & (tcold <= 20), 6),
        ((ty > 18) & (ty <= 22) & (tcold > 20), 7),
        ((ty > 22) & (ty <= 25), 7),
        (ty > 25, 8),
    ]:
        thermal[cond] = c
    if thermal_table == 'literal':
        thermal[tcold > 20] = 7

    moisture = np.zeros(shape, dtype=np.uint8)
    for cond, c in [
        ((DEF > 800) & (SUR == 0), 1),
        ((DEF > 150) & (DEF <= 800) & (SUR == 0), 2),
        ((DEF > 150) & (SUR > 0) & (SUR <= 200), 3),
        ((DEF > 150) & (SUR > 200), 4),
        ((DEF <= 150) & (SUR <= 200), 5),
        ((DEF > 0) & (DEF <= 150) & (SUR > 200), 6),
        ((DEF == 0) & (SUR > 200) & (SUR <= 1000), 7),
        ((DEF == 0) & (SUR > 1000), 8),
    ]:
        moisture[cond] = c

    def_seasons = np.stack([_season_sum(wb['DEF'], e, north)
                             for e in ('summer', 'autumn', 'winter', 'spring')])
    with_letter = np.isin(moisture, [3, 4, 5, 6]) & (def_seasons.max(0) >= _ZERO_MM)
    dry_season = np.where(with_letter, def_seasons.argmax(0) + 1, 0).astype(np.uint8)

    for c in (thermal, moisture, dry_season):
        c[~valid] = 0
    klass = (_labels(thermal, CAMARGO_THERMAL) + '-' + _labels(moisture, CAMARGO_MOISTURE)
              + _labels(dry_season, CAMARGO_DRY_SEASON))
    klass = np.where(valid, klass, '').astype(str)

    nan = np.where(valid, 1.0, np.nan)
    return {
        'class': klass, 'thermal': thermal, 'moisture': moisture, 'dry_season': dry_season,
        'T_annual': ty * nan, 'T_coldest': tcold * nan,
        'DEF_annual': def_annual * nan, 'SUR_annual': sur_annual * nan,
    }


def camargo_classification(T, P, PET, lat, AWC=100.0, thermal_table='consistent'):
    """
    Camargo (1991) mod. Maluf (2000) climate classification of one
    location. See `camargo_classification_grid`.

    Parameters
    ----------
    T : list or array of 12 values
        Monthly mean air temperature (Jan-Dec), in °C.
    P, PET : list or array of 12 values
        Monthly precipitation and potential evapotranspiration, in mm/month.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    AWC : float, optional
        Available water capacity of the soil, in mm. Default 100.
    thermal_table : {'consistent', 'literal'}, optional
        Reading of Table 6 of Aparecido et al. (2016). Default 'consistent'.

    Returns
    -------
    dict
        'class' (e.g. 'ST-UMi'), 'description', codes and annual variables.
    """
    _single_site(T, P, PET)
    out = _scalar(camargo_classification_grid(T, P, PET, lat, AWC, thermal_table))
    if out['class']:
        out['description'] = (f"{CAMARGO_THERMAL_NAMES[out['thermal']]} "
                              f"{CAMARGO_MOISTURE_NAMES[out['moisture']].lower()}")
    else:
        out['description'] = ''
    return out


# ---------------------------------------------------------------------------
# Holdridge life zones (38 zones; Leemans, 1990; Jungkunst et al., 2021)
# ---------------------------------------------------------------------------
# Original names (Jungkunst et al., 2021, Table 1).
HOLDRIDGE_ZONES = {
    1: "Polar ice", 2: "Polar desert",
    3: "Subpolar dry tundra", 4: "Subpolar moist tundra", 5: "Subpolar wet tundra",
    6: "Subpolar rain tundra",
    7: "Boreal desert", 8: "Boreal dry bush", 9: "Boreal moist forest", 10: "Boreal wet forest",
    11: "Boreal rain forest",
    12: "Cool temperate desert", 13: "Cool temperate desert bush", 14: "Cool temperate steppe",
    15: "Cool temperate moist forest", 16: "Cool temperate wet forest",
    17: "Cool temperate rain forest",
    18: "Warm temperate desert", 19: "Warm temperate desert bush", 20: "Warm temperate thorn steppe",
    21: "Warm temperate dry forest", 22: "Warm temperate moist forest",
    23: "Warm temperate wet forest", 24: "Warm temperate rain forest",
    25: "Subtropical desert", 26: "Subtropical desert bush", 27: "Subtropical thorn steppe",
    28: "Subtropical dry forest", 29: "Subtropical moist forest", 30: "Subtropical wet forest",
    31: "Subtropical rain forest",
    32: "Tropical desert", 33: "Tropical desert bush", 34: "Tropical thorn steppe",
    35: "Tropical very dry forest", 36: "Tropical dry forest", 37: "Tropical moist forest",
    38: "Tropical wet forest",
}

# (biotemperature class, minimum and maximum PET/P class, zone). Each thermal band starts at
# a different PET/P ratio (Jungkunst et al., 2021, Fig. 1); values beyond the ends fall in the
# extreme zone of the band.
_HOLDRIDGE_TABLE = [
    (2, 1, 6, 3), (2, 7, 7, 4), (2, 8, 8, 5), (2, 9, 10, 6),
    (3, 1, 5, 7), (3, 6, 6, 8), (3, 7, 7, 9), (3, 8, 8, 10), (3, 9, 10, 11),
    (4, 1, 4, 12), (4, 5, 5, 13), (4, 6, 6, 14), (4, 7, 7, 15), (4, 8, 8, 16), (4, 9, 10, 17),
    (5, 1, 3, 18), (5, 4, 4, 19), (5, 5, 5, 20), (5, 6, 6, 21), (5, 7, 7, 22), (5, 8, 8, 23),
    (5, 9, 10, 24),
    (6, 1, 3, 25), (6, 4, 4, 26), (6, 5, 5, 27), (6, 6, 6, 28), (6, 7, 7, 29), (6, 8, 8, 30),
    (6, 9, 10, 31),
    (7, 1, 2, 32), (7, 3, 3, 33), (7, 4, 4, 34), (7, 5, 5, 35), (7, 6, 6, 36), (7, 7, 7, 37),
    (7, 8, 10, 38),
]

HOLDRIDGE_PET_FACTOR = 58.93  # mm/year per °C of biotemperature (Holdridge, 1967)


def holdridge_classification_grid(T, P, lat, PET=None, correction_threshold=24.0):
    """
    Holdridge life zones (38 zones) for many locations, from the
    biotemperature and the potential evapotranspiration ratio (PET/P).

    Biotemperature: mean of the monthly temperatures corrected for
    latitude, t - 3|lat|/100 (t - 24)^2 in the months with
    t > `correction_threshold`, and clipped to [0, 30] °C.

    Parameters
    ----------
    T : array_like, shape (12, ...)
        Monthly mean air temperature, in °C.
    P : array_like, shape (12, ...)
        Monthly precipitation, in mm/month.
    lat : float or array_like
        Latitude, in degrees.
    PET : array_like, shape (12, ...), optional
        Monthly PET, in mm/month (e.g. Penman-Monteith). If None (default),
        uses Holdridge's PET: 58.93 x biotemperature (mm/year).
    correction_threshold : float or None, optional
        Monthly temperature (°C) above which the latitude correction is
        applied. Default 24 (Holdridge, 1967). None applies it to every month.

    Returns
    -------
    dict of numpy.ndarray
        'id' (1-38, see HOLDRIDGE_ZONES; 0 = no data), 'class' (zone name,
        Jungkunst et al., 2021), 'biotemperature' (°C), 'P_annual' and
        'PET_annual' (mm) and 'PET_ratio' (PET/P).
    """
    T = _monthly(T, "T")
    P = _monthly(P, "P")
    shape = T.shape[1:]
    alat = np.abs(_spatial(lat, shape))
    valid = np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)

    correction = 3 * alat / 100 * (T - 24) ** 2
    if correction_threshold is not None:
        correction = np.where(T > correction_threshold, correction, 0.0)
    biotemp = np.clip(T - correction, 0, 30).mean(0)

    p_annual = P.sum(0)
    if PET is None:
        pet_annual = HOLDRIDGE_PET_FACTOR * biotemp
    else:
        PET = _monthly(PET, "ETP")
        pet_annual = PET.sum(0)
        valid &= np.all(np.isfinite(PET), axis=0)
    ratio = pet_annual / np.maximum(p_annual, 1.0)

    t_cls = _bands(biotemp, [1.5, 3, 6, 12, 18, 24], [1, 2, 3, 4, 5, 6, 7])
    u_cls = _bands(ratio, [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32], [10, 9, 8, 7, 6, 5, 4, 3, 2, 1])

    zone = np.zeros(shape, dtype=np.uint8)
    zone[(t_cls == 1) & (biotemp < 1)] = 1
    zone[(t_cls == 1) & (biotemp >= 1)] = 2
    for tc, umin, umax, z in _HOLDRIDGE_TABLE:
        zone[(t_cls == tc) & (u_cls >= umin) & (u_cls <= umax)] = z
    zone[~valid] = 0

    nan = np.where(valid, 1.0, np.nan)
    return {
        'id': zone, 'class': _labels(zone, HOLDRIDGE_ZONES).astype(str),
        'biotemperature': biotemp * nan, 'P_annual': p_annual * nan,
        'PET_annual': pet_annual * nan, 'PET_ratio': ratio * nan,
    }


def holdridge_classification(T, P, lat, PET=None, correction_threshold=24.0):
    """
    Holdridge life zone of one location. See `holdridge_classification_grid`.

    Parameters
    ----------
    T : list or array of 12 values
        Monthly mean air temperature (Jan-Dec), in °C.
    P : list or array of 12 values
        Monthly precipitation, in mm/month.
    lat : float
        Latitude, in degrees.
    PET : list or array of 12 values, optional
        Monthly PET, in mm/month. If None, uses 58.93 x biotemperature.
    correction_threshold : float or None, optional
        Threshold (°C) of the latitude correction. Default 24.

    Returns
    -------
    dict
        'id', 'class' (zone name), 'biotemperature', 'P_annual', 'PET_annual'
        and 'PET_ratio'.
    """
    _single_site(T, P, PET)
    return _scalar(holdridge_classification_grid(T, P, lat, PET, correction_threshold))
