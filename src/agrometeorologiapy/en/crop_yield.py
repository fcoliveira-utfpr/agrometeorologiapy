"""Potential and attainable yield (FAO Agro-Ecological Zone model).

Model of Doorenbos & Kassam (1979), with the temperature corrections of
Barbieri & Tuon (1992) (in Pereira et al., 2002) and the relative sunshine
duration from the inverted Angstrom formula.

Names follow the yield gap literature (van Ittersum et al., 2013): Yp is the
potential yield and Yw the attainable (water-limited) yield; CE = Yw/Yp is
the climatic efficiency (Battisti et al., 2013).
"""

from calendar import monthrange

import numpy as np
import pandas as pd

from .radiation import (
    extraterrestrial_radiation,
    inverse_relative_distance,
    relative_sunshine_duration,
    solar_declination,
    sunset_hour_angle,
)

__all__ = [
    "temperature_correction_overcast",
    "temperature_correction_clear",
    "gross_potential_productivity",
    "lai_correction",
    "respiration_correction",
    "potential_yield",
    "attainable_yield",
]


def _like_input(value, ref):
    """Returns `value` with the type of `ref`: float for a scalar, a Series
    with the same index for a Series, numpy.ndarray otherwise."""
    if isinstance(ref, pd.Series):
        return pd.Series(np.asarray(value, dtype=float), index=ref.index)
    if np.ndim(ref) == 0:
        return float(value)
    return np.asarray(value, dtype=float)


def _pathway(pathway):
    r = str(pathway).upper()
    if r not in ("C3", "C4"):
        raise ValueError("pathway must be 'C3' or 'C4'.")
    return r


def temperature_correction_overcast(T, pathway, T_threshold=16.5):
    """
    Temperature correction of photosynthesis for overcast sky (cTo), from
    the polynomials of Barbieri & Tuon (1992), clipped at zero.

    C3: cTo = -0.0425 + 0.035 T + 0.00325 T² - 0.0000925 T³
    C4: cTo = -1.064 + 0.173 T - 0.0029 T²       (T >= 16.5 °C)
        cTo = -4.16 + 0.4325 T - 0.00725 T²      (T < 16.5 °C)

    Parameters
    ----------
    T : float or array_like
        Mean air temperature of the period, in °C.
    pathway : str
        Photosynthetic pathway of the crop: 'C3' (soybean, wheat, bean,
        sunflower, rice...) or 'C4' (maize, sorghum, sugarcane...).
    T_threshold : float, optional
        Temperature that separates the two C4 polynomials, in °C. Default
        16.5. Ignored for C3.

    Returns
    -------
    cTo : float or array_like
        Temperature correction for overcast sky (dimensionless, >= 0).

    Examples
    --------
    >>> round(temperature_correction_overcast(25, 'C3'), 3), round(temperature_correction_overcast(25, 'C4'), 3)
    (1.418, 1.448)
    """
    pathway = _pathway(pathway)
    t = np.asarray(T, dtype=float)
    if pathway == "C4":
        c = np.where(t >= T_threshold,
                     -1.064 + 0.173 * t - 0.0029 * t ** 2,
                     -4.16 + 0.4325 * t - 0.00725 * t ** 2)
    else:
        c = -0.0425 + 0.035 * t + 0.00325 * t ** 2 - 0.0000925 * t ** 3
    return _like_input(np.clip(c, 0.0, None), T)


def temperature_correction_clear(T, pathway, T_threshold=16.5):
    """
    Temperature correction of photosynthesis for clear sky (cTc), from the
    polynomials of Barbieri & Tuon (1992), clipped at zero.

    C3: cTc = 0.583 + 0.014 T + 0.0013 T² - 0.000037 T³
    C4: cTc = -4.16 + 0.4325 T - 0.00725 T²      (T >= 16.5 °C)
        cTc = -9.32 + 0.865 T - 0.0145 T²        (T < 16.5 °C)

    Parameters
    ----------
    T : float or array_like
        Mean air temperature of the period, in °C.
    pathway : str
        Photosynthetic pathway of the crop: 'C3' or 'C4'.
    T_threshold : float, optional
        Temperature that separates the two C4 polynomials, in °C. Default
        16.5. Ignored for C3.

    Returns
    -------
    cTc : float or array_like
        Temperature correction for clear sky (dimensionless, >= 0).

    Examples
    --------
    >>> round(temperature_correction_clear(25, 'C3'), 3), round(temperature_correction_clear(25, 'C4'), 3)
    (1.167, 2.121)
    """
    pathway = _pathway(pathway)
    t = np.asarray(T, dtype=float)
    if pathway == "C4":
        c = np.where(t >= T_threshold,
                     -4.16 + 0.4325 * t - 0.00725 * t ** 2,
                     -9.32 + 0.865 * t - 0.0145 * t ** 2)
    else:
        c = 0.583 + 0.014 * t + 0.0013 * t ** 2 - 0.000037 * t ** 3
    return _like_input(np.clip(c, 0.0, None), T)


def gross_potential_productivity(Ra, nN, cTc, cTo, a_c=107.2, b_c=8.604, a_o=31.7, b_o=5.234):
    """
    Gross potential productivity (GP) of the hypothetical standard crop
    (LAI = 5, no water or nutrient limitation): the sum of the clear-sky
    and overcast fractions of the day.

    GPc = (a_c + b_c Ra) cTc (n/N)
    GPo = (a_o + b_o Ra) cTo (1 - n/N)
    GP = GPc + GPo

    The default coefficients are those of Doorenbos & Kassam converted to
    Ra in MJ m-2 day-1 (0.36 and 0.219 with Ra in cal cm-2 day-1).

    Parameters
    ----------
    Ra : float or array_like
        Extraterrestrial radiation, in MJ m-2 day-1.
    nN : float or array_like
        Relative sunshine duration n/N (0 to 1). See
        `relative_sunshine_duration`.
    cTc, cTo : float or array_like
        Temperature corrections for clear and overcast sky. See
        `temperature_correction_clear` and `temperature_correction_overcast`.
    a_c, b_c : float, optional
        Coefficients of the clear-sky fraction. Default 107.2 and 8.604.
    a_o, b_o : float, optional
        Coefficients of the overcast fraction. Default 31.7 and 5.234.

    Returns
    -------
    GP : float or array_like
        Gross potential productivity, in kg DM ha-1 day-1.

    Examples
    --------
    >>> round(gross_potential_productivity(42.96, 0.501, cTc=2.210, cTo=1.484), 0)
    718.0
    """
    GPc = (a_c + b_c * np.asarray(Ra, dtype=float)) * cTc * nN
    GPo = (a_o + b_o * np.asarray(Ra, dtype=float)) * cTo * (1 - np.asarray(nN, dtype=float))
    return _like_input(GPc + GPo, Ra)


def lai_correction(LAI, LAI_standard=5.0):
    """
    Leaf area index correction (LAI_corr): adjusts the standard crop
    (LAI = 5) to the maximum LAI of the crop.

    LAI_corr = 0.0093 + 0.185 LAI - 0.0175 LAI²   (LAI < 5)
    LAI_corr = 0.5                                (LAI >= 5)

    Parameters
    ----------
    LAI : float
        Maximum leaf area index of the crop (m² m-2).
    LAI_standard : float, optional
        LAI from which LAI_corr = 0.5. Default 5.0.

    Returns
    -------
    LAI_corr : float
        Leaf area index correction (dimensionless).

    Examples
    --------
    >>> round(lai_correction(4.5), 3), lai_correction(5)
    (0.487, 0.5)
    """
    i = np.asarray(LAI, dtype=float)
    c = np.where(i >= LAI_standard, 0.5, 0.0093 + 0.185 * i - 0.0175 * i ** 2)
    return _like_input(c, LAI)


def respiration_correction(T, T_threshold=20.0, CR_cool=0.6, CR_warm=0.5):
    """
    Maintenance respiration correction (CR), from the mean temperature of
    the period.

    CR = 0.6 if T < 20 °C;  CR = 0.5 if T >= 20 °C.

    Parameters
    ----------
    T : float or array_like
        Mean air temperature of the period, in °C.
    T_threshold : float, optional
        Temperature at which CR changes, in °C. Default 20.0.
    CR_cool, CR_warm : float, optional
        CR below and from `T_threshold` on. Default 0.6 and 0.5.

    Returns
    -------
    CR : float or array_like
        Respiration correction (dimensionless).

    Examples
    --------
    >>> respiration_correction(18.0), respiration_correction(26.5)
    (0.6, 0.5)
    """
    c = np.where(np.asarray(T, dtype=float) < T_threshold, CR_cool, CR_warm)
    return _like_input(c, T)


def _period_starts(df, interval, year):
    """Start date of each period (pandas.DatetimeIndex), or None if the df
    has no dates. Accepts the 'date' column, a DatetimeIndex or the 'day'
    and 'month' columns (+ optional 'year'; without it, the year starts at
    `year` and advances when the month decreases)."""
    if 'date' in df.columns:
        dates = pd.DatetimeIndex(pd.to_datetime(df['date']))
    elif isinstance(df.index, pd.DatetimeIndex):
        dates = df.index
    elif {'day', 'month'} <= set(df.columns):
        months = df['month'].astype(int).to_numpy()
        if 'year' in df.columns:
            years = df['year'].astype(int).to_numpy()
        else:
            years = year + np.concatenate([[0], np.cumsum(np.diff(months) < 0)])
        dates = pd.DatetimeIndex([pd.Timestamp(int(a), int(m), int(d))
                                  for a, m, d in zip(years, months, df['day'].astype(int))])
    else:
        return None

    # Normalizes to the start of the calendar period.
    if interval == 'dekad':
        day0 = np.minimum((dates.day - 1) // 10, 2) * 10 + 1
        dates = pd.DatetimeIndex([pd.Timestamp(d.year, d.month, int(di)) for d, di in zip(dates, day0)])
    elif interval == 'M':
        dates = pd.DatetimeIndex([pd.Timestamp(d.year, d.month, 1) for d in dates])
    return dates.normalize()


def _period_days(start, interval):
    """Number of days of the calendar period that starts at `start`."""
    if interval == 'd':
        return 1
    days_in_month = monthrange(start.year, start.month)[1]
    if interval == 'M':
        return days_in_month
    return 10 if start.day <= 11 else days_in_month - 20


def potential_yield(df, lat, pathway, LAI=5.0, HI=0.35, moisture=13.0, interval='dekad',
                            sowing=None, cycle=None, year=2023, a=None, b=0.52,
                            T_threshold_CR=20.0, CR_cool=0.6, CR_warm=0.5):
    """
    Potential yield (Yp) by the FAO Agro-Ecological Zone model (Doorenbos
    & Kassam, 1979), accumulated period by period over the crop cycle.

    In each period (day, dekad or month):

    1. n/N from the inverted Angstrom formula (if not given);
    2. cTo and cTc of the C3 or C4 pathway at the mean temperature;
    3. GP = GPc + GPo (kg DM ha-1 day-1);
    4. CR from the mean temperature of the period;
    5. Yp_period = GP x ND x LAI_corr x CR x HI / (1 - 0.01 moisture),
       with ND = days of the period with the crop in the field.

    Yp = sum of Yp_period over the cycle (column 'Yp_cum').

    Parameters
    ----------
    df : pandas.DataFrame
        One row per period, in chronological order. Columns:

        - 'Tmean' (°C), required: mean temperature of the period.
        - Radiation: 'nN' (relative sunshine duration), or 'n' and 'N'
          (hours), or 'Rs' (MJ m-2 day-1, daily mean of the period), in
          this order of priority.
        - 'Ra' (MJ m-2 day-1, daily mean of the period), optional: if
          missing, it is computed from `lat` and the dates.
        - Dates: 'date' column, DatetimeIndex or 'day' and 'month' columns
          (+ optional 'year'), with the start of each period. Needed to
          compute Ra and ND and to use `sowing`.
        - 'ND', optional: number of days of the period to consider. If
          missing, it comes from the calendar (dekad: 10, 10 and 8 to 11
          days).
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    pathway : str
        Photosynthetic pathway of the crop: 'C3' or 'C4'.
    LAI : float, optional
        Maximum leaf area index of the crop. Default 5.0.
    HI : float, optional
        Harvest index. Default 0.35.
    moisture : float, optional
        Moisture content of the harvested product, in %. Default 13.0.
    interval : str, optional
        Time step of the rows of df: 'd' (daily), 'dekad' (calendar
        dekads: 1-10, 11-20, 21-end of month) or 'M' (monthly). Default
        'dekad'.
    sowing : str or datetime, optional
        Sowing date ('YYYY-MM-DD'). With `cycle`, only the days from sowing
        to sowing + cycle - 1 count in ND (the first and last periods
        become partial and periods outside the cycle are dropped).
    cycle : int, optional
        Length of the cycle, in days, counting the sowing day.
    year : int, optional
        Year of the first period when df has 'day' and 'month' but no
        'year'. Default 2023.
    a, b : float, optional
        Angstrom coefficients for n/N. Default a = 0.29 cos(lat), b = 0.52.
    T_threshold_CR, CR_cool, CR_warm : float, optional
        Parameters of the respiration correction (see
        `respiration_correction`). Default 20.0, 0.6 and 0.5.

    Returns
    -------
    result : pandas.DataFrame
        Copy of df (only the periods in the cycle), with the added columns
        'ND', 'Ra', 'nN', 'cTo', 'cTc', 'GPc', 'GPo', 'GP'
        (kg DM ha-1 day-1), 'CR', 'Yp_period' and 'Yp_cum' (kg/ha, at the
        given moisture). The last row of 'Yp_cum' is Yp.
        `result.attrs` holds 'Yp' and 'LAI_corr'.

    Examples
    --------
    First dekad of January in Santa Helena-PR, Brazil, maize (C4):

    >>> df = pd.DataFrame({'day': [1], 'month': [1], 'Tmean': [26.5], 'Rs': [22.5]})
    >>> yp = potential_yield(df, lat=-24.86, pathway='C4', LAI=5.0, HI=0.35, moisture=13)
    >>> round(float(yp['GP'].iloc[0]), 0), round(yp.attrs['Yp'], 0)
    (718.0, 722.0)
    """
    if interval not in ('d', 'dekad', 'M'):
        raise ValueError("interval must be 'd', 'dekad' or 'M'")
    if (sowing is None) != (cycle is None):
        raise ValueError("Give sowing and cycle together.")
    pathway = _pathway(pathway)

    res = df.copy()
    start = _period_starts(res, interval, year)

    # ---------- Days of each period in the field (ND) ----------
    if sowing is not None:
        if start is None:
            raise ValueError("To use sowing, df needs dates ('date', a DatetimeIndex or 'day' and 'month').")
        sow = pd.Timestamp(sowing).normalize()
        cycle_end = sow + pd.Timedelta(days=int(cycle) - 1)
        nd_cal = (res['ND'].to_numpy(dtype=float) if 'ND' in res.columns
                  else np.array([_period_days(d, interval) for d in start], dtype=float))
        end = start + pd.to_timedelta(nd_cal - 1, unit='D')
        field_start = start.where(start > sow, sow)
        field_end = end.where(end < cycle_end, cycle_end)
        nd = np.maximum(np.asarray((field_end - field_start).days) + 1, 0)
        in_cycle = nd > 0
        if start[0] > sow or end[-1] < cycle_end:
            raise ValueError("df does not cover the whole cycle (from sowing to sowing + cycle - 1).")
        res = res.loc[in_cycle].copy()
        res['ND'] = nd[in_cycle]
        start = field_start[in_cycle]
    elif 'ND' not in res.columns:
        if interval == 'd':
            res['ND'] = 1
        elif start is None:
            raise ValueError("For interval 'dekad' or 'M', df needs dates or the 'ND' column.")
        else:
            res['ND'] = [_period_days(d, interval) for d in start]

    # ---------- Ra: daily mean over the days of the period in the field ----------
    if 'Ra' not in res.columns:
        if start is None:
            raise ValueError("Without the 'Ra' column, df needs dates to compute it.")
        ra = []
        for d0, nd in zip(start, res['ND'].astype(int)):
            doy = np.array([(d0 + pd.Timedelta(days=k)).day_of_year for k in range(nd)])
            decl = solar_declination(doy)
            ws = sunset_hour_angle(lat, decl)
            ra.append(extraterrestrial_radiation(lat, decl, ws, inverse_relative_distance(doy)).mean())
        res['Ra'] = ra

    # ---------- Relative sunshine duration ----------
    if 'nN' in res.columns:
        res['nN'] = res['nN'].clip(0, 1)
    elif {'n', 'N'} <= set(res.columns):
        res['nN'] = (res['n'] / res['N']).clip(0, 1)
    elif 'Rs' in res.columns:
        res['nN'] = relative_sunshine_duration(res['Rs'], res['Ra'], lat, b=b, a=a)
    else:
        raise ValueError("df needs 'nN', 'n' and 'N', or 'Rs'.")

    # ---------- GP and Yp ----------
    res['cTo'] = temperature_correction_overcast(res['Tmean'], pathway)
    res['cTc'] = temperature_correction_clear(res['Tmean'], pathway)
    res['GPc'] = (107.2 + 8.604 * res['Ra']) * res['cTc'] * res['nN']
    res['GPo'] = (31.7 + 5.234 * res['Ra']) * res['cTo'] * (1 - res['nN'])
    res['GP'] = res['GPc'] + res['GPo']
    res['CR'] = respiration_correction(res['Tmean'], T_threshold=T_threshold_CR, CR_cool=CR_cool, CR_warm=CR_warm)

    lai_corr = lai_correction(LAI)
    res['Yp_period'] = res['GP'] * res['ND'] * lai_corr * res['CR'] * HI / (1 - 0.01 * moisture)
    res['Yp_cum'] = res['Yp_period'].cumsum()

    res.attrs['Yp'] = float(res['Yp_cum'].iloc[-1]) if len(res) else 0.0
    res.attrs['LAI_corr'] = lai_corr
    return res


def attainable_yield(Yp, ETa, ETc, ky, phase=None, method='product'):
    """
    Attainable (water-limited) yield (Yw): the potential yield penalized by
    the water deficit, with the relation of Doorenbos & Kassam (1979):

    1 - Yw/Yp = ky (1 - ETa/ETc)

    - 'product' (Battisti et al., 2013): applied phase by phase, in
      sequence; the yield at the end of one phase is the "potential" of
      the next.  Yw = Yp x prod[1 - ky_i (1 - ETa_i/ETc_i)]
    - 'single': whole-cycle ky and sums of ETa and ETc over the cycle.
      Yw = Yp x [1 - ky_cycle (1 - sum ETa / sum ETc)]

    Each factor is clipped at zero (Yw is never negative). Phases with
    ETc = 0 are not penalized.

    Parameters
    ----------
    Yp : float
        Potential yield, in kg/ha (see `potential_yield`).
    ETa, ETc : array_like
        Actual and crop evapotranspiration (mm), per phase or per period
        (e.g. the 'ETa' and 'ETc' columns of `crop_water_balance`).
    ky : float, sequence or dict
        Yield response factor. With 'product', one ky per phase: a sequence
        in phase order or a dict {phase: ky}. With 'single', the
        whole-cycle ky (scalar).
    phase : array_like, optional
        Phase of each period (same length as ETa and ETc). ETa and ETc are
        summed per phase, in the order in which the phases appear. If None,
        each element of ETa and ETc is already a phase.
    method : str, optional
        'product' (default) or 'single'.

    Returns
    -------
    result : pandas.DataFrame
        One row per phase (or a single 'cycle' row with 'single'): 'ETc',
        'ETa', 'ETa/ETc', 'ky', 'factor', 'Yw_start', 'Yw_end', 'loss'
        (kg/ha) and 'CE' (climatic efficiency Yw_end / Yp, cumulative). The
        last row of 'Yw_end' is Yw. `result.attrs` holds 'Yp', 'Yw' and 'CE'.

    Examples
    --------
    Maize with a deficit at flowering:

    >>> etc, eta = [60, 110, 120, 40], [57, 77, 102, 36]
    >>> yw = attainable_yield(1000, eta, etc, ky=[0.4, 1.5, 0.5, 0.2])
    >>> round(yw.attrs['CE'], 3)
    0.489
    >>> yw1 = attainable_yield(1000, eta, etc, ky=1.25, method='single')
    >>> round(yw1.attrs['CE'], 3)
    0.78
    """
    ETa = np.asarray(ETa, dtype=float)
    ETc = np.asarray(ETc, dtype=float)
    if ETa.shape != ETc.shape:
        raise ValueError("ETa and ETc must have the same length.")

    if method == 'single':
        if np.ndim(ky) != 0:
            raise ValueError("With method='single', ky must be the whole-cycle ky (scalar).")
        tab = pd.DataFrame({'ETc': [ETc.sum()], 'ETa': [ETa.sum()]}, index=pd.Index(['cycle'], name='phase'))
        ky_phases = [float(ky)]
    elif method == 'product':
        if phase is None:
            tab = pd.DataFrame({'ETc': ETc, 'ETa': ETa}, index=pd.RangeIndex(len(ETc), name='phase'))
        else:
            phase = np.asarray(phase)
            if phase.shape != ETa.shape:
                raise ValueError("phase must have the same length as ETa and ETc.")
            tab = (pd.DataFrame({'phase': phase, 'ETc': ETc, 'ETa': ETa})
                   .groupby('phase', sort=False)[['ETc', 'ETa']].sum())
        if isinstance(ky, dict):
            ky_phases = [float(ky[f]) for f in tab.index]
        else:
            ky_phases = list(np.atleast_1d(np.asarray(ky, dtype=float)))
            if len(ky_phases) != len(tab):
                raise ValueError(f"ky tem {len(ky_phases)} valores para {len(tab)} fases.")
    else:
        raise ValueError("method must be 'product' or 'single'")

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(tab['ETc'] > 0, tab['ETa'] / tab['ETc'], 1.0)
    tab['ETa/ETc'] = ratio
    tab['ky'] = ky_phases
    tab['factor'] = np.clip(1 - tab['ky'] * (1 - tab['ETa/ETc']), 0.0, None)

    yw = Yp * np.cumprod(tab['factor'].to_numpy())
    tab['Yw_start'] = np.concatenate([[Yp], yw[:-1]])
    tab['Yw_end'] = yw
    tab['loss'] = tab['Yw_start'] - tab['Yw_end']
    tab['CE'] = yw / Yp if Yp else np.nan

    tab.attrs['Yp'] = float(Yp)
    tab.attrs['Yw'] = float(yw[-1])
    tab.attrs['CE'] = float(tab['CE'].iloc[-1])
    return tab
