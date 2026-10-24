"""Evapotranspiration.

Potential evapotranspiration (PET) by Thornthwaite, Camargo, Hargreaves-
Samani and Priestley-Taylor, and reference evapotranspiration (ETo) by FAO
Penman-Monteith (Allen et al., 1998).
"""

import numpy as np

from .._trig import acosd, cosd, sind
from .humidity import saturation_vapour_pressure

__all__ = [
    "pet_thornthwaite",
    "pet_camargo",
    "pet_hargreaves_samani",
    "slope_vapour_pressure_curve",
    "pet_priestley_taylor",
    "wind_speed_2m",
    "eto_penman_monteith",
]

_DAYS_IN_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
_MID_MONTH_DOY = np.array([17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344])


def pet_thornthwaite(df, col_T='T_mean_C', lat=None):
    """
    Monthly potential evapotranspiration (PET) by Thornthwaite (1948).

    PET from the monthly mean air temperature, corrected by the daylight
    hours (as a function of latitude) and by the number of days of each
    month, following the classic formulation of the method.

    Parameters
    ----------
    df : pandas.DataFrame
        12 rows (January to December, in this order) with at least the
        monthly mean temperature column.
    col_T : str, optional
        Name of the monthly mean temperature column (°C). Default 'T_mean_C'.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).

    Returns
    -------
    df : pandas.DataFrame
        Copy of the input with the added column 'PET_mm_month' (corrected
        potential evapotranspiration, mm/month).
    """
    df = df.copy()
    T = df[col_T].to_numpy(dtype=float)

    i = np.where(T > 0, (T / 5) ** 1.514, 0)
    I = i.sum()
    a = 6.75e-7 * I ** 3 - 7.71e-5 * I ** 2 + 1.792e-2 * I + 0.49239
    PET_unadj = np.where(T > 0, 16 * (10 * T / I) ** a, 0)
    # Above 26.5 °C the original method uses a table, not the exponential equation.
    PET_unadj = np.where(T >= 26.5, -415.85 + 32.24 * T - 0.43 * T ** 2, PET_unadj)

    declination = 23.45 * sind(360 * (284 + _MID_MONTH_DOY) / 365)
    tan_lat = sind(lat) / cosd(lat)
    tan_dec = sind(declination) / cosd(declination)
    ws = acosd(-tan_lat * tan_dec)
    N = (2 / 15) * ws
    K = (N / 12) * (_DAYS_IN_MONTH / 30)

    df['PET_mm_month'] = np.round(PET_unadj * K, 2)
    return df


def pet_camargo(df, col_T='T_mean_C', lat=None):
    """
    Monthly potential evapotranspiration (PET) by Camargo (1971), with the
    coefficient F modified by Maluf (Camargo et al., 1999).

    PET from the monthly mean air temperature and the extraterrestrial
    radiation (Ra, computed internally from the solar declination, the
    sunset hour angle and the Earth-Sun distance), with an adjustment
    coefficient F that depends on the annual mean temperature.

    Parameters
    ----------
    df : pandas.DataFrame
        12 rows (January to December, in this order) with at least the
        monthly mean temperature column.
    col_T : str, optional
        Name of the monthly mean temperature column (°C). Default 'T_mean_C'.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).

    Returns
    -------
    df : pandas.DataFrame
        Copy of the input with the added columns 'Ra_MJ_m2day'
        (extraterrestrial radiation, MJ m-2 day-1), 'Ra_mm_day' (its
        evaporation equivalent, mm/day) and 'PET_mm_month' (mm/month).
    """
    df = df.copy()
    T = df[col_T].to_numpy(dtype=float)
    J = _MID_MONTH_DOY

    declination = 23.45 * np.sin(np.radians(360 * (284 + J) / 365))

    # Sunset hour angle (degrees)
    tan_lat = np.sin(np.radians(lat)) / np.cos(np.radians(lat))
    tan_dec = np.sin(np.radians(declination)) / np.cos(np.radians(declination))
    ws = np.degrees(np.arccos(np.clip(-tan_lat * tan_dec, -1, 1)))

    # Inverse relative distance Earth-Sun
    dr = 1 + 0.033 * np.cos(np.radians(360 * J / 365))

    # Daily extraterrestrial radiation, MJ m-2 day-1
    ws_rad = np.radians(ws)
    Ra_MJ = 37.6 * dr * (
        ws_rad * np.sin(np.radians(lat)) * np.sin(np.radians(declination))
        + np.cos(np.radians(lat)) * np.cos(np.radians(declination)) * np.sin(np.radians(ws))
    )
    Ra_mm = 0.408 * Ra_MJ  # evaporation equivalent, mm/day

    # Coefficient F as a function of the annual mean temperature (Maluf)
    T_annual = T.mean()
    if T_annual < 23:
        F = 0.01
    elif T_annual < 24:
        F = 0.0105
    else:
        F = 0.011

    PET = F * Ra_mm * T * _DAYS_IN_MONTH

    df['Ra_MJ_m2day'] = np.round(Ra_MJ, 2)
    df['Ra_mm_day'] = np.round(Ra_mm, 2)
    df['PET_mm_month'] = np.round(PET, 2)
    return df


def pet_hargreaves_samani(Ra, Tmax, Tmin, Tmean):
    """
    Potential evapotranspiration by Hargreaves & Samani (1985) (FAO-56
    eq. 52), from the daily temperature range and the extraterrestrial
    radiation.

    Parameters
    ----------
    Ra : float
        Extraterrestrial radiation, in MJ m-2 day-1.
    Tmax, Tmin, Tmean : float
        Daily maximum, minimum and mean air temperature, in °C.

    Returns
    -------
    PET : float
        Potential evapotranspiration, in mm/day.
    """
    Ra_mm = 0.408 * Ra  # evaporation equivalent, mm/day
    return 0.0023 * Ra_mm * (Tmax - Tmin) ** 0.5 * (Tmean + 17.8)


def slope_vapour_pressure_curve(T):
    """
    Slope of the saturation vapour pressure curve (Δ, FAO-56 eq. 13), the
    derivative of the Tetens equation with respect to temperature.

    Parameters
    ----------
    T : float
        Air temperature, in °C.

    Returns
    -------
    Delta : float
        Slope of the saturation vapour pressure curve, in kPa °C-1.
    """
    es = saturation_vapour_pressure(T)
    return 4098 * es / (T + 237.3) ** 2


def pet_priestley_taylor(Rn, G, Delta, gamma, alpha=1.26):
    """
    Potential evapotranspiration by Priestley & Taylor (1972).

    Simplification of the aerodynamic term of the Penman equation, valid
    for well-watered surfaces. The empirical coefficient alpha (1.26)
    accounts for the advection left out by ignoring the aerodynamic term.

    Parameters
    ----------
    Rn : float
        Net radiation, in MJ m-2 day-1.
    G : float
        Soil heat flux, in MJ m-2 day-1 (usually G = 0 at the daily scale).
    Delta : float
        Slope of the saturation vapour pressure curve, in kPa °C-1.
    gamma : float
        Psychrometric constant, in kPa °C-1.
    alpha : float, optional
        Priestley-Taylor coefficient. Default 1.26.

    Returns
    -------
    PET : float
        Potential evapotranspiration, in mm/day.
    """
    lambda_v = 2.45  # latent heat of vaporization, MJ/kg
    return alpha * (Delta / (Delta + gamma)) * (Rn - G) / lambda_v


def wind_speed_2m(uz, z):
    """
    Converts the wind speed measured at any height to the standard height
    of 2 m, with the logarithmic wind profile (FAO-56 eq. 47).

    Required by the FAO Penman-Monteith equation, which uses the wind speed
    at 2 m above a grass surface.

    Parameters
    ----------
    uz : float or array
        Wind speed measured at height z, in m/s.
    z : float or array
        Measurement height above the ground surface, in m.

    Returns
    -------
    u2 : float or array
        Wind speed at 2 m, in m/s.
    """
    return uz * 4.87 / np.log(67.8 * z - 5.42)


def eto_penman_monteith(Rn, G, Tmean, u2, es, ea, Delta, gamma):
    """
    Reference evapotranspiration (ETo) by the FAO Penman-Monteith equation
    (FAO-56 eq. 6), for the hypothetical grass reference crop.

    Parameters
    ----------
    Rn : float
        Net radiation at the crop surface, in MJ m-2 day-1.
    G : float
        Soil heat flux density, in MJ m-2 day-1.
    Tmean : float
        Mean daily air temperature at 2 m, in °C.
    u2 : float
        Wind speed at 2 m, in m/s.
    es, ea : float
        Saturation and actual vapour pressure, in kPa.
    Delta : float
        Slope of the saturation vapour pressure curve, in kPa °C-1.
    gamma : float
        Psychrometric constant, in kPa °C-1.

    Returns
    -------
    ETo : float
        Reference evapotranspiration, in mm/day.
    """
    numerator = 0.408 * Delta * (Rn - G) + gamma * (900 / (Tmean + 273)) * u2 * (es - ea)
    denominator = Delta + gamma * (1 + 0.34 * u2)
    return numerator / denominator
