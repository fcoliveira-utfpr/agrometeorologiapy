"""Solar radiation.

Terminology and symbols follow FAO Irrigation and Drainage Paper 56
(Allen et al., 1998, ch. 3). Angles are in degrees.
"""

from datetime import datetime

import numpy as np

from .._trig import acosd, cosd, sind, tand

__all__ = [
    "day_of_year",
    "solar_declination",
    "hour_angle",
    "solar_zenith_angle",
    "solar_azimuth",
    "shadow_length",
    "daylight_hours",
    "sunset_hour_angle",
    "inverse_relative_distance",
    "extraterrestrial_radiation",
    "sunshine_duration",
    "solar_radiation_angstrom",
    "relative_sunshine_duration",
    "solar_radiation_hargreaves",
]


def day_of_year(day, month, year=2023):
    """
    Day of the year (J, FAO-56).

    Parameters
    ----------
    day : int
        Day of the month (1–31).
    month : int
        Month (1–12).
    year : int, optional
        Year. Default 2023 (not a leap year).

    Returns
    -------
    int
        Day of the year (1 to 365, or 366 in leap years).

    Examples
    --------
    >>> day_of_year(1, 1)
    1
    >>> day_of_year(25, 12)
    359
    >>> day_of_year(29, 2, 2024)
    60
    """
    return datetime(year, month, day).timetuple().tm_yday


def solar_declination(doy):
    """
    Solar declination (δ) for a given day of the year.

    Sine approximation of the annual variation of the declination caused
    by the 23.45° tilt of the Earth's axis. Positive when the Sun is over
    the northern hemisphere, negative over the southern hemisphere.

    Parameters
    ----------
    doy : int
        Day of the year (1 Jan = 1, ..., 31 Dec = 365).

    Returns
    -------
    delta : float
        Solar declination, in degrees.
    """
    return 23.45 * sind(360 * (doy - 80) / 365)


def hour_angle(hour, minute=0):
    """
    Hour angle from the local (true solar) time.

    The Earth rotates 15° per hour. The hour angle is zero at local solar
    noon, negative in the morning and positive in the afternoon.

    Parameters
    ----------
    hour : int
        Hour of the day (0–23).
    minute : int, optional
        Minutes (0–59). Default 0.

    Returns
    -------
    h : float
        Hour angle, in degrees.
    """
    return (hour + minute / 60 - 12) * 15


def solar_zenith_angle(lat, declination):
    """
    Solar zenith angle at local solar noon (h = 0°), the smallest zenith
    angle (highest Sun position) of the day.

    Parameters
    ----------
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    declination : float
        Solar declination of the day, in degrees.

    Returns
    -------
    Z : float
        Zenith angle, in degrees.
    """
    cos_Z = sind(lat) * sind(declination) + cosd(lat) * cosd(declination)
    cos_Z = max(-1.0, min(1.0, cos_Z))  # avoids a domain error in acosd
    return acosd(cos_Z)


def solar_azimuth(lat, declination, zenith):
    """
    Solar azimuth (α): horizontal direction of the Sun relative to the
    north-south line (referenced to the south in the southern hemisphere).

    Parameters
    ----------
    lat : float
        Latitude, in degrees.
    declination : float
        Solar declination of the day, in degrees.
    zenith : float
        Solar zenith angle, in degrees.

    Returns
    -------
    alpha : float
        Solar azimuth, in degrees.
    """
    num = sind(lat) * cosd(zenith) - sind(declination)
    den = cosd(lat) * sind(zenith)
    return acosd(num / den)


def shadow_length(height, zenith):
    """
    Length of the shadow cast by an object of a given height.

    Parameters
    ----------
    height : float
        Height of the object (e.g. m).
    zenith : float
        Solar zenith angle at the time considered, in degrees.

    Returns
    -------
    S : float
        Shadow length, in the unit of `height`.
    """
    return height * tand(zenith)


def daylight_hours(ws):
    """
    Daylight hours (N, FAO-56 eq. 34), the maximum possible duration of
    sunshine, from the sunset hour angle.

    Parameters
    ----------
    ws : float
        Sunset hour angle, in degrees.

    Returns
    -------
    N : float
        Daylight hours, in hours.
    """
    return 2 * ws / 15


def sunset_hour_angle(lat, declination):
    """
    Sunset hour angle (ωs, FAO-56 eq. 25), obtained by setting Z = 90° in
    the zenith angle equation. By symmetry it equals the sunrise hour
    angle in absolute value.

    Parameters
    ----------
    lat : float
        Latitude, in degrees.
    declination : float
        Solar declination of the day, in degrees.

    Returns
    -------
    ws : float
        Sunset hour angle, in degrees.
    """
    return acosd(-tand(lat) * tand(declination))


def inverse_relative_distance(doy):
    """
    Inverse relative distance Earth-Sun (dr, FAO-56 eq. 23), the
    eccentricity correction (d/D)² of the Earth's orbit.

    Parameters
    ----------
    doy : int
        Day of the year.

    Returns
    -------
    dr : float
        Inverse relative distance Earth-Sun, dimensionless.
    """
    return 1 + 0.033 * cosd(doy * 360 / 365)


def extraterrestrial_radiation(lat, declination, ws, dr):
    """
    Daily extraterrestrial radiation (Ra, FAO-56 eq. 21) on a horizontal
    surface at the top of the atmosphere.

    Parameters
    ----------
    lat : float
        Latitude, in degrees.
    declination : float
        Solar declination of the day, in degrees.
    ws : float
        Sunset hour angle, in degrees.
    dr : float
        Inverse relative distance Earth-Sun.

    Returns
    -------
    Ra : float
        Extraterrestrial radiation, in MJ m-2 day-1.
    """
    ws_rad = np.radians(ws)
    return 37.6 * dr * (ws_rad * sind(lat) * sind(declination) + cosd(lat) * cosd(declination) * sind(ws))


def sunshine_duration(N, Tmax, Tmin, lat, k=0.19):
    """
    Estimated actual duration of sunshine (n) from the daily temperature
    range (Tmax - Tmin) and the daylight hours, corrected by latitude.

    Parameters
    ----------
    N : float
        Daylight hours, in hours.
    Tmax, Tmin : float
        Daily maximum and minimum temperature, in °C.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    k : float, optional
        Empirical coefficient (0.16 for inland and 0.19 for coastal
        regions). Default 0.19.

    Returns
    -------
    n : float
        Estimated actual duration of sunshine, in hours.
    """
    return (N / 0.52) * (k * (Tmax - Tmin) ** 0.5 - 0.29 * cosd(lat))


def solar_radiation_angstrom(n, N, Ra, lat, b=0.52):
    """
    Solar radiation (Rs) by the Angstrom formula (FAO-56 eq. 35), in the
    Glover-McCulloch variant with a = 0.29 cos(lat).

    Parameters
    ----------
    n : float
        Actual duration of sunshine, in hours.
    N : float
        Daylight hours, in hours.
    Ra : float
        Extraterrestrial radiation, in MJ m-2 day-1.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    b : float, optional
        Empirical regression coefficient. Default 0.52.

    Returns
    -------
    Rs : float
        Solar radiation, in MJ m-2 day-1.
    """
    a = 0.29 * cosd(lat)
    return Ra * (a + b * n / N)


def relative_sunshine_duration(Rs, Ra, lat, b=0.52, a=None, clip=True):
    """
    Relative sunshine duration (n/N) from the Angstrom formula solved for
    n/N (the inverse of `solar_radiation_angstrom`), for when solar
    radiation is available (pyranometer or gridded data such as BR-DWGD)
    but sunshine duration is not.

    n/N = (Rs/Ra - a) / b

    Parameters
    ----------
    Rs : float or array_like
        Solar radiation, in MJ m-2 day-1.
    Ra : float or array_like
        Extraterrestrial radiation, in MJ m-2 day-1.
    lat : float
        Latitude, in degrees (negative in the southern hemisphere).
    b : float, optional
        Empirical regression coefficient. Default 0.52.
    a : float, optional
        Intercept. Default None, which uses a = 0.29 cos(lat) as in
        `solar_radiation_angstrom`.
    clip : bool, optional
        If True (default), clips n/N to 0–1. Values outside that range mean
        Rs/Ra beyond the reach of the formula for the coefficients used.

    Returns
    -------
    nN : float or array_like
        Relative sunshine duration, dimensionless.

    Examples
    --------
    >>> round(relative_sunshine_duration(22.5, 42.96, -24.86), 3)
    0.501
    >>> Rs = solar_radiation_angstrom(9.0, 12.0, 30.0, -22.0)
    >>> round(relative_sunshine_duration(Rs, 30.0, -22.0), 3)
    0.75
    """
    if a is None:
        a = 0.29 * cosd(lat)
    nN = (Rs / Ra - a) / b
    if clip:
        nN = np.clip(nN, 0.0, 1.0)
    return float(nN) if np.ndim(nN) == 0 else nN


def solar_radiation_hargreaves(Tmax, Tmin, Ra, k=0.16):
    """
    Solar radiation (Rs) by the Hargreaves radiation formula (FAO-56
    eq. 50), from the daily temperature range, for when sunshine data are
    not available.

    Parameters
    ----------
    Tmax, Tmin : float
        Daily maximum and minimum temperature, in °C.
    Ra : float
        Extraterrestrial radiation, in MJ m-2 day-1.
    k : float, optional
        Adjustment coefficient kRs (0.16 for inland and 0.19 for coastal
        locations). Default 0.16.

    Returns
    -------
    Rs : float
        Solar radiation, in MJ m-2 day-1.
    """
    return k * (Tmax - Tmin) ** 0.5 * Ra
