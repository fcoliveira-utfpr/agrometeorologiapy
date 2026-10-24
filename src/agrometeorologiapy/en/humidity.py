"""Air humidity.

Terminology and symbols follow FAO-56 (Allen et al., 1998, ch. 3).
"""

import math

__all__ = [
    "saturation_vapour_pressure",
    "actual_vapour_pressure",
    "vapour_pressure_deficit",
    "atmospheric_pressure",
    "absolute_humidity",
    "saturation_humidity",
    "relative_humidity",
    "dewpoint_temperature",
    "psychrometric_constant",
]


def saturation_vapour_pressure(T):
    """
    Saturation vapour pressure, e°(T), by the Tetens equation (FAO-56
    eq. 11).

    Parameters
    ----------
    T : float
        Air temperature, in °C.

    Returns
    -------
    es : float
        Saturation vapour pressure, in kPa.
    """
    return 0.6108 * 10 ** ((7.5 * T) / (237.3 + T))


def actual_vapour_pressure(es, RH):
    """
    Actual vapour pressure (ea) from the saturation vapour pressure and the
    relative humidity (FAO-56 eq. 10 solved for ea).

    Parameters
    ----------
    es : float
        Saturation vapour pressure, in kPa.
    RH : float
        Relative humidity, in % (e.g. 65).

    Returns
    -------
    ea : float
        Actual vapour pressure, in kPa.
    """
    return es * (RH / 100)


def vapour_pressure_deficit(es, ea):
    """
    Vapour pressure deficit (es - ea): how far the air is from saturation,
    proportional to its evaporative demand.

    Parameters
    ----------
    es : float
        Saturation vapour pressure, in kPa.
    ea : float
        Actual vapour pressure, in kPa.

    Returns
    -------
    VPD : float
        Vapour pressure deficit, in kPa.
    """
    return es - ea


def atmospheric_pressure(z):
    """
    Atmospheric pressure (P) from the elevation (FAO-56 eq. 7), a
    simplification of the ideal gas law for a standard atmosphere at 20 °C.

    Parameters
    ----------
    z : float
        Elevation above sea level, in m.

    Returns
    -------
    Patm : float
        Atmospheric pressure, in kPa.
    """
    return 101.3 * ((293 - 0.0065 * z) / 293) ** 5.26


def absolute_humidity(ea, T):
    """
    Absolute humidity: mass of water vapour per unit volume of air.

    Parameters
    ----------
    ea : float
        Actual vapour pressure, in kPa.
    T : float
        Air temperature, in °C (converted to kelvin internally).

    Returns
    -------
    AH : float
        Absolute humidity, in g of H2O per m³ of air.
    """
    T_K = T + 273.15
    return 2168 * ea / T_K


def saturation_humidity(es, T):
    """
    Saturation (absolute) humidity: the maximum mass of water vapour per
    unit volume of air at a given temperature.

    Parameters
    ----------
    es : float
        Saturation vapour pressure, in kPa.
    T : float
        Air temperature, in °C.

    Returns
    -------
    SH : float
        Saturation humidity, in g of H2O per m³ of air.
    """
    T_K = T + 273.15
    return 2168 * es / T_K


def relative_humidity(ea, es):
    """
    Relative humidity (RH, FAO-56 eq. 10): ratio of the actual to the
    saturation vapour pressure.

    Parameters
    ----------
    ea : float
        Actual vapour pressure, in kPa.
    es : float
        Saturation vapour pressure, in kPa.

    Returns
    -------
    RH : float
        Relative humidity, in %.
    """
    return 100 * (ea / es)


def dewpoint_temperature(ea):
    """
    Dewpoint temperature (Tdew): the temperature to which the air must be
    cooled, at constant vapour content, to become saturated (inverse of the
    Tetens equation).

    Parameters
    ----------
    ea : float
        Actual vapour pressure, in kPa.

    Returns
    -------
    Tdew : float
        Dewpoint temperature, in °C.
    """
    x = math.log10(ea / 0.6108)
    return (237.3 * x) / (7.5 - x)


def psychrometric_constant(Patm):
    """
    Psychrometric constant (γ, FAO-56 eq. 8).

    Parameters
    ----------
    Patm : float
        Atmospheric pressure, in kPa (see `atmospheric_pressure`).

    Returns
    -------
    gamma : float
        Psychrometric constant, in kPa °C-1.
    """
    return 0.665e-3 * Patm
