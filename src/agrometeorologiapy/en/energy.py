"""Radiation (energy) balance.

Sign convention: the net longwave radiation is returned as a negative value
(a net loss) and the net radiation is Rn = Rns + Rnl. FAO-56 writes Rnl as a
positive outgoing flux and Rn = Rns - Rnl; the result is the same.
"""

__all__ = ["net_radiation", "net_shortwave_radiation", "net_longwave_radiation"]


def net_radiation(Rns, Rnl):
    """
    Net radiation (Rn, FAO-56 eq. 40), the sum of the net shortwave and
    the (negative) net longwave radiation: the radiant energy available at
    the surface for physical and biological processes.

    Parameters
    ----------
    Rns : float
        Net shortwave radiation, in the unit of Rn.
    Rnl : float
        Net longwave radiation, negative (net loss), in the unit of Rn.

    Returns
    -------
    Rn : float
        Net radiation, in MJ m-2 day-1 (or W m-2).
    """
    return Rns + Rnl


def net_shortwave_radiation(Rs, albedo=0.25):
    """
    Net shortwave radiation (Rns, FAO-56 eq. 38): incoming solar radiation
    minus the part reflected by the surface.

    Parameters
    ----------
    Rs : float
        Incoming solar radiation (MJ m-2 day-1 or W m-2).
    albedo : float, optional
        Albedo (reflection coefficient) of the surface, 0–1. Default 0.25
        (grass; 0.23 for the FAO-56 hypothetical reference crop).

    Returns
    -------
    Rns : float
        Net shortwave radiation, in the unit of `Rs`.
    """
    return Rs * (1 - albedo)


def net_longwave_radiation(Tmax, Tmin, ea, Rs, Rso):
    """
    Net longwave radiation (Rnl, FAO-56 eq. 39): Stefan-Boltzmann law
    corrected for air humidity (via ea) and cloudiness (via Rs/Rso).

    Temperatures are taken in °C, as in the rest of the package, and
    converted to kelvin internally.

    Parameters
    ----------
    Tmax, Tmin : float
        Daily maximum and minimum temperature, in °C.
    ea : float
        Actual vapour pressure, in kPa.
    Rs : float
        Measured or estimated solar radiation, in MJ m-2 day-1.
    Rso : float
        Clear-sky solar radiation, in MJ m-2 day-1.

    Returns
    -------
    Rnl : float
        Net longwave radiation, in MJ m-2 day-1, negative (net loss of
        energy by terrestrial emission).
    """
    Tmax_K = Tmax + 273.15
    Tmin_K = Tmin + 273.15
    temperature_term = 4.903e-9 * ((Tmax_K ** 4 + Tmin_K ** 4) / 2)
    humidity_term = 0.34 - 0.14 * ea ** 0.5
    cloudiness_term = 1.35 * (Rs / Rso) - 0.35
    return -(temperature_term * humidity_term * cloudiness_term)
