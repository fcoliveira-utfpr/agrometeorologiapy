"""Agroclimatic risks: frost, heat stress and leaf wetness.

Frost risk from the distribution of the annual minimum temperatures; heat
stress of animals by the temperature-humidity index (NRC, 1971; LCI, 1970)
and of crops by days and degree days above a critical temperature; leaf
wetness duration by the number of hours with relative humidity above a
threshold (Sentelhas et al., 2008).
"""

from statistics import NormalDist

import numpy as np

__all__ = [
    "frost_probability",
    "return_period",
    "temperature_humidity_index",
    "thi_category",
    "heat_stress_days",
    "heat_stress_degree_days",
    "leaf_wetness_duration",
]

# Livestock Conservation Institute (1970): (lower bound, class).
THI_CATEGORIES = [
    (84.0, "Emergency"),
    (79.0, "Danger"),
    (75.0, "Alert"),
    (-np.inf, "Normal"),
]


def frost_probability(Tmin, threshold=3.0, method='normal'):
    """
    Probability that the minimum temperature of a year falls to `threshold`
    or below, from a series of annual absolute minimum temperatures.

    - 'normal': the annual minima are taken as normally distributed,
      P = Phi[(threshold - mean) / s], with the sample standard deviation s;
    - 'empirical': fraction of the years with Tmin <= threshold.

    Parameters
    ----------
    Tmin : array_like
        Absolute minimum air temperature of each year, in °C (NaN is
        ignored). Use the minima of the period of interest (e.g. of the
        winter, or of the flowering months) for a seasonal risk.
    threshold : float, optional
        Critical temperature, in °C. Default 3.0: temperatures in the
        meteorological shelter are a few degrees higher than at the surface,
        and a shelter minimum of about 3 °C is a usual criterion for frost
        at the grass level. Use the critical temperature of the crop for
        damage risk.
    method : str, optional
        'normal' (default) or 'empirical'.

    Returns
    -------
    probability : float
        Annual probability, 0-1.

    Examples
    --------
    >>> Tmin = [1.2, 4.5, -0.8, 3.1, 2.4, 5.6, 0.3, 6.2, 2.9, 4.0]
    >>> round(frost_probability(Tmin), 3)
    0.511
    >>> frost_probability(Tmin, method='empirical')
    0.5
    """
    T = np.asarray(Tmin, dtype=float)
    T = T[~np.isnan(T)]
    if method == 'normal':
        return NormalDist(T.mean(), T.std(ddof=1)).cdf(threshold)
    if method == 'empirical':
        return float(np.mean(T <= threshold))
    raise ValueError("method must be 'normal' or 'empirical'")


def return_period(probability):
    """
    Return period of an event with an annual probability of occurrence:
    TR = 1 / probability.

    Parameters
    ----------
    probability : float or array_like
        Annual probability, 0-1.

    Returns
    -------
    TR : float or numpy.ndarray
        Return period, in years (inf for probability 0).

    Examples
    --------
    >>> return_period(0.2)
    5.0
    """
    with np.errstate(divide='ignore'):
        TR = 1 / np.asarray(probability, dtype=float)
    return TR.item() if TR.ndim == 0 else TR


def temperature_humidity_index(T, RH):
    """
    Temperature-humidity index, THI (NRC, 1971), a measure of heat stress
    of livestock:

    THI = (1.8 T + 32) - (0.55 - 0.0055 RH)(1.8 T - 26)

    Parameters
    ----------
    T : float or array_like
        Air temperature, in °C.
    RH : float or array_like
        Relative humidity, in % (e.g. 65).

    Returns
    -------
    THI : float or numpy.ndarray
        Temperature-humidity index (dimensionless); see `thi_category`.

    Examples
    --------
    >>> round(temperature_humidity_index(30, 70), 1)
    81.4
    """
    return (1.8 * T + 32) - (0.55 - 0.0055 * RH) * (1.8 * T - 26)


def thi_category(THI):
    """
    Heat stress class of a temperature-humidity index (Livestock
    Conservation Institute, 1970): 'Normal' (THI < 75), 'Alert' (75 to
    78), 'Danger' (79 to 83) and 'Emergency' (>= 84).

    Parameters
    ----------
    THI : float or array_like
        Temperature-humidity index.

    Returns
    -------
    category : str or numpy.ndarray of str
        Class name ('' for NaN).

    Examples
    --------
    >>> thi_category(81.4)
    'Danger'
    """
    v = np.asarray(THI, dtype=float)
    out = np.select([v >= lower for lower, _ in THI_CATEGORIES],
                    [name for _, name in THI_CATEGORIES], default='')
    return str(out) if out.ndim == 0 else out


def heat_stress_days(Tmax, Tcrit):
    """
    Number of days with maximum temperature above the critical temperature
    of the crop (e.g. for flower abortion).

    Parameters
    ----------
    Tmax : array_like
        Daily maximum air temperature, in °C, days on the first axis (NaN
        is not counted).
    Tcrit : float
        Critical temperature, in °C.

    Returns
    -------
    days : int or numpy.ndarray
        Number of days with Tmax > Tcrit.

    Examples
    --------
    >>> heat_stress_days([31.0, 34.5, 36.0, 33.0], Tcrit=34)
    2
    """
    n = np.sum(np.asarray(Tmax, dtype=float) > Tcrit, axis=0)
    return int(n) if np.ndim(n) == 0 else n


def heat_stress_degree_days(Tmax, Tcrit):
    """
    Heat stress degree days: accumulated excess of the maximum temperature
    over the critical temperature of the crop,

    HDD = sum of max(Tmax - Tcrit, 0)

    Parameters
    ----------
    Tmax : array_like
        Daily maximum air temperature, in °C, days on the first axis (NaN
        is ignored).
    Tcrit : float
        Critical temperature, in °C.

    Returns
    -------
    HDD : float or numpy.ndarray
        Heat stress degree days, in °C·day.

    Examples
    --------
    >>> heat_stress_degree_days([31.0, 34.5, 36.0, 33.0], Tcrit=34)
    2.5
    """
    excess = np.clip(np.asarray(Tmax, dtype=float) - Tcrit, 0, None)
    HDD = np.nansum(excess, axis=0)
    return float(HDD) if np.ndim(HDD) == 0 else HDD


def leaf_wetness_duration(RH, threshold=90.0, interval=1.0):
    """
    Leaf wetness duration estimated by the number of hours with relative
    humidity at or above a threshold (NHRH; Sentelhas et al., 2008, who
    found 90 % a good threshold for several crops and climates).

    Parameters
    ----------
    RH : array_like
        Relative humidity, in %, one value per observation (e.g. the 24
        hourly values of a day), observations on the first axis. NaN is not
        counted as wet.
    threshold : float, optional
        Relative humidity from which the leaf is taken as wet, in %.
        Default 90.
    interval : float, optional
        Time between observations, in hours. Default 1 (hourly data).

    Returns
    -------
    LWD : float or numpy.ndarray
        Leaf wetness duration, in hours.

    Examples
    --------
    >>> RH = [95, 96, 97, 93, 88, 75, 62, 55, 58, 70, 84, 91]
    >>> leaf_wetness_duration(RH, interval=2)
    10.0
    """
    LWD = np.sum(np.asarray(RH, dtype=float) >= threshold, axis=0) * float(interval)
    return float(LWD) if np.ndim(LWD) == 0 else LWD
