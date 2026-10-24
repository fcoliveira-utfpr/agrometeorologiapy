"""Air temperature."""

__all__ = ["mean_temperature", "mean_temperature_observations"]


def mean_temperature(Tmax, Tmin):
    """
    Mean daily air temperature (Tmean, FAO-56 eq. 9) as the mean of the
    daily maximum and minimum.

    The simplest and most widely used method. It tends to overestimate the
    "true" mean, because the daily temperature course is not symmetric
    around the mean (it stays longer near the nighttime minimum than near
    the daytime maximum).

    Parameters
    ----------
    Tmax, Tmin : float
        Daily maximum and minimum temperature, in °C.

    Returns
    -------
    Tmean : float
        Mean daily air temperature, in °C.
    """
    return (Tmax + Tmin) / 2


def mean_temperature_observations(temperatures):
    """
    Mean temperature from automatic weather station observations (simple
    arithmetic mean of the observations in the period).

    The more observations, the closer the estimate is to the "true" value
    (the continuous integral of temperature over the day).

    Parameters
    ----------
    temperatures : list of float
        Temperature of each observation, in °C.

    Returns
    -------
    Tmean : float
        Mean temperature of the period, in °C.
    """
    return sum(temperatures) / len(temperatures)
