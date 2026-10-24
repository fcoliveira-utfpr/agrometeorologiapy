"""Growing degree days."""

from calendar import monthrange
from datetime import date, timedelta

import numpy as np
import pandas as pd

__all__ = ["maturity_date", "sowing_date"]

_MONTHS = {
    1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June',
    7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 12: 'December',
}
_MONTHS_ABBR = {m: name[:3] for m, name in _MONTHS.items()}


def _record(month, day_date, Tmean, GDD, GDD_period, GDD_cycle):
    return {'month': _MONTHS_ABBR[month], 'date': day_date, 'Tmean': Tmean, 'GDD': round(GDD, 2),
            'GDD_period': round(GDD_period, 2), 'GDD_cycle': round(GDD_cycle, 2)}


def maturity_date(df, Tb, TT, sowing_day, sowing_month, interval='d', year=2023, verbose=True):
    """
    Physiological maturity date of a crop from its sowing date, by
    accumulating growing degree days (GDD) up to the thermal time
    requirement of the cycle (TT).

    Daily growing degree days (GDD):

    - if Tb < Tmin:  GDD = Tmean - Tb
    - if Tb >= Tmin: GDD = (Tmax - Tb)^2 / [2 (Tmax - Tmin)]

    Parameters
    ----------
    df : pandas.DataFrame
        Columns 'day', 'month', 'Tmean', 'Tmax', 'Tmin' (one row per period,
        in chronological order). Monthly: 'day' is only a marker (e.g. 1);
        dekadal: 'day' is the first day of the dekad (1, 11 or 21); daily:
        'day' is the actual day of the month.
    Tb : float
        Base temperature of the crop, in °C.
    TT : float
        Thermal time requirement of the cycle (thermal constant), in °C·day.
    sowing_day, sowing_month : int
        Day and month of sowing.
    interval : str, optional
        'd' (daily), 'dekad' (10-day period, always 10 days) or 'M'
        (monthly). Default 'd'.
    year : int, optional
        Reference year (default 2023, not a leap year).
    verbose : bool, optional
        If True (default), prints the sowing and maturity dates.

    Returns
    -------
    result : pandas.DataFrame
        One row per period: 'month' (abbreviated), 'date', 'Tmean' (°C),
        'GDD' (daily growing degree days, °C·day), 'GDD_period' (GDD x days
        counted in the period, °C·day) and 'GDD_cycle' (accumulated in the
        cycle, °C·day). From sowing to maturity; the last row is the
        maturity date and the GDD accumulated up to it. The year advances
        when the cycle crosses December -> January.
    """
    df = df.reset_index(drop=True)

    if interval == 'M':
        idx_ref = df.index[df['month'] == sowing_month]
    else:
        idx_ref = df.index[(df['month'] == sowing_month) & (df['day'] == sowing_day)]
    if len(idx_ref) == 0:
        raise ValueError(f"No row with day={sowing_day}, month={sowing_month} in df.")
    idx_ref = idx_ref[0]

    n_rows = len(df)
    order = [(idx_ref + i) % n_rows for i in range(n_rows)]

    records = []
    total = 0.0
    maturity = None
    row_year = year
    previous_month = None

    for k, i in enumerate(order):
        row = df.loc[i]
        Tmean, Tmax, Tmin = row['Tmean'], row['Tmax'], row['Tmin']
        row_month = int(row['month'])
        if previous_month is not None and row_month < previous_month:
            row_year += 1
        previous_month = row_month

        if Tb < Tmin:
            GDD = Tmean - Tb
        else:
            GDD = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if interval == 'd':
            n_period = 1
        elif interval == 'dekad':
            n_period = 10
        elif interval == 'M':
            n_period = monthrange(row_year, row_month)[1]
        else:
            raise ValueError("interval must be 'd', 'dekad' or 'M'")

        if k == 0 and interval == 'M':
            n_effective = n_period - sowing_day
        else:
            n_effective = n_period

        GDD_period = GDD * n_effective
        previous_total = total
        total += GDD_period

        if interval in ('d', 'dekad'):
            period_date = date(row_year, row_month, int(row['day']))
        else:
            period_date = date(row_year, row_month, sowing_day if k == 0 else 1)

        if total >= TT:
            if interval == 'd':
                maturity = period_date
            else:
                # Days of the period needed to complete TT; the last row
                # becomes the maturity date itself.
                missing = TT - previous_total
                days_needed = int(np.ceil(missing / GDD))
                days_needed = min(max(days_needed, 1), n_effective)
                # In the first month (monthly) counting starts the day after sowing.
                offset = days_needed if (k == 0 and interval == 'M') else days_needed - 1
                maturity = period_date + timedelta(days=offset)
                GDD_period = GDD * days_needed
                total = previous_total + GDD_period
            records.append(_record(row_month, maturity, Tmean, GDD, GDD_period, total))
            break

        records.append(_record(row_month, period_date, Tmean, GDD, GDD_period, total))

    if maturity is None:
        raise ValueError("The thermal sum of df (one full cycle) does not reach TT.")


    result = pd.DataFrame(records)
    if verbose:
        print(f"Sowing date: {_MONTHS[sowing_month]} {sowing_day}")
        print(f"Physiological maturity date: {_MONTHS[maturity.month]} {maturity.day}, {maturity.year}")
    return result


def sowing_date(df, Tb, TT, maturity_day, maturity_month, interval='d', year=2023, verbose=True):
    """
    Sowing date needed for a crop to reach physiological maturity on a
    given date (e.g. a target harvest date), by accumulating growing degree
    days (GDD) backwards up to the thermal time requirement (TT).

    Daily growing degree days (GDD):

    - if Tb < Tmin:  GDD = Tmean - Tb
    - if Tb >= Tmin: GDD = (Tmax - Tb)^2 / [2 (Tmax - Tmin)]

    Parameters
    ----------
    df : pandas.DataFrame
        Columns 'day', 'month', 'Tmean', 'Tmax', 'Tmin' (one row per period,
        in chronological order). Monthly: 'day' is only a marker (e.g. 1);
        dekadal: 'day' is the first day of the dekad (1, 11 or 21); daily:
        'day' is the actual day of the month.
    Tb : float
        Base temperature of the crop, in °C.
    TT : float
        Thermal time requirement of the cycle (thermal constant), in °C·day.
    maturity_day, maturity_month : int
        Day and month of maturity (known reference date).
    interval : str, optional
        'd' (daily), 'dekad' (10-day period, always 10 days) or 'M'
        (monthly). Default 'd'.
    year : int, optional
        Reference year (default 2023, not a leap year).
    verbose : bool, optional
        If True (default), prints the sowing and maturity dates.

    Returns
    -------
    result : pandas.DataFrame
        One row per period, with the columns of `maturity_date`, from
        maturity (reference) back to sowing; the last row is the sowing
        date and the GDD accumulated up to it. `year` is the year of
        maturity; it goes back when the cycle crosses January -> December.
    """
    df = df.reset_index(drop=True)

    if interval == 'M':
        idx_ref = df.index[df['month'] == maturity_month]
    else:
        idx_ref = df.index[(df['month'] == maturity_month) & (df['day'] == maturity_day)]
    if len(idx_ref) == 0:
        raise ValueError(f"No row with day={maturity_day}, month={maturity_month} in df.")
    idx_ref = idx_ref[0]

    n_rows = len(df)
    order = [(idx_ref - i) % n_rows for i in range(n_rows)]

    records = []
    total = 0.0
    sowing = None
    row_year = year
    previous_month = None

    for k, i in enumerate(order):
        row = df.loc[i]
        Tmean, Tmax, Tmin = row['Tmean'], row['Tmax'], row['Tmin']
        row_month = int(row['month'])
        if previous_month is not None and row_month > previous_month:
            row_year -= 1
        previous_month = row_month

        if Tb < Tmin:
            GDD = Tmean - Tb
        else:
            GDD = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if interval == 'd':
            n_period = 1
        elif interval == 'dekad':
            n_period = 10
        elif interval == 'M':
            n_period = monthrange(row_year, row_month)[1]
        else:
            raise ValueError("interval must be 'd', 'dekad' or 'M'")

        if k == 0 and interval == 'M':
            n_effective = maturity_day
        else:
            n_effective = n_period

        GDD_period = GDD * n_effective
        previous_total = total
        total += GDD_period

        if interval in ('d', 'dekad'):
            period_date = date(row_year, row_month, int(row['day']))
        else:
            period_date = date(row_year, row_month, maturity_day if k == 0 else 1)

        if total >= TT:
            if interval == 'd':
                sowing = period_date
            else:
                # Counts back from the last day of the period.
                missing = TT - previous_total
                days_needed = int(np.ceil(missing / GDD))
                days_needed = min(max(days_needed, 1), n_effective)
                if interval == 'M':
                    period_end = date(row_year, row_month, maturity_day if k == 0 else n_period)
                else:
                    period_end = period_date + timedelta(days=n_period - 1)
                # Monthly: as in maturity_date, the sowing day is not summed;
                # sowing is the day before the first day counted.
                back = days_needed if interval == 'M' else days_needed - 1
                sowing = period_end - timedelta(days=back)
                GDD_period = GDD * days_needed
                total = previous_total + GDD_period
            records.append(_record(row_month, sowing, Tmean, GDD, GDD_period, total))
            break

        records.append(_record(row_month, period_date, Tmean, GDD, GDD_period, total))

    if sowing is None:
        raise ValueError("The thermal sum of df (one full cycle) does not reach TT.")


    result = pd.DataFrame(records)
    if verbose:
        print(f"Maturity date (reference): {_MONTHS[maturity_month]} {maturity_day}")
        print(f"Required sowing date: {_MONTHS[sowing.month]} {sowing.day}, {sowing.year}")
    return result
