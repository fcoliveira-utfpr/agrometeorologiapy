"""Soil water balance by Thornthwaite & Mather (1955).

Terms: STO (soil water storage), APWL (accumulated potential water loss),
dSTO (change in storage), ETa (actual evapotranspiration), DEF (water
deficit), SUR (water surplus), AWC (available water capacity) and WRSI
(water requirement satisfaction index, ETa/ETc; Frère & Popov, 1979).
"""

import numpy as np

__all__ = [
    "climatic_water_balance",
    "climatic_water_balance_grid",
    "crop_water_balance",
]


def _storage_cycle(p_pet, AWC, sto_initial):
    """Runs the 12 months once from `sto_initial` (storage at the end of the
    month before January). Returns the storage of each month."""
    sto = np.empty_like(p_pet)
    prev = sto_initial
    for m in range(p_pet.shape[0]):
        pm = p_pet[m]
        # Drying: STO = AWC * exp(APWL / AWC), equivalent to STO_prev * exp((P-PET) / AWC).
        # Recharge: STO = min(STO_prev + (P-PET), AWC).
        prev = np.where(pm < 0, prev * np.exp(pm / AWC), np.minimum(prev + pm, AWC))
        sto[m] = prev
    return sto


def climatic_water_balance_grid(P, PET, AWC=100.0, cyclic=True, tol=0.01, max_iter=100):
    """
    Climatic water balance of Thornthwaite & Mather (1955), vectorized for
    many locations at once (e.g. every pixel of a raster with the monthly
    climate normal).

    Parameters
    ----------
    P, PET : array_like, shape (12, ...)
        Monthly precipitation and potential evapotranspiration (mm/month),
        January to December on the first axis. The other dimensions are
        free (e.g. (12, rows, columns)).
    AWC : float or array_like, optional
        Available water capacity of the soil (mm), > 0. Scalar or array
        with the shape of the spatial dimensions (broadcast). Default 100.0.
    cyclic : bool, optional
        If True (default), the initial storage is the equilibrium of the
        annual cycle: the 12 months are repeated, with December's storage
        as the starting point of January, until convergence. This is the
        right choice for a climate normal, whatever month the dry season
        starts. If False, starts from a full soil (STO = AWC) before January.
    tol : float, optional
        Convergence tolerance of December's storage (mm). Default 0.01.
    max_iter : int, optional
        Maximum number of annual cycles. Default 100. Locations without any
        month with P > PET dry out asymptotically and may not reach `tol`;
        their storage is then already negligible (close to zero).

    Returns
    -------
    dict of numpy.ndarray, shape (12, ...)
        'P-PET', 'STO', 'APWL', 'dSTO', 'ETa', 'DEF', 'SUR' (mm/month; STO
        and APWL in mm).

    Examples
    --------
    >>> import numpy as np
    >>> P = np.array([200, 190, 190, 120, 15, 5, 2, 5, 20, 100, 220, 230.0])
    >>> PET = np.full(12, 110.0)
    >>> wb = climatic_water_balance_grid(P, PET, AWC=100.0)
    >>> bool(np.isclose(P.sum(), wb['ETa'].sum() + wb['SUR'].sum()))
    True
    """
    P = np.asarray(P, dtype=float)
    PET = np.asarray(PET, dtype=float)
    if P.shape != PET.shape or P.shape[0] != 12:
        raise ValueError("P and PET must have the same shape, with 12 months on the first axis.")
    AWC = np.broadcast_to(np.asarray(AWC, dtype=float), P.shape[1:])
    if np.any(AWC <= 0):
        raise ValueError("AWC must be > 0 (use NaN for locations without data).")

    p_pet = P - PET

    sto_initial = AWC.copy()
    sto = _storage_cycle(p_pet, AWC, sto_initial)
    if cyclic:
        for _ in range(max_iter):
            if not np.any(np.abs(sto[-1] - sto_initial) > tol):  # NaN counts as converged
                break
            sto_initial = sto[-1]
            sto = _storage_cycle(p_pet, AWC, sto_initial)

    sto_prev = np.concatenate([sto_initial[np.newaxis], sto[:-1]])
    dsto = sto - sto_prev
    with np.errstate(divide="ignore"):
        apwl = np.where(sto < AWC, AWC * np.log(sto / AWC), 0.0)
    eta = np.where(p_pet < 0, P + np.abs(dsto), PET)
    sur = np.where(p_pet > 0, np.maximum(p_pet - dsto, 0.0), 0.0)

    return {
        "P-PET": p_pet,
        "STO": sto,
        "APWL": apwl,
        "dSTO": dsto,
        "ETa": eta,
        "DEF": PET - eta,
        "SUR": sur,
    }


def climatic_water_balance(df, AWC=100.0, cyclic=True):
    """
    Climatic water balance by Thornthwaite & Mather (1955), from monthly
    precipitation (P) and potential evapotranspiration (PET).

    Parameters
    ----------
    df : pandas.DataFrame
        12 rows (January to December). Required columns: 'P (mm/month)' and
        'PET (mm/month)' (a 'month' column, if present, is kept).
    AWC : float, optional
        Available water capacity of the soil, in mm. Default 100.0.
    cyclic : bool, optional
        If True (default), uses the equilibrium storage of the annual cycle
        (December's storage feeds January). If False, starts from a full
        soil before January. See `climatic_water_balance_grid`.

    Returns
    -------
    df_wb : pandas.DataFrame
        Copy of the input with the added columns 'AWC', 'P-PET',
        'STO (mm/month)', 'APWL (mm)', 'dSTO (mm/month)', 'ETa (mm/month)',
        'DEF (mm/month)' and 'SUR (mm/month)'.
    """
    df_wb = df.copy()
    wb = climatic_water_balance_grid(
        df_wb['P (mm/month)'].to_numpy(), df_wb['PET (mm/month)'].to_numpy(), AWC, cyclic=cyclic
    )

    df_wb['AWC'] = AWC
    df_wb['P-PET'] = wb['P-PET']
    df_wb['STO (mm/month)'] = wb['STO']
    df_wb['APWL (mm)'] = wb['APWL']
    df_wb['dSTO (mm/month)'] = wb['dSTO']
    df_wb['ETa (mm/month)'] = wb['ETa']
    df_wb['DEF (mm/month)'] = wb['DEF']
    df_wb['SUR (mm/month)'] = wb['SUR']

    return df_wb


def crop_water_balance(df):
    """
    Crop water balance by Thornthwaite & Mather, from a DataFrame with
    precipitation, crop evapotranspiration and available water capacity
    per period.

    Independent of the time step (daily, dekadal, monthly...): the user
    computes 'ETc' (= Kc x ETo, for the crop phenological phase) and 'AWC'
    (= root depth x available water, following root growth) at the chosen
    time step; this function only does the period-by-period accounting.

    Parameters
    ----------
    df : pandas.DataFrame
        Required columns: 'P' (precipitation, mm/period), 'ETc' (mm/period)
        and 'AWC' (available water capacity in the period, mm), in
        chronological order.

    Returns
    -------
    df_wb : pandas.DataFrame
        Copy of the input with the added columns 'P-ETc', 'STO', 'dSTO',
        'ETa', 'DEF', 'SUR' and 'WRSI' (water requirement satisfaction
        index, ETa/ETc).
    """
    df_wb = df.copy()

    df_wb['P-ETc'] = df_wb['P'] - df_wb['ETc']

    p_etc = df_wb['P-ETc'].to_numpy()
    AWC = df_wb['AWC'].to_numpy()

    STO = [AWC[0]]  # soil full at the AWC of the first period
    for p, awc in zip(p_etc, AWC):
        prev = STO[-1]
        if p < 0:
            STO.append(prev * np.exp(p / awc))
        elif p + prev >= awc:
            STO.append(awc)
        else:
            STO.append(prev + p)
    STO = STO[1:]
    df_wb['STO'] = STO

    dSTO = [0] + list(np.array(STO[1:]) - np.array(STO[:-1]))
    df_wb['dSTO'] = dSTO

    df_wb['ETa'] = np.where(
        df_wb['P-ETc'] < 0,
        df_wb['P'] + df_wb['dSTO'].abs(),
        df_wb['ETc']
    )
    df_wb['DEF'] = df_wb['ETc'] - df_wb['ETa']
    df_wb['SUR'] = np.where(
        df_wb['STO'] < df_wb['AWC'],
        0,
        df_wb['P-ETc'] - df_wb['dSTO']
    )
    df_wb['WRSI'] = df_wb['ETa'] / df_wb['ETc']

    return df_wb
