"""Model evaluation: agreement between observed and predicted values.

Errors (MBE, MAE, RMSE), Pearson's correlation, Willmott's (1981) index of
agreement d, the confidence index c of Camargo & Sentelhas (1997), the
Nash-Sutcliffe efficiency (Nash & Sutcliffe, 1970) and the Kling-Gupta
efficiency (Gupta et al., 2009; Kling et al., 2012).

Every function takes the observed and the predicted (estimated, simulated)
values as `obs` and `pred`, of the same length; pairs with NaN in either
are left out.
"""

import numpy as np

__all__ = [
    "mbe",
    "mae",
    "rmse",
    "pearson_r",
    "willmott_d",
    "confidence_index",
    "confidence_index_category",
    "nse",
    "kge",
    "evaluate_model",
]

# Camargo & Sentelhas (1997): (exclusive lower bound, class).
CONFIDENCE_CATEGORIES = [
    (0.85, "Excellent"),
    (0.75, "Very good"),
    (0.65, "Good"),
    (0.60, "Fair"),
    (0.50, "Poor"),
    (0.40, "Bad"),
    (-np.inf, "Very bad"),
]


def _pairs(obs, pred):
    O = np.asarray(obs, dtype=float).ravel()
    P = np.asarray(pred, dtype=float).ravel()
    if O.shape != P.shape:
        raise ValueError("obs and pred must have the same length.")
    ok = ~(np.isnan(O) | np.isnan(P))
    if not ok.any():
        raise ValueError("No pair of obs and pred without NaN.")
    return O[ok], P[ok]


def mbe(obs, pred):
    """
    Mean bias error: MBE = mean(P - O). Positive when the model
    overestimates.

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    MBE : float
        In the unit of the data.

    Examples
    --------
    >>> round(mbe([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 3)
    0.0
    """
    O, P = _pairs(obs, pred)
    return float(np.mean(P - O))


def mae(obs, pred):
    """
    Mean absolute error: MAE = mean(|P - O|).

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    MAE : float
        In the unit of the data.

    Examples
    --------
    >>> mae([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0])
    0.5
    """
    O, P = _pairs(obs, pred)
    return float(np.mean(np.abs(P - O)))


def rmse(obs, pred):
    """
    Root mean square error: RMSE = sqrt[mean((P - O)^2)].

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    RMSE : float
        In the unit of the data.

    Examples
    --------
    >>> round(rmse([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.6124
    """
    O, P = _pairs(obs, pred)
    return float(np.sqrt(np.mean((P - O) ** 2)))


def pearson_r(obs, pred):
    """
    Pearson's correlation coefficient, r: the precision of the model (how
    close the points are to a straight line). r² is the coefficient of
    determination of the linear regression of P on O.

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    r : float
        Between -1 and 1.

    Examples
    --------
    >>> round(pearson_r([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.8485
    """
    O, P = _pairs(obs, pred)
    return float(np.corrcoef(O, P)[0, 1])


def willmott_d(obs, pred):
    """
    Index of agreement of Willmott (1981): the accuracy of the model (how
    close the points are to the 1:1 line),

    d = 1 - sum[(P - O)^2] / sum[(|P - mean(O)| + |O - mean(O)|)^2]

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    d : float
        Between 0 (no agreement) and 1 (perfect agreement).

    Examples
    --------
    >>> round(willmott_d([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.8889
    """
    O, P = _pairs(obs, pred)
    Om = O.mean()
    return float(1 - np.sum((P - O) ** 2) / np.sum((np.abs(P - Om) + np.abs(O - Om)) ** 2))


def confidence_index(obs, pred):
    """
    Confidence (performance) index of Camargo & Sentelhas (1997), c = r d:
    precision (Pearson's r) times accuracy (Willmott's d). See
    `confidence_index_category`.

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    c : float

    Examples
    --------
    >>> round(confidence_index([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.7542
    """
    return pearson_r(obs, pred) * willmott_d(obs, pred)


def confidence_index_category(c):
    """
    Performance class of the confidence index c (Camargo & Sentelhas,
    1997): 'Excellent' (c > 0.85), 'Very good' (0.76 to 0.85), 'Good' (0.66
    to 0.75), 'Fair' (0.61 to 0.65), 'Poor' (0.51 to 0.60), 'Bad' (0.41 to
    0.50) and 'Very bad' (<= 0.40). As in the table, c is first rounded
    to two decimals.

    Parameters
    ----------
    c : float or array_like
        Confidence index.

    Returns
    -------
    category : str or numpy.ndarray of str
        Class name ('' for NaN).

    Examples
    --------
    >>> confidence_index_category(0.754)
    'Good'
    """
    v = np.round(np.asarray(c, dtype=float), 2)
    out = np.select([v > lower for lower, _ in CONFIDENCE_CATEGORIES],
                    [name for _, name in CONFIDENCE_CATEGORIES], default='')
    return str(out) if out.ndim == 0 else out


def nse(obs, pred):
    """
    Nash-Sutcliffe efficiency (Nash & Sutcliffe, 1970):

    NSE = 1 - sum[(O - P)^2] / sum[(O - mean(O))^2]

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    NSE : float
        1 is a perfect model; 0, a model as good as the observed mean; < 0,
        worse than the mean.

    Examples
    --------
    >>> round(nse([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.7
    """
    O, P = _pairs(obs, pred)
    return float(1 - np.sum((O - P) ** 2) / np.sum((O - O.mean()) ** 2))


def kge(obs, pred, version=2009):
    """
    Kling-Gupta efficiency, which decomposes the error into correlation,
    variability and bias:

    KGE = 1 - sqrt[(r - 1)^2 + (a - 1)^2 + (b - 1)^2]

    with r Pearson's correlation, b = mean(P) / mean(O) and a the
    variability ratio: std(P) / std(O) in the original version (Gupta et
    al., 2009) or the ratio of the coefficients of variation, CV(P) / CV(O),
    in the modified version (Kling et al., 2012), which keeps bias and
    variability independent.

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.
    version : int, optional
        2009 (default) or 2012.

    Returns
    -------
    KGE : float
        1 is a perfect model; KGE > -0.41 is better than the observed mean
        (Knoben et al., 2019).

    Examples
    --------
    >>> round(kge([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.6703
    """
    O, P = _pairs(obs, pred)
    r = np.corrcoef(O, P)[0, 1]
    b = P.mean() / O.mean()
    if version == 2009:
        a = P.std() / O.std()
    elif version == 2012:
        a = (P.std() / P.mean()) / (O.std() / O.mean())
    else:
        raise ValueError("version must be 2009 or 2012")
    return float(1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + (b - 1) ** 2))


def evaluate_model(obs, pred):
    """
    All the statistics of this module for one model.

    Parameters
    ----------
    obs, pred : array_like
        Observed (O) and predicted (P) values.

    Returns
    -------
    dict
        'n' (pairs used), 'MBE', 'MAE', 'RMSE', 'r', 'd', 'c', 'NSE' and
        'KGE' (2009).

    Examples
    --------
    >>> s = evaluate_model([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0])
    >>> s['n'], round(s['RMSE'], 4), round(s['d'], 4)
    (4, 0.6124, 0.8889)
    """
    O, P = _pairs(obs, pred)
    return {'n': int(O.size), 'MBE': mbe(O, P), 'MAE': mae(O, P), 'RMSE': rmse(O, P),
            'r': pearson_r(O, P), 'd': willmott_d(O, P), 'c': confidence_index(O, P),
            'NSE': nse(O, P), 'KGE': kge(O, P)}
