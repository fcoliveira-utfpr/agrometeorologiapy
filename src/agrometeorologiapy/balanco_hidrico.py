"""Capítulo 11 — Balanço hídrico."""

import numpy as np

__all__ = [
    "balanco_hidrico_climatologico",
    "balanco_hidrico_climatologico_grade",
    "balanco_hidrico_cultura",
]


def _ciclo_armazenamento(p_etp, CAD, arm_inicial):
    """Percorre os 12 meses uma vez a partir de `arm_inicial` (armazenamento
    no fim do mês anterior a janeiro). Retorna o ARM de cada mês."""
    arm = np.empty_like(p_etp)
    prev = arm_inicial
    for m in range(p_etp.shape[0]):
        pm = p_etp[m]
        # Secagem: ARM = CAD * exp(NEG.ACUM / CAD), equivalente a ARM_ant * exp((P-ETP) / CAD).
        # Reposição: ARM = min(ARM_ant + (P-ETP), CAD).
        prev = np.where(pm < 0, prev * np.exp(pm / CAD), np.minimum(prev + pm, CAD))
        arm[m] = prev
    return arm


def balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0, ciclico=True, tol=0.01, max_iter=100):
    """
    Balanço Hídrico Climatológico (BHC) de Thornthwaite & Mather (1955),
    vetorizado para muitos locais de uma vez (ex.: todos os pixels de um
    raster com a normal climatológica mensal).

    Parâmetros
    ----------
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais (mm/mês),
        janeiro a dezembro no primeiro eixo. As demais dimensões são livres
        (ex.: (12, linhas, colunas)).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo (mm), > 0. Escalar ou array
        com o formato das dimensões espaciais (broadcast). Padrão 100.0.
    ciclico : bool, opcional
        Se True (padrão), o armazenamento inicial é o de equilíbrio do ciclo
        anual: repete os 12 meses, usando o ARM de dezembro como ponto de
        partida de janeiro, até convergir. É o correto para uma normal
        climatológica, qualquer que seja o mês de início da estação seca.
        Se False, parte do solo cheio (ARM = CAD) antes de janeiro.
    tol : float, opcional
        Tolerância de convergência do ARM de dezembro (mm). Padrão 0.01.
    max_iter : int, opcional
        Número máximo de ciclos anuais. Padrão 100. Locais sem nenhum mês
        com P > ETP secam assintoticamente e podem não atingir `tol`; nesse
        caso o ARM já é desprezível (próximo de zero).

    Retorna
    -------
    dict de numpy.ndarray, formato (12, ...)
        'P-ETP', 'ARM', 'NEG.ACUM', 'ALT', 'ETR', 'DEF', 'EXC' (mm/mês;
        ARM e NEG.ACUM em mm).

    Exemplos
    --------
    >>> import numpy as np
    >>> P = np.array([200, 190, 190, 120, 15, 5, 2, 5, 20, 100, 220, 230.0])
    >>> ETP = np.full(12, 110.0)
    >>> bh = balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0)
    >>> bool(np.isclose(P.sum(), bh['ETR'].sum() + bh['EXC'].sum()))
    True
    """
    P = np.asarray(P, dtype=float)
    ETP = np.asarray(ETP, dtype=float)
    if P.shape != ETP.shape or P.shape[0] != 12:
        raise ValueError("P e ETP devem ter o mesmo formato, com 12 meses no primeiro eixo.")
    CAD = np.broadcast_to(np.asarray(CAD, dtype=float), P.shape[1:])
    if np.any(CAD <= 0):
        raise ValueError("CAD deve ser > 0 (use NaN para locais sem dado).")

    p_etp = P - ETP

    arm_inicial = CAD.copy()
    arm = _ciclo_armazenamento(p_etp, CAD, arm_inicial)
    if ciclico:
        for _ in range(max_iter):
            if not np.any(np.abs(arm[-1] - arm_inicial) > tol):  # NaN conta como convergido
                break
            arm_inicial = arm[-1]
            arm = _ciclo_armazenamento(p_etp, CAD, arm_inicial)

    arm_ant = np.concatenate([arm_inicial[np.newaxis], arm[:-1]])
    alt = arm - arm_ant
    with np.errstate(divide="ignore"):
        neg_acum = np.where(arm < CAD, CAD * np.log(arm / CAD), 0.0)
    etr = np.where(p_etp < 0, P + np.abs(alt), ETP)
    exc = np.where(p_etp > 0, np.maximum(p_etp - alt, 0.0), 0.0)

    return {
        "P-ETP": p_etp,
        "ARM": arm,
        "NEG.ACUM": neg_acum,
        "ALT": alt,
        "ETR": etr,
        "DEF": ETP - etr,
        "EXC": exc,
    }


def balanco_hidrico_climatologico(df, CAD=100.0, ciclico=True):
    """
    Calcula o Balanço Hídrico Climatológico (BHC), pelo método de
    Thornthwaite & Mather (1955), a partir de uma série de precipitação
    (P) e evapotranspiração potencial (ETP).

    Parâmetros
    ----------
    df : pandas.DataFrame
        12 linhas (janeiro a dezembro). Colunas obrigatórias: 'Meses',
        'P (mm/mês)', 'ETP (mm/mês)'.
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.0.
    ciclico : bool, opcional
        Se True (padrão), usa o armazenamento de equilíbrio do ciclo anual
        (o ARM de dezembro alimenta janeiro). Se False, parte do solo cheio
        antes de janeiro. Veja `balanco_hidrico_climatologico_grade`.

    Retorna
    -------
    df_bh : pandas.DataFrame
        Cópia do df de entrada, acrescido das colunas:
        'CAD', 'P-ETP', 'ARM (mm/mês)', 'NEG.ACUM (mm)', 'ALT (mm/mês)',
        'ETR (mm/mês)', 'DEF (mm/mês)', 'EXC (mm/mês)'.
    """
    df_bh = df.copy()
    bh = balanco_hidrico_climatologico_grade(
        df_bh['P (mm/mês)'].to_numpy(), df_bh['ETP (mm/mês)'].to_numpy(), CAD, ciclico=ciclico
    )

    df_bh['CAD'] = CAD
    df_bh['P-ETP'] = bh['P-ETP']
    df_bh['ARM (mm/mês)'] = bh['ARM']
    df_bh['NEG.ACUM (mm)'] = bh['NEG.ACUM']
    df_bh['ALT (mm/mês)'] = bh['ALT']
    df_bh['ETR (mm/mês)'] = bh['ETR']
    df_bh['DEF (mm/mês)'] = bh['DEF']
    df_bh['EXC (mm/mês)'] = bh['EXC']

    return df_bh


def balanco_hidrico_cultura(df):
    """
    Calcula o Balanço Hídrico de Cultura (BHc), pelo método de
    Thornthwaite & Mather, a partir de um df já estruturado com Chuva,
    ETc e CAD por período.

    Independente da escala temporal (diária, decendial, mensal etc.) —
    o usuário é responsável por pré-calcular 'ETc' (= Kc x ETo, já
    considerando a fase fenológica da cultura) e 'CAD' (= z x DTA, já
    considerando o avanço da profundidade radicular) na escala desejada;
    esta função só executa a contabilidade hídrica período a período.

    Parâmetros
    ----------
    df : pandas.DataFrame
        Colunas obrigatórias: 'Chuva' (mm/período), 'ETc' (mm/período) e
        'CAD' (mm, capacidade de água disponível no período), em ordem
        cronológica.

    Retorna
    -------
    df_bhc : pandas.DataFrame
        Cópia do df de entrada, acrescido de: 'P-ETc', 'ARM', 'ALT',
        'ETR', 'DEF', 'EXC', 'ISNA'.
    """
    df_bhc = df.copy()

    df_bhc['P-ETc'] = df_bhc['Chuva'] - df_bhc['ETc']

    PETc = df_bhc['P-ETc'].to_numpy()
    CAD = df_bhc['CAD'].to_numpy()

    ARM = [CAD[0]]  # solo cheio na CAD do primeiro período
    for p, cad in zip(PETc, CAD):
        prev = ARM[-1]
        if p < 0:
            ARM.append(prev * np.exp(p / cad))
        elif p + prev >= cad:
            ARM.append(cad)
        else:
            ARM.append(prev + p)
    ARM = ARM[1:]
    df_bhc['ARM'] = ARM

    ALT = [0] + list(np.array(ARM[1:]) - np.array(ARM[:-1]))
    df_bhc['ALT'] = ALT

    df_bhc['ETR'] = np.where(
        df_bhc['P-ETc'] < 0,
        df_bhc['Chuva'] + df_bhc['ALT'].abs(),
        df_bhc['ETc']
    )
    df_bhc['DEF'] = df_bhc['ETc'] - df_bhc['ETR']
    df_bhc['EXC'] = np.where(
        df_bhc['ARM'] < df_bhc['CAD'],
        0,
        df_bhc['P-ETc'] - df_bhc['ALT']
    )
    df_bhc['ISNA'] = df_bhc['ETR'] / df_bhc['ETc']

    return df_bhc
