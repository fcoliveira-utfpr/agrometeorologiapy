"""Capítulo 11 — Balanço hídrico."""

import numpy as np

__all__ = ["balanco_hidrico_climatologico", "balanco_hidrico_cultura"]


def balanco_hidrico_climatologico(df, CAD=100.0):
    """
    Calcula o Balanço Hídrico Climatológico (BHC), pelo método de
    Thornthwaite & Mather (1955), a partir de uma série de precipitação
    (P) e evapotranspiração potencial (ETP).

    Parâmetros
    ----------
    df : pandas.DataFrame
        Colunas obrigatórias: 'Meses', 'P (mm/mês)', 'ETP (mm/mês)'.
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.0.

    Retorna
    -------
    df_bh : pandas.DataFrame
        Cópia do df de entrada, acrescido das colunas:
        'CAD', 'P-ETP', 'ARM (mm/mês)', 'NEG.ACUM (mm)', 'ALT (mm/mês)',
        'ETR (mm/mês)', 'DEF (mm/mês)', 'EXC (mm/mês)'.
    """
    df_bh = df.copy()

    df_bh['CAD'] = CAD

    df_bh['P-ETP'] = df_bh['P (mm/mês)'] - df_bh['ETP (mm/mês)']

    ARM = [df_bh['CAD'].iloc[0]]      # solo cheio no início
    NEG_ACUM = [0.0]                  # sem déficit acumulado inicial

    for i in range(len(df_bh)):
        p_etp = df_bh['P-ETP'].iloc[i]
        cad = df_bh['CAD'].iloc[i]
        arm_prev = ARM[-1]
        neg_prev = NEG_ACUM[-1]

        if p_etp < 0:
            # Acumula déficit
            neg = neg_prev + p_etp
            arm = cad * np.exp(neg / cad)
        else:
            # Reposição hídrica
            arm = min(arm_prev + p_etp, cad)

            # Recalcula NEG.ACUM pela inversão
            if arm < cad:
                neg = cad * np.log(arm / cad)
            else:
                neg = 0.0

        ARM.append(arm)
        NEG_ACUM.append(neg)

    # Remove o valor inicial extra
    df_bh['ARM (mm/mês)'] = ARM[1:]
    df_bh['NEG.ACUM (mm)'] = NEG_ACUM[1:]

    df_bh['ALT (mm/mês)'] = df_bh['ARM (mm/mês)'].diff().fillna(0)

    df_bh['ETR (mm/mês)'] = np.where(
        df_bh['P-ETP'] < 0,
        df_bh['P (mm/mês)'] + df_bh['ALT (mm/mês)'].abs(),
        df_bh['ETP (mm/mês)']
    )

    df_bh['DEF (mm/mês)'] = df_bh['ETP (mm/mês)'] - df_bh['ETR (mm/mês)']

    df_bh['EXC (mm/mês)'] = np.where(
        (df_bh['P-ETP'] > 0) & (df_bh['ARM (mm/mês)'] == df_bh['CAD']),
        df_bh['P-ETP'] - df_bh['ALT (mm/mês)'],
        0
    )

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
