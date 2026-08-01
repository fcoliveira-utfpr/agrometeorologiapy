"""Capítulo 10 — Grau-dias."""

from calendar import monthrange
from datetime import date, timedelta

import numpy as np
import pandas as pd

__all__ = ["data_maturacao_fisiologica", "data_semeadura"]

_MESES_PT = {
    1: 'janeiro', 2: 'fevereiro', 3: 'março', 4: 'abril', 5: 'maio', 6: 'junho',
    7: 'julho', 8: 'agosto', 9: 'setembro', 10: 'outubro', 11: 'novembro', 12: 'dezembro',
}


def data_maturacao_fisiologica(df, Tb, CT, dia_semeadura, mes_semeadura, intervalo='d', ano=2023):
    """
    Calcula a data de maturação fisiológica de uma cultura, a partir da
    data de semeadura, por acúmulo de graus-dia (GDA) até atingir a
    constante térmica do ciclo (CT).

    Regra de cálculo do GD diário (GDi):
    - Se Tb < Tmin:  GDi = Tmed - Tb
    - Se Tb >= Tmin: GDi = (Tmax - Tb)^2 / [2*(Tmax - Tmin)]

    Parâmetros
    ----------
    df : pandas.DataFrame
        Colunas: 'dia', 'mes', 'Tmed', 'Tmax', 'Tmin' (uma linha por período,
        em ordem cronológica). No mensal, 'dia' é só um marcador (ex.: 1);
        no decendial, 'dia' é o dia de início do decêndio (1, 11 ou 21); no
        diário, 'dia' é o dia real do mês.
    Tb : float
        Temperatura base da cultura, em °C.
    CT : float
        Soma térmica total do ciclo (constante térmica), em °C·dia.
    dia_semeadura, mes_semeadura : int
        Dia e mês da semeadura.
    intervalo : str, opcional
        'd' (diário), 'M' (mensal) ou 'dec' (decendial, sempre 10 dias). Padrão 'd'.
    ano : int, opcional
        Ano de referência (padrão 2023, não-bissexto).

    Retorna
    -------
    resultado : pandas.DataFrame
        Colunas 'data' e 'GD_ciclo' (GDA acumulado), da semeadura até a maturação.
    """
    df = df.reset_index(drop=True)

    if intervalo == 'M':
        idx_ref = df.index[df['mes'] == mes_semeadura]
    else:
        idx_ref = df.index[(df['mes'] == mes_semeadura) & (df['dia'] == dia_semeadura)]
    if len(idx_ref) == 0:
        raise ValueError(f"Não encontrei a linha com dia={dia_semeadura}, mes={mes_semeadura} no df.")
    idx_ref = idx_ref[0]

    n_linhas = len(df)
    ordem = [(idx_ref + i) % n_linhas for i in range(n_linhas)]

    registros = []
    acumulado = 0.0
    data_final = None

    for pos, i in enumerate(ordem):
        row = df.loc[i]
        Tmed, Tmax, Tmin = row['Tmed'], row['Tmax'], row['Tmin']
        mes_row = int(row['mes'])

        if Tb < Tmin:
            GDi = Tmed - Tb
        else:
            GDi = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if intervalo == 'd':
            n_periodo = 1
        elif intervalo == 'dec':
            n_periodo = 10
        elif intervalo == 'M':
            n_periodo = monthrange(ano, mes_row)[1]
        else:
            raise ValueError("intervalo deve ser 'd', 'M' ou 'dec'")

        if pos == 0 and intervalo == 'M':
            n_efetivo = n_periodo - dia_semeadura
        else:
            n_efetivo = n_periodo

        GD_periodo = GDi * n_efetivo
        acumulado_anterior = acumulado
        acumulado += GD_periodo

        if intervalo in ('d', 'dec'):
            data_periodo = date(ano, mes_row, int(row['dia']))
        else:
            data_periodo = date(ano, mes_row, dia_semeadura if pos == 0 else 1)

        registros.append({'data': data_periodo, 'GD_ciclo': round(acumulado, 2)})

        if acumulado >= CT:
            if intervalo == 'd':
                data_final = data_periodo
            else:
                faltante = CT - acumulado_anterior
                dias_necessarios = int(np.ceil(faltante / GDi))
                dias_necessarios = min(dias_necessarios, n_efetivo)
                data_final = data_periodo + timedelta(days=dias_necessarios - 1)
            break

    resultado = pd.DataFrame(registros)
    print(f"Data de semeadura: {dia_semeadura:02d} de {_MESES_PT[mes_semeadura]}")
    print(f"Data de maturação fisiológica: {data_final.day:02d} de {_MESES_PT[data_final.month]}")
    return resultado


def data_semeadura(df, Tb, CT, dia_maturacao, mes_maturacao, intervalo='d', ano=2023):
    """
    Calcula a data de semeadura necessária para que uma cultura atinja a
    maturação fisiológica em uma data de referência conhecida (ex.: data
    de colheita desejada), por acúmulo retroativo de graus-dia (GDA) até
    a constante térmica do ciclo (CT).

    Regra de cálculo do GD diário (GDi):
    - Se Tb < Tmin:  GDi = Tmed - Tb
    - Se Tb >= Tmin: GDi = (Tmax - Tb)^2 / [2*(Tmax - Tmin)]

    Parâmetros
    ----------
    df : pandas.DataFrame
        Colunas: 'dia', 'mes', 'Tmed', 'Tmax', 'Tmin' (uma linha por período,
        em ordem cronológica). No mensal, 'dia' é só um marcador (ex.: 1);
        no decendial, 'dia' é o dia de início do decêndio (1, 11 ou 21); no
        diário, 'dia' é o dia real do mês.
    Tb : float
        Temperatura base da cultura, em °C.
    CT : float
        Soma térmica total do ciclo, em °C·dia.
    dia_maturacao, mes_maturacao : int
        Dia e mês da maturação (data de referência conhecida).
    intervalo : str, opcional
        'd' (diário), 'M' (mensal) ou 'dec' (decendial, sempre 10 dias). Padrão 'd'.
    ano : int, opcional
        Ano de referência (padrão 2023, não-bissexto).

    Retorna
    -------
    resultado : pandas.DataFrame
        Colunas 'data' e 'GD_ciclo' (GDA acumulado), da maturação (referência)
        até a semeadura.
    """
    df = df.reset_index(drop=True)

    if intervalo == 'M':
        idx_ref = df.index[df['mes'] == mes_maturacao]
    else:
        idx_ref = df.index[(df['mes'] == mes_maturacao) & (df['dia'] == dia_maturacao)]
    if len(idx_ref) == 0:
        raise ValueError(f"Não encontrei a linha com dia={dia_maturacao}, mes={mes_maturacao} no df.")
    idx_ref = idx_ref[0]

    n_linhas = len(df)
    ordem = [(idx_ref - i) % n_linhas for i in range(n_linhas)]

    registros = []
    acumulado = 0.0
    data_sem = None

    for pos, i in enumerate(ordem):
        row = df.loc[i]
        Tmed, Tmax, Tmin = row['Tmed'], row['Tmax'], row['Tmin']
        mes_row = int(row['mes'])

        if Tb < Tmin:
            GDi = Tmed - Tb
        else:
            GDi = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if intervalo == 'd':
            n_periodo = 1
        elif intervalo == 'dec':
            n_periodo = 10
        elif intervalo == 'M':
            n_periodo = monthrange(ano, mes_row)[1]
        else:
            raise ValueError("intervalo deve ser 'd', 'M' ou 'dec'")

        if pos == 0 and intervalo == 'M':
            n_efetivo = dia_maturacao
        else:
            n_efetivo = n_periodo

        GD_periodo = GDi * n_efetivo
        acumulado_anterior = acumulado
        acumulado += GD_periodo

        if intervalo in ('d', 'dec'):
            data_periodo = date(ano, mes_row, int(row['dia']))
        else:
            data_periodo = date(ano, mes_row, 1)

        registros.append({'data': data_periodo, 'GD_ciclo': round(acumulado, 2)})

        if acumulado >= CT:
            if intervalo == 'd':
                data_sem = data_periodo
            else:
                faltante = CT - acumulado_anterior
                dias_necessarios = int(np.ceil(faltante / GDi))
                dias_necessarios = min(dias_necessarios, n_efetivo)
                data_sem = data_periodo + timedelta(days=dias_necessarios - 1)
            break

    resultado = pd.DataFrame(registros)
    print(f"Data de maturação (referência): {dia_maturacao:02d} de {_MESES_PT[mes_maturacao]}")
    print(f"Data de semeadura necessária: {data_sem.day:02d} de {_MESES_PT[data_sem.month]}")
    return resultado
