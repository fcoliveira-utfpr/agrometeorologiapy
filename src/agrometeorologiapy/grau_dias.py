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
_MESES_ABREV = {m: nome[:3] for m, nome in _MESES_PT.items()}


def _registro(mes, data, Tmed, GDi, GDA_mes, GDA_ciclo):
    return {'Meses': _MESES_ABREV[mes], 'data': data, 'Tmed': Tmed, 'GDi': round(GDi, 2),
            'GDA_mes': round(GDA_mes, 2), 'GDA_ciclo': round(GDA_ciclo, 2)}


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
        Uma linha por período: 'Meses' (mês abreviado), 'data', 'Tmed' (°C),
        'GDi' (grau-dia diário, °C·dia), 'GDA_mes' (GDi x dias contados no
        período, °C·dia) e 'GDA_ciclo' (GDA acumulado no ciclo, °C·dia).
        Da semeadura até a maturação; a última linha é a data de maturação e o GDA acumulado até ela.
        O ano avança quando o ciclo cruza dezembro -> janeiro.
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
    ano_row = ano
    mes_anterior = None

    for pos, i in enumerate(ordem):
        row = df.loc[i]
        Tmed, Tmax, Tmin = row['Tmed'], row['Tmax'], row['Tmin']
        mes_row = int(row['mes'])
        if mes_anterior is not None and mes_row < mes_anterior:
            ano_row += 1
        mes_anterior = mes_row

        if Tb < Tmin:
            GDi = Tmed - Tb
        else:
            GDi = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if intervalo == 'd':
            n_periodo = 1
        elif intervalo == 'dec':
            n_periodo = 10
        elif intervalo == 'M':
            n_periodo = monthrange(ano_row, mes_row)[1]
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
            data_periodo = date(ano_row, mes_row, int(row['dia']))
        else:
            data_periodo = date(ano_row, mes_row, dia_semeadura if pos == 0 else 1)

        if acumulado >= CT:
            if intervalo == 'd':
                data_final = data_periodo
            else:
                # Dias do período necessários para completar CT; a última
                # linha passa a ser a própria data de maturação.
                faltante = CT - acumulado_anterior
                dias_necessarios = int(np.ceil(faltante / GDi))
                dias_necessarios = min(max(dias_necessarios, 1), n_efetivo)
                # No 1º mês (mensal) a contagem começa no dia seguinte à semeadura.
                deslocamento = dias_necessarios if (pos == 0 and intervalo == 'M') else dias_necessarios - 1
                data_final = data_periodo + timedelta(days=deslocamento)
                GD_periodo = GDi * dias_necessarios
                acumulado = acumulado_anterior + GD_periodo
            registros.append(_registro(mes_row, data_final, Tmed, GDi, GD_periodo, acumulado))
            break

        registros.append(_registro(mes_row, data_periodo, Tmed, GDi, GD_periodo, acumulado))

    if data_final is None:
        raise ValueError("A soma térmica do df (um ciclo completo) não atinge CT.")

    resultado = pd.DataFrame(registros)
    print(f"Data de semeadura: {dia_semeadura:02d} de {_MESES_PT[mes_semeadura]}")
    print(f"Data de maturação fisiológica: {data_final.day:02d} de {_MESES_PT[data_final.month]} de {data_final.year}")
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
        Uma linha por período: 'Meses' (mês abreviado), 'data', 'Tmed' (°C),
        'GDi' (grau-dia diário, °C·dia), 'GDA_mes' (GDi x dias contados no
        período, °C·dia) e 'GDA_ciclo' (GDA acumulado no ciclo, °C·dia).
        Da maturação (referência) até a semeadura; a última linha é a data de semeadura e o GDA
        acumulado até ela. `ano` é o ano da maturação; o ano recua quando o
        ciclo cruza janeiro -> dezembro.
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
    ano_row = ano
    mes_anterior = None

    for pos, i in enumerate(ordem):
        row = df.loc[i]
        Tmed, Tmax, Tmin = row['Tmed'], row['Tmax'], row['Tmin']
        mes_row = int(row['mes'])
        if mes_anterior is not None and mes_row > mes_anterior:
            ano_row -= 1
        mes_anterior = mes_row

        if Tb < Tmin:
            GDi = Tmed - Tb
        else:
            GDi = (Tmax - Tb) ** 2 / (2 * (Tmax - Tmin))

        if intervalo == 'd':
            n_periodo = 1
        elif intervalo == 'dec':
            n_periodo = 10
        elif intervalo == 'M':
            n_periodo = monthrange(ano_row, mes_row)[1]
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
            data_periodo = date(ano_row, mes_row, int(row['dia']))
        else:
            data_periodo = date(ano_row, mes_row, dia_maturacao if pos == 0 else 1)

        if acumulado >= CT:
            if intervalo == 'd':
                data_sem = data_periodo
            else:
                # Retrocede a partir do último dia do período.
                faltante = CT - acumulado_anterior
                dias_necessarios = int(np.ceil(faltante / GDi))
                dias_necessarios = min(max(dias_necessarios, 1), n_efetivo)
                if intervalo == 'M':
                    fim_periodo = date(ano_row, mes_row, dia_maturacao if pos == 0 else n_periodo)
                else:
                    fim_periodo = data_periodo + timedelta(days=n_periodo - 1)
                # No mensal, como em data_maturacao_fisiologica, o dia da semeadura
                # não entra na soma: a semeadura é a véspera do 1º dia contado.
                recuo = dias_necessarios if intervalo == 'M' else dias_necessarios - 1
                data_sem = fim_periodo - timedelta(days=recuo)
                GD_periodo = GDi * dias_necessarios
                acumulado = acumulado_anterior + GD_periodo
            registros.append(_registro(mes_row, data_sem, Tmed, GDi, GD_periodo, acumulado))
            break

        registros.append(_registro(mes_row, data_periodo, Tmed, GDi, GD_periodo, acumulado))

    if data_sem is None:
        raise ValueError("A soma térmica do df (um ciclo completo) não atinge CT.")

    resultado = pd.DataFrame(registros)
    print(f"Data de maturação (referência): {dia_maturacao:02d} de {_MESES_PT[mes_maturacao]}")
    print(f"Data de semeadura necessária: {data_sem.day:02d} de {_MESES_PT[data_sem.month]} de {data_sem.year}")
    return resultado
