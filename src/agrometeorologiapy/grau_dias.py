"""Capítulo 10 — Grau-dias.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.degree_days`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import degree_days as _en

_MESES_PT = {
    1: 'janeiro', 2: 'fevereiro', 3: 'março', 4: 'abril', 5: 'maio', 6: 'junho',
    7: 'julho', 8: 'agosto', 9: 'setembro', 10: 'outubro', 11: 'novembro', 12: 'dezembro',
}
_ABREV_EN_PT = {nome_en[:3]: nome_pt[:3] for nome_en, nome_pt in zip(
    ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
     'September', 'October', 'November', 'December'], _MESES_PT.values())}

_COLUNAS_ENTRADA = {'dia': 'day', 'mes': 'month', 'Tmed': 'Tmean'}
_COLUNAS_SAIDA = {'month': 'Meses', 'date': 'data', 'Tmean': 'Tmed', 'GDD': 'GDi',
                  'GDD_period': 'GDA_mes', 'GDD_cycle': 'GDA_ciclo'}
_VALORES = {'interval': {'dec': 'dekad'}}


def _meses_em_portugues(resultado):
    resultado['Meses'] = resultado['Meses'].map(_ABREV_EN_PT)
    return resultado


def _pos_maturacao(resultado, kw):
    resultado = _meses_em_portugues(resultado)
    fim = resultado['data'].iloc[-1]
    print(f"Data de semeadura: {kw['sowing_day']:02d} de {_MESES_PT[kw['sowing_month']]}")
    print(f"Data de maturação fisiológica: {fim.day:02d} de {_MESES_PT[fim.month]} de {fim.year}")
    return resultado


def _pos_semeadura(resultado, kw):
    resultado = _meses_em_portugues(resultado)
    sem = resultado['data'].iloc[-1]
    print(f"Data de maturação (referência): {kw['maturity_day']:02d} de {_MESES_PT[kw['maturity_month']]}")
    print(f"Data de semeadura necessária: {sem.day:02d} de {_MESES_PT[sem.month]} de {sem.year}")
    return resultado

_traduzir = tradutor(__name__)

__all__ = [
    "data_maturacao_fisiologica",
    "data_semeadura",
]


data_maturacao_fisiologica = _traduzir(
    _en.maturity_date, "data_maturacao_fisiologica",
    {'CT': 'TT', 'dia_semeadura': 'sowing_day', 'mes_semeadura': 'sowing_month',
     'intervalo': 'interval', 'ano': 'year'},
    valores=_VALORES,
    colunas_entrada=_COLUNAS_ENTRADA,
    colunas_saida=_COLUNAS_SAIDA,
    fixos={'verbose': False},
    pos=_pos_maturacao,
    doc="""
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
    """)


data_semeadura = _traduzir(
    _en.sowing_date, "data_semeadura",
    {'CT': 'TT', 'dia_maturacao': 'maturity_day', 'mes_maturacao': 'maturity_month',
     'intervalo': 'interval', 'ano': 'year'},
    valores=_VALORES,
    colunas_entrada=_COLUNAS_ENTRADA,
    colunas_saida=_COLUNAS_SAIDA,
    fixos={'verbose': False},
    pos=_pos_semeadura,
    doc="""
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
    """)
