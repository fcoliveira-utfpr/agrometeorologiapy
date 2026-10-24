"""Capítulo 5 — Radiação solar.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.radiation`; aqui ficam só os nomes, a tradução dos
parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import radiation as _en

_traduzir = tradutor(__name__)

__all__ = [
    "nda",
    "declinacao_solar",
    "angulo_horario",
    "angulo_zenital",
    "azimute_solar",
    "comprimento_sombra",
    "fotoperiodo",
    "angulo_horario_nascer",
    "fator_correcao_distancia",
    "irradiancia_extraterrestre",
    "insolacao",
    "Qg_angstrom",
    "relacao_n_N",
    "Qg_hargreaves",
]


nda = _traduzir(_en.day_of_year, "nda", {'dia': 'day', 'mes': 'month', 'ano': 'year'}, doc="""
    Retorna o Número do Dia do Ano (NDA).

    Parâmetros
    ----------
    dia : int
        Dia do mês (1–31).
    mes : int
        Mês (1–12).
    ano : int, opcional
        Ano. Padrão é 2023 (não bissexto).

    Retorna
    -------
    int
        Número do dia no ano (1 a 365).

    Exemplos
    --------
    >>> nda(1, 1)
    1
    >>> nda(25, 12)
    359
    >>> nda(29, 2, 2024)
    60
    """)


declinacao_solar = _traduzir(_en.solar_declination, "declinacao_solar", {'NDA': 'doy'}, doc="""
    Declinação solar (delta) para um dado dia do ano.

    Aproximação senoidal da variação anual da declinação solar, decorrente
    da inclinação do eixo de rotação da Terra (23,45°). Positiva quando o
    Sol está aparentemente no hemisfério norte, negativa no hemisfério sul.

    Parâmetros
    ----------
    NDA : int
        Número do dia do ano (1/jan = 1, 1/fev = 32, ..., 31/dez = 365).

    Retorna
    -------
    delta : float
        Declinação solar, em graus.
    """)


angulo_horario = _traduzir(
    _en.hour_angle, "angulo_horario", {'hora': 'hour', 'minuto': 'minute'},
    doc="""
    Converte a hora local (hora solar verdadeira) em ângulo horário.

    A Terra gira 360° em 24h, ou seja, 15° por hora. O ângulo horário é nulo
    (h = 0°) exatamente ao meio-dia local, quando o Sol cruza o meridiano do
    observador; é negativo pela manhã e positivo à tarde.

    Parâmetros
    ----------
    hora : int
        Hora do dia (0–23).
    minuto : int, opcional
        Minutos (0–59). Padrão é 0.

    Retorna
    -------
    h : float
        Ângulo horário, em graus.
    """)


angulo_zenital = _traduzir(
    _en.solar_zenith_angle, "angulo_zenital", {'declinacao': 'declination'},
    doc="""
    Ângulo zenital do Sol.

    Obtida da equação geral do ângulo zenital (eq. 5.4) fazendo-se h = 0°
    (meio-dia local). Representa a menor altura zenital (posição mais alta
    do Sol) atingida no dia.

    Parâmetros
    ----------
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    declinacao : float
        Declinação solar do dia, em graus.

    Retorna
    -------
    Z : float
        Ângulo zenital, em graus.
    """)


azimute_solar = _traduzir(
    _en.solar_azimuth, "azimute_solar", {'declinacao': 'declination', 'Z': 'zenith'},
    doc="""
    Azimute solar (alpha): direção horizontal do Sol em
    relação à linha Norte-Sul (referência Sul, no hemisfério sul).

    Parâmetros
    ----------
    lat : float
        Latitude do local, em graus.
    declinacao : float
        Declinação solar do dia, em graus.
    Z : float
        Ângulo zenital (eq. 5.4), em graus.

    Retorna
    -------
    alpha : float
        Azimute solar, em graus.
    """)


comprimento_sombra = _traduzir(
    _en.shadow_length, "comprimento_sombra", {'d': 'height', 'Z': 'zenith'},
    doc="""
    Comprimento da sombra (S) projetada por um objeto de altura d.

    Quanto maior o ângulo zenital (Sol mais baixo no horizonte), maior a
    tangente de Z e, portanto, mais longa a sombra.

    Parâmetros
    ----------
    d : float
        Altura do objeto (m).
    Z : float
        Ângulo zenital no instante considerado, em graus.

    Retorna
    -------
    S : float
        Comprimento da sombra, na mesma unidade de d.
    """)


fotoperiodo = _traduzir(_en.daylight_hours, "fotoperiodo", {'Hn': 'ws'}, doc="""
    Fotoperíodo (N), ou duração do dia, a partir do ângulo
    horário do nascer do Sol.

    Decorre da simetria da trajetória solar em relação ao meio-dia: o
    fotoperíodo é o dobro do ângulo horário do nascer, convertido de graus
    para horas (15°/hora).

    Parâmetros
    ----------
    Hn : float
        Ângulo horário no nascer do Sol, em graus.

    Retorna
    -------
    N : float
        Fotoperíodo, em horas.
    """)


angulo_horario_nascer = _traduzir(
    _en.sunset_hour_angle, "angulo_horario_nascer", {'declinacao': 'declination'},
    doc="""
    Ângulo horário no nascer do Sol (Hn).

    Obtida impondo-se Z = 90° (cos Z = 0) na equação geral do ângulo
    zenital, já que no nascer/pôr do Sol o astro está exatamente no
    horizonte.

    Parâmetros
    ----------
    lat : float
        Latitude do local, em graus.
    declinacao : float
        Declinação solar do dia, em graus.

    Retorna
    -------
    Hn : float
        Ângulo horário no nascer do Sol, em graus.
    """)


fator_correcao_distancia = _traduzir(
    _en.inverse_relative_distance, "fator_correcao_distancia", {'NDA': 'doy'},
    doc="""
    Fator de correção (d/D)^2 da excentricidade da órbita
    terrestre, para o dia do ano considerado.

    Corrige a constante solar em função da variação da distância real
    Terra-Sol (D) em torno da distância média (d = 1 UA) ao longo da
    órbita elíptica terrestre.

    Parâmetros
    ----------
    NDA : int
        Número do dia do ano.

    Retorna
    -------
    (d/D)^2 : float
        Fator de correção adimensional.
    """)


irradiancia_extraterrestre = _traduzir(
    _en.extraterrestrial_radiation, "irradiancia_extraterrestre",
    {'declinacao': 'declination', 'Hn': 'ws', 'dD2': 'dr'},
    doc="""
    Irradiância solar global extraterrestre diária (Qo), no
    topo da atmosfera, para uma superfície horizontal.

    Representa o total diário máximo de energia solar teórica que
    incidiria sobre uma superfície horizontal caso não houvesse atenuação
    atmosférica. Depende apenas da latitude, da declinação solar e da
    correção pela distância Terra-Sol do dia.

    Parâmetros
    ----------
    lat : float
        Latitude do local, em graus.
    declinacao : float
        Declinação solar do dia, em graus.
    Hn : float
        Ângulo horário no nascer do Sol, em graus (eq. 5.20).
    dD2 : float
        Fator de correção (d/D)^2 da distância Terra-Sol (eq. 5.30).

    Retorna
    -------
    Qo : float
        Irradiância solar global extraterrestre diária, em MJ/m² dia.
    """)


insolacao = _traduzir(_en.sunshine_duration, "insolacao", {}, doc="""
    Estimativa do número de horas de insolação (brilho solar), insol.

    Fórmula empírica baseada na amplitude térmica diária (Tmax - Tmin) e na
    duração astronômica do dia (N), com correção pela latitude local.

    Parâmetros
    ----------
    N : float
        Duração astronômica do dia (fotoperíodo), em horas.
    k : float, opcional
        Coeficiente empírico de ajuste (padrão 0,16 para regiões interioranas;
        0,19 é comumente usado para regiões costeiras).
    Tmax : float
        Temperatura máxima diária, em °C.
    Tmin : float
        Temperatura mínima diária, em °C.
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).

    Retorna
    -------
    insol : float
        Número de horas de insolação (brilho solar) estimado, em horas.
    """)


Qg_angstrom = _traduzir(
    _en.solar_radiation_angstrom, "Qg_angstrom", {'insolacao': 'n', 'Qo': 'Ra'},
    doc="""
    Equação de Angström-Prescott (variante de Glover-McCulloch): estimativa
    da radiação solar global a partir da razão de insolação n/N.

    Parâmetros
    ----------
    insolacao : float
        Número de horas de brilho solar (insolação) medido no dia, em horas.
    N : float
        Fotoperíodo do dia (número máximo de horas de brilho solar), em horas.
    Qo : float
        Irradiância solar no topo da atmosfera (radiação extraterrestre), em MJ/m².
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    b : float, opcional
        Coeficiente empírico de regressão (padrão 0,52).

    Retorna
    -------
    Qg : float
        Irradiância solar global, em MJ/m².
    """)


relacao_n_N = _traduzir(
    _en.relative_sunshine_duration, "relacao_n_N", {'Qg': 'Rs', 'Qo': 'Ra', 'limitar': 'clip'},
    doc="""
    Razão de insolação n/N pela equação de Angström-Prescott (variante de
    Glover-McCulloch) invertida: a mesma fórmula de `Qg_angstrom`, resolvida
    para n/N a partir da radiação solar global.

    Útil quando há Qg (piranômetro ou bases em grade, como o BR-DWGD), mas
    não há heliógrafo.

    Parâmetros
    ----------
    Qg : float ou array_like
        Irradiância solar global, em MJ/m² dia.
    Qo : float ou array_like
        Irradiância solar no topo da atmosfera (radiação extraterrestre), em MJ/m² dia.
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    b : float, opcional
        Coeficiente empírico de regressão (padrão 0,52).
    a : float, opcional
        Coeficiente linear. Padrão None, que usa a = 0,29 cos(lat), como em
        `Qg_angstrom`.
    limitar : bool, opcional
        Se True (padrão), limita n/N entre 0 e 1. Valores fora dessa faixa
        indicam Qg/Qo fora do alcance da equação (dia muito nublado ou
        muito limpo para os coeficientes usados).

    Retorna
    -------
    nN : float ou array_like
        Razão de insolação n/N, adimensional.

    Exemplos
    --------
    >>> round(relacao_n_N(22.5, 42.96, -24.86), 3)
    0.501
    >>> Qg = Qg_angstrom(9.0, 12.0, 30.0, -22.0)
    >>> round(relacao_n_N(Qg, 30.0, -22.0), 3)
    0.75
    """)


Qg_hargreaves = _traduzir(_en.solar_radiation_hargreaves, "Qg_hargreaves", {'Qo': 'Ra'}, doc="""
    Equação de Hargreaves-Samani: estimativa da radiação solar global a partir da
    amplitude térmica diária (Tmax - Tmin) e da radiação extraterrestre (Qo).

    Método útil quando não há dados de insolação (n) disponíveis, exigindo
    apenas temperaturas máxima e mínima diárias.

    Parâmetros
    ----------
    Tmax : float
        Temperatura máxima diária, em °C.
    Tmin : float
        Temperatura mínima diária, em °C.
    Qo : float
        Irradiância solar no topo da atmosfera (radiação extraterrestre), em MJ/m².
    k : float, opcional
        Coeficiente empírico de ajuste (padrão 0,16 para regiões interioranas;
        0,19 é comumente usado para regiões costeiras).

    Retorna
    -------
    Qg : float
        Irradiância solar global, em MJ/m².
    """)
