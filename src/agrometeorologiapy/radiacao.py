"""Capítulo 5 — Radiação solar."""

from datetime import datetime

import numpy as np

from ._trig import acosd, cosd, sind, tand

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
    "Qg_hargreaves",
]


def nda(dia, mes, ano=2023):
    """
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
    """
    dt = datetime(ano, mes, dia)
    return dt.timetuple().tm_yday


def declinacao_solar(NDA):
    """
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
    """
    return 23.45 * sind(360 * (NDA - 80) / 365)


def angulo_horario(hora, minuto=0):
    """
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
    """
    hora_decimal = hora + minuto / 60
    return (hora_decimal - 12) * 15


def angulo_zenital(lat, declinacao):
    """
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
    """
    cos_Z = sind(lat) * sind(declinacao) + cosd(lat) * cosd(declinacao)
    cos_Z = max(-1.0, min(1.0, cos_Z))  # evita erro de domínio no acosd
    return acosd(cos_Z)


def azimute_solar(lat, declinacao, Z):
    """
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
    """
    num = sind(lat) * cosd(Z) - sind(declinacao)
    den = cosd(lat) * sind(Z)
    return acosd(num / den)


def comprimento_sombra(d, Z):
    """
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
    """
    return d * tand(Z)


def fotoperiodo(Hn):
    """
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
    """
    return 2 * Hn / 15


def angulo_horario_nascer(lat, declinacao):
    """
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
    """
    return acosd(-tand(lat) * tand(declinacao))


def fator_correcao_distancia(NDA):
    """
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
    """
    return 1 + 0.033 * cosd(NDA * 360 / 365)


def irradiancia_extraterrestre(lat, declinacao, Hn, dD2):
    """
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
    """
    hn_rad = np.radians(Hn)
    return 37.6 * dD2 * (hn_rad * sind(lat) * sind(declinacao) + cosd(lat) * cosd(declinacao) * sind(Hn))


def insolacao(N, Tmax, Tmin, lat, k=0.19):
    """
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
    """
    insol = (N / 0.52) * (k * (Tmax - Tmin) ** 0.5 - 0.29 * cosd(lat))
    return insol


def Qg_angstrom(insolacao, N, Qo, lat, b=0.52):
    """
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
    """
    a = 0.29 * cosd(lat)
    razao_insolacao = insolacao / N
    Qg = Qo * (a + b * razao_insolacao)
    return Qg


def Qg_hargreaves(Tmax, Tmin, Qo, k=0.16):
    """
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
    """
    Qg = k * (Tmax - Tmin) ** 0.5 * Qo
    return Qg
