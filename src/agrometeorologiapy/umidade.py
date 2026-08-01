"""Capítulo 7 — Umidade do ar."""

import math

__all__ = [
    "es_tetens",
    "ea_umidade",
    "deficit_saturacao",
    "patm_altitude",
    "umidade_absoluta",
    "umidade_saturacao",
    "umidade_relativa",
    "ponto_orvalho",
    "constante_psicrometrica",
]


def es_tetens(T_ar):
    """
    Equação de Tetens: pressão de saturação de vapor d'água
    (es) em função da temperatura do ar.

    Expressão empírica que fornece a pressão máxima que o vapor d'água
    pode exercer no ar numa dada temperatura, antes de condensar.

    Parâmetros
    ----------
    T_ar : float
        Temperatura do ar, em °C.

    Retorna
    -------
    es : float
        Pressão de saturação de vapor, em kPa.
    """
    return 0.6108 * 10 ** ((7.5 * T_ar) / (237.3 + T_ar))


def ea_umidade(es, UR):
    """
    Pressão parcial (atual) de vapor d'água (ea) a partir da pressão de
    saturação (es) e da umidade relativa do ar (UR).

    Representa a pressão de vapor d'água realmente exercida na atmosfera,
    proporcional à fração de saturação indicada pela umidade relativa.

    Parâmetros
    ----------
    es : float
        Pressão de saturação de vapor d'água, em kPa.
    UR : float
        Umidade relativa do ar, em porcentagem (ex. 65%).

    Retorna
    -------
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.
    """
    return es * (UR / 100)


def deficit_saturacao(es, ea):
    """
    Déficit de saturação de vapor do ar (delta e).

    Mede o "quanto falta" para o ar atingir a saturação, sendo
    proporcional ao poder evaporante da atmosfera (equivalente ao VPD,
    vapor pressure deficit).

    Parâmetros
    ----------
    es : float
        Pressão de saturação de vapor, em kPa (eq. 7.2).
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.

    Retorna
    -------
    delta_e : float
        Déficit de saturação de vapor, em kPa.
    """
    return es - ea


def patm_altitude(A):
    """
    Pressão atmosférica (Patm) em função da altitude do local.

    Equação derivada da lei dos gases ideais, considerando o gradiente
    térmico padrão da atmosfera, para estimar a pressão atmosférica local
    a partir apenas da elevação (altitude) acima do nível do mar.

    Parâmetros
    ----------
    A : float
        Altitude do local, em metros.

    Retorna
    -------
    Patm : float
        Pressão atmosférica local, em kPa.
    """
    return 101.3 * ((293 - 0.0065 * A) / 293) ** 5.26


def umidade_absoluta(ea, T_ar_C):
    """
    Umidade absoluta do ar (UA): massa de vapor d'água por
    unidade de volume de ar.

    Deriva da equação de estado dos gases ideais aplicada ao vapor
    d'água; a constante 2168 vem da razão entre a massa molar da água e
    a constante universal dos gases.

    Parâmetros
    ----------
    ea : float
        Pressão parcial de vapor d'água, em kPa.
    T_ar_C : float
        Temperatura do ar, em °C (internamente convertida para Kelvin).

    Retorna
    -------
    UA : float
        Umidade absoluta, em g de H2O por m³ de ar.
    """
    T_K = T_ar_C + 273.15
    return 2168 * ea / T_K


def umidade_saturacao(es, T_ar_C):
    """
    Umidade de saturação do ar (US): massa máxima de vapor
    d'água que o ar pode reter por unidade de volume, na temperatura T.

    Calculada do mesmo modo que a umidade absoluta (eq. 7.8), porém
    usando a pressão de saturação (es) no lugar da pressão parcial (ea).

    Parâmetros
    ----------
    es : float
        Pressão de saturação de vapor, em kPa (eq. 7.2).
    T_ar_C : float
        Temperatura do ar, em °C.

    Retorna
    -------
    US : float
        Umidade de saturação, em g de H2O por m³ de ar.
    """
    T_K = T_ar_C + 273.15
    return 2168 * es / T_K


def umidade_relativa(ea, es):
    """
    Umidade relativa do ar (UR%), razão entre a pressão
    parcial e a pressão de saturação de vapor (equivalente à razão entre
    UA e US).

    Parâmetros
    ----------
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.
    es : float
        Pressão de saturação de vapor, em kPa.

    Retorna
    -------
    UR : float
        Umidade relativa do ar, em %.
    """
    return 100 * (ea / es)


def ponto_orvalho(ea):
    """
    Temperatura do ponto de orvalho (To): temperatura à qual
    o ar, mantendo o mesmo teor de vapor d'água, atingiria a saturação.

    Obtida invertendo-se algebricamente a equação de Tetens (eq. 7.2)
    para se obter To a partir de um valor conhecido de ea.

    Parâmetros
    ----------
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.

    Retorna
    -------
    To : float
        Temperatura do ponto de orvalho, em °C.
    """
    log_termo = math.log10(ea / 0.6108)
    return (237.3 * log_termo) / (7.5 - log_termo)


def constante_psicrometrica(Patm):
    """
    Constante psicrométrica (gamma).

    Relaciona a pressão parcial de vapor d'água com a temperatura do ar
    numa dada pressão atmosférica, sendo essencial para converter energia
    disponível em déficit de pressão de vapor equivalente.

    Parâmetros
    ----------
    Patm : float
        Pressão atmosférica local, em kPa (eq. 7.6, patm_altitude).

    Retorna
    -------
    gamma : float
        Constante psicrométrica, em kPa/°C.
    """
    return 0.665e-3 * Patm
