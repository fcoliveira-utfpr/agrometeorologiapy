"""Capítulo 7 — Umidade do ar.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.humidity`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import humidity as _en

_traduzir = tradutor(__name__)

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


es_tetens = _traduzir(_en.saturation_vapour_pressure, "es_tetens", {'T_ar': 'T'}, doc="""
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
    """)


ea_umidade = _traduzir(_en.actual_vapour_pressure, "ea_umidade", {'UR': 'RH'}, doc="""
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
    """)


deficit_saturacao = _traduzir(_en.vapour_pressure_deficit, "deficit_saturacao", {}, doc="""
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
    """)


patm_altitude = _traduzir(_en.atmospheric_pressure, "patm_altitude", {'A': 'z'}, doc="""
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
    """)


umidade_absoluta = _traduzir(_en.absolute_humidity, "umidade_absoluta", {'T_ar_C': 'T'}, doc="""
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
    """)


umidade_saturacao = _traduzir(
    _en.saturation_humidity, "umidade_saturacao", {'T_ar_C': 'T'},
    doc="""
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
    """)


umidade_relativa = _traduzir(_en.relative_humidity, "umidade_relativa", {}, doc="""
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
    """)


ponto_orvalho = _traduzir(_en.dewpoint_temperature, "ponto_orvalho", {}, doc="""
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
    """)


constante_psicrometrica = _traduzir(
    _en.psychrometric_constant, "constante_psicrometrica", {},
    doc="""
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
    """)
