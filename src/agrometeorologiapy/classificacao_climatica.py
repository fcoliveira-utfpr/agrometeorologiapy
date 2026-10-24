"""Classificação climática: Köppen-Geiger, Thornthwaite (1948), Camargo (1991)
mod. Maluf (2000) e zonas de vida de Holdridge.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.climate_classification`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import climate_classification as _en

import numpy as np

from .en.climate_classification import _labels

# Códigos das classes (iguais nas duas APIs).
KOPPEN_CLASSES = _en.KOPPEN_CLASSES
TH_UMIDADE = _en.TH_MOISTURE
TH_SUBTIPO = _en.TH_SUBTYPE
TH_TERMICA = _en.TH_THERMAL
TH_CONCENTRACAO = _en.TH_CONCENTRATION
CAMARGO_TERMICA = _en.CAMARGO_THERMAL
CAMARGO_HIDRICA = _en.CAMARGO_MOISTURE
CAMARGO_ESTACAO_SECA = _en.CAMARGO_DRY_SEASON  # verão, outono, inverno, primavera
HOLDRIDGE_ZONAS_EN = _en.HOLDRIDGE_ZONES  # nomes originais (Jungkunst et al., 2021, Tab. 1)
FATOR_ETP_HOLDRIDGE = _en.HOLDRIDGE_PET_FACTOR  # mm/ano por °C de biotemperatura (Holdridge, 1967)

# Nomes em português.
TH_UMIDADE_NOMES = {1: "Superúmido", 2: "Úmido", 3: "Úmido", 4: "Úmido", 5: "Úmido",
                    6: "Subúmido", 7: "Subúmido seco", 8: "Semiárido", 9: "Árido"}
TH_TERMICA_NOMES = {1: "Megatérmico", 2: "Mesotérmico", 3: "Mesotérmico", 4: "Mesotérmico",
                    5: "Mesotérmico", 6: "Microtérmico", 7: "Microtérmico", 8: "Tundra",
                    9: "Gelo perpétuo"}
CAMARGO_TERMICA_NOMES = {1: 'Glacial', 2: 'Frio', 3: 'Frio moderado', 4: 'Temperado',
                         5: 'Subtemperado', 6: 'Subtropical', 7: 'Tropical', 8: 'Equatorial'}
CAMARGO_HIDRICA_NOMES = {1: 'Desértico', 2: 'Árido', 3: 'Seco', 4: 'Monçônico', 5: 'Subúmido',
                         6: 'Úmido', 7: 'Superúmido', 8: 'Extremamente úmido'}
HOLDRIDGE_ZONAS = {
    1: "Gelo polar", 2: "Deserto polar",
    3: "Tundra seca subpolar", 4: "Tundra úmida subpolar", 5: "Tundra muito úmida subpolar",
    6: "Tundra pluvial subpolar",
    7: "Deserto boreal", 8: "Arbustal seco boreal", 9: "Floresta úmida boreal",
    10: "Floresta muito úmida boreal", 11: "Floresta pluvial boreal",
    12: "Deserto temperado frio", 13: "Arbustal desértico temperado frio",
    14: "Estepe temperada fria", 15: "Floresta úmida temperada fria",
    16: "Floresta muito úmida temperada fria", 17: "Floresta pluvial temperada fria",
    18: "Deserto temperado quente", 19: "Arbustal desértico temperado quente",
    20: "Estepe espinhosa temperada quente", 21: "Floresta seca temperada quente",
    22: "Floresta úmida temperada quente", 23: "Floresta muito úmida temperada quente",
    24: "Floresta pluvial temperada quente",
    25: "Deserto subtropical", 26: "Arbustal desértico subtropical",
    27: "Estepe espinhosa subtropical", 28: "Floresta seca subtropical",
    29: "Floresta úmida subtropical", 30: "Floresta muito úmida subtropical",
    31: "Floresta pluvial subtropical",
    32: "Deserto tropical", 33: "Arbustal desértico tropical", 34: "Estepe espinhosa tropical",
    35: "Floresta muito seca tropical", 36: "Floresta seca tropical", 37: "Floresta úmida tropical",
    38: "Floresta muito úmida tropical",
}

_CHAVES = {'class': 'classe', 'description': 'descricao', 'moisture': 'umidade', 'subtype': 'subtipo',
           'thermal': 'termica', 'concentration': 'concentracao', 'PET_annual': 'ETP_anual',
           'DEF_annual': 'DEF_anual', 'SUR_annual': 'EXC_anual', 'PET_summer_pct': 'ETP_verao_pct'}
_CHAVES_CAMARGO = {'class': 'classe', 'description': 'descricao', 'thermal': 'termica',
                   'moisture': 'hidrica', 'dry_season': 'estacao_seca', 'T_annual': 'T_anual',
                   'T_coldest': 'T_mes_frio', 'DEF_annual': 'DEF_anual', 'SUR_annual': 'EXC_anual'}
_CHAVES_HOLDRIDGE = {'class': 'classe_en', 'biotemperature': 'biotemperatura', 'P_annual': 'P_anual',
                     'PET_annual': 'ETP_anual', 'PET_ratio': 'razao_ETP'}
_ETP_CAD = {'ETP': 'PET', 'CAD': 'AWC'}
_CAMARGO = {**_ETP_CAD, 'tabela_termica': 'thermal_table'}
_HOLDRIDGE = {'ETP': 'PET', 'limiar_correcao': 'correction_threshold'}
_TABELA_TERMICA = {'thermal_table': {'coerente': 'consistent'}}


def _descricao_thornthwaite(saida, kw):
    if saida['classe']:
        saida['descricao'] = (f"{TH_UMIDADE_NOMES[saida['umidade']]}, "
                              f"{TH_TERMICA_NOMES[saida['termica']].lower()}")
    return saida


def _descricao_camargo(saida, kw):
    if saida['classe']:
        saida['descricao'] = (f"{CAMARGO_TERMICA_NOMES[saida['termica']]} "
                              f"{CAMARGO_HIDRICA_NOMES[saida['hidrica']].lower()}")
    return saida


def _zona_holdridge(saida, kw):
    """Põe 'classe' (nome em português) antes de 'classe_en', como na API original."""
    zona = saida['id']
    if isinstance(zona, np.ndarray):
        classe = _labels(zona, HOLDRIDGE_ZONAS).astype(str)
    else:
        classe = HOLDRIDGE_ZONAS.get(zona, '')
    ordem = {'id': zona, 'classe': classe}
    ordem.update((k, v) for k, v in saida.items() if k != 'id')
    return ordem

_traduzir = tradutor(__name__)

__all__ = [
    "classificacao_koppen_grade",
    "classificacao_koppen",
    "classificacao_thornthwaite_grade",
    "classificacao_thornthwaite",
    "classificacao_camargo_grade",
    "classificacao_camargo",
    "classificacao_holdridge_grade",
    "classificacao_holdridge",
]


classificacao_koppen_grade = _traduzir(
    _en.koppen_classification_grid, "classificacao_koppen_grade", {},
    colunas_saida=_CHAVES,
    doc="""
    Classificação climática de Köppen-Geiger para muitos locais, pela
    metodologia de Alvares et al. (2013), com a sazonalidade f/s/w dos
    grupos C e D de Kottek et al. (2006).

    O verão é o semestre out-mar no hemisfério sul e abr-set no norte.

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C (janeiro a dezembro).
    P : array_like, formato (12, ...)
        Precipitação mensal, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (negativa no hemisfério sul); escalar ou com o
        formato das dimensões espaciais.

    Retorna
    -------
    dict de numpy.ndarray
        'id' (1-31, ver KOPPEN_CLASSES; 0 = sem dado) e 'classe' (ex.: 'Cfa').
    """)


classificacao_koppen = _traduzir(
    _en.koppen_classification, "classificacao_koppen", {},
    colunas_saida=_CHAVES,
    doc="""
    Classificação climática de Köppen-Geiger de um local (Alvares et al.,
    2013; Kottek et al., 2006). Ver `classificacao_koppen_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P : lista ou array de 12 valores
        Precipitação mensal (jan-dez), em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).

    Retorna
    -------
    dict
        'id' (int) e 'classe' (str, ex.: 'Cfa').
    """)


classificacao_thornthwaite_grade = _traduzir(
    _en.thornthwaite_classification_grid, "classificacao_thornthwaite_grade", _ETP_CAD,
    colunas_saida=_CHAVES,
    doc="""
    Classificação climática de Thornthwaite (1948) para muitos locais, a
    partir do balanço hídrico climatológico de Thornthwaite & Mather (1955).

    Tabelas de Aparecido et al. (2016), com três correções pelo original de
    Thornthwaite (1948): limite B'3/B'2 em 855 mm (o artigo traz 885); nos
    climas secos, s2/w2 a partir de Ih >= 20 (o artigo traz 33,3) e excedente
    maior no verão = w, no inverno = s (o artigo traz invertido).

    Índices (%): Ih = 100 EXC/ETP, Ia = 100 DEF/ETP, Im = Ih - 0,6 Ia.
    Concentração estival: ETP do verão (1/3 dez + jan + fev + 2/3 mar, no
    hemisfério sul) em % da ETP anual.

    Parâmetros
    ----------
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (define o hemisfério para verão e inverno).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.

    Retorna
    -------
    dict de numpy.ndarray
        'classe' (ex.: "B1rA'a'"), 'umidade', 'subtipo', 'termica',
        'concentracao' (códigos inteiros, 0 = sem dado), 'Ih', 'Ia', 'Im',
        'ETP_anual', 'DEF_anual', 'EXC_anual' (mm) e 'ETP_verao_pct'.
    """)


classificacao_thornthwaite = _traduzir(
    _en.thornthwaite_classification, "classificacao_thornthwaite", _ETP_CAD,
    colunas_saida=_CHAVES,
    pos=_descricao_thornthwaite,
    doc="""
    Classificação climática de Thornthwaite (1948) de um local. Ver
    `classificacao_thornthwaite_grade`.

    Parâmetros
    ----------
    P, ETP : lista ou array de 12 valores
        Precipitação e evapotranspiração potencial mensais (jan-dez), em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.

    Retorna
    -------
    dict
        'classe' (ex.: "B1rA'a'"), 'descricao', códigos e índices (Ih, Ia, Im...).
    """)


classificacao_camargo_grade = _traduzir(
    _en.camargo_classification_grid, "classificacao_camargo_grade", _CAMARGO,
    valores=_TABELA_TERMICA,
    colunas_saida=_CHAVES_CAMARGO,
    doc="""
    Classificação climática de Camargo (1991) modificada por Maluf (2000)
    para muitos locais (Tabelas 6-8 de Aparecido et al., 2016).

    Classe térmica pela temperatura média anual (Ty) e do mês mais frio
    (Tcold); classe hídrica pela deficiência (DEF) e pelo excedente (EXC)
    anuais do balanço hídrico climatológico; letra da estação com maior
    deficiência (v, o, i, p) nas classes SE, MO, SB e UM.

    A Tabela 6 do artigo é ambígua para Tcold > 20 °C. Com
    `tabela_termica='coerente'` (padrão), Tcold só subdivide a faixa
    18 < Ty <= 22 (STE <= 13 < ST <= 20 < TR) e Ty > 25 é sempre EQ. Com
    'literal', Tcold > 20 °C leva a TR em qualquer caso (como no script do
    GEE do repositório climas_brasil), e a classe EQ praticamente não ocorre.

    A Tabela 7 não cobre 0 < DEF <= 150 mm com EXC = 0 nem DEF = 0 com
    EXC <= 200 mm; esses casos entram em SB (0 <= DEF <= 150 e 0 <= EXC <= 200).

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C.
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (define o hemisfério das estações).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.
    tabela_termica : {'coerente', 'literal'}, opcional
        Leitura da Tabela 6 (ver acima). Padrão 'coerente'.

    Retorna
    -------
    dict de numpy.ndarray
        'classe' (ex.: 'ST-UMi'), 'termica', 'hidrica', 'estacao_seca'
        (códigos inteiros, 0 = sem dado/sem letra), 'T_anual', 'T_mes_frio',
        'DEF_anual' e 'EXC_anual'.
    """)


classificacao_camargo = _traduzir(
    _en.camargo_classification, "classificacao_camargo", _CAMARGO,
    valores=_TABELA_TERMICA,
    colunas_saida=_CHAVES_CAMARGO,
    pos=_descricao_camargo,
    doc="""
    Classificação climática de Camargo (1991) mod. Maluf (2000) de um local.
    Ver `classificacao_camargo_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P, ETP : lista ou array de 12 valores
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.
    tabela_termica : {'coerente', 'literal'}, opcional
        Leitura da Tabela 6 de Aparecido et al. (2016). Padrão 'coerente'.

    Retorna
    -------
    dict
        'classe' (ex.: 'ST-UMi'), 'descricao', códigos e variáveis anuais.
    """)


classificacao_holdridge_grade = _traduzir(
    _en.holdridge_classification_grid, "classificacao_holdridge_grade", _HOLDRIDGE,
    colunas_saida=_CHAVES_HOLDRIDGE,
    pos=_zona_holdridge,
    doc="""
    Zonas de vida de Holdridge (38 zonas) para muitos locais, pela
    biotemperatura e pela razão de evapotranspiração potencial (ETP/P).

    Biotemperatura: média das temperaturas mensais corrigidas pela
    latitude, t - 3|lat|/100 (t - 24)^2 nos meses com t > `limiar_correcao`,
    e limitadas a [0, 30] °C.

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C.
    P : array_like, formato (12, ...)
        Precipitação mensal, em mm/mês.
    lat : float ou array_like
        Latitude, em graus.
    ETP : array_like, formato (12, ...), opcional
        ETP mensal, em mm/mês (ex.: Penman-Monteith). Se None (padrão), usa
        a ETP de Holdridge: 58,93 x biotemperatura (mm/ano).
    limiar_correcao : float ou None, opcional
        Temperatura mensal (°C) acima da qual se aplica a correção por
        latitude. Padrão 24 (Holdridge, 1967). None aplica a todos os meses.

    Retorna
    -------
    dict de numpy.ndarray
        'id' (1-38, ver HOLDRIDGE_ZONAS; 0 = sem dado), 'classe' (nome da
        zona em português), 'classe_en' (nome original em inglês),
        'biotemperatura' (°C), 'P_anual' e 'ETP_anual' (mm) e
        'razao_ETP' (ETP/P).
    """)


classificacao_holdridge = _traduzir(
    _en.holdridge_classification, "classificacao_holdridge", _HOLDRIDGE,
    colunas_saida=_CHAVES_HOLDRIDGE,
    pos=_zona_holdridge,
    doc="""
    Zona de vida de Holdridge de um local. Ver `classificacao_holdridge_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P : lista ou array de 12 valores
        Precipitação mensal, em mm/mês.
    lat : float
        Latitude, em graus.
    ETP : lista ou array de 12 valores, opcional
        ETP mensal, em mm/mês. Se None, usa 58,93 x biotemperatura.
    limiar_correcao : float ou None, opcional
        Limiar (°C) da correção por latitude. Padrão 24.

    Retorna
    -------
    dict
        'id', 'classe' (nome da zona), 'classe_en', 'biotemperatura', 'P_anual',
        'ETP_anual' e 'razao_ETP'.
    """)
