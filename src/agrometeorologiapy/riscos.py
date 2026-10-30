"""Riscos agroclimáticos: geada, estresse térmico e molhamento foliar.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.climate_risk`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import rotulos, tradutor
from .en import climate_risk as _en

CATEGORIAS_ITU = {"Normal": "Normal", "Alert": "Alerta", "Danger": "Perigo", "Emergency": "Emergência"}

_traduzir = tradutor(__name__)

__all__ = [
    "probabilidade_geada",
    "periodo_retorno",
    "indice_temperatura_umidade",
    "categoria_itu",
    "dias_estresse_calor",
    "graus_dia_estresse_calor",
    "duracao_molhamento_foliar",
]


probabilidade_geada = _traduzir(
    _en.frost_probability, "probabilidade_geada", {'limiar': 'threshold', 'metodo': 'method'},
    valores={'method': {'empirico': 'empirical'}},
    doc="""
    Probabilidade de a temperatura mínima de um ano chegar a `limiar` ou
    menos, a partir de uma série de temperaturas mínimas absolutas anuais.

    - 'normal': as mínimas anuais seguem a distribuição normal,
      P = Phi[(limiar - média) / s], com o desvio padrão amostral s;
    - 'empirico': fração dos anos com Tmin <= limiar.

    Parâmetros
    ----------
    Tmin : array_like
        Temperatura mínima absoluta do ar de cada ano, em °C (NaN é
        ignorado). Para um risco sazonal, use as mínimas do período de
        interesse (ex.: do inverno, ou dos meses da florada).
    limiar : float, opcional
        Temperatura crítica, em °C. Padrão 3,0: a temperatura no abrigo
        meteorológico é alguns graus mais alta que junto à superfície, e
        uma mínima de cerca de 3 °C no abrigo é um critério usual de geada
        na relva. Para risco de dano, use a temperatura crítica da cultura.
    metodo : str, opcional
        'normal' (padrão) ou 'empirico'.

    Retorna
    -------
    probabilidade : float
        Probabilidade anual, 0-1.

    Exemplos
    --------
    >>> Tmin = [1.2, 4.5, -0.8, 3.1, 2.4, 5.6, 0.3, 6.2, 2.9, 4.0]
    >>> round(probabilidade_geada(Tmin), 3)
    0.511
    >>> probabilidade_geada(Tmin, metodo='empirico')
    0.5
    """)


periodo_retorno = _traduzir(
    _en.return_period, "periodo_retorno", {'probabilidade': 'probability'},
    doc="""
    Período de retorno de um evento com uma dada probabilidade anual de
    ocorrência: TR = 1 / probabilidade.

    Parâmetros
    ----------
    probabilidade : float ou array_like
        Probabilidade anual, 0-1.

    Retorna
    -------
    TR : float ou numpy.ndarray
        Período de retorno, em anos (inf para probabilidade 0).

    Exemplos
    --------
    >>> periodo_retorno(0.2)
    5.0
    """)


indice_temperatura_umidade = _traduzir(
    _en.temperature_humidity_index, "indice_temperatura_umidade", {'UR': 'RH'},
    doc="""
    Índice de temperatura e umidade, ITU (NRC, 1971), medida do estresse
    térmico de animais:

    ITU = (1,8 T + 32) - (0,55 - 0,0055 UR)(1,8 T - 26)

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura do ar, em °C.
    UR : float ou array_like
        Umidade relativa do ar, em % (ex.: 65).

    Retorna
    -------
    ITU : float ou numpy.ndarray
        Índice de temperatura e umidade (adimensional); ver `categoria_itu`.

    Exemplos
    --------
    >>> round(indice_temperatura_umidade(30, 70), 1)
    81.4
    """)


categoria_itu = _traduzir(
    _en.thi_category, "categoria_itu", {'ITU': 'THI'},
    pos=rotulos(CATEGORIAS_ITU),
    doc="""
    Classe de estresse térmico de um índice de temperatura e umidade
    (Livestock Conservation Institute, 1970): 'Normal' (ITU < 75), 'Alerta'
    (75 a 78), 'Perigo' (79 a 83) e 'Emergência' (>= 84).

    Parâmetros
    ----------
    ITU : float ou array_like
        Índice de temperatura e umidade.

    Retorna
    -------
    categoria : str ou numpy.ndarray de str
        Nome da classe ('' para NaN).

    Exemplos
    --------
    >>> categoria_itu(81.4)
    'Perigo'
    """)


dias_estresse_calor = _traduzir(
    _en.heat_stress_days, "dias_estresse_calor", {},
    doc="""
    Número de dias com temperatura máxima acima da temperatura crítica da
    cultura (ex.: para abortamento de flores).

    Parâmetros
    ----------
    Tmax : array_like
        Temperatura máxima diária do ar, em °C, com os dias no primeiro
        eixo (NaN não é contado).
    Tcrit : float
        Temperatura crítica, em °C.

    Retorna
    -------
    dias : int ou numpy.ndarray
        Número de dias com Tmax > Tcrit.

    Exemplos
    --------
    >>> dias_estresse_calor([31.0, 34.5, 36.0, 33.0], Tcrit=34)
    2
    """)


graus_dia_estresse_calor = _traduzir(
    _en.heat_stress_degree_days, "graus_dia_estresse_calor", {},
    doc="""
    Graus-dia de estresse por calor: soma do excesso da temperatura máxima
    sobre a temperatura crítica da cultura,

    GDE = soma de max(Tmax - Tcrit, 0)

    Parâmetros
    ----------
    Tmax : array_like
        Temperatura máxima diária do ar, em °C, com os dias no primeiro
        eixo (NaN é ignorado).
    Tcrit : float
        Temperatura crítica, em °C.

    Retorna
    -------
    GDE : float ou numpy.ndarray
        Graus-dia de estresse por calor, em °C·dia.

    Exemplos
    --------
    >>> graus_dia_estresse_calor([31.0, 34.5, 36.0, 33.0], Tcrit=34)
    2.5
    """)


duracao_molhamento_foliar = _traduzir(
    _en.leaf_wetness_duration, "duracao_molhamento_foliar",
    {'UR': 'RH', 'limiar': 'threshold', 'intervalo': 'interval'},
    doc="""
    Duração do período de molhamento foliar estimada pelo número de horas
    com umidade relativa igual ou acima de um limiar (NHUR; Sentelhas et
    al., 2008, que encontraram em 90 % um bom limiar para várias culturas
    e climas).

    Parâmetros
    ----------
    UR : array_like
        Umidade relativa do ar, em %, um valor por observação (ex.: os 24
        valores horários de um dia), com as observações no primeiro eixo.
        NaN não conta como molhado.
    limiar : float, opcional
        Umidade relativa a partir da qual a folha é considerada molhada, em
        %. Padrão 90.
    intervalo : float, opcional
        Tempo entre observações, em horas. Padrão 1 (dados horários).

    Retorna
    -------
    DPM : float ou numpy.ndarray
        Duração do período de molhamento, em horas.

    Exemplos
    --------
    >>> UR = [95, 96, 97, 93, 88, 75, 62, 55, 58, 70, 84, 91]
    >>> duracao_molhamento_foliar(UR, intervalo=2)
    10.0
    """)
