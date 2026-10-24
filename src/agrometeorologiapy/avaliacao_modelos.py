"""Avaliação de modelos: concordância entre valores observados e estimados.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.model_evaluation`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.

Todas as funções recebem os valores observados (`obs`) e os estimados
(`est`, previstos ou simulados), de mesmo tamanho; pares com NaN em
qualquer um deles são descartados.
"""

from ._i18n import rotulos, tradutor
from .en import model_evaluation as _en

CATEGORIAS_INDICE_CONFIANCA = {
    "Excellent": "Ótimo",
    "Very good": "Muito bom",
    "Good": "Bom",
    "Fair": "Mediano",
    "Poor": "Sofrível",
    "Bad": "Mau",
    "Very bad": "Péssimo",
}

_EST = {'est': 'pred'}

_traduzir = tradutor(__name__)

__all__ = [
    "mbe",
    "mae",
    "rmse",
    "r_pearson",
    "d_willmott",
    "indice_confianca",
    "categoria_indice_confianca",
    "nse",
    "kge",
    "avaliar_modelo",
]

_PARAMETROS = """
    Parâmetros
    ----------
    obs, est : array_like
        Valores observados (O) e estimados (E)."""


mbe = _traduzir(_en.mbe, "mbe", _EST, doc=f"""
    Erro médio (viés), MBE = média(E - O). Positivo quando o modelo
    superestima.
{_PARAMETROS}

    Retorna
    -------
    MBE : float
        Na unidade dos dados.

    Exemplos
    --------
    >>> round(mbe([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 3)
    0.0
    """)


mae = _traduzir(_en.mae, "mae", _EST, doc=f"""
    Erro absoluto médio, MAE = média(|E - O|).
{_PARAMETROS}

    Retorna
    -------
    MAE : float
        Na unidade dos dados.

    Exemplos
    --------
    >>> mae([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0])
    0.5
    """)


rmse = _traduzir(_en.rmse, "rmse", _EST, doc=f"""
    Raiz do erro quadrático médio, RMSE = raiz[média((E - O)^2)].
{_PARAMETROS}

    Retorna
    -------
    RMSE : float
        Na unidade dos dados.

    Exemplos
    --------
    >>> round(rmse([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.6124
    """)


r_pearson = _traduzir(_en.pearson_r, "r_pearson", _EST, doc=f"""
    Coeficiente de correlação de Pearson, r: a precisão do modelo (quão
    próximos de uma reta estão os pontos). r² é o coeficiente de
    determinação da regressão linear de E em O.
{_PARAMETROS}

    Retorna
    -------
    r : float
        Entre -1 e 1.

    Exemplos
    --------
    >>> round(r_pearson([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.8485
    """)


d_willmott = _traduzir(_en.willmott_d, "d_willmott", _EST, doc=f"""
    Índice de concordância de Willmott (1981): a exatidão do modelo (quão
    próximos da reta 1:1 estão os pontos),

    d = 1 - soma[(E - O)^2] / soma[(|E - média(O)| + |O - média(O)|)^2]
{_PARAMETROS}

    Retorna
    -------
    d : float
        Entre 0 (nenhuma concordância) e 1 (concordância perfeita).

    Exemplos
    --------
    >>> round(d_willmott([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.8889
    """)


indice_confianca = _traduzir(_en.confidence_index, "indice_confianca", _EST, doc=f"""
    Índice de confiança (desempenho) de Camargo & Sentelhas (1997), c = r d:
    precisão (r de Pearson) vezes exatidão (d de Willmott). Ver
    `categoria_indice_confianca`.
{_PARAMETROS}

    Retorna
    -------
    c : float

    Exemplos
    --------
    >>> round(indice_confianca([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.7542
    """)


categoria_indice_confianca = _traduzir(
    _en.confidence_index_category, "categoria_indice_confianca", {},
    pos=rotulos(CATEGORIAS_INDICE_CONFIANCA),
    doc="""
    Desempenho segundo o índice de confiança c (Camargo & Sentelhas, 1997):
    'Ótimo' (c > 0,85), 'Muito bom' (0,76 a 0,85), 'Bom' (0,66 a 0,75),
    'Mediano' (0,61 a 0,65), 'Sofrível' (0,51 a 0,60), 'Mau' (0,41 a 0,50)
    e 'Péssimo' (<= 0,40). Como na tabela, c é antes arredondado a duas
    casas.

    Parâmetros
    ----------
    c : float ou array_like
        Índice de confiança.

    Retorna
    -------
    categoria : str ou numpy.ndarray de str
        Nome da classe ('' para NaN).

    Exemplos
    --------
    >>> categoria_indice_confianca(0.754)
    'Bom'
    """)


nse = _traduzir(_en.nse, "nse", _EST, doc=f"""
    Eficiência de Nash-Sutcliffe (Nash & Sutcliffe, 1970):

    NSE = 1 - soma[(O - E)^2] / soma[(O - média(O))^2]
{_PARAMETROS}

    Retorna
    -------
    NSE : float
        1 é um modelo perfeito; 0, um modelo tão bom quanto a média
        observada; < 0, pior que a média.

    Exemplos
    --------
    >>> round(nse([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.7
    """)


kge = _traduzir(_en.kge, "kge", {**_EST, 'versao': 'version'}, doc=f"""
    Eficiência de Kling-Gupta, que decompõe o erro em correlação,
    variabilidade e viés:

    KGE = 1 - raiz[(r - 1)^2 + (a - 1)^2 + (b - 1)^2]

    com r a correlação de Pearson, b = média(E) / média(O) e a a razão de
    variabilidade: dp(E) / dp(O) na versão original (Gupta et al., 2009) ou
    a razão dos coeficientes de variação, CV(E) / CV(O), na versão
    modificada (Kling et al., 2012), que mantém viés e variabilidade
    independentes.
{_PARAMETROS}
    versao : int, opcional
        2009 (padrão) ou 2012.

    Retorna
    -------
    KGE : float
        1 é um modelo perfeito; KGE > -0,41 é melhor que a média observada
        (Knoben et al., 2019).

    Exemplos
    --------
    >>> round(kge([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0]), 4)
    0.6703
    """)


avaliar_modelo = _traduzir(_en.evaluate_model, "avaliar_modelo", _EST, doc=f"""
    Todas as estatísticas deste módulo para um modelo.
{_PARAMETROS}

    Retorna
    -------
    dict
        'n' (pares usados), 'MBE', 'MAE', 'RMSE', 'r', 'd', 'c', 'NSE' e
        'KGE' (2009).

    Exemplos
    --------
    >>> s = avaliar_modelo([1, 2, 3, 4], [1.5, 2.0, 3.5, 3.0])
    >>> s['n'], round(s['RMSE'], 4), round(s['d'], 4)
    (4, 0.6124, 0.8889)
    """)
