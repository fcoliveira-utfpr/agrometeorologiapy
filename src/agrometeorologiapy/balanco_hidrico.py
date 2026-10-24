"""Capítulo 11 — Balanço hídrico.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.water_balance`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import water_balance as _en

_CHAVES_GRADE = {'P-PET': 'P-ETP', 'STO': 'ARM', 'APWL': 'NEG.ACUM', 'dSTO': 'ALT', 'ETa': 'ETR',
                 'SUR': 'EXC'}
_ENTRADA_BHC = {'Meses': 'month', 'P (mm/mês)': 'P (mm/month)', 'ETP (mm/mês)': 'PET (mm/month)'}
_SAIDA_BHC = {'month': 'Meses', 'P (mm/month)': 'P (mm/mês)', 'PET (mm/month)': 'ETP (mm/mês)',
              'AWC': 'CAD',
              'P-PET': 'P-ETP', 'STO (mm/month)': 'ARM (mm/mês)', 'APWL (mm)': 'NEG.ACUM (mm)',
              'dSTO (mm/month)': 'ALT (mm/mês)', 'ETa (mm/month)': 'ETR (mm/mês)',
              'DEF (mm/month)': 'DEF (mm/mês)', 'SUR (mm/month)': 'EXC (mm/mês)'}
_ENTRADA_CULTURA = {'Chuva': 'P', 'CAD': 'AWC'}
_SAIDA_CULTURA = {'P': 'Chuva', 'AWC': 'CAD', 'STO': 'ARM', 'dSTO': 'ALT', 'ETa': 'ETR', 'SUR': 'EXC',
                  'WRSI': 'ISNA'}

_traduzir = tradutor(__name__)

__all__ = [
    "balanco_hidrico_climatologico_grade",
    "balanco_hidrico_climatologico",
    "balanco_hidrico_cultura",
]


balanco_hidrico_climatologico_grade = _traduzir(
    _en.climatic_water_balance_grid, "balanco_hidrico_climatologico_grade",
    {'ETP': 'PET', 'CAD': 'AWC', 'ciclico': 'cyclic'},
    colunas_saida=_CHAVES_GRADE,
    doc="""
    Balanço Hídrico Climatológico (BHC) de Thornthwaite & Mather (1955),
    vetorizado para muitos locais de uma vez (ex.: todos os pixels de um
    raster com a normal climatológica mensal).

    Parâmetros
    ----------
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais (mm/mês),
        janeiro a dezembro no primeiro eixo. As demais dimensões são livres
        (ex.: (12, linhas, colunas)).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo (mm), > 0. Escalar ou array
        com o formato das dimensões espaciais (broadcast). Padrão 100.0.
    ciclico : bool, opcional
        Se True (padrão), o armazenamento inicial é o de equilíbrio do ciclo
        anual: repete os 12 meses, usando o ARM de dezembro como ponto de
        partida de janeiro, até convergir. É o correto para uma normal
        climatológica, qualquer que seja o mês de início da estação seca.
        Se False, parte do solo cheio (ARM = CAD) antes de janeiro.
    tol : float, opcional
        Tolerância de convergência do ARM de dezembro (mm). Padrão 0.01.
    max_iter : int, opcional
        Número máximo de ciclos anuais. Padrão 100. Locais sem nenhum mês
        com P > ETP secam assintoticamente e podem não atingir `tol`; nesse
        caso o ARM já é desprezível (próximo de zero).

    Retorna
    -------
    dict de numpy.ndarray, formato (12, ...)
        'P-ETP', 'ARM', 'NEG.ACUM', 'ALT', 'ETR', 'DEF', 'EXC' (mm/mês;
        ARM e NEG.ACUM em mm).

    Exemplos
    --------
    >>> import numpy as np
    >>> P = np.array([200, 190, 190, 120, 15, 5, 2, 5, 20, 100, 220, 230.0])
    >>> ETP = np.full(12, 110.0)
    >>> bh = balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0)
    >>> bool(np.isclose(P.sum(), bh['ETR'].sum() + bh['EXC'].sum()))
    True
    """)


balanco_hidrico_climatologico = _traduzir(
    _en.climatic_water_balance, "balanco_hidrico_climatologico",
    {'CAD': 'AWC', 'ciclico': 'cyclic'},
    colunas_entrada=_ENTRADA_BHC,
    colunas_saida=_SAIDA_BHC,
    doc="""
    Calcula o Balanço Hídrico Climatológico (BHC), pelo método de
    Thornthwaite & Mather (1955), a partir de uma série de precipitação
    (P) e evapotranspiração potencial (ETP).

    Parâmetros
    ----------
    df : pandas.DataFrame
        12 linhas (janeiro a dezembro). Colunas obrigatórias: 'Meses',
        'P (mm/mês)', 'ETP (mm/mês)'.
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.0.
    ciclico : bool, opcional
        Se True (padrão), usa o armazenamento de equilíbrio do ciclo anual
        (o ARM de dezembro alimenta janeiro). Se False, parte do solo cheio
        antes de janeiro. Veja `balanco_hidrico_climatologico_grade`.

    Retorna
    -------
    df_bh : pandas.DataFrame
        Cópia do df de entrada, acrescido das colunas:
        'CAD', 'P-ETP', 'ARM (mm/mês)', 'NEG.ACUM (mm)', 'ALT (mm/mês)',
        'ETR (mm/mês)', 'DEF (mm/mês)', 'EXC (mm/mês)'.
    """)


balanco_hidrico_cultura = _traduzir(
    _en.crop_water_balance, "balanco_hidrico_cultura", {},
    colunas_entrada=_ENTRADA_CULTURA,
    colunas_saida=_SAIDA_CULTURA,
    doc="""
    Calcula o Balanço Hídrico de Cultura (BHc), pelo método de
    Thornthwaite & Mather, a partir de um df já estruturado com Chuva,
    ETc e CAD por período.

    Independente da escala temporal (diária, decendial, mensal etc.) —
    o usuário é responsável por pré-calcular 'ETc' (= Kc x ETo, já
    considerando a fase fenológica da cultura) e 'CAD' (= z x DTA, já
    considerando o avanço da profundidade radicular) na escala desejada;
    esta função só executa a contabilidade hídrica período a período.

    Parâmetros
    ----------
    df : pandas.DataFrame
        Colunas obrigatórias: 'Chuva' (mm/período), 'ETc' (mm/período) e
        'CAD' (mm, capacidade de água disponível no período), em ordem
        cronológica.

    Retorna
    -------
    df_bhc : pandas.DataFrame
        Cópia do df de entrada, acrescido de: 'P-ETc', 'ARM', 'ALT',
        'ETR', 'DEF', 'EXC', 'ISNA'.
    """)
