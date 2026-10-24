"""Capítulo 12 — Produtividade potencial e atingível (Zona Agroecológica da FAO).

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.crop_yield`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

import pandas as pd  # usado nos exemplos (doctests)

from ._i18n import tradutor
from .en import crop_yield as _en

_ENTRADA_PP = {'data': 'date', 'dia': 'day', 'mes': 'month', 'ano': 'year', 'Tmed': 'Tmean',
               'Qg': 'Rs', 'Qo': 'Ra'}
_SAIDA_PP = {**{en: pt for pt, en in _ENTRADA_PP.items()},
             'cTo': 'cTn', 'GPc': 'PPBc', 'GPo': 'PPBn', 'GP': 'PPBp', 'Yp_period': 'PPf_periodo',
             'Yp_cum': 'PPf_acum', 'Yp': 'PPf', 'LAI_corr': 'CIAF'}
_SAIDA_PA = {'ETa': 'ETR', 'ETa/ETc': 'ETR/ETc', 'factor': 'fator', 'Yw_start': 'PA_inicial',
             'Yw_end': 'PA_final', 'loss': 'perda', 'CE': 'EC', 'Yp': 'PPf', 'Yw': 'PA'}


def _indice_em_portugues(tab, kw):
    """Índice da tabela da PA: 'phase' -> 'fase' e o rótulo 'cycle' -> 'ciclo'."""
    if kw.get('method') == 'single':
        tab.index = tab.index.map({'cycle': 'ciclo'})
    tab.index.name = 'fase'
    return tab

_traduzir = tradutor(__name__)

__all__ = [
    "cTn",
    "cTc",
    "PPBp",
    "CIAF",
    "CR",
    "produtividade_potencial",
    "produtividade_atingivel",
]


cTn = _traduzir(
    _en.temperature_correction_overcast, "cTn", {'rota': 'pathway', 'T_limiar': 'T_threshold'},
    doc="""
    Correção de temperatura da fotossíntese para o céu nublado (cTn),
    pelos polinômios de Barbieri & Tuon (1992), limitada a zero.

    C3: cTn = -0,0425 + 0,035 T + 0,00325 T² - 0,0000925 T³
    C4: cTn = -1,064 + 0,173 T - 0,0029 T²       (T >= 16,5 °C)
        cTn = -4,16 + 0,4325 T - 0,00725 T²      (T < 16,5 °C)

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    rota : str
        Rota fotossintética da cultura: 'C3' (soja, trigo, feijão,
        girassol, arroz...) ou 'C4' (milho, sorgo, cana-de-açúcar...).
    T_limiar : float, opcional
        Temperatura que separa os dois polinômios das culturas C4, em °C.
        Padrão 16,5. Ignorado para C3.

    Retorna
    -------
    cTn : float ou array_like
        Correção de temperatura para céu nublado (adimensional, >= 0).

    Exemplos
    --------
    >>> round(cTn(25, 'C3'), 3), round(cTn(25, 'C4'), 3)
    (1.418, 1.448)
    """)


cTc = _traduzir(
    _en.temperature_correction_clear, "cTc", {'rota': 'pathway', 'T_limiar': 'T_threshold'},
    doc="""
    Correção de temperatura da fotossíntese para o céu claro (cTc),
    pelos polinômios de Barbieri & Tuon (1992), limitada a zero.

    C3: cTc = 0,583 + 0,014 T + 0,0013 T² - 0,000037 T³
    C4: cTc = -4,16 + 0,4325 T - 0,00725 T²      (T >= 16,5 °C)
        cTc = -9,32 + 0,865 T - 0,0145 T²        (T < 16,5 °C)

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    rota : str
        Rota fotossintética da cultura: 'C3' ou 'C4'.
    T_limiar : float, opcional
        Temperatura que separa os dois polinômios das culturas C4, em °C.
        Padrão 16,5. Ignorado para C3.

    Retorna
    -------
    cTc : float ou array_like
        Correção de temperatura para céu claro (adimensional, >= 0).

    Exemplos
    --------
    >>> round(cTc(25, 'C3'), 3), round(cTc(25, 'C4'), 3)
    (1.167, 2.121)
    """)


PPBp = _traduzir(
    _en.gross_potential_productivity, "PPBp",
    {'Qo': 'Ra', 'cTn': 'cTo', 'a_n': 'a_o', 'b_n': 'b_o'},
    doc="""
    Produtividade potencial bruta padrão (PPBp) da cultura padrão
    hipotética (IAF = 5, sem restrição de água e nutrientes), soma das
    frações de céu claro e de céu nublado do dia.

    PPBc = (a_c + b_c Qo) cTc (n/N)
    PPBn = (a_n + b_n Qo) cTn (1 - n/N)
    PPBp = PPBc + PPBn

    Os coeficientes padrão são os de Doorenbos & Kassam convertidos para
    Qo em MJ/m² dia (0,36 e 0,219 com Qo em cal/cm² dia).

    Parâmetros
    ----------
    Qo : float ou array_like
        Irradiância solar extraterrestre, em MJ/m² dia.
    nN : float ou array_like
        Razão de insolação n/N (0 a 1). Ver `relacao_n_N`.
    cTc, cTn : float ou array_like
        Correções de temperatura para céu claro e nublado. Ver `cTc` e `cTn`.
    a_c, b_c : float, opcional
        Coeficientes da fração de céu claro. Padrão 107,2 e 8,604.
    a_n, b_n : float, opcional
        Coeficientes da fração de céu nublado. Padrão 31,7 e 5,234.

    Retorna
    -------
    PPBp : float ou array_like
        Produtividade potencial bruta padrão, em kg MS/ha dia.

    Exemplos
    --------
    >>> round(PPBp(42.96, 0.501, cTc=2.210, cTn=1.484), 0)
    718.0
    """)


CIAF = _traduzir(_en.lai_correction, "CIAF", {'IAF': 'LAI', 'IAF_padrao': 'LAI_standard'}, doc="""
    Correção do índice de área foliar (CIAF): ajusta a cultura padrão
    (IAF = 5) ao IAF máximo da cultura.

    CIAF = 0,0093 + 0,185 IAF - 0,0175 IAF²   (IAF < 5)
    CIAF = 0,5                                (IAF >= 5)

    Parâmetros
    ----------
    IAF : float
        Índice de área foliar máximo da cultura (m²/m²).
    IAF_padrao : float, opcional
        IAF a partir do qual CIAF = 0,5. Padrão 5,0.

    Retorna
    -------
    CIAF : float
        Correção do índice de área foliar (adimensional).

    Exemplos
    --------
    >>> round(CIAF(4.5), 3), CIAF(5)
    (0.487, 0.5)
    """)


CR = _traduzir(
    _en.respiration_correction, "CR",
    {'T_limiar': 'T_threshold', 'CR_frio': 'CR_cool', 'CR_quente': 'CR_warm'},
    doc="""
    Correção para a respiração de manutenção (CR), conforme a temperatura
    média do período.

    CR = 0,6 se T < 20 °C;  CR = 0,5 se T >= 20 °C.

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    T_limiar : float, opcional
        Temperatura de troca do CR, em °C. Padrão 20,0.
    CR_frio, CR_quente : float, opcional
        CR abaixo e a partir de `T_limiar`. Padrão 0,6 e 0,5.

    Retorna
    -------
    CR : float ou array_like
        Correção para a respiração (adimensional).

    Exemplos
    --------
    >>> CR(18.0), CR(26.5)
    (0.6, 0.5)
    """)


produtividade_potencial = _traduzir(
    _en.potential_yield, "produtividade_potencial",
    {'rota': 'pathway', 'IAF': 'LAI', 'Cc': 'HI', 'U': 'moisture', 'intervalo': 'interval',
     'semeadura': 'sowing', 'ciclo': 'cycle', 'ano': 'year', 'T_limiar_CR': 'T_threshold_CR',
     'CR_frio': 'CR_cool', 'CR_quente': 'CR_warm'},
    valores={'interval': {'dec': 'dekad'}},
    colunas_entrada=_ENTRADA_PP,
    colunas_saida=_SAIDA_PP,
    doc="""
    Produtividade potencial (PPf) pelo Modelo da Zona Agroecológica da FAO
    (Doorenbos & Kassam, 1979), acumulada período a período ao longo do
    ciclo da cultura.

    Em cada período (dia, decêndio ou mês):

    1. n/N pela inversão de Angström-Prescott (se não for informada);
    2. cTn e cTc da rota C3 ou C4 com a temperatura média do período;
    3. PPBp = PPBc + PPBn (kg MS/ha dia);
    4. CR conforme a temperatura do período;
    5. PPf_periodo = PPBp x ND x CIAF x CR x Cc / (1 - 0,01 U),
       com ND = dias do período em que a cultura estava no campo.

    PPf = soma de PPf_periodo no ciclo (coluna 'PPf_acum').

    Parâmetros
    ----------
    df : pandas.DataFrame
        Uma linha por período, em ordem cronológica. Colunas:

        - 'Tmed' (°C), obrigatória: temperatura média do período.
        - Radiação: 'nN' (razão de insolação), ou 'n' e 'N' (horas), ou
          'Qg' (MJ/m² dia, média diária do período), nessa prioridade.
        - 'Qo' (MJ/m² dia, média diária do período), opcional: se faltar,
          é calculada a partir de `lat` e das datas.
        - Datas: coluna 'data', DatetimeIndex ou colunas 'dia' e 'mes'
          (+ 'ano' opcional), com o início de cada período. Necessárias
          para calcular Qo, ND e para usar `semeadura`.
        - 'ND', opcional: número de dias do período a considerar. Se
          faltar, vem do calendário (decêndio: 10, 10 e 8 a 11 dias).
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    rota : str
        Rota fotossintética da cultura: 'C3' ou 'C4'.
    IAF : float, opcional
        Índice de área foliar máximo da cultura. Padrão 5,0.
    Cc : float, opcional
        Coeficiente (índice) de colheita. Padrão 0,35.
    U : float, opcional
        Umidade do produto colhido, em %. Padrão 13,0.
    intervalo : str, opcional
        Escala de tempo das linhas do df: 'd' (diário), 'dec' (decendial
        do calendário: 1-10, 11-20, 21-fim do mês) ou 'M' (mensal).
        Padrão 'dec'.
    semeadura : str ou datetime, opcional
        Data de semeadura ('AAAA-MM-DD'). Com `ciclo`, só os dias de
        semeadura a semeadura + ciclo - 1 entram em ND (os períodos das
        pontas ficam parciais e os de fora do ciclo são descartados).
    ciclo : int, opcional
        Duração do ciclo, em dias, contando o dia da semeadura.
    ano : int, opcional
        Ano do primeiro período quando o df tem 'dia' e 'mes' sem 'ano'.
        Padrão 2023.
    a, b : float, opcional
        Coeficientes de Angström-Prescott para n/N. Padrão a = 0,29 cos(lat)
        e b = 0,52.
    T_limiar_CR, CR_frio, CR_quente : float, opcional
        Parâmetros da correção de respiração (ver `CR`). Padrão 20,0, 0,6 e 0,5.

    Retorna
    -------
    resultado : pandas.DataFrame
        Cópia do df (só os períodos do ciclo), acrescido de: 'ND', 'Qo',
        'nN', 'cTn', 'cTc', 'PPBc', 'PPBn', 'PPBp' (kg MS/ha dia), 'CR',
        'PPf_periodo' e 'PPf_acum' (kg/ha, na umidade U). A última linha
        de 'PPf_acum' é a PPf. `resultado.attrs` guarda 'PPf' e 'CIAF'.

    Exemplos
    --------
    1º decêndio de janeiro em Santa Helena-PR, milho (C4):

    >>> df = pd.DataFrame({'dia': [1], 'mes': [1], 'Tmed': [26.5], 'Qg': [22.5]})
    >>> pp = produtividade_potencial(df, lat=-24.86, rota='C4', IAF=5.0, Cc=0.35, U=13)
    >>> round(float(pp['PPBp'].iloc[0]), 0), round(pp.attrs['PPf'], 0)
    (718.0, 722.0)
    """)


produtividade_atingivel = _traduzir(
    _en.attainable_yield, "produtividade_atingivel",
    {'PPf': 'Yp', 'ETR': 'ETa', 'fase': 'phase', 'metodo': 'method'},
    valores={'method': {'produtorio': 'product', 'etapa_unica': 'single'}},
    colunas_saida=_SAIDA_PA,
    pos=_indice_em_portugues,
    doc="""
    Produtividade atingível (PA): a PPf penalizada pelo déficit hídrico,
    pela relação de Doorenbos & Kassam (1979):

    1 - PA/PPf = ky (1 - ETR/ETc)

    - 'produtorio' (Battisti et al., 2013): aplicada fase a fase, em
      sequência; a PA ao fim de uma fase é a "PPf" da seguinte.
      PA = PPf x prod[1 - ky_i (1 - ETR_i/ETc_i)]
    - 'etapa_unica': ky do ciclo e somas de ETR e ETc no ciclo.
      PA = PPf x [1 - ky_ciclo (1 - soma ETR / soma ETc)]

    Cada fator é limitado a zero (PA nunca negativa). Fases com ETc = 0
    não penalizam.

    Parâmetros
    ----------
    PPf : float
        Produtividade potencial, em kg/ha (ver `produtividade_potencial`).
    ETR, ETc : array_like
        Evapotranspiração real e da cultura (mm), por fase ou por período
        (ex.: as colunas 'ETR' e 'ETc' de `balanco_hidrico_cultura`).
    ky : float, sequência ou dict
        Fator de resposta da produtividade ao déficit hídrico. No
        produtório, um ky por fase: sequência na ordem das fases ou dict
        {fase: ky}. Na etapa única, o ky do ciclo (escalar).
    fase : array_like, opcional
        Fase de cada período (mesmo tamanho de ETR e ETc). ETR e ETc são
        somados por fase, na ordem em que as fases aparecem. Se None, cada
        elemento de ETR e ETc já é uma fase.
    metodo : str, opcional
        'produtorio' (padrão) ou 'etapa_unica'.

    Retorna
    -------
    resultado : pandas.DataFrame
        Uma linha por fase (ou uma linha 'ciclo' na etapa única):
        'ETc', 'ETR', 'ETR/ETc', 'ky', 'fator', 'PA_inicial', 'PA_final',
        'perda' (kg/ha) e 'EC' (PA_final / PPf, acumulada). A última linha
        de 'PA_final' é a PA. `resultado.attrs` guarda 'PPf', 'PA' e 'EC'.

    Exemplos
    --------
    Milho com déficit na floração:

    >>> etc, etr = [60, 110, 120, 40], [57, 77, 102, 36]
    >>> pa = produtividade_atingivel(1000, etr, etc, ky=[0.4, 1.5, 0.5, 0.2])
    >>> round(pa.attrs['EC'], 3)
    0.489
    >>> pa1 = produtividade_atingivel(1000, etr, etc, ky=1.25, metodo='etapa_unica')
    >>> round(pa1.attrs['EC'], 3)
    0.78
    """)
