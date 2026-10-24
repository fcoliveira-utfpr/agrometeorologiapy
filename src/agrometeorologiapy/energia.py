"""Capítulo 8 — Balanço de energia.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.energy`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import energy as _en

_traduzir = tradutor(__name__)

__all__ = [
    "saldo_radiacao",
    "boc_saldo",
    "bol_saldo",
]


saldo_radiacao = _traduzir(
    _en.net_radiation, "saldo_radiacao", {'BOC': 'Rns', 'BOL': 'Rnl'},
    doc="""
    Saldo de radiação (Rn), soma do balanço de ondas curtas
    (BOC) com o balanço de ondas longas (BOL).

    Representa a energia radiante líquida disponível numa superfície para
    os processos físicos e biológicos que nela ocorrem.

    Parâmetros
    ----------
    BOC : float
        Balanço de ondas curtas, mesma unidade de Rn.
    BOL : float
        Balanço de ondas longas, mesma unidade de Rn.

    Retorna
    -------
    Rn : float
        Saldo de radiação, em MJ/m² dia (ou W/m²).
    """)


boc_saldo = _traduzir(
    _en.net_shortwave_radiation, "boc_saldo", {'Qg': 'Rs', 'r': 'albedo'},
    doc="""
    Balanço de ondas curtas (BOC): irradiância solar global
    incidente menos a parcela refletida pela superfície (albedo).

    Parâmetros
    ----------
    Qg : float
        Irradiância solar global incidente na superfície (MJ/m² dia ou W/m²).
    r : float
        Coeficiente de reflexão da superfície (albedo), adimensional (0-1).
        r = 0.25 média para gramado
    Retorna
    -------
    BOC : float
        Balanço de ondas curtas, mesma unidade de Qg.
    """)


bol_saldo = _traduzir(
    _en.net_longwave_radiation, "bol_saldo", {'Qg': 'Rs', 'Qg_cs': 'Rso'},
    doc="""
    Balanço de radiação de onda longa (BOL), segundo a equação de
    Stefan-Boltzmann corrigida pela FAO-56.

    Estima o saldo líquido de radiação de onda longa emitida pela
    superfície, corrigido pela nebulosidade (razão Qg/Qg_cs) e pela
    umidade do ar (via pressão de vapor atual, ea). O termo de
    Stefan-Boltzmann exige temperatura absoluta; por consistência com as
    demais funções do pacote (que recebem Tmax/Tmin em °C), a conversão
    para Kelvin é feita internamente.

    Parâmetros
    ----------
    Tmax : float
        Temperatura máxima diária, em °C.
    Tmin : float
        Temperatura mínima diária, em °C.
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.
    Qg : float
        Radiação solar global medida/estimada no dia, em MJ/m².
    Qg_cs : float
        Radiação solar de céu claro (clear-sky), em MJ/m².

    Retorna
    -------
    BOL : float
        Balanço de radiação de onda longa, em MJ/m² (valor negativo,
        representando perda líquida de energia por emissão terrestre).
    """)
