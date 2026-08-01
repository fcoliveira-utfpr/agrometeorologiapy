"""Capítulo 8 — Balanço de energia."""

__all__ = ["saldo_radiacao", "boc_saldo", "bol_saldo"]


def saldo_radiacao(BOC, BOL):
    """
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
    """
    return BOC + BOL


def boc_saldo(Qg, r=0.25):
    """
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
    """
    return Qg * (1 - r)


def bol_saldo(Tmax, Tmin, ea, Qg, Qg_cs):
    """
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
    """
    Tmax_K = Tmax + 273.15
    Tmin_K = Tmin + 273.15
    termo_temp = 4.903e-9 * ((Tmax_K ** 4 + Tmin_K ** 4) / 2)
    termo_umidade = 0.34 - 0.14 * ea ** 0.5
    termo_nebulosidade = 1.35 * (Qg / Qg_cs) - 0.35
    BOL = -(termo_temp * termo_umidade * termo_nebulosidade)
    return BOL
