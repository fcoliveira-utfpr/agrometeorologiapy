"""Fórmulas de Agrometeorologia em Python.

Funções para radiação solar, temperatura, umidade do ar, balanço de
energia, evapotranspiração (Thornthwaite, Camargo-Maluf, Hargreaves-Samani,
Priestley-Taylor, Penman-Monteith FAO-56), grau-dias e balanço hídrico
(climatológico e de cultura).
"""

from ._trig import acosd, cosd, sind, tand
from .balanco_hidrico import balanco_hidrico_climatologico, balanco_hidrico_cultura
from .energia import boc_saldo, bol_saldo, saldo_radiacao
from .evapotranspiracao import (
    camargo_maluf_mensal,
    declive_pressao_vapor,
    eto_penman_monteith_fao56,
    etp_hargreaves_samani,
    etp_priestley_taylor,
    thornthwaite_mensal,
)
from .grau_dias import data_maturacao_fisiologica, data_semeadura
from .radiacao import (
    Qg_angstrom,
    Qg_hargreaves,
    angulo_horario,
    angulo_horario_nascer,
    angulo_zenital,
    azimute_solar,
    comprimento_sombra,
    declinacao_solar,
    fator_correcao_distancia,
    fotoperiodo,
    insolacao,
    irradiancia_extraterrestre,
    nda,
)
from .temperatura import temp_media_estacao_automatica, temp_media_extremos
from .umidade import (
    constante_psicrometrica,
    deficit_saturacao,
    ea_umidade,
    es_tetens,
    patm_altitude,
    ponto_orvalho,
    umidade_absoluta,
    umidade_relativa,
    umidade_saturacao,
)

__version__ = "0.1.0"

__all__ = [
    "sind", "cosd", "tand", "acosd",
    "nda", "declinacao_solar", "angulo_horario", "angulo_zenital",
    "azimute_solar", "comprimento_sombra", "fotoperiodo",
    "angulo_horario_nascer", "fator_correcao_distancia",
    "irradiancia_extraterrestre", "insolacao", "Qg_angstrom", "Qg_hargreaves",
    "temp_media_extremos", "temp_media_estacao_automatica",
    "es_tetens", "ea_umidade", "deficit_saturacao", "patm_altitude",
    "umidade_absoluta", "umidade_saturacao", "umidade_relativa",
    "ponto_orvalho", "constante_psicrometrica",
    "saldo_radiacao", "boc_saldo", "bol_saldo",
    "thornthwaite_mensal", "camargo_maluf_mensal", "etp_hargreaves_samani",
    "declive_pressao_vapor", "etp_priestley_taylor", "eto_penman_monteith_fao56",
    "data_maturacao_fisiologica", "data_semeadura",
    "balanco_hidrico_climatologico", "balanco_hidrico_cultura",
]
