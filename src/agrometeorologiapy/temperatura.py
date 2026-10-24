"""Capítulo 6 — Temperatura.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.temperature`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import temperature as _en

_traduzir = tradutor(__name__)

__all__ = [
    "temp_media_extremos",
    "temp_media_estacao_automatica",
]


temp_media_extremos = _traduzir(_en.mean_temperature, "temp_media_extremos", {}, doc="""
    Temperatura média diária estimada pela média dos valores
    extremos (máxima e mínima).

    Método mais simples e mais usado, porém tende a superestimar a
    temperatura média "real", pois o ritmo diário da temperatura não é
    simétrico em torno da média (fica mais tempo próximo da mínima
    noturna do que da máxima diurna).

    Parâmetros
    ----------
    Tmax, Tmin : float
        Temperaturas máxima e mínima do dia, em °C.

    Retorna
    -------
    Tmed : float
        Temperatura média diária estimada, em °C.
    """)


temp_media_estacao_automatica = _traduzir(
    _en.mean_temperature_observations, "temp_media_estacao_automatica",
    {'temperaturas': 'temperatures'},
    doc="""
    Temperatura média a partir de observações de estações
    automáticas (média aritmética simples de N observações no período).

    Quanto maior o número de observações, mais próxima a estimativa fica
    do valor "real" (integral contínua da temperatura ao longo do dia).

    Parâmetros
    ----------
    temperaturas : list[float]
        Lista com as temperaturas (°C) de cada observação (No valores).

    Retorna
    -------
    Tmed : float
        Temperatura média do período, em °C.
    """)
