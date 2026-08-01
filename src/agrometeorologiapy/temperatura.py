"""Capítulo 6 — Temperatura."""

__all__ = ["temp_media_extremos", "temp_media_estacao_automatica"]


def temp_media_extremos(Tmax, Tmin):
    """
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
    """
    return (Tmax + Tmin) / 2


def temp_media_estacao_automatica(temperaturas):
    """
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
    """
    return sum(temperaturas) / len(temperaturas)
