"""Funções auxiliares de trigonometria em graus."""

import numpy as np

__all__ = ["sind", "cosd", "tand", "acosd"]


def sind(x):
    """Seno de um ângulo x expresso em GRAUS."""
    return np.sin(np.radians(x))


def cosd(x):
    """Cosseno de um ângulo x expresso em GRAUS."""
    return np.cos(np.radians(x))


def tand(x):
    """Tangente de um ângulo x expresso em GRAUS."""
    return np.tan(np.radians(x))


def acosd(x):
    """Arco-cosseno que retorna o ângulo em GRAUS (em vez de radianos)."""
    return np.degrees(np.arccos(x))
