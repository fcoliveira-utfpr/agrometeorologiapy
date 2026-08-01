"""Capítulo 9 — Evapo(transpi)ração."""

import numpy as np

from ._trig import acosd, cosd, sind
from .umidade import es_tetens

__all__ = [
    "thornthwaite_mensal",
    "camargo_maluf_mensal",
    "etp_hargreaves_samani",
    "declive_pressao_vapor",
    "etp_priestley_taylor",
    "eto_penman_monteith_fao56",
]


def thornthwaite_mensal(df, col_T='T_media_C', lat=None):
    """
    Evapotranspiração potencial mensal (ETP) pelo método de Thornthwaite (1948).

    Estima a ETP a partir da temperatura média mensal do ar, corrigida pelo
    fotoperíodo (número de horas de brilho solar) em função da latitude e
    do número de dias de cada mês, conforme a formulação clássica do método.

    Parâmetros
    ----------
    df : pandas.DataFrame
        DataFrame com 12 linhas (Janeiro a Dezembro, nessa ordem) contendo
        ao menos a coluna de temperatura média mensal.
    col_T : str, opcional
        Nome da coluna de temperatura média mensal (°C) em `df`. Padrão 'T_media_C'.
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).

    Retorna
    -------
    df : pandas.DataFrame
        Cópia do DataFrame de entrada com a coluna adicional 'ETP_mm_mes'
        (evapotranspiração potencial corrigida, em mm/mês).
    """
    df = df.copy()
    T_mensal = df[col_T].to_numpy(dtype=float)
    dias_mes = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    dia_juliano_medio = np.array([17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344])

    i_mensal = np.where(T_mensal > 0, (T_mensal / 5) ** 1.514, 0)
    I = i_mensal.sum()
    a = 6.75e-7 * I ** 3 - 7.71e-5 * I ** 2 + 1.792e-2 * I + 0.49239
    ETP_nc = np.where(T_mensal > 0, 16 * (10 * T_mensal / I) ** a, 0)

    declinacao = 23.45 * sind(360 * (284 + dia_juliano_medio) / 365)
    tand_lat = sind(lat) / cosd(lat)
    tand_dec = sind(declinacao) / cosd(declinacao)
    ws = acosd(-tand_lat * tand_dec)
    N = (2 / 15) * ws
    K = (N / 12) * (dias_mes / 30)

    df['ETP_mm_mes'] = np.round(ETP_nc * K, 2)
    return df


def camargo_maluf_mensal(df, col_T='T_media_C', lat=None):
    """
    Evapotranspiração potencial mensal (ETP) pelo método de Camargo (1971),
    com o coeficiente F modificado por Maluf (Camargo et al., 1999).

    Estima a ETP a partir da temperatura média mensal do ar e da radiação
    solar extraterrestre (Qo), calculada internamente a partir da declinação
    solar, do ângulo horário do nascer do Sol e da correção pela distância
    Terra-Sol, com coeficiente de ajuste F dependente da temperatura média
    anual do local. Função autocontida — não depende de funções auxiliares
    externas (sind, cosd, tand, irradiancia_extraterrestre, etc.).

    Parâmetros
    ----------
    df : pandas.DataFrame
        DataFrame com 12 linhas (Janeiro a Dezembro, nessa ordem) contendo
        ao menos a coluna de temperatura média mensal.
    col_T : str, opcional
        Nome da coluna de temperatura média mensal (°C) em `df`. Padrão 'T_media_C'.
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).

    Retorna
    -------
    df : pandas.DataFrame
        Cópia do DataFrame de entrada com as colunas adicionais 'Qo_MJ_m2dia'
        (radiação extraterrestre, MJ/m² dia), 'Qo_mm_dia' (equivalente de
        evaporação, mm/dia) e 'ETP_mm_mes' (ETP mensal, em mm/mês).
    """
    df = df.copy()
    T_mensal = df[col_T].to_numpy(dtype=float)
    dias_mes = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    J = np.array([17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344])  # dia juliano médio

    # Declinação solar (graus)
    declinacao = 23.45 * np.sin(np.radians(360 * (284 + J) / 365))

    # Ângulo horário do nascer do Sol (graus) - eq. 5.20
    tan_lat = np.sin(np.radians(lat)) / np.cos(np.radians(lat))
    tan_dec = np.sin(np.radians(declinacao)) / np.cos(np.radians(declinacao))
    Hn = np.degrees(np.arccos(np.clip(-tan_lat * tan_dec, -1, 1)))

    # Correção da distância Terra-Sol (d/D)^2 - eq. 5.30
    dD2 = 1 + 0.033 * np.cos(np.radians(360 * J / 365))

    # Irradiância solar extraterrestre diária (Qo), em MJ/m² dia
    hn_rad = np.radians(Hn)
    Qo_MJ = 37.6 * dD2 * (
        hn_rad * np.sin(np.radians(lat)) * np.sin(np.radians(declinacao))
        + np.cos(np.radians(lat)) * np.cos(np.radians(declinacao)) * np.sin(np.radians(Hn))
    )
    Qo_mm = 0.408 * Qo_MJ  # equivalente de evaporação, mm/dia

    # Coeficiente F em função da temperatura média anual (Camargo modif. Maluf)
    Tmed_anual = T_mensal.mean()
    if Tmed_anual < 23:
        F = 0.01
    elif Tmed_anual < 24:
        F = 0.0105
    else:
        F = 0.011

    ETP = F * Qo_mm * T_mensal * dias_mes

    df['Qo_MJ_m2dia'] = np.round(Qo_MJ, 2)
    df['Qo_mm_dia'] = np.round(Qo_mm, 2)
    df['ETP_mm_mes'] = np.round(ETP, 2)
    return df


def etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed):
    """
    Evapotranspiração potencial pelo método de
    Hargreaves & Samani (1985), que usa a amplitude térmica diária como
    substituto indireto da nebulosidade/radiação solar efetiva.

    Parâmetros
    ----------
    Qo : float
        Irradiância solar extraterrestre, em MJ/m² dia.
    Tmax : float
        Temperatura máxima do ar, em °C.
    Tmin : float
        Temperatura mínima do ar, em °C.
    Tmed : float
        Temperatura média do ar, em °C.

    Retorna
    -------
    ETP : float
        Evapotranspiração potencial, em mm/dia.
    """
    Qo_mm = 0.408 * Qo  # equivalente de evaporação, mm/dia
    return 0.0023 * Qo_mm * (Tmax - Tmin) ** 0.5 * (Tmed + 17.8)


def declive_pressao_vapor(T_ar):
    """
    Declive da curva de pressão de saturação de vapor (Delta), na
    temperatura do ar considerada.

    Obtida derivando-se a equação de Tetens (eq. 7.2) em relação à
    temperatura. Necessária nos métodos combinados (energia + aerodinâmico),
    como Priestley-Taylor e Penman-Monteith.

    Parâmetros
    ----------
    T_ar : float
        Temperatura do ar, em °C (em geral, a temperatura média diária).

    Retorna
    -------
    Delta : float
        Declive da curva de pressão de saturação de vapor, em kPa/°C.
    """
    es = es_tetens(T_ar)
    return 4098 * es / (T_ar + 237.3) ** 2


def etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26):
    """
    Evapotranspiração potencial pelo método de Priestley & Taylor (1972).

    Simplificação do termo aerodinâmico da equação de Penman, válida para
    superfícies bem supridas de água. O coeficiente empírico alfa (1,26)
    compensa a advecção não contabilizada quando se ignora o termo
    aerodinâmico.

    Parâmetros
    ----------
    Rn : float
        Saldo de radiação, em MJ/m² dia.
    G : float
        Fluxo de calor no solo, em MJ/m² dia (geralmente desprezado,
        G = 0, em escala diária).
    Delta : float
        Declive da curva de pressão de saturação de vapor, em kPa/°C.
    gamma : float
        Constante psicrométrica, em kPa/°C.
    alfa : float, opcional
        Coeficiente de Priestley-Taylor (padrão 1,26).

    Retorna
    -------
    ETP : float
        Evapotranspiração potencial, em mm/dia.
    """
    lambda_v = 2.45  # calor latente de vaporização, MJ/kg
    return alfa * (Delta / (Delta + gamma)) * (Rn - G) / lambda_v


def eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma):
    """
    Evapotranspiração de referência (ETo) pela equação de
    Penman-Monteith, padronizada pelo boletim FAO-56 (Allen et al., 1998).

    Combina os termos de energia (radiativo) e aerodinâmico (advectivo),
    referenciados a uma cultura hipotética (grama, altura 0,12 m, albedo
    0,23, resistência de superfície fixa), sendo o método-padrão
    internacional para estimativa da evapotranspiração de referência.

    Parâmetros
    ----------
    Rn : float
        Saldo de radiação, em MJ/m² dia.
    G : float
        Fluxo de calor no solo, em MJ/m² dia (G = 0 para escala diária).
    Tmed : float
        Temperatura média diária do ar, em °C.
    u2 : float
        Velocidade do vento a 2 m de altura, em m/s.
    es : float
        Pressão de saturação de vapor, em kPa.
    ea : float
        Pressão parcial (atual) de vapor d'água, em kPa.
    Delta : float
        Declive da curva de pressão de saturação de vapor, em kPa/°C.
    gamma : float
        Constante psicrométrica, em kPa/°C.

    Retorna
    -------
    ETo : float
        Evapotranspiração de referência, em mm/dia.
    """
    numerador = 0.408 * Delta * (Rn - G) + gamma * (900 / (Tmed + 273)) * u2 * (es - ea)
    denominador = Delta + gamma * (1 + 0.34 * u2)
    return numerador / denominador
