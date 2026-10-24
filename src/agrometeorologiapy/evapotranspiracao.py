"""Capítulo 9 — Evapotranspiração.

API em português. A implementação está no núcleo em inglês,
`agrometeorologiapy.en.evapotranspiration`; aqui ficam só os nomes, a tradução
dos parâmetros e a documentação em português.
"""

from ._i18n import tradutor
from .en import evapotranspiration as _en

_COLUNAS_SAIDA = {'T_mean_C': 'T_media_C', 'PET_mm_month': 'ETP_mm_mes',
                  'Ra_MJ_m2day': 'Qo_MJ_m2dia', 'Ra_mm_day': 'Qo_mm_dia'}

_traduzir = tradutor(__name__)

__all__ = [
    "thornthwaite_mensal",
    "camargo_maluf_mensal",
    "etp_hargreaves_samani",
    "declive_pressao_vapor",
    "etp_priestley_taylor",
    "vento_2m",
    "eto_penman_monteith_fao56",
]


thornthwaite_mensal = _traduzir(
    _en.pet_thornthwaite, "thornthwaite_mensal", {},
    valores={'col_T': {'T_media_C': 'T_mean_C'}},
    colunas_entrada={'T_media_C': 'T_mean_C'},
    colunas_saida=_COLUNAS_SAIDA,
    doc="""
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
    """)


camargo_maluf_mensal = _traduzir(
    _en.pet_camargo, "camargo_maluf_mensal", {},
    valores={'col_T': {'T_media_C': 'T_mean_C'}},
    colunas_entrada={'T_media_C': 'T_mean_C'},
    colunas_saida=_COLUNAS_SAIDA,
    doc="""
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
    """)


etp_hargreaves_samani = _traduzir(
    _en.pet_hargreaves_samani, "etp_hargreaves_samani", {'Qo': 'Ra', 'Tmed': 'Tmean'},
    doc="""
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
    """)


declive_pressao_vapor = _traduzir(
    _en.slope_vapour_pressure_curve, "declive_pressao_vapor", {'T_ar': 'T'},
    doc="""
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
    """)


etp_priestley_taylor = _traduzir(
    _en.pet_priestley_taylor, "etp_priestley_taylor", {'alfa': 'alpha'},
    doc="""
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
    """)


vento_2m = _traduzir(_en.wind_speed_2m, "vento_2m", {}, doc="""
    Converte a velocidade do vento medida a uma altura qualquer para a
    altura-padrão de 2 m, pelo perfil logarítmico do vento (FAO-56, eq. 47).

    Necessária para a equação de Penman-Monteith FAO-56, que exige o vento
    a 2 m acima de uma superfície gramada.

    Parâmetros
    ----------
    uz : float ou array
        Velocidade do vento medida à altura z, em m/s.
    z : float ou array
        Altura de medição acima da superfície do solo, em m.

    Retorna
    -------
    u2 : float ou array
        Velocidade do vento a 2 m de altura, em m/s.
    """)


eto_penman_monteith_fao56 = _traduzir(
    _en.eto_penman_monteith, "eto_penman_monteith_fao56", {'Tmed': 'Tmean'},
    doc="""
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
    """)
