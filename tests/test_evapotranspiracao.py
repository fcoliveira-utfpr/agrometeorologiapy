import pandas as pd
import pytest

import agrometeorologiapy as amp

MESES = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
T_MENSAL = [24.5, 24.8, 23.9, 21.3, 18.2, 16.5, 16.0, 17.4, 18.9, 21.0, 22.6, 23.8]


def test_aplicacao_3_thornthwaite_mensal():
    df_in = pd.DataFrame({'Mes': MESES, 'T_media_C': T_MENSAL})
    df_out = amp.thornthwaite_mensal(df_in, lat=-24.85)
    esperado = [129.94, 115.77, 111.65, 77.97, 53.15, 39.91, 38.95, 49.55, 61.64, 86.1, 103.42, 123.01]
    assert df_out['ETP_mm_mes'].tolist() == pytest.approx(esperado)


def test_aplicacao_4_camargo_maluf_mensal():
    df_in = pd.DataFrame({'Mes': MESES, 'T_media_C': T_MENSAL})
    df_out = amp.camargo_maluf_mensal(df_in, lat=-24.85)
    esperado_etp = [131.81, 113.02, 107.2, 76.83, 55.59, 43.58, 45.76, 59.34, 75.98, 101.43, 115.23, 129.48]
    assert df_out['ETP_mm_mes'].tolist() == pytest.approx(esperado_etp)


def test_aplicacao_5_etp_hargreaves_samani():
    Tmax, Tmin, Qo = 21.2, 7.4, 22.8
    Tmed = (Tmax + Tmin) / 2
    ETP = amp.etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)
    assert Tmed == pytest.approx(14.3)
    assert ETP == pytest.approx(2.551334617209329)


def test_aplicacao_6_priestley_taylor():
    Delta = 0.10551678957120231
    gamma = 0.06553693475060951
    Rn = 6.192321952102397
    G = 0
    ETP_PT = amp.etp_priestley_taylor(Rn, G, Delta, gamma)
    assert ETP_PT == pytest.approx(1.9644773392383408)


def test_aplicacao_7_penman_monteith_fao56():
    Rn, G, Tmed, u2 = 6.192321952102397, 0, 14.3, 2.0
    es, ea = 1.773644115842008, 1.3407926169150235
    Delta, gamma = 0.10551678957120231, 0.06553693475060951
    ETo_PM = amp.eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)
    assert ETo_PM == pytest.approx(2.06065115224452)


def test_declive_pressao_vapor():
    assert amp.declive_pressao_vapor(14.3) == pytest.approx(0.10551678957120231)
