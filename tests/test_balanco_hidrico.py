import numpy as np
import pandas as pd
import pytest

import agrometeorologiapy as amp


def test_aplicacao_10_balanco_hidrico_climatologico():
    dados_bhc = {
        'Meses': ['jan', 'fev', 'mar', 'abr', 'maio', 'jun', 'jul', 'ago', 'set', 'out', 'nov', 'dez'],
        'P (mm/mês)': [203.400000, 189.533333, 196.266667, 116.933333, 14.033333, 3.333333,
                       1.300000, 4.366667, 18.033333, 99.466667, 223.733333, 234.966667],
        'ETP (mm/mês)': [128.730000, 106.286667, 117.300000, 109.263333, 106.800000, 102.033333,
                          114.986667, 134.523333, 141.426667, 149.146667, 118.616667, 119.960000],
    }
    df_ini = pd.DataFrame(dados_bhc)
    df_bh = amp.balanco_hidrico_climatologico(df_ini, CAD=100.0)

    esperado_arm = [100.0, 100.0, 100.0, 100.0, 39.548, 14.739, 4.729, 1.287, 0.375, 0.228, 100.0, 100.0]
    esperado_etr = [128.73, 106.287, 117.3, 109.263, 74.486, 28.142, 11.31, 7.809, 18.945, 99.613, 118.617, 119.96]
    esperado_def = [0.0, 0.0, 0.0, 0.0, 32.314, 73.892, 103.676, 126.715, 122.481, 49.533, 0.0, 0.0]
    esperado_exc = [74.67, 83.247, 78.967, 7.67, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.345, 115.007]

    assert df_bh['ARM (mm/mês)'].round(3).tolist() == pytest.approx(esperado_arm, abs=1e-2)
    assert df_bh['ETR (mm/mês)'].round(3).tolist() == pytest.approx(esperado_etr, abs=1e-2)
    assert df_bh['DEF (mm/mês)'].round(3).tolist() == pytest.approx(esperado_def, abs=1e-2)
    assert df_bh['EXC (mm/mês)'].round(3).tolist() == pytest.approx(esperado_exc, abs=1e-2)


# Semiárido nordestino (tipo Petrolina): janeiro é seco e a estação chuvosa
# é fev-abr, então partir do solo cheio em janeiro não serve.
P_SEMIARIDO = np.array([40.0, 80, 110, 70, 20, 10, 5, 2, 3, 8, 20, 35])
ETP_SEMIARIDO = np.array([170.0, 145, 150, 135, 125, 110, 120, 150, 170, 190, 185, 180])
P_UMIDO = np.array([203.4, 189.53, 196.27, 116.93, 14.03, 3.33, 1.3, 4.37, 18.03, 99.47, 223.73, 234.97])
ETP_UMIDO = np.array([128.73, 106.29, 117.3, 109.26, 106.8, 102.03, 114.99, 134.52, 141.43, 149.15, 118.62, 119.96])


def test_bhc_grade_igual_ao_dataframe():
    df_ini = pd.DataFrame({'Meses': range(1, 13), 'P (mm/mês)': P_UMIDO, 'ETP (mm/mês)': ETP_UMIDO})
    df_bh = amp.balanco_hidrico_climatologico(df_ini, CAD=100.0)
    bh = amp.balanco_hidrico_climatologico_grade(P_UMIDO, ETP_UMIDO, CAD=100.0)
    for chave, coluna in [('ARM', 'ARM (mm/mês)'), ('ETR', 'ETR (mm/mês)'),
                          ('DEF', 'DEF (mm/mês)'), ('EXC', 'EXC (mm/mês)')]:
        assert bh[chave] == pytest.approx(df_bh[coluna].to_numpy())


def test_bhc_grade_varios_pixels_e_cad_por_pixel():
    # 2 x 3 pixels: cada coluna é um clima; cada linha, uma CAD diferente.
    P = np.stack([P_UMIDO, P_SEMIARIDO, P_UMIDO], axis=1)[:, np.newaxis, :].repeat(2, axis=1)
    ETP = np.stack([ETP_UMIDO, ETP_SEMIARIDO, ETP_UMIDO], axis=1)[:, np.newaxis, :].repeat(2, axis=1)
    CAD = np.array([[100.0, 100.0, 100.0], [50.0, 150.0, 200.0]])
    bh = amp.balanco_hidrico_climatologico_grade(P, ETP, CAD)
    assert bh['ARM'].shape == (12, 2, 3)
    for i in range(2):
        for j in range(3):
            um = amp.balanco_hidrico_climatologico_grade(P[:, i, j], ETP[:, i, j], CAD[i, j])
            for chave in ('ARM', 'ETR', 'DEF', 'EXC'):
                # Com vários pixels o laço segue até todos convergirem: diferenças < tol.
                assert bh[chave][:, i, j] == pytest.approx(um[chave], abs=0.02)


def test_bhc_ciclico_fecha_o_ano_e_conserva_agua():
    CAD = 100.0
    bh = amp.balanco_hidrico_climatologico_grade(P_SEMIARIDO, ETP_SEMIARIDO, CAD)
    # O ARM de dezembro, aplicado a janeiro, reproduz o ARM de janeiro.
    p_etp_jan = P_SEMIARIDO[0] - ETP_SEMIARIDO[0]
    assert bh['ARM'][0] == pytest.approx(bh['ARM'][-1] * np.exp(p_etp_jan / CAD), abs=0.05)
    # Em equilíbrio a soma das variações de armazenamento é ~0: P = ETR + EXC e ETP = ETR + DEF.
    assert bh['ALT'].sum() == pytest.approx(0.0, abs=0.05)
    assert P_SEMIARIDO.sum() == pytest.approx(bh['ETR'].sum() + bh['EXC'].sum(), abs=0.05)
    assert ETP_SEMIARIDO.sum() == pytest.approx(bh['ETR'].sum() + bh['DEF'].sum())
    # Nesse clima o solo nunca enche: sem excedente e armazenamento bem abaixo da CAD.
    assert bh['EXC'].sum() == 0.0
    assert bh['ARM'].max() < 0.5 * CAD


def test_bhc_nao_ciclico_parte_do_solo_cheio():
    CAD = 100.0
    bh = amp.balanco_hidrico_climatologico_grade(P_SEMIARIDO, ETP_SEMIARIDO, CAD, ciclico=False)
    arm_jan = CAD * np.exp((P_SEMIARIDO[0] - ETP_SEMIARIDO[0]) / CAD)
    assert bh['ARM'][0] == pytest.approx(arm_jan)
    # ALT de janeiro conta a saída de água a partir do solo cheio (não é zerado).
    assert bh['ALT'][0] == pytest.approx(arm_jan - CAD)
    assert bh['ETR'][0] == pytest.approx(P_SEMIARIDO[0] + CAD - arm_jan)


def test_bhc_grade_valida_entradas():
    with pytest.raises(ValueError):
        amp.balanco_hidrico_climatologico_grade(np.ones(11), np.ones(11))
    with pytest.raises(ValueError):
        amp.balanco_hidrico_climatologico_grade(np.ones(12), np.ones(12), CAD=0.0)


def test_aplicacao_11_balanco_hidrico_cultura():
    dados_bhc_cultura = {
        'Chuva': [45.2, 38.7, 22.4, 15.1, 8.3, 5.6, 4.2, 9.8, 18.5, 35.9, 48.3, 52.1],
        'ETc':   [17.6, 23.93, 33.9, 44.46, 54.16, 53.47, 52.78, 49.28, 41.04, 34.0, 28.42, 24.54],
        'CAD':   [22.5, 37.5, 52.5, 67.5, 82.5, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0],
    }
    df_ini = pd.DataFrame(dados_bhc_cultura)
    df_bhc = amp.balanco_hidrico_cultura(df_ini)

    esperado_arm = [22.5, 37.27, 29.938, 19.379, 11.115, 6.53, 3.806, 2.455, 1.911, 3.811, 23.691, 51.251]
    esperado_etr = [17.6, 23.93, 29.732, 25.66, 16.564, 10.185, 6.924, 11.152, 19.044, 34.0, 28.42, 24.54]
    esperado_def = [0.0, 0.0, 4.168, 18.8, 37.596, 43.285, 45.856, 38.128, 21.996, 0.0, 0.0, 0.0]
    esperado_exc = [27.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    esperado_isna = [1.0, 1.0, 0.877, 0.577, 0.306, 0.19, 0.131, 0.226, 0.464, 1.0, 1.0, 1.0]

    assert df_bhc['ARM'].round(3).tolist() == pytest.approx(esperado_arm, abs=1e-2)
    assert df_bhc['ETR'].round(3).tolist() == pytest.approx(esperado_etr, abs=1e-2)
    assert df_bhc['DEF'].round(3).tolist() == pytest.approx(esperado_def, abs=1e-2)
    assert df_bhc['EXC'].round(3).tolist() == pytest.approx(esperado_exc, abs=1e-2)
    assert df_bhc['ISNA'].round(3).tolist() == pytest.approx(esperado_isna, abs=1e-2)
