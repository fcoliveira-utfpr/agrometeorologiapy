import numpy as np
import pandas as pd
import pytest

import agrometeorologiapy as amp

# Tabela de conferência da Aula 03 (Relações Físicas do Ambiente Agrícola).
T = [15.0, 17.5, 20.0, 25.0, 30.0, 35.0]
C3_CTN = [0.902, 1.070, 1.218, 1.418, 1.435, 1.198]
C3_CTC = [0.961, 1.028, 1.087, 1.167, 1.174, 1.079]
C4_CTN = [0.696, 1.075, 1.236, 1.448, 1.516, 1.438]
C4_CTC = [0.392, 1.188, 1.590, 2.121, 2.290, 2.096]


def test_correcoes_temperatura_tabela():
    assert amp.cTn(np.array(T), 'C3') == pytest.approx(C3_CTN, abs=1e-3)
    assert amp.cTc(np.array(T), 'C3') == pytest.approx(C3_CTC, abs=1e-3)
    assert amp.cTn(np.array(T), 'c4') == pytest.approx(C4_CTN, abs=1e-3)
    assert amp.cTc(np.array(T), 'C4') == pytest.approx(C4_CTC, abs=1e-3)


def test_correcoes_c4_limitadas_a_zero():
    assert amp.cTn(10, 'C4') == 0.0
    assert amp.cTc(12.5, 'C4') == 0.0


def test_rota_invalida():
    with pytest.raises(ValueError):
        amp.cTn(20, 'CAM')


def test_tipo_de_saida_segue_entrada():
    s = pd.Series([20.0, 25.0], index=['a', 'b'])
    assert isinstance(amp.cTc(25.0, 'C3'), float)
    assert amp.cTc(s, 'C3').index.tolist() == ['a', 'b']


def test_ciaf_e_cr():
    assert [amp.CIAF(i) for i in (3.0, 3.5, 4.0, 4.5, 5.0, 6.0)] == pytest.approx(
        [0.407, 0.442, 0.469, 0.487, 0.5, 0.5], abs=1e-3)
    assert amp.CR(np.array([19.9, 20.0])).tolist() == [0.6, 0.5]
    assert amp.CR(15, CR_frio=0.65) == 0.65


def test_exemplo_decendio_santa_helena():
    # 1º decêndio de janeiro, φ = -24,86°, T = 26,5 °C, Qg = 22,5 MJ/m² dia.
    df = pd.DataFrame({'dia': [1], 'mes': [1], 'Tmed': [26.5], 'Qg': [22.5]})
    milho = amp.produtividade_potencial(df, lat=-24.86, rota='C4', IAF=5.0, Cc=0.35, U=13)
    soja = amp.produtividade_potencial(df, lat=-24.86, rota='C3', IAF=4.5, Cc=0.35, U=13)

    assert milho['Qo'].iloc[0] == pytest.approx(42.96, abs=0.01)
    assert milho['nN'].iloc[0] == pytest.approx(0.501, abs=1e-3)
    assert milho['ND'].iloc[0] == 10
    assert milho['PPBc'].iloc[0] == pytest.approx(528.2, abs=0.2)
    assert milho['PPBn'].iloc[0] == pytest.approx(189.9, abs=0.2)
    assert milho['PPBp'].iloc[0] == pytest.approx(718.0, abs=0.2)
    assert soja['PPBp'].iloc[0] == pytest.approx(466.6, abs=0.2)
    assert milho.attrs['PPf'] == pytest.approx(722, abs=0.5)
    assert soja.attrs['PPf'] == pytest.approx(458, abs=0.5)


def test_ppbp_funcao_isolada():
    qo, nN = 42.96, 0.501
    v = amp.PPBp(qo, nN, cTc=amp.cTc(26.5, 'C4'), cTn=amp.cTn(26.5, 'C4'))
    assert v == pytest.approx(718.0, abs=0.2)


def test_dias_dos_decendios_do_calendario():
    df = pd.DataFrame({'dia': [1, 11, 21, 1, 11, 21], 'mes': [1, 1, 1, 2, 2, 2],
                       'Tmed': 25.0, 'nN': 0.5, 'Qo': 40.0})
    pp = amp.produtividade_potencial(df, lat=-25, rota='C3', ano=2024)
    assert pp['ND'].tolist() == [10, 10, 11, 10, 10, 9]  # 2024 é bissexto
    assert pp['PPf_acum'].iloc[-1] == pytest.approx(pp['PPf_periodo'].sum())


def test_mensal_e_diario():
    mensal = pd.DataFrame({'dia': [1, 1], 'mes': [2, 3], 'Tmed': 25.0, 'nN': 0.5, 'Qo': 40.0})
    pp = amp.produtividade_potencial(mensal, lat=-25, rota='C3', intervalo='M')
    assert pp['ND'].tolist() == [28, 31]

    diario = pd.DataFrame({'Tmed': [25.0] * 3, 'nN': 0.5, 'Qo': 40.0})
    pp = amp.produtividade_potencial(diario, lat=-25, rota='C3', intervalo='d')
    assert pp['ND'].tolist() == [1, 1, 1]
    # O mesmo dia repetido: PPf diária x 3.
    assert pp.attrs['PPf'] == pytest.approx(3 * pp['PPf_periodo'].iloc[0])


def test_semeadura_e_ciclo_cortam_decendios_das_pontas():
    datas = pd.date_range('2023-10-01', '2024-02-21', freq='D')
    datas = datas[datas.day.isin([1, 11, 21])]
    df = pd.DataFrame({'data': datas, 'Tmed': 25.0, 'nN': 0.5, 'Qo': 40.0})
    # Semeadura em 05/10, 125 dias: até 06/02/2024.
    pp = amp.produtividade_potencial(df, lat=-25, rota='C3', semeadura='2023-10-05', ciclo=125)
    assert pp['ND'].iloc[0] == 6
    assert pp['ND'].iloc[-1] == 6
    assert pp['ND'].sum() == 125
    assert pp['data'].iloc[-1] == pd.Timestamp('2024-02-01')


def test_semeadura_fora_do_df():
    df = pd.DataFrame({'dia': [1, 11], 'mes': [1, 1], 'Tmed': 25.0, 'nN': 0.5, 'Qo': 40.0})
    with pytest.raises(ValueError):
        amp.produtividade_potencial(df, lat=-25, rota='C3', semeadura='2023-01-01', ciclo=60)


def test_cruza_o_ano_sem_coluna_ano():
    df = pd.DataFrame({'dia': [21, 1], 'mes': [12, 1], 'Tmed': 25.0, 'nN': 0.5, 'Qo': 40.0})
    pp = amp.produtividade_potencial(df, lat=-25, rota='C3', semeadura='2023-12-25', ciclo=10)
    assert pp['ND'].tolist() == [7, 3]


def test_pa_produtorio_e_etapa_unica():
    etc, etr = [60, 110, 120, 40], [57, 77, 102, 36]
    pa = amp.produtividade_atingivel(1000, etr, etc, ky=[0.4, 1.5, 0.5, 0.2])
    assert pa['fator'].tolist() == pytest.approx([0.980, 0.550, 0.925, 0.980], abs=1e-3)
    assert pa.attrs['EC'] == pytest.approx(0.489, abs=1e-3)
    assert pa['PA_final'].iloc[-1] == pytest.approx(pa.attrs['PA'])
    assert pa['perda'].sum() == pytest.approx(1000 - pa.attrs['PA'])

    unica = amp.produtividade_atingivel(1000, etr, etc, ky=1.25, metodo='etapa_unica')
    assert unica.attrs['EC'] == pytest.approx(0.780, abs=1e-3)


def test_pa_agrupa_periodos_por_fase():
    fase = ['veg', 'veg', 'flor', 'flor', 'form']
    etc = [30, 30, 55, 55, 120]
    etr = [28.5, 28.5, 38.5, 38.5, 102]
    pa = amp.produtividade_atingivel(1000, etr, etc, ky={'veg': 0.4, 'flor': 1.5, 'form': 0.5}, fase=fase)
    assert pa.index.tolist() == ['veg', 'flor', 'form']
    assert pa['ETc'].tolist() == [60, 110, 120]
    assert pa.attrs['EC'] == pytest.approx(0.98 * 0.55 * 0.925, abs=1e-3)


def test_pa_fator_nao_negativo_e_etc_zero():
    pa = amp.produtividade_atingivel(1000, [0, 0], [100, 0], ky=[1.5, 1.0])
    assert pa['fator'].tolist() == [0.0, 1.0]
    assert pa.attrs['PA'] == 0.0


def test_pa_ky_incompativel():
    with pytest.raises(ValueError):
        amp.produtividade_atingivel(1000, [1, 2], [2, 2], ky=[0.5])
