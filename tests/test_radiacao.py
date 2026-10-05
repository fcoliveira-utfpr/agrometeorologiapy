import pytest

import agrometeorologiapy as amp


def test_nda_doctest_examples():
    assert amp.nda(1, 1) == 1
    assert amp.nda(25, 12) == 359
    assert amp.nda(29, 2, 2024) == 60


def test_aplicacao_1_sombra_e_fotoperiodo():
    """Comprimento de sombra e fotoperíodo de um poste (Aplicação 1 do material)."""
    d, lat, dia, mes, hora, minuto = 10, -22, 2, 3, 10, 27

    NDA = amp.nda(dia, mes)
    delta = amp.declinacao_solar(NDA)
    h = amp.angulo_horario(hora, minuto)
    Z = amp.angulo_zenital(lat, delta)
    alfa = amp.azimute_solar(lat, delta, Z)
    S = amp.comprimento_sombra(d, Z)
    Hn = amp.angulo_horario_nascer(lat, delta)
    N = amp.fotoperiodo(Hn)

    assert NDA == 61
    assert delta == pytest.approx(-7.533773566685945)
    assert h == pytest.approx(-23.25000000000001)
    assert Z == pytest.approx(14.466226433314066)
    assert alfa == pytest.approx(179.9999973001307)
    assert S == pytest.approx(2.579887953183402)
    assert Hn == pytest.approx(93.06296504278654)
    assert N == pytest.approx(12.408395339038206)


def test_aplicacao_2_radiacao_global():
    """Radiação solar global por Angström-Prescott e Hargreaves-Samani (Aplicação 2)."""
    Tmax, Tmin, lat, dia, mes = 21.2, 7.4, -25.6, 21, 5

    NDA = amp.nda(dia, mes)
    delta = amp.declinacao_solar(NDA)
    Hn = amp.angulo_horario_nascer(lat, delta)
    N = amp.fotoperiodo(Hn)
    dD2 = amp.fator_correcao_distancia(NDA)
    Qo = amp.irradiancia_extraterrestre(lat, delta, Hn, dD2)
    insol = amp.insolacao(N, Tmax, Tmin, lat)
    Qg_AP = amp.Qg_angstrom(insol, N, Qo, lat, b=0.52)
    Qg_HS = amp.Qg_hargreaves(Tmax, Tmin, Qo)

    assert NDA == 141
    assert Qo == pytest.approx(22.84185036361831)
    assert insol == pytest.approx(9.087185925749393)
    assert Qg_AP == pytest.approx(16.122204526178894)
    assert Qg_HS == pytest.approx(13.576593285203279)
    # Inversa de Qg_angstrom: devolve a n/N usada (insol / N)
    assert amp.relacao_n_N(Qg_AP, Qo, lat) == pytest.approx(insol / N)


def test_relacao_n_N():
    # 1º decêndio de janeiro em Santa Helena-PR (Aula 03)
    assert amp.relacao_n_N(22.5, 42.96, -24.86) == pytest.approx(0.501, abs=1e-3)
    # Limitada entre 0 e 1, a menos que limitar=False
    assert amp.relacao_n_N(5.0, 40.0, -25.0) == 0.0
    assert amp.relacao_n_N(5.0, 40.0, -25.0, limitar=False) < 0
    assert amp.relacao_n_N(36.0, 40.0, -25.0) == 1.0
    # Coeficientes como entrada
    assert amp.relacao_n_N(20.0, 40.0, -25.0, b=0.5, a=0.25) == pytest.approx(0.5)
    # Arrays e Series
    import pandas as pd
    s = amp.relacao_n_N(pd.Series([15.0, 20.0], index=['a', 'b']), 40.0, -25.0)
    assert s.index.tolist() == ['a', 'b']
