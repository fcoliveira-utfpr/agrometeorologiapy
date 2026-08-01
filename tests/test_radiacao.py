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
