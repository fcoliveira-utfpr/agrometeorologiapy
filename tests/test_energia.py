import pytest

import agrometeorologiapy as amp


def test_boc_saldo():
    assert amp.boc_saldo(13.576593285203279, r=0.23) == pytest.approx(10.453976829606525)


def test_bol_saldo():
    """Tmax/Tmin em °C: a conversão para Kelvin é feita internamente (padrão FAO-56)."""
    ea = 1.3407926169150235
    Qg = 13.576593285203279
    Qg_cs = 17.23878558475339
    BOL = amp.bol_saldo(21.2, 7.4, ea, Qg, Qg_cs)
    assert BOL == pytest.approx(-4.261654877504128)


def test_saldo_radiacao():
    assert amp.saldo_radiacao(10.453976829606525, -4.261654877504128) == pytest.approx(6.192321952102397)
