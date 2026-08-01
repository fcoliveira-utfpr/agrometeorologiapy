import pytest

import agrometeorologiapy as amp


def test_temp_media_extremos():
    assert amp.temp_media_extremos(21.2, 7.4) == pytest.approx(14.3)


def test_temp_media_estacao_automatica():
    assert amp.temp_media_estacao_automatica([20, 22, 24]) == pytest.approx(22.0)
