import pytest

import agrometeorologiapy as amp


def test_es_tetens():
    assert amp.es_tetens(21.2) == pytest.approx(2.5175961687628154)


def test_ea_umidade_e_deficit():
    es = amp.es_tetens(21.2)
    ea = amp.ea_umidade(es, 65)
    assert ea == pytest.approx(es * 0.65)
    assert amp.deficit_saturacao(es, ea) == pytest.approx(es - ea)


def test_patm_altitude():
    assert amp.patm_altitude(235.09) == pytest.approx(98.55178157986391)


def test_umidade_absoluta_saturacao_relativa():
    es = amp.es_tetens(21.2)
    ea = amp.ea_umidade(es, 65)
    UA = amp.umidade_absoluta(ea, 21.2)
    US = amp.umidade_saturacao(es, 21.2)
    UR = amp.umidade_relativa(ea, es)
    assert UR == pytest.approx(65.0)
    assert UA == pytest.approx(US * 0.65)


def test_ponto_orvalho_inverte_tetens():
    es = amp.es_tetens(21.2)
    ea = amp.ea_umidade(es, 65)
    To = amp.ponto_orvalho(ea)
    assert amp.es_tetens(To) == pytest.approx(ea)


def test_constante_psicrometrica():
    Patm = amp.patm_altitude(235.09)
    assert amp.constante_psicrometrica(Patm) == pytest.approx(0.06553693475060951)
