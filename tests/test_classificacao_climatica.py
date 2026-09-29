import numpy as np
import pytest

import agrometeorologiapy as amp

# Normais ilustrativas (T em °C, P em mm/mês, jan-dez).
PONTA_GROSSA = ([21.6, 21.7, 20.8, 18.6, 15.7, 14.2, 13.9, 15.2, 16.3, 18.2, 19.6, 21.0],
                [180, 160, 130, 100, 110, 110, 95, 80, 140, 150, 130, 160], -25.1)
LONDRINA = ([24.7, 24.7, 23.4, 22.5, 18.9, 17.6, 17.6, 19.6, 21.4, 23.1, 23.9, 24.8],
            [210, 180, 140, 110, 110, 90, 60, 50, 120, 150, 160, 200], -23.3)
MANAUS = ([26.1, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 27.1, 27.5, 27.4, 27.1, 26.6],
          [264, 295, 315, 300, 255, 114, 87, 58, 83, 126, 183, 217], -3.1)
PETROLINA = ([27.5, 27.2, 27.0, 26.5, 25.6, 24.6, 24.2, 24.9, 26.4, 27.8, 28.3, 28.0],
             [95, 90, 130, 80, 25, 12, 8, 3, 5, 12, 50, 80], -9.4)
LISBOA = ([11.6, 12.8, 14.8, 15.9, 18.3, 21.4, 23.2, 23.6, 22.2, 19.1, 15.1, 12.4],
          [99, 87, 57, 68, 55, 15, 5, 6, 33, 95, 114, 120], 38.7)

# Dados do BHC da Aplicação 10 (test_balanco_hidrico.py).
P_UMIDO = [203.4, 189.53, 196.27, 116.93, 14.03, 3.33, 1.3, 4.37, 18.03, 99.47, 223.73, 234.97]
ETP_UMIDO = [128.73, 106.29, 117.3, 109.26, 106.8, 102.03, 114.99, 134.52, 141.43, 149.15, 118.62, 119.96]


@pytest.mark.parametrize("local, esperado", [
    (PONTA_GROSSA, 'Cfb'), (LONDRINA, 'Cfa'), (MANAUS, 'Am'), (PETROLINA, 'BSh'), (LISBOA, 'Csa'),
])
def test_koppen(local, esperado):
    T, P, lat = local
    assert amp.classificacao_koppen(T, P, lat)['classe'] == esperado


def test_koppen_troca_estacoes_no_hemisferio_sul():
    # Mesmo regime de Lisboa, mas com a seca em out-mar (inverno no sul): Cwa.
    T, P, _ = LISBOA
    assert amp.classificacao_koppen(T, P, -38.7)['classe'] == 'Cwa'


def test_koppen_grade_igual_ao_local_e_nan():
    locais = [PONTA_GROSSA, LONDRINA, MANAUS, PETROLINA, LISBOA]
    T = np.array([l[0] for l in locais], dtype=float).T
    P = np.array([l[1] for l in locais], dtype=float).T
    lat = np.array([l[2] for l in locais])
    T[5, 0] = np.nan
    saida = amp.classificacao_koppen_grade(T, P, lat)
    assert saida['classe'].tolist() == ['', 'Cfa', 'Am', 'BSh', 'Csa']
    assert saida['id'][0] == 0


def test_thornthwaite():
    # Ih = 25,2; Ia = 35,1; Im = 4,1 -> C2; Ia >= 33,3 com DEF no inverno -> w2;
    # ETP anual 1449 mm -> A'; ETP do verão = 24,4% -> a'.
    r = amp.classificacao_thornthwaite(P_UMIDO, ETP_UMIDO, lat=-23)
    assert r['classe'] == "C2w2A'a'"
    assert r['Im'] == pytest.approx(4.12, abs=0.01)
    assert r['ETP_verao_pct'] == pytest.approx(24.37, abs=0.01)


def test_thornthwaite_limite_855_mm():
    P = np.full(12, 200.0)
    for etp_anual, termica in [(860, "B'3"), (850, "B'2")]:
        r = amp.classificacao_thornthwaite(P, np.full(12, etp_anual / 12), lat=-23)
        assert r['classe'][2:5] == termica


def test_camargo():
    T, _, _ = LONDRINA
    r = amp.classificacao_camargo(T, P_UMIDO, ETP_UMIDO, lat=-23)
    # Ty = 21,9 e Tcold = 17,6 -> ST; DEF 509 e EXC 365 mm -> MO; seca no inverno -> i
    assert r['classe'] == 'ST-MOi'


def test_camargo_equatorial_coerente_x_literal():
    T, P, lat = MANAUS
    ETP = np.full(12, 130.0)
    assert amp.classificacao_camargo(T, P, ETP, lat)['classe'].startswith('EQ-')
    assert amp.classificacao_camargo(T, P, ETP, lat, tabela_termica='literal')['classe'].startswith('TR-')


def test_camargo_casos_fora_da_tabela_7_viram_subumido():
    T = np.full(12, 20.0)
    # DEF pequena e sem excedente
    P = np.array([100.0] * 6 + [85.0] * 6)
    r = amp.classificacao_camargo(T, P, np.full(12, 100.0), lat=-23)
    assert r['EXC_anual'] == pytest.approx(0.0, abs=1e-6)
    assert 0 < r['DEF_anual'] <= 150
    assert r['hidrica'] == 5
    # sem DEF e excedente pequeno
    r = amp.classificacao_camargo(T, np.full(12, 110.0), np.full(12, 100.0), lat=-23)
    assert r['DEF_anual'] == pytest.approx(0.0, abs=1e-6)
    assert r['classe'] == 'ST-SB'


def test_holdridge():
    r = amp.classificacao_holdridge(np.full(12, 26.0), np.full(12, 2000 / 12), lat=0)
    assert r['biotemperatura'] == pytest.approx(26.0)
    assert r['ETP_anual'] == pytest.approx(58.93 * 26)
    assert r['classe'] == 'Floresta úmida tropical'
    assert r['classe_en'] == 'Tropical moist forest'
    r = amp.classificacao_holdridge(np.full(12, 26.0), np.full(12, 1000 / 12), lat=0)
    assert r['classe'] == 'Floresta seca tropical'


def test_holdridge_correcao_por_latitude():
    # Com limiar 24 °C, só os meses acima de 24 são corrigidos: 28 - 0,6 * 16 = 18,4.
    T = np.array([28.0] * 6 + [20.0] * 6)
    r = amp.classificacao_holdridge(T, np.full(12, 100.0), lat=-20)
    assert r['biotemperatura'] == pytest.approx((18.4 * 6 + 20 * 6) / 12)
    r = amp.classificacao_holdridge(T, np.full(12, 100.0), lat=-20, limiar_correcao=None)
    assert r['biotemperatura'] == pytest.approx((18.4 * 6 + 10.4 * 6) / 12)


def test_versao_local_rejeita_varios_locais():
    with pytest.raises(ValueError):
        amp.classificacao_koppen(np.zeros((12, 2)), np.zeros((12, 2)), -20)
