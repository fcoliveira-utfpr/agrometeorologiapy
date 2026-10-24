import numpy as np
import pytest

import agrometeorologiapy as amp
from agrometeorologiapy import en as agm

OBS = np.array([1.0, 2.0, 3.0, 4.0])
EST = np.array([1.5, 2.0, 3.5, 3.0])


def test_erros():
    assert amp.mbe(OBS, EST) == pytest.approx(0.0)
    assert amp.mbe(OBS, EST + 1) == pytest.approx(1.0)       # superestima
    assert amp.mae(OBS, EST) == pytest.approx(0.5)
    assert amp.rmse(OBS, EST) == pytest.approx(np.sqrt(1.5 / 4))


def test_indices_de_concordancia():
    # r: scipy.stats.pearsonr; d e NSE pelas definições
    assert amp.r_pearson(OBS, EST) == pytest.approx(0.8485281374238571)
    assert amp.d_willmott(OBS, EST) == pytest.approx(1 - 1.5 / 13.5)
    assert amp.indice_confianca(OBS, EST) == pytest.approx(0.8485281374238571 * (1 - 1.5 / 13.5))
    assert amp.nse(OBS, EST) == pytest.approx(1 - 1.5 / 5)


def test_modelo_perfeito_e_media():
    assert amp.d_willmott(OBS, OBS) == amp.nse(OBS, OBS) == amp.kge(OBS, OBS) == 1.0
    assert amp.rmse(OBS, OBS) == 0.0
    media = np.full(4, OBS.mean())
    assert amp.nse(OBS, media) == pytest.approx(0.0)


def test_kge():
    r = 0.8485281374238571
    a, b = EST.std() / OBS.std(), EST.mean() / OBS.mean()
    assert amp.kge(OBS, EST) == pytest.approx(1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + (b - 1) ** 2))
    # 2012: razão dos CVs; com viés multiplicativo puro, a = 1 e só b pesa.
    assert amp.kge(OBS, 1.2 * OBS, versao=2012) == pytest.approx(0.8)
    assert amp.kge(OBS, 1.2 * OBS, versao=2009) == pytest.approx(1 - np.sqrt(0.08))
    with pytest.raises(ValueError, match="version must be"):
        agm.kge(OBS, EST, version=2021)


def test_categoria_indice_confianca():
    valores = [0.86, 0.85, 0.76, 0.754, 0.66, 0.61, 0.51, 0.41, 0.40, np.nan]
    assert amp.categoria_indice_confianca(valores).tolist() == [
        'Ótimo', 'Muito bom', 'Muito bom', 'Bom', 'Bom', 'Mediano', 'Sofrível', 'Mau', 'Péssimo', '']
    assert agm.confidence_index_category(0.9) == 'Excellent'


def test_nan_e_tamanhos():
    obs = [1.0, np.nan, 3.0, 4.0, 2.0]
    est = [1.5, 2.0, np.nan, 3.0, 2.0]
    assert amp.avaliar_modelo(obs, est)['n'] == 3
    assert amp.rmse(obs, est) == pytest.approx(np.sqrt((0.25 + 1.0 + 0.0) / 3))
    with pytest.raises(ValueError, match="same length"):
        amp.rmse([1, 2], [1, 2, 3])
    with pytest.raises(ValueError, match="without NaN"):
        amp.rmse([np.nan], [1.0])


def test_avaliar_modelo_igual_nas_duas_apis():
    pt = amp.avaliar_modelo(OBS, est=EST)
    en = agm.evaluate_model(OBS, pred=EST)
    assert list(pt) == ['n', 'MBE', 'MAE', 'RMSE', 'r', 'd', 'c', 'NSE', 'KGE']
    assert pt == en
