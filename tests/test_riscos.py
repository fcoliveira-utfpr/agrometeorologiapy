import numpy as np
import pytest

import agrometeorologiapy as amp
from agrometeorologiapy import en as agm

TMIN = [1.2, 4.5, -0.8, 3.1, 2.4, 5.6, 0.3, 6.2, 2.9, 4.0, np.nan]


def test_probabilidade_geada():
    # scipy.stats.norm.cdf(3, média, desvio padrão amostral)
    assert amp.probabilidade_geada(TMIN) == pytest.approx(0.5106370250384482, abs=1e-12)
    assert amp.probabilidade_geada(TMIN, limiar=0.0) < amp.probabilidade_geada(TMIN)
    assert amp.probabilidade_geada(TMIN, metodo='empirico') == 0.5
    assert agm.frost_probability(TMIN, method='empirical') == 0.5
    with pytest.raises(ValueError, match="method must be"):
        amp.probabilidade_geada(TMIN, metodo='gumbel')


def test_periodo_retorno():
    assert amp.periodo_retorno(0.25) == 4.0
    assert amp.periodo_retorno([0.5, 0.0]).tolist() == [2.0, np.inf]


def test_indice_temperatura_umidade_e_categoria():
    # (1,8*30 + 32) - (0,55 - 0,0055*70)(1,8*30 - 26) = 86 - 0,165*28
    assert amp.indice_temperatura_umidade(30, 70) == pytest.approx(81.38)
    # Com UR = 100 %, o ITU é a temperatura em °F.
    assert amp.indice_temperatura_umidade(25, 100) == pytest.approx(77.0)
    assert [amp.categoria_itu(v) for v in (74.9, 75, 78.9, 79, 84)] == \
        ['Normal', 'Alerta', 'Alerta', 'Perigo', 'Emergência']
    assert agm.thi_category(np.array([70.0, 90.0])).tolist() == ['Normal', 'Emergency']


def test_estresse_calor():
    Tmax = np.array([[31.0, 35.0], [34.5, 33.0], [36.0, np.nan], [33.0, 37.0]])  # 4 dias x 2 locais
    assert amp.dias_estresse_calor(Tmax, Tcrit=34).tolist() == [2, 2]
    assert amp.graus_dia_estresse_calor(Tmax, Tcrit=34) == pytest.approx([2.5, 4.0])
    assert amp.dias_estresse_calor([30, 31], 34) == 0


def test_duracao_molhamento_foliar():
    UR = [95, 96, 97, 93, 88, 75, 62, 55, 58, 70, 84, 91, 90, 89, 92, np.nan]
    assert amp.duracao_molhamento_foliar(UR) == 7.0
    assert amp.duracao_molhamento_foliar(UR, limiar=88) == 9.0
    assert agm.leaf_wetness_duration(UR, interval=0.5) == 3.5
    dias = np.array(UR[:12] * 2).reshape(2, 12).T  # 12 observações x 2 dias
    assert amp.duracao_molhamento_foliar(dias).tolist() == [5.0, 5.0]
