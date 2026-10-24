"""Camada em português sobre o núcleo em inglês (agrometeorologiapy.en)."""

import inspect

import numpy as np
import pandas as pd
import pytest

import agrometeorologiapy as amp
from agrometeorologiapy import en as agm
from agrometeorologiapy._i18n import traduzir

MODULOS_PT = [amp.radiacao, amp.temperatura, amp.umidade, amp.energia, amp.evapotranspiracao,
              amp.grau_dias, amp.balanco_hidrico, amp.classificacao_climatica, amp.produtividade]


def _funcoes_traduzidas():
    for mod in MODULOS_PT:
        for nome in mod.__all__:
            yield mod, nome, getattr(mod, nome)


@pytest.mark.parametrize("mod,nome,f", list(_funcoes_traduzidas()), ids=lambda x: getattr(x, '__name__', str(x)))
def test_wrapper_consistente_com_o_nucleo(mod, nome, f):
    en = f.__en__
    assert en.__name__ in agm.__all__ and getattr(agm, en.__name__) is en
    assert f.__name__ == nome and f.__module__ == mod.__name__
    assert f.__doc__ and f.__doc__.strip()          # docstring em português
    sig_pt, sig_en = inspect.signature(f), inspect.signature(en)
    params_en = [p for p in sig_en.parameters.values() if p.name not in f.__fixos__]
    assert len(sig_pt.parameters) == len(params_en)
    assert [p.kind for p in sig_pt.parameters.values()] == [p.kind for p in params_en]
    assert getattr(amp, nome) is f                  # exportada no nível do pacote


def test_correspondencia_um_para_um():
    pt = {getattr(amp, n).__en__.__name__ for n in amp.__all__ if hasattr(getattr(amp, n), '__en__')}
    en = set(agm.__all__) - {'sind', 'cosd', 'tand', 'acosd'}
    assert pt == en
    assert len(amp.__all__) == len(agm.__all__)


def test_resultados_iguais_nas_duas_apis():
    lat, nda = -25.09, amp.nda(15, 9)
    assert nda == agm.day_of_year(15, 9)
    dec = amp.declinacao_solar(nda)
    assert dec == agm.solar_declination(nda)
    Hn = amp.angulo_horario_nascer(lat, dec)
    assert Hn == agm.sunset_hour_angle(lat, dec)
    Qo = amp.irradiancia_extraterrestre(lat, dec, Hn, amp.fator_correcao_distancia(nda))
    Ra = agm.extraterrestrial_radiation(lat, dec, Hn, agm.inverse_relative_distance(nda))
    assert Qo == Ra
    assert amp.Qg_angstrom(insolacao=8, N=12, Qo=Qo, lat=lat) == agm.solar_radiation_angstrom(n=8, N=12, Ra=Ra, lat=lat)
    assert amp.relacao_n_N(Qg=20, Qo=Qo, lat=lat, limitar=False) == agm.relative_sunshine_duration(Rs=20, Ra=Ra, lat=lat, clip=False)


def test_argumentos_por_posicao_nome_e_arrays():
    assert amp.nda(29, 2, 2024) == amp.nda(dia=29, mes=2, ano=2024) == 60
    nda = np.array([1, 100, 200])
    assert np.allclose(amp.declinacao_solar(NDA=nda), agm.solar_declination(doy=nda))
    s = amp.relacao_n_N(pd.Series([15.0, 20.0], index=['a', 'b']), 40.0, -25.0)
    assert s.index.tolist() == ['a', 'b']


def test_cada_api_aceita_so_os_seus_nomes():
    with pytest.raises(TypeError):
        amp.nda(day=1, month=1)
    with pytest.raises(TypeError):
        agm.day_of_year(dia=1, mes=1)


def test_assinatura_em_portugues():
    assert list(inspect.signature(amp.irradiancia_extraterrestre).parameters) == ['lat', 'declinacao', 'Hn', 'dD2']
    assert list(inspect.signature(agm.extraterrestrial_radiation).parameters) == ['lat', 'declination', 'ws', 'dr']
    assert inspect.signature(amp.relacao_n_N).parameters['limitar'].default is True


# ---------------------------------------------------------------- tradutor
def _soma(df, interval='dekad', scale=1.0):
    out = pd.DataFrame({'total': (df['rain'] + df['irrigation']) * scale})
    out.attrs['label'] = interval
    return out


def test_traduz_colunas_valores_e_attrs():
    soma = traduzir(_soma, 'soma', {'intervalo': 'interval', 'escala': 'scale'}, modulo=__name__,
                    valores={'interval': {'dec': 'dekad', 'M': 'month'}},
                    colunas_entrada={'chuva': 'rain', 'irrigacao': 'irrigation'},
                    colunas_saida={'total': 'soma', 'label': 'rotulo'})
    assert inspect.signature(soma).parameters['intervalo'].default == 'dec'
    df = pd.DataFrame({'chuva': [1.0, 2.0], 'irrigacao': [0.5, 0.0]})
    r = soma(df, intervalo='M', escala=2)
    assert r['soma'].tolist() == [3.0, 4.0]
    assert r.attrs == {'rotulo': 'month'}
    assert df.columns.tolist() == ['chuva', 'irrigacao']   # entrada não é alterada
    assert soma(df).attrs['rotulo'] == 'dekad'


def test_parametro_inexistente_no_nucleo():
    with pytest.raises(ValueError):
        traduzir(_soma, 'soma', {'x': 'nao_existe'}, modulo=__name__)
