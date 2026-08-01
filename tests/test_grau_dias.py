import datetime

import pandas as pd
import pytest

import agrometeorologiapy as amp

MESES = list(range(1, 13))
TMED = [24.7, 24.7, 23.4, 22.5, 18.9, 17.6, 17.6, 19.6, 21.4, 23.1, 23.9, 24.8]
TMIN = [20.4, 20.4, 18.6, 17.4, 13.9, 12.7, 12.0, 13.4, 15.5, 17.6, 18.5, 19.9]
TMAX = [29.7, 30.0, 29.3, 28.0, 24.9, 23.5, 23.8, 25.9, 26.9, 28.0, 28.9, 29.6]


def test_aplicacao_8_data_maturacao_fisiologica():
    dia_semeadura, mes_semeadura = 12, 11
    df = pd.DataFrame({
        'dia': [1] * dia_semeadura,
        'mes': MESES,
        'Tmed': TMED,
        'Tmax': TMAX,
        'Tmin': TMIN,
    })
    resultado = amp.data_maturacao_fisiologica(df, Tb=14, CT=1030,
                                                dia_semeadura=dia_semeadura,
                                                mes_semeadura=mes_semeadura,
                                                intervalo='M')
    assert resultado['data'].tolist()[-1] == datetime.date(2023, 2, 1)
    assert resultado['GD_ciclo'].tolist() == pytest.approx([178.2, 513.0, 844.7, 1144.3])


def test_aplicacao_9_data_semeadura():
    dia_maturacao, mes_maturacao = 22, 6
    df = pd.DataFrame({
        'dia': [1] * 12,
        'mes': MESES,
        'Tmed': TMED,
        'Tmax': TMAX,
        'Tmin': TMIN,
    })
    resultado = amp.data_semeadura(df, Tb=10, CT=800,
                                    dia_maturacao=dia_maturacao,
                                    mes_maturacao=mes_maturacao,
                                    intervalo='M')
    assert resultado['data'].tolist()[-1] == datetime.date(2023, 4, 1)
    assert resultado['GD_ciclo'].tolist() == pytest.approx([167.2, 443.1, 818.1])
