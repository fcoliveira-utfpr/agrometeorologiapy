"""English API (agrometeorologiapy.en): names, columns, keys and messages."""

import datetime

import numpy as np
import pandas as pd
import pytest

from agrometeorologiapy import en as agm

MONTHS = list(range(1, 13))
TMEAN = [24.7, 24.7, 23.4, 22.5, 18.9, 17.6, 17.6, 19.6, 21.4, 23.1, 23.9, 24.8]
TMIN = [20.4, 20.4, 18.6, 17.4, 13.9, 12.7, 12.0, 13.4, 15.5, 17.6, 18.5, 19.9]
TMAX = [29.7, 30.0, 29.3, 28.0, 24.9, 23.5, 23.8, 25.9, 26.9, 28.0, 28.9, 29.6]
P = [190, 160, 130, 110, 120, 100, 80, 70, 130, 180, 170, 180]
PET = [135, 115, 110, 75, 50, 38, 40, 55, 70, 100, 115, 130]


def test_radiation_and_humidity():
    J = agm.day_of_year(21, 5)
    decl = agm.solar_declination(J)
    ws = agm.sunset_hour_angle(-25.6, decl)
    Ra = agm.extraterrestrial_radiation(-25.6, decl, ws, agm.inverse_relative_distance(J))
    assert J == 141
    assert Ra == pytest.approx(22.84185036361831)
    assert agm.daylight_hours(ws) == pytest.approx(10.63576948863949)
    es = agm.saturation_vapour_pressure(25)
    assert agm.relative_humidity(agm.actual_vapour_pressure(es, RH=60), es) == pytest.approx(60)
    assert agm.dewpoint_temperature(agm.saturation_vapour_pressure(18.0)) == pytest.approx(18.0)


def test_energy_sign_convention():
    Rnl = agm.net_longwave_radiation(Tmax=28, Tmin=15, ea=2, Rs=20, Rso=30)
    assert Rnl < 0
    assert agm.net_radiation(agm.net_shortwave_radiation(20, albedo=0.23), Rnl) == pytest.approx(20 * 0.77 + Rnl)


def test_pet_monthly_columns():
    df = pd.DataFrame({'month': MONTHS, 'T_mean_C': TMEAN})
    th = agm.pet_thornthwaite(df, lat=-23.3)
    assert th.columns.tolist() == ['month', 'T_mean_C', 'PET_mm_month']
    ca = agm.pet_camargo(df, lat=-23.3)
    assert ca.columns.tolist() == ['month', 'T_mean_C', 'Ra_MJ_m2day', 'Ra_mm_day', 'PET_mm_month']


def test_degree_days_columns_and_messages(capsys):
    df = pd.DataFrame({'day': [1] * 12, 'month': MONTHS, 'Tmean': TMEAN, 'Tmax': TMAX, 'Tmin': TMIN})
    r = agm.maturity_date(df, Tb=14, TT=1030, sowing_day=12, sowing_month=11, interval='M')
    assert r.columns.tolist() == ['month', 'date', 'Tmean', 'GDD', 'GDD_period', 'GDD_cycle']
    assert r['month'].tolist() == ['Nov', 'Dec', 'Jan', 'Feb']
    assert r['date'].iloc[-1] == datetime.date(2024, 2, 18)
    out = capsys.readouterr().out
    assert out == "Sowing date: November 12\nPhysiological maturity date: February 18, 2024\n"
    r = agm.sowing_date(df, Tb=10, TT=800, maturity_day=22, maturity_month=6, interval='M', verbose=False)
    assert r['date'].iloc[-1] == datetime.date(2023, 4, 1)
    assert capsys.readouterr().out == ''


def test_degree_days_dekad_value():
    dates = pd.date_range('2023-01-01', '2023-12-31')
    dates = dates[dates.day.isin([1, 11, 21])]
    df = pd.DataFrame({'day': dates.day, 'month': dates.month, 'Tmean': 22.0, 'Tmax': 28.0, 'Tmin': 16.0})
    r = agm.maturity_date(df, 10, 1200, 11, 10, interval='dekad', verbose=False)
    assert r['GDD_cycle'].iloc[-1] >= 1200
    with pytest.raises(ValueError, match="interval must be"):
        agm.maturity_date(df, 10, 1200, 11, 10, interval='dec', verbose=False)


def test_water_balance():
    df = pd.DataFrame({'month': MONTHS, 'P (mm/month)': P, 'PET (mm/month)': PET})
    wb = agm.climatic_water_balance(df, AWC=100)
    assert wb.columns.tolist()[3:] == ['AWC', 'P-PET', 'STO (mm/month)', 'APWL (mm)', 'dSTO (mm/month)',
                                       'ETa (mm/month)', 'DEF (mm/month)', 'SUR (mm/month)']
    grid = agm.climatic_water_balance_grid(P, PET, AWC=100)
    assert list(grid) == ['P-PET', 'STO', 'APWL', 'dSTO', 'ETa', 'DEF', 'SUR']
    crop = agm.crop_water_balance(pd.DataFrame({'P': [52.3, 10.0], 'ETc': [12.5, 30.0], 'AWC': [18.0, 30.0]}))
    assert crop.columns.tolist() == ['P', 'ETc', 'AWC', 'P-ETc', 'STO', 'dSTO', 'ETa', 'DEF', 'SUR', 'WRSI']
    assert crop['WRSI'].iloc[0] == 1.0


def test_climate_classification():
    assert agm.koppen_classification(TMEAN, P, -23.3) == {'id': 9, 'class': 'Cfa'}
    th = agm.thornthwaite_classification(P, PET, -23.3)
    assert th['class'] == "B2rB'4a'" and th['description'] == 'Humid, mesothermal'
    assert 'SUR_annual' in th and 'PET_summer_pct' in th
    cam = agm.camargo_classification(TMEAN, P, PET, -23.3, thermal_table='consistent')
    assert cam['class'] == 'ST-PU' and cam['description'] == 'Subtropical perhumid'
    hd = agm.holdridge_classification(TMEAN, P, -23.3)
    assert hd['class'] == 'Subtropical moist forest'
    assert list(hd) == ['id', 'class', 'biotemperature', 'P_annual', 'PET_annual', 'PET_ratio']
    g = agm.holdridge_classification_grid(np.array([TMEAN, TMEAN]).T, np.array([P, P]).T, -23.3)
    assert g['class'].tolist() == ['Subtropical moist forest'] * 2


def test_crop_yield():
    df = pd.DataFrame({'day': [1], 'month': [1], 'Tmean': [26.5], 'Rs': [22.5]})
    yp = agm.potential_yield(df, lat=-24.86, pathway='C4', LAI=5.0, HI=0.35, moisture=13)
    assert yp.columns.tolist() == ['day', 'month', 'Tmean', 'Rs', 'ND', 'Ra', 'nN', 'cTo', 'cTc',
                                   'GPc', 'GPo', 'GP', 'CR', 'Yp_period', 'Yp_cum']
    assert yp.attrs['Yp'] == pytest.approx(722, abs=0.5) and yp.attrs['LAI_corr'] == 0.5
    etc, eta = [60, 110, 120, 40], [57, 77, 102, 36]
    yw = agm.attainable_yield(1000, eta, etc, ky=[0.4, 1.5, 0.5, 0.2])
    assert yw.columns.tolist() == ['ETc', 'ETa', 'ETa/ETc', 'ky', 'factor', 'Yw_start', 'Yw_end', 'loss', 'CE']
    assert yw.index.name == 'phase' and yw.attrs['CE'] == pytest.approx(0.489, abs=1e-3)
    single = agm.attainable_yield(1000, eta, etc, ky=1.25, method='single')
    assert single.index.tolist() == ['cycle'] and list(single.attrs) == ['Yp', 'Yw', 'CE']
