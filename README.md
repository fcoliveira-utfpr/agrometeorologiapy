# agrometeorologiapy

[![PyPI](https://img.shields.io/pypi/v/agrometeorologiapy.svg)](https://pypi.org/project/agrometeorologiapy/)
[![CI](https://github.com/fcoliveira-utfpr/agrometeorologiapy/actions/workflows/ci.yml/badge.svg)](https://github.com/fcoliveira-utfpr/agrometeorologiapy/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/LICENSE)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab_en.ipynb)

*[🇧🇷 Português](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/LEIAME.md)*

Agrometeorology formulas in Python: solar radiation and relative sunshine
duration, air temperature and humidity, radiation balance, evapotranspiration
(Thornthwaite, Camargo-Maluf, Hargreaves-Samani, Priestley-Taylor, FAO
Penman-Monteith), growing degree days, soil water balance (climatic and crop),
climate classification (Köppen-Geiger, Thornthwaite, Camargo and Holdridge)
and potential and attainable crop yield (FAO Agro-Ecological Zone model).

Terms and symbols follow FAO Irrigation and Drainage Paper 56 (Allen et al.,
1998) wherever they exist. Every function is available in two languages over
a single implementation: an English API, `agrometeorologiapy.en`, and a
Portuguese API at the top level of the package (see
[Portuguese API](#portuguese-api)).

## Installation

```bash
pip install agrometeorologiapy
```

Requires Python 3.9+, NumPy and pandas.

## Interactive tutorial

Want to see the package in action, function by function, without installing
anything? Open the tutorial in Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab_en.ipynb)

## Quick start

```python
from agrometeorologiapy import en as agm
```

### Solar radiation and relative sunshine duration

```python
lat = -24.86                                   # Santa Helena, Paraná, Brazil
J = agm.day_of_year(5, 1)                      # January 5 -> 5
decl = agm.solar_declination(J)                # -22.54°
ws = agm.sunset_hour_angle(lat, decl)
N = agm.daylight_hours(ws)                     # 13.48 h
Ra = agm.extraterrestrial_radiation(lat, decl, ws, agm.inverse_relative_distance(J))  # 42.98 MJ m-2 day-1

# n/N without a sunshine recorder: the Angstrom formula solved for n/N from measured Rs
nN = agm.relative_sunshine_duration(Rs=22.5, Ra=Ra, lat=lat)   # 0.501
```

### Evapotranspiration and water balance

```python
Tmax, Tmin = 31.2, 20.4
Tmean = agm.mean_temperature(Tmax, Tmin)
pet = agm.pet_hargreaves_samani(Ra, Tmax, Tmin, Tmean)   # 5.78 mm/day

import pandas as pd
df = pd.DataFrame({
    'month': ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
    'P (mm/month)':   [190, 160, 130, 110, 120, 100, 80, 70, 130, 180, 170, 180],
    'PET (mm/month)': [135, 115, 110, 75, 50, 38, 40, 55, 70, 100, 115, 130],
})
wb = agm.climatic_water_balance(df, AWC=100)   # STO, APWL, ETa, DEF and SUR per month
```

### Potential and attainable yield (FAO Agro-Ecological Zone model)

Soybean (C3) sown on 2023-10-15, 125-day cycle, with **monthly** data
(`interval` accepts `'d'`, `'dekad'` or `'M'`). The sowing and harvest months
only count the days with the crop in the field, and Ra is computed from the
dates:

```python
season = pd.DataFrame({
    'day':   [1, 1, 1, 1, 1],
    'month': [10, 11, 12, 1, 2],
    'Tmean': [22.8, 24.1, 25.6, 26.3, 25.9],     # °C
    'Rs':    [19.8, 22.4, 23.1, 22.6, 21.0],     # MJ m-2 day-1
})
yp = agm.potential_yield(season, lat=-24.86, pathway='C3', LAI=4.5, HI=0.35, moisture=13,
                         interval='M', sowing='2023-10-15', cycle=125)
yp[['month', 'ND', 'nN', 'cTc', 'cTo', 'GP', 'Yp_cum']]   # potential yield accumulated month by month
Yp = yp.attrs['Yp']                                        # 5531 kg/ha

# Water-deficit penalty, with ETc and ETa of each phase
# (e.g. summed from agm.crop_water_balance)
ETc, ETa = [190, 110, 175, 60], [185, 80, 140, 55]

# By phase (product): one ky per phase, applied in sequence
yw = agm.attainable_yield(Yp, ETa, ETc, ky=[0.2, 0.8, 1.0, 0.0])
yw.attrs['Yw'], yw.attrs['CE']                             # 3441 kg/ha, 0.622

# Single step: whole-cycle ky
yw1 = agm.attainable_yield(Yp, ETa, ETc, ky=0.85, method='single')
yw1.attrs['Yw'], yw1.attrs['CE']                           # 4872 kg/ha, 0.881
```

Names follow the yield gap literature (van Ittersum et al., 2013): Yp is the
potential yield, Yw the attainable (water-limited) yield and CE = Yw/Yp the
climatic efficiency. Each formula of the model also has its own function
(`temperature_correction_overcast` and `temperature_correction_clear`, with
`pathway='C3'` or `'C4'`, `gross_potential_productivity`, `lai_correction`,
`respiration_correction`), and every constant or tabulated value is a
parameter with the usual value as default (e.g. `moisture=13`, `b=0.52`).

Functions can also be reached by submodule (`agm.radiation`, `agm.humidity`,
`agm.evapotranspiration`, etc.).

## Portuguese API

The same functions are available with Portuguese names (following Pereira,
Angelocci & Sentelhas, 2002, the reference textbook in Brazil) at the top
level of the package. Both APIs share one implementation and give identical
results; the Portuguese API is permanent.

```python
import agrometeorologiapy as amp

amp.relacao_n_N(Qg=22.5, Qo=42.98, lat=-24.86)                 # 0.501
agm.relative_sunshine_duration(Rs=22.5, Ra=42.98, lat=-24.86)  # 0.501
```

Full name mapping: [`docs/API_PT_EN.md`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/docs/API_PT_EN.md).

## Modules

| Module (en · pt) | Contents |
| --- | --- |
| `en.radiation` · `radiacao` | Day of the year, solar declination, hour angle, zenith angle, solar azimuth, daylight hours, extraterrestrial radiation, solar radiation (Angstrom, Hargreaves) and relative sunshine duration n/N (inverted Angstrom formula) |
| `en.temperature` · `temperatura` | Mean daily temperature (from extremes and from station observations) |
| `en.humidity` · `umidade` | Saturation and actual vapour pressure, vapour pressure deficit, atmospheric pressure, absolute, saturation and relative humidity, dewpoint temperature, psychrometric constant |
| `en.energy` · `energia` | Net shortwave, net longwave and net radiation |
| `en.evapotranspiration` · `evapotranspiracao` | PET by Thornthwaite, Camargo-Maluf, Hargreaves-Samani and Priestley-Taylor; ETo by FAO Penman-Monteith; wind speed at 2 m |
| `en.degree_days` · `grau_dias` | Physiological maturity date / sowing date by accumulated growing degree days |
| `en.water_balance` · `balanco_hidrico` | Climatic and crop water balance (Thornthwaite & Mather), with the WRSI |
| `en.climate_classification` · `classificacao_climatica` | Köppen-Geiger, Thornthwaite, Camargo and Holdridge, for one location or on a grid (rasters) |
| `en.crop_yield` · `produtividade` | FAO Agro-Ecological Zone model (Doorenbos & Kassam): temperature corrections for C3/C4 crops, gross potential productivity, LAI and respiration corrections, potential yield (daily, dekadal or monthly) and attainable yield (by phase or in a single step) with climatic efficiency |

- [`docs/FORMULAS.en.md`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/docs/FORMULAS.en.md) ([Português](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/docs/FORMULAS.md)) — mathematical documentation of every function: the original formula, variables and units.
- [`examples/tutorial_colab_en.ipynb`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab_en.ipynb) ([Português](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab.ipynb)) — guided tour, ready for Colab, with a runnable example for every function.
- [`examples/formulas_agrometeorologia.ipynb`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/formulas_agrometeorologia.ipynb) — original development notebook (in Portuguese), with the derivation of each formula.

## Development

```bash
git clone https://github.com/fcoliveira-utfpr/agrometeorologiapy.git
cd agrometeorologiapy
pip install -e ".[dev]"
pytest
```

## License

BSD-3-Clause. See [LICENSE](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/LICENSE).
