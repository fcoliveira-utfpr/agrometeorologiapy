# agrometeorologiapy

[![PyPI](https://img.shields.io/pypi/v/agrometeorologiapy.svg)](https://pypi.org/project/agrometeorologiapy/)
[![CI](https://github.com/fcoliveira-utfpr/agrometeorologiapy/actions/workflows/ci.yml/badge.svg)](https://github.com/fcoliveira-utfpr/agrometeorologiapy/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/LICENSE)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab.ipynb)

Fórmulas de agrometeorologia em Python: radiação solar, temperatura, umidade
do ar, balanço de energia, evapotranspiração (Thornthwaite, Camargo-Maluf,
Hargreaves-Samani, Priestley-Taylor, Penman-Monteith FAO-56), grau-dias,
balanço hídrico (climatológico e de cultura), classificação climática
(Köppen-Geiger, Thornthwaite, Camargo e Holdridge) e produtividade potencial
e atingível (Zona Agroecológica da FAO).

## Instalação

```bash
pip install agrometeorologiapy
```

## Tutorial interativo

Quer ver o pacote em ação, função por função, sem instalar nada? Abra o
tutorial direto no Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab.ipynb)

## Uso rápido

```python
import agrometeorologiapy as amp
```

### Radiação solar e razão de insolação

```python
lat = -24.86                                   # Santa Helena-PR
nda = amp.nda(5, 1)                            # 5 de janeiro -> 5
delta = amp.declinacao_solar(nda)              # -22,54°
Hn = amp.angulo_horario_nascer(lat, delta)
N = amp.fotoperiodo(Hn)                        # 13,48 h
Qo = amp.irradiancia_extraterrestre(lat, delta, Hn, amp.fator_correcao_distancia(nda))  # 42,98 MJ/m² dia

# n/N sem heliógrafo: Angström-Prescott invertida a partir da radiação global medida
nN = amp.relacao_n_N(Qg=22.5, Qo=Qo, lat=lat)  # 0,501
```

### Evapotranspiração e balanço hídrico

```python
Tmax, Tmin = 31.2, 20.4
Tmed = amp.temp_media_extremos(Tmax, Tmin)
etp = amp.etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)   # 5,78 mm/dia

import pandas as pd
df = pd.DataFrame({
    'Meses': ['jan', 'fev', 'mar', 'abr', 'mai', 'jun', 'jul', 'ago', 'set', 'out', 'nov', 'dez'],
    'P (mm/mês)':   [190, 160, 130, 110, 120, 100, 80, 70, 130, 180, 170, 180],
    'ETP (mm/mês)': [135, 115, 110, 75, 50, 38, 40, 55, 70, 100, 115, 130],
})
bhc = amp.balanco_hidrico_climatologico(df, CAD=100)    # ARM, ETR, DEF e EXC por mês
```

### Produtividade potencial e atingível (Zona Agroecológica da FAO)

Soja (C3) semeada em 15/10/2023, ciclo de 125 dias, com dados **mensais**
(`intervalo` aceita `'d'`, `'dec'` ou `'M'`). Os meses da semeadura e da
colheita contam só os dias com a cultura no campo, e Qo é calculada a partir
das datas:

```python
safra = pd.DataFrame({
    'dia':  [1, 1, 1, 1, 1],
    'mes':  [10, 11, 12, 1, 2],
    'Tmed': [22.8, 24.1, 25.6, 26.3, 25.9],      # °C
    'Qg':   [19.8, 22.4, 23.1, 22.6, 21.0],      # MJ/m² dia
})
pp = amp.produtividade_potencial(safra, lat=-24.86, rota='C3', IAF=4.5, Cc=0.35, U=13,
                                 intervalo='M', semeadura='2023-10-15', ciclo=125)
pp[['mes', 'ND', 'nN', 'cTc', 'cTn', 'PPBp', 'PPf_acum']]   # PPf acumulada mês a mês
PPf = pp.attrs['PPf']                                      # 5531 kg/ha

# Penalização pelo déficit hídrico, com ETc e ETR de cada fase
# (ex.: somadas a partir de amp.balanco_hidrico_cultura)
ETc, ETR = [190, 110, 175, 60], [185, 80, 140, 55]

# Por fase (produtório): um ky por fase, aplicado em sequência
pa = amp.produtividade_atingivel(PPf, ETR, ETc, ky=[0.2, 0.8, 1.0, 0.0])
pa.attrs['PA'], pa.attrs['EC']                             # 3441 kg/ha, 0,622

# Etapa única: ky do ciclo
pa1 = amp.produtividade_atingivel(PPf, ETR, ETc, ky=0.85, metodo='etapa_unica')
pa1.attrs['PA'], pa1.attrs['EC']                           # 4872 kg/ha, 0,881
```

Cada fórmula do modelo também tem sua função (`cTn`, `cTc`, com `rota='C3'`
ou `'C4'`, `PPBp`, `CIAF`, `CR`), e todas as constantes e valores tabelados
são parâmetros com o valor usual como padrão (ex.: `U=13`, `b=0.52`).

As funções também podem ser acessadas por submódulo (`amp.radiacao`,
`amp.umidade`, `amp.evapotranspiracao`, etc.).

## Módulos

| Módulo | Conteúdo |
| --- | --- |
| `radiacao` | NDA, declinação solar, ângulo horário, ângulo zenital, azimute solar, fotoperíodo, irradiância extraterrestre, radiação global (Angström-Prescott, Hargreaves-Samani) e razão de insolação n/N (`relacao_n_N`, Angström-Prescott invertida) |
| `temperatura` | Temperatura média diária (extremos e estação automática) |
| `umidade` | Pressão de saturação/parcial de vapor, déficit de saturação, pressão atmosférica, umidade absoluta/relativa/de saturação, ponto de orvalho, constante psicrométrica |
| `energia` | Saldo de radiação, balanço de ondas curtas e longas |
| `evapotranspiracao` | ETP por Thornthwaite, Camargo-Maluf, Hargreaves-Samani, Priestley-Taylor; ETo por Penman-Monteith FAO-56; conversão do vento para 2 m |
| `grau_dias` | Data de maturação fisiológica / data de semeadura por acúmulo de graus-dia |
| `balanco_hidrico` | Balanço hídrico climatológico e de cultura (Thornthwaite & Mather) |
| `classificacao_climatica` | Köppen-Geiger, Thornthwaite, Camargo e Holdridge, para um local ou em grade (rasters) |
| `produtividade` | Zona Agroecológica da FAO (Doorenbos & Kassam): correções de temperatura cTn e cTc (C3/C4), PPBp, CIAF, CR, produtividade potencial (PPf, em escala diária, decendial ou mensal) e atingível (PA, por fase ou em etapa única) com eficiência climática |

- [`docs/FORMULAS.md`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/docs/FORMULAS.md) ([English](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/docs/FORMULAS.en.md)) — documentação matemática de cada função: a fórmula original, variáveis e unidades.
- [`examples/tutorial_colab.ipynb`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/tutorial_colab.ipynb) — tutorial guiado, pronto para o Colab, com um exemplo executável para cada função.
- [`examples/formulas_agrometeorologia.ipynb`](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/examples/formulas_agrometeorologia.ipynb) — notebook original de desenvolvimento, com a dedução teórica de cada fórmula.

## Desenvolvimento

```bash
git clone https://github.com/fcoliveira-utfpr/agrometeorologiapy.git
cd agrometeorologiapy
pip install -e ".[dev]"
pytest
```

## Licença

BSD-3-Clause. Veja [LICENSE](https://github.com/fcoliveira-utfpr/agrometeorologiapy/blob/main/LICENSE).
