# agrometeorologiapy

[![PyPI](https://img.shields.io/pypi/v/agrometeorologiapy.svg)](https://pypi.org/project/agrometeorologiapy/)
[![CI](https://github.com/fcoliveira/agrometeorologiapy/actions/workflows/ci.yml/badge.svg)](https://github.com/fcoliveira/agrometeorologiapy/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

Fórmulas de agrometeorologia em Python: radiação solar, temperatura, umidade
do ar, balanço de energia, evapotranspiração (Thornthwaite, Camargo-Maluf,
Hargreaves-Samani, Priestley-Taylor, Penman-Monteith FAO-56), grau-dias e
balanço hídrico (climatológico e de cultura).

## Instalação

```bash
pip install agrometeorologiapy
```

## Uso rápido

```python
import agrometeorologiapy as amp

nda = amp.nda(2, 3)                       # Número do Dia do Ano
delta = amp.declinacao_solar(nda)         # Declinação solar (graus)
Hn = amp.angulo_horario_nascer(lat=-22, declinacao=delta)
N = amp.fotoperiodo(Hn)                   # Fotoperíodo (horas)
```

As funções também podem ser acessadas por submódulo (`amp.radiacao`,
`amp.umidade`, `amp.evapotranspiracao`, etc.).

## Módulos

| Módulo | Conteúdo |
| --- | --- |
| `radiacao` | NDA, declinação solar, ângulo horário, ângulo zenital, azimute solar, fotoperíodo, irradiância extraterrestre, radiação global (Angström-Prescott, Hargreaves-Samani) |
| `temperatura` | Temperatura média diária (extremos e estação automática) |
| `umidade` | Pressão de saturação/parcial de vapor, déficit de saturação, pressão atmosférica, umidade absoluta/relativa/de saturação, ponto de orvalho, constante psicrométrica |
| `energia` | Saldo de radiação, balanço de ondas curtas e longas |
| `evapotranspiracao` | ETP por Thornthwaite, Camargo-Maluf, Hargreaves-Samani, Priestley-Taylor; ETo por Penman-Monteith FAO-56 |
| `grau_dias` | Data de maturação fisiológica / data de semeadura por acúmulo de graus-dia |
| `balanco_hidrico` | Balanço hídrico climatológico e de cultura (Thornthwaite & Mather) |

Um notebook com exemplos numéricos resolvidos para cada método está em
[`examples/formulas_agrometeorologia.ipynb`](examples/formulas_agrometeorologia.ipynb).

## Desenvolvimento

```bash
git clone https://github.com/fcoliveira/agrometeorologiapy.git
cd agrometeorologiapy
pip install -e ".[dev]"
pytest
```

## Licença

BSD-3-Clause. Veja [LICENSE](LICENSE).
