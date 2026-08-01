# Changelog

## [0.1.1] - 2026-08-01

Sem mudanças de código. Só documentação e metadados:

- README: links relativos (`LICENSE`, `docs/FORMULAS.md`, notebooks em
  `examples/`) trocados por URLs absolutas do GitHub — no PyPI, a página do
  projeto é renderizada isolada do repositório, então links relativos
  quebravam.
- Adiciona `docs/FORMULAS.md`: documentação matemática de cada função (a
  fórmula original, variáveis e unidades).
- Adiciona `examples/tutorial_colab.ipynb`: tutorial interativo pronto para
  o Google Colab, com um exemplo por função.
- `pyproject.toml`: adiciona os links "Documentation", "Tutorial" e
  "Changelog" em `[project.urls]`, exibidos na barra lateral do PyPI.

## [0.1.0] - 2026-08-01

Primeira versão publicável do pacote, extraída do notebook
`formulas_agrometeorologia.ipynb`.

- Radiação solar, temperatura, umidade do ar, balanço de energia.
- Evapotranspiração: Thornthwaite, Camargo-Maluf, Hargreaves-Samani,
  Priestley-Taylor, Penman-Monteith FAO-56.
- Grau-dias (data de maturação fisiológica / data de semeadura).
- Balanço hídrico climatológico e de cultura.

### Corrigido em relação ao notebook original
- `bol_saldo`: o notebook original aplicava o termo de Stefan-Boltzmann
  diretamente sobre Tmax/Tmin em °C (sem converter para Kelvin), o que
  subestimava drasticamente o balanço de ondas longas (BOL ≈ -6e-5 em vez
  de ≈ -4.3 MJ/m² dia num exemplo típico). A função agora converte
  Tmax/Tmin de °C para Kelvin internamente, mantendo a assinatura em °C
  (consistente com o resto do pacote) mas com o cálculo fisicamente
  correto — isso também corrige `saldo_radiacao`, `etp_priestley_taylor`
  e `eto_penman_monteith_fao56` quando alimentados a partir de `bol_saldo`.
