# Changelog

## [0.2.0] - não publicada

### Adicionado
- `balanco_hidrico_climatologico_grade`: BHC de Thornthwaite & Mather
  vetorizado com numpy, para muitos locais de uma vez (ex.: todos os pixels
  de um raster com a normal mensal), com CAD escalar ou por local.

### Corrigido
- `balanco_hidrico_climatologico`: o armazenamento inicial passa a ser o de
  equilíbrio do ciclo anual (`ciclico=True`, padrão), e não mais o solo
  cheio em janeiro. Faz diferença onde janeiro é seco (ex.: semiárido,
  hemisfério norte). `ciclico=False` mantém o início com solo cheio.
- `balanco_hidrico_climatologico`: o ALT de janeiro era zerado
  (`diff().fillna(0)`), ignorando a variação a partir do armazenamento
  inicial; isso distorcia ETR e DEF de janeiro quando ele é seco.
- `thornthwaite_mensal`: para T >= 26,5 °C usa a equação da tabela do método
  original (-415,85 + 32,24 T - 0,43 T²) no lugar da exponencial.

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
