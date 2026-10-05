# Changelog

## [0.3.0] - 2026-10-04

### Adicionado
- `relacao_n_N` (módulo `radiacao`): razão de insolação n/N pela mesma
  equação de `Qg_angstrom`, resolvida para n/N a partir de Qg (para quando
  não há heliógrafo). Limitada entre 0 e 1 (`limitar=True`).
- Módulo `produtividade`: Modelo da Zona Agroecológica da FAO (Doorenbos &
  Kassam, 1979). `cTn` e `cTc` (correções de temperatura C3/C4 de Barbieri & Tuon, 1992),
  `PPBp`, `CIAF`, `CR`, `produtividade_potencial` (PPf acumulada período a
  período, em escala diária, decendial ou mensal, com decêndios parciais na
  semeadura e na colheita) e `produtividade_atingivel` (PA e EC pelo
  produtório por fase ou em etapa única). Constantes e valores tabelados são
  parâmetros com os valores usuais como padrão.

### Documentação
- `docs/FORMULAS.md` e `.en.md`: seção 10 (produtividade) e `relacao_n_N`
  na seção 2. As fórmulas em bloco passam para ```` ```math ````, porque o
  GitHub não renderizava `$$ ... $$` colado a texto e engolia os escapes
  (`\,`, `\;`, `\_`, `*`) dentro das fórmulas.
- Tutorial do Colab: seções 2.14 (`relacao_n_N`) e 10 (produtividade).
  Notebook de fórmulas: Capítulo 12 e Aplicações 12 a 14.
- README: exemplos de uso de radiação, n/N, evapotranspiração, balanço
  hídrico e produtividade.

## [0.2.1] - 2026-09-29

### Corrigido
- `data_semeadura` (intervalo mensal): o dia da semeadura deixa de entrar na
  soma de graus-dia, como já ocorre em `data_maturacao_fisiologica`. A data
  calculada passa a ser um dia antes (Aplicação 9: 01/04 em vez de 02/04), e
  as duas funções ficam inversas uma da outra.

## [0.2.0] - 2026-09-29

### Adicionado
- `balanco_hidrico_climatologico_grade`: BHC de Thornthwaite & Mather
  vetorizado com numpy, para muitos locais de uma vez (ex.: todos os pixels
  de um raster com a normal mensal), com CAD escalar ou por local.
- Módulo `classificacao_climatica`: Köppen-Geiger (Alvares et al., 2013;
  Kottek et al., 2006), Thornthwaite (1948), Camargo (1991) mod. Maluf (2000)
  (Aparecido et al., 2016) e zonas de vida de Holdridge (38 zonas), cada um
  para um local (`classificacao_koppen`, `classificacao_thornthwaite`,
  `classificacao_camargo`, `classificacao_holdridge`) e em grade para
  rasters (sufixo `_grade`).
- `vento_2m`: converte o vento medido a qualquer altura para 2 m (FAO-56,
  eq. 47), para uso em `eto_penman_monteith_fao56`.

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
- `data_maturacao_fisiologica` e `data_semeadura`: a última linha do
  DataFrame passa a ser a data calculada (maturação ou semeadura) com o GDA
  acumulado até ela, e não o início do último período com o GD do período
  inteiro. O ano agora vira quando o ciclo cruza dezembro/janeiro.
- `data_maturacao_fisiologica` e `data_semeadura`: o DataFrame passa a ter
  as colunas `Meses`, `data`, `Tmed`, `GDi`, `GDA_mes` e `GDA_ciclo`. A antiga
  `GD_ciclo` foi renomeada para `GDA_ciclo` (**quebra compatibilidade**).
- `data_semeadura`: no último período, os dias passam a ser contados para
  trás a partir do fim do período (antes contava a partir do início).

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
