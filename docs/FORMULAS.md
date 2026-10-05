# Fórmulas de Agrometeorologia — Referência das Funções

*[🇺🇸 English version](FORMULAS.en.md)*

Documentação matemática de cada função pública do `agrometeorologiapy`: a
fórmula original, o significado de cada variável e a unidade esperada.
Para instalação e exemplos de uso, veja o [README](../README.md); para um
tutorial executável, veja
[`examples/tutorial_colab.ipynb`](../examples/tutorial_colab.ipynb).

**Fonte principal:** Pereira, Angelocci & Sentelhas (2002) — *Agrometeorologia:
fundamentos e aplicações práticas*, ESALQ/USP, complementado por Allen et al.
(1998, FAO-56) nos métodos de evapotranspiração.

> ⚠️ **Convenção de ângulos:** todas as funções trigonométricas do pacote
> (`sind`, `cosd`, `tand`, `acosd`, e por extensão toda função que calcula
> ângulos) trabalham em **graus**, não em radianos.

## Sumário

- [1. Funções auxiliares de trigonometria](#1-funções-auxiliares-de-trigonometria)
- [2. Radiação Solar](#2-radiação-solar)
- [3. Temperatura](#3-temperatura)
- [4. Umidade do Ar](#4-umidade-do-ar)
- [5. Balanço de Energia](#5-balanço-de-energia)
- [6. Evapotranspiração](#6-evapotranspiração)
- [7. Grau-Dias](#7-grau-dias)
- [8. Balanço Hídrico](#8-balanço-hídrico)
- [9. Classificação Climática](#9-classificação-climática)
- [10. Produtividade Potencial e Atingível](#10-produtividade-potencial-e-atingível)

---

## 1. Funções auxiliares de trigonometria

Módulo `agrometeorologiapy._trig`. Atalhos usados internamente por quase
todas as demais funções, para evitar conversões manuais grau↔radiano.

### `sind(x)`

```math
\text{sind}(x) = \sin\left(x \cdot \frac{\pi}{180}\right)
```

**Onde:** `x` — ângulo, em graus.

### `cosd(x)`

```math
\text{cosd}(x) = \cos\left(x \cdot \frac{\pi}{180}\right)
```

**Onde:** `x` — ângulo, em graus.

### `tand(x)`

```math
\text{tand}(x) = \tan\left(x \cdot \frac{\pi}{180}\right)
```

**Onde:** `x` — ângulo, em graus.

### `acosd(x)`

```math
\text{acosd}(x) = \arccos(x) \cdot \frac{180}{\pi}
```

**Onde:** `x` — cosseno de um ângulo, adimensional (entre -1 e 1). Retorna o
ângulo correspondente em graus.

---

## 2. Radiação Solar

Módulo `agrometeorologiapy.radiacao`.

### `nda(dia, mes, ano=2023)`
Número do Dia do Ano (NDA): posição ordinal do dia dentro do ano
(1 de janeiro = 1, 31 de dezembro = 365 ou 366). Calculado via
`datetime.timetuple().tm_yday`, sem fórmula fechada.

**Onde:**
- `dia` — dia do mês, inteiro (1–31)
- `mes` — mês, inteiro (1–12)
- `ano` — ano, inteiro (padrão 2023, não bissexto)
- `NDA` (retorno) — número do dia do ano, inteiro (1–366)

### `declinacao_solar(NDA)`
Aproximação senoidal da declinação solar δ, decorrente da inclinação do
eixo terrestre (23,45°):

```math
\delta = 23{,}45 \cdot \sin\left(\frac{360 (NDA - 80)}{365}\right) \quad [\text{graus}]
```

**Onde:**
- $\delta$ (retorno) — declinação solar do dia, em graus
- $NDA$ — número do dia do ano (1–366)

### `angulo_horario(hora, minuto=0)`

```math
h = (\text{hora} + \tfrac{\text{minuto}}{60} - 12) \times 15 \quad [\text{graus}]
```

(15°/hora, nulo ao meio-dia solar, negativo pela manhã, positivo à tarde)

**Onde:**
- $h$ (retorno) — ângulo horário, em graus
- `hora` — hora do dia, inteiro (0–23)
- `minuto` — minuto da hora, inteiro (0–59)

### `angulo_zenital(lat, declinacao)`
Ângulo zenital do Sol ao meio-dia local (h = 0°):

```math
\cos Z = \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)
```

```math
Z = \arccos(\cos Z)
```

(o valor de $\cos Z$ é limitado ao intervalo $[-1, 1]$ antes do arco-cosseno,
para evitar erro de domínio por arredondamento)

**Onde:**
- $Z$ (retorno) — ângulo zenital, em graus
- $\varphi$ (`lat`) — latitude do local, em graus (negativa no hemisfério sul)
- $\delta$ (`declinacao`) — declinação solar do dia, em graus

### `azimute_solar(lat, declinacao, Z)`
Direção horizontal do Sol em relação à linha Norte-Sul:

```math
\cos \alpha = \frac{\sin(\varphi)\cos(Z) - \sin(\delta)}{\cos(\varphi)\sin(Z)}
```

```math
\alpha = \arccos(\cos \alpha)
```

**Onde:**
- $\alpha$ (retorno) — azimute solar, em graus
- $\varphi$ (`lat`) — latitude do local, em graus
- $\delta$ (`declinacao`) — declinação solar do dia, em graus
- $Z$ — ângulo zenital, em graus (eq. anterior)

### `comprimento_sombra(d, Z)`
Comprimento da sombra projetada por um objeto de altura `d`:

```math
S = d \cdot \tan(Z)
```

**Onde:**
- $S$ (retorno) — comprimento da sombra, na mesma unidade de `d`
- `d` — altura do objeto, em metros
- $Z$ — ângulo zenital no instante considerado, em graus

### `fotoperiodo(Hn)`
Duração do dia (dobro do ângulo horário do nascer, convertido para horas):

```math
N = \frac{2 \, H_n}{15} \quad [\text{horas}]
```

**Onde:**
- $N$ (retorno) — fotoperíodo (duração do dia), em horas
- $H_n$ (`Hn`) — ângulo horário no nascer do Sol, em graus

### `angulo_horario_nascer(lat, declinacao)`
Obtida impondo-se $Z = 90°$ (cos Z = 0) na equação do ângulo zenital:

```math
H_n = \arccos(-\tan(\varphi)\tan(\delta))
```

**Onde:**
- $H_n$ (retorno) — ângulo horário no nascer do Sol, em graus
- $\varphi$ (`lat`) — latitude do local, em graus
- $\delta$ (`declinacao`) — declinação solar do dia, em graus

### `fator_correcao_distancia(NDA)`
Correção $(d/D)^2$ da excentricidade da órbita terrestre:

```math
\left(\frac{d}{D}\right)^2 = 1 + 0{,}033 \cdot \cos\left(\frac{360 \cdot NDA}{365}\right)
```

**Onde:**
- $(d/D)^2$ (retorno) — fator de correção da distância Terra-Sol, adimensional
- $NDA$ — número do dia do ano (1–366)

### `irradiancia_extraterrestre(lat, declinacao, Hn, dD2)`
Irradiância solar global no topo da atmosfera (Qo):

```math
Q_o = 37{,}6 \cdot \left(\frac{d}{D}\right)^2 \cdot \left[ H_{n,\text{rad}} \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)\sin(H_n) \right] \quad [\text{MJ/m}^2\text{dia}]
```

**Onde:**
- $Q_o$ (retorno) — irradiância solar extraterrestre diária, em MJ/m² dia
- $\varphi$ (`lat`) — latitude do local, em graus
- $\delta$ (`declinacao`) — declinação solar do dia, em graus
- $H_n$ (`Hn`) — ângulo horário no nascer do Sol, em graus; $H_{n,\text{rad}}$
  é o mesmo valor convertido para radianos, usado apenas nesse termo
- $(d/D)^2$ (`dD2`) — fator de correção da distância Terra-Sol, adimensional

### `insolacao(N, Tmax, Tmin, lat, k=0.19)`
Estimativa do número de horas de brilho solar (n) a partir da amplitude
térmica diária:

```math
n = \frac{N}{0{,}52} \cdot \left[ k \sqrt{T_{max} - T_{min}} - 0{,}29 \cos(\varphi) \right] \quad [\text{horas}]
```

**Onde:**
- $n$ (retorno) — número de horas de insolação (brilho solar) estimado, em horas
- $N$ — fotoperíodo do dia, em horas
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — temperaturas máxima e mínima diárias, em °C
- $\varphi$ (`lat`) — latitude do local, em graus
- $k$ — coeficiente empírico, adimensional (0,16 para regiões interioranas,
  0,19 — padrão da função — para regiões costeiras)

### `Qg_angstrom(insolacao, N, Qo, lat, b=0.52)`
Equação de Angström-Prescott (variante Glover-McCulloch):

```math
a = 0{,}29 \cos(\varphi)
```

```math
Q_g = Q_o \cdot \left( a + b \cdot \frac{n}{N} \right) \quad [\text{MJ/m}^2]
```

**Onde:**
- $Q_g$ (retorno) — irradiância solar global, em MJ/m²
- $n$ (`insolacao`) — horas de brilho solar medidas no dia, em horas
- $N$ — fotoperíodo do dia, em horas
- $Q_o$ (`Qo`) — irradiância solar extraterrestre, em MJ/m²
- $\varphi$ (`lat`) — latitude do local, em graus
- $a$ — coeficiente de regressão derivado da latitude, adimensional
- $b$ — coeficiente empírico de regressão, adimensional (padrão 0,52)

### `relacao_n_N(Qg, Qo, lat, b=0.52, a=None, limitar=True)`
A mesma equação de Angström-Prescott de `Qg_angstrom`, resolvida para a razão
de insolação — para quando há $Q_g$ (piranômetro ou bases em grade, como o
BR-DWGD), mas não há heliógrafo:

```math
\frac{n}{N} = \frac{Q_g / Q_o - a}{b}, \qquad a = 0{,}29 \cos(\varphi)
```

**Onde:**
- $n/N$ (retorno) — razão de insolação, adimensional; com `limitar=True`
  (padrão), limitada entre 0 e 1
- $Q_g$ (`Qg`) — irradiância solar global, em MJ/m² dia
- $Q_o$ (`Qo`) — irradiância solar extraterrestre, em MJ/m² dia
- $\varphi$ (`lat`) — latitude do local, em graus
- $b$ — coeficiente empírico de regressão, adimensional (padrão 0,52)
- $a$ — coeficiente linear, adimensional (padrão `None`: $0{,}29 \cos\varphi$)

### `Qg_hargreaves(Tmax, Tmin, Qo, k=0.16)`
Equação de Hargreaves-Samani, sem depender de dados de insolação:

```math
Q_g = k \sqrt{T_{max} - T_{min}} \cdot Q_o \quad [\text{MJ/m}^2]
```

**Onde:**
- $Q_g$ (retorno) — irradiância solar global, em MJ/m²
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — temperaturas máxima e mínima diárias, em °C
- $Q_o$ (`Qo`) — irradiância solar extraterrestre, em MJ/m²
- $k$ — coeficiente empírico, adimensional (0,16 — padrão da função — para
  regiões interioranas, 0,19 para regiões costeiras)

---

## 3. Temperatura

Módulo `agrometeorologiapy.temperatura`.

### `temp_media_extremos(Tmax, Tmin)`

```math
T_{med} = \frac{T_{max} + T_{min}}{2}
```

**Onde:**
- $T_{med}$ (retorno) — temperatura média diária estimada, em °C
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — temperaturas máxima e mínima do dia, em °C

### `temp_media_estacao_automatica(temperaturas)`
Média aritmética simples de $n$ observações no período:

```math
T_{med} = \frac{1}{n}\sum_{i=1}^{n} T_i
```

**Onde:**
- $T_{med}$ (retorno) — temperatura média do período, em °C
- $T_i$ (`temperaturas`) — cada observação de temperatura, em °C
- $n$ — número de observações na lista

---

## 4. Umidade do Ar

Módulo `agrometeorologiapy.umidade`.

### `es_tetens(T_ar)`
Equação de Tetens — pressão de saturação de vapor d'água:

```math
e_s = 0{,}6108 \cdot 10^{\frac{7{,}5 \, T}{237{,}3 + T}} \quad [\text{kPa}]
```

**Onde:**
- $e_s$ (retorno) — pressão de saturação de vapor, em kPa
- $T$ (`T_ar`) — temperatura do ar, em °C

### `ea_umidade(es, UR)`
Pressão parcial (atual) de vapor:

```math
e_a = e_s \cdot \frac{UR}{100}
```

**Onde:**
- $e_a$ (retorno) — pressão parcial (atual) de vapor d'água, em kPa
- $e_s$ (`es`) — pressão de saturação de vapor, em kPa
- $UR$ — umidade relativa do ar, em % (ex.: 65)

### `deficit_saturacao(es, ea)`

```math
\Delta e = e_s - e_a
```

**Onde:**
- $\Delta e$ (retorno) — déficit de saturação de vapor, em kPa
- $e_s$ (`es`) — pressão de saturação de vapor, em kPa
- $e_a$ (`ea`) — pressão parcial (atual) de vapor d'água, em kPa

### `patm_altitude(A)`
Pressão atmosférica local em função da altitude:

```math
P_{atm} = 101{,}3 \cdot \left( \frac{293 - 0{,}0065 A}{293} \right)^{5{,}26} \quad [\text{kPa}]
```

**Onde:**
- $P_{atm}$ (retorno) — pressão atmosférica local, em kPa
- $A$ — altitude do local, em metros

### `umidade_absoluta(ea, T_ar_C)`

```math
UA = \frac{2168 \cdot e_a}{T + 273{,}15} \quad [\text{g de H}_2\text{O / m}^3]
```

**Onde:**
- $UA$ (retorno) — umidade absoluta do ar, em g de H₂O por m³ de ar
- $e_a$ (`ea`) — pressão parcial de vapor d'água, em kPa
- $T$ (`T_ar_C`) — temperatura do ar, em °C (convertida para Kelvin
  internamente, $T + 273{,}15$)

### `umidade_saturacao(es, T_ar_C)`

```math
US = \frac{2168 \cdot e_s}{T + 273{,}15} \quad [\text{g de H}_2\text{O / m}^3]
```

**Onde:**
- $US$ (retorno) — umidade de saturação do ar, em g de H₂O por m³ de ar
- $e_s$ (`es`) — pressão de saturação de vapor, em kPa
- $T$ (`T_ar_C`) — temperatura do ar, em °C

### `umidade_relativa(ea, es)`

```math
UR = 100 \cdot \frac{e_a}{e_s} \quad [\%]
```

**Onde:**
- $UR$ (retorno) — umidade relativa do ar, em %
- $e_a$ (`ea`) — pressão parcial (atual) de vapor d'água, em kPa
- $e_s$ (`es`) — pressão de saturação de vapor, em kPa

### `ponto_orvalho(ea)`
Inversão algébrica da equação de Tetens:

```math
x = \log_{10}\left(\frac{e_a}{0{,}6108}\right)
```

```math
T_o = \frac{237{,}3 \, x}{7{,}5 - x} \quad [°C]
```

**Onde:**
- $T_o$ (retorno) — temperatura do ponto de orvalho, em °C
- $e_a$ (`ea`) — pressão parcial (atual) de vapor d'água, em kPa
- $x$ — variável auxiliar (logaritmo da razão $e_a/0{,}6108$), adimensional

### `constante_psicrometrica(Patm)`

```math
\gamma = 0{,}665 \times 10^{-3} \cdot P_{atm} \quad [\text{kPa}/°C]
```

**Onde:**
- $\gamma$ (retorno) — constante psicrométrica, em kPa/°C
- $P_{atm}$ (`Patm`) — pressão atmosférica local, em kPa

---

## 5. Balanço de Energia

Módulo `agrometeorologiapy.energia`.

### `boc_saldo(Qg, r=0.25)`
Balanço de ondas curtas — radiação global incidente menos a fração refletida
(albedo `r`; 0,25 é a média para gramado):

```math
BOC = Q_g (1 - r)
```

**Onde:**
- $BOC$ (retorno) — balanço de ondas curtas, na mesma unidade de $Q_g$
- $Q_g$ (`Qg`) — irradiância solar global incidente na superfície, em
  MJ/m² dia (ou W/m²)
- $r$ — coeficiente de reflexão da superfície (albedo), adimensional (0–1;
  padrão 0,25, típico de gramado)

### `bol_saldo(Tmax, Tmin, ea, Qg, Qg_cs)`
Balanço de ondas longas, equação de Stefan-Boltzmann corrigida pela FAO-56.
`Tmax`/`Tmin` são recebidos em °C (consistente com o resto do pacote) e
convertidos internamente para Kelvin, porque o termo de Stefan-Boltzmann
exige temperatura absoluta:

```math
T_{max,K} = T_{max} + 273{,}15, \qquad T_{min,K} = T_{min} + 273{,}15
```

```math
BOL = -\left[ 4{,}903 \times 10^{-9} \cdot \frac{T_{max,K}^4 + T_{min,K}^4}{2} \right] \cdot \left[ 0{,}34 - 0{,}14\sqrt{e_a} \right] \cdot \left[ 1{,}35 \frac{Q_g}{Q_{g,cs}} - 0{,}35 \right]
```

**Onde:**
- $BOL$ (retorno) — balanço de radiação de onda longa, em MJ/m² (negativo:
  perda líquida de energia por emissão terrestre)
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — temperaturas máxima e mínima
  diárias, em °C; $T_{max,K}$, $T_{min,K}$ são as mesmas em Kelvin
- $e_a$ (`ea`) — pressão parcial (atual) de vapor d'água, em kPa
- $Q_g$ (`Qg`) — radiação solar global medida/estimada no dia, em MJ/m²
- $Q_{g,cs}$ (`Qg_cs`) — radiação solar de céu claro (*clear-sky*), em MJ/m²
- $4{,}903 \times 10^{-9}$ — constante de Stefan-Boltzmann, em MJ K⁻⁴ m⁻² dia⁻¹

> 🐛 **Nota:** o notebook original de origem aplicava esse termo diretamente
> sobre Tmax/Tmin em °C (sem converter para Kelvin), o que subestimava o BOL
> em várias ordens de grandeza. Foi corrigido nesta biblioteca — veja o
> [CHANGELOG](../CHANGELOG.md).

### `saldo_radiacao(BOC, BOL)`

```math
R_n = BOC + BOL \quad [\text{MJ/m}^2\text{dia}]
```

**Onde:**
- $R_n$ (retorno) — saldo de radiação, em MJ/m² dia (ou W/m²)
- $BOC$ — balanço de ondas curtas, mesma unidade de $R_n$
- $BOL$ — balanço de ondas longas, mesma unidade de $R_n$

---

## 6. Evapotranspiração

Módulo `agrometeorologiapy.evapotranspiracao`.

### `thornthwaite_mensal(df, col_T='T_media_C', lat=None)`
Método de Thornthwaite (1948), mensal. Para cada um dos 12 meses:

```math
i = \left(\frac{T}{5}\right)^{1{,}514} \quad \text{se } T > 0, \text{ senão } 0
```

```math
I = \sum_{m=1}^{12} i_m \qquad \text{(índice de calor anual)}
```

```math
a = 6{,}75 \times 10^{-7} I^3 - 7{,}71 \times 10^{-5} I^2 + 1{,}792 \times 10^{-2} I + 0{,}49239
```

```math
ETP_{nc} = 16 \left( \frac{10 T}{I} \right)^{a} \quad \text{(não corrigida, mm/mês de 30 dias e 12h de sol)}
```

Para $T \ge 26{,}5$ °C o método original usa uma tabela no lugar da equação acima:

```math
ETP_{nc} = -415{,}85 + 32{,}24\,T - 0{,}43\,T^2
```

Correção pelo fotoperíodo real do mês (via declinação solar do dia juliano
médio de cada mês) e pelo número real de dias do mês:

```math
K = \frac{N}{12} \cdot \frac{\text{dias}_{\text{mês}}}{30}
```

```math
ETP = ETP_{nc} \cdot K \quad [\text{mm/mês}]
```

**Onde:**
- $ETP$ (retorno, coluna `ETP_mm_mes`) — evapotranspiração potencial
  corrigida, em mm/mês
- $T$ (coluna `col_T` de `df`) — temperatura média de cada mês, em °C
- $i$ — índice de calor mensal, adimensional; $I$ — soma dos 12 índices
  mensais (índice de calor anual), adimensional
- $a$ — expoente empírico, função de $I$, adimensional
- $ETP_{nc}$ — evapotranspiração potencial não corrigida, em mm/mês
- $N$ — fotoperíodo do mês (calculado internamente a partir de `lat` e da
  declinação solar do dia juliano médio do mês), em horas
- `dias_mes` — número real de dias do mês (28–31)

### `camargo_maluf_mensal(df, col_T='T_media_C', lat=None)`
Método de Camargo (1971), coeficiente F modificado por Maluf (Camargo et
al., 1999). Função autocontida: recalcula internamente a radiação
extraterrestre mensal ($Q_o$, a partir da declinação solar do dia juliano
médio de cada mês) e converte para equivalente de evaporação:

```math
Q_{o,mm} = 0{,}408 \cdot Q_{o}
```

O coeficiente F depende da temperatura média **anual** ($\bar{T}_{anual}$):

- $F = 0{,}0100$ se $\bar{T}_{anual} < 23°C$
- $F = 0{,}0105$ se $23°C \le \bar{T}_{anual} < 24°C$
- $F = 0{,}0110$ se $\bar{T}_{anual} \ge 24°C$

```math
ETP = F \cdot Q_{o,mm} \cdot T \cdot \text{dias}_{\text{mês}} \quad [\text{mm/mês}]
```

**Onde:**
- $ETP$ (retorno, coluna `ETP_mm_mes`) — evapotranspiração potencial mensal,
  em mm/mês
- $Q_o$ (coluna `Qo_MJ_m2dia`) — radiação extraterrestre mensal, em
  MJ/m² dia; $Q_{o,mm}$ (coluna `Qo_mm_dia`) — a mesma, em equivalente de
  evaporação, mm/dia
- $F$ — coeficiente empírico, adimensional, escolhido pela temperatura
  média anual
- $\bar{T}_{anual}$ — temperatura média anual (média das 12 temperaturas
  mensais de entrada), em °C
- $T$ (coluna `col_T` de `df`) — temperatura média de cada mês, em °C
- `dias_mes` — número real de dias do mês (28–31)

### `etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)`

```math
Q_{o,mm} = 0{,}408 \cdot Q_o
```

```math
ETP = 0{,}0023 \cdot Q_{o,mm} \cdot \sqrt{T_{max} - T_{min}} \cdot (T_{med} + 17{,}8) \quad [\text{mm/dia}]
```

**Onde:**
- $ETP$ (retorno) — evapotranspiração potencial, em mm/dia
- $Q_o$ (`Qo`) — irradiância solar extraterrestre, em MJ/m² dia;
  $Q_{o,mm}$ — a mesma, em equivalente de evaporação, mm/dia
- $T_{max}$, $T_{min}$, $T_{med}$ (`Tmax`, `Tmin`, `Tmed`) — temperaturas
  máxima, mínima e média do ar, em °C

### `declive_pressao_vapor(T_ar)`
Derivada da equação de Tetens em relação à temperatura:

```math
\Delta = \frac{4098 \cdot e_s(T)}{(T + 237{,}3)^2} \quad [\text{kPa}/°C]
```

**Onde:**
- $\Delta$ (retorno) — declive da curva de pressão de saturação de vapor,
  em kPa/°C
- $T$ (`T_ar`) — temperatura do ar, em °C (em geral, a temperatura média diária)
- $e_s(T)$ — pressão de saturação de vapor na temperatura $T$ (`es_tetens`), em kPa

### `etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26)`

```math
ETP = \alpha \cdot \frac{\Delta}{\Delta + \gamma} \cdot \frac{R_n - G}{\lambda} \quad [\text{mm/dia}]
```

**Onde:**
- $ETP$ (retorno) — evapotranspiração potencial, em mm/dia
- $R_n$ (`Rn`) — saldo de radiação, em MJ/m² dia
- $G$ — fluxo de calor no solo, em MJ/m² dia (geralmente $G = 0$ em escala diária)
- $\Delta$ (`Delta`) — declive da curva de pressão de saturação de vapor, em kPa/°C
- $\gamma$ (`gamma`) — constante psicrométrica, em kPa/°C
- $\lambda = 2{,}45$ — calor latente de vaporização, em MJ/kg (constante fixa no código)
- $\alpha$ (`alfa`) — coeficiente de Priestley-Taylor, adimensional (padrão 1,26)

### `vento_2m(uz, z)`
Converte o vento medido a uma altura qualquer para 2 m, pelo perfil
logarítmico do vento (FAO-56, eq. 47):

```math
u_2 = u_z \cdot \frac{4{,}87}{\ln(67{,}8 \, z - 5{,}42)} \quad [\text{m/s}]
```

**Onde:**
- $u_2$ (retorno) — velocidade do vento a 2 m de altura, em m/s
- $u_z$ (`uz`) — velocidade do vento medida à altura $z$, em m/s
- $z$ (`z`) — altura de medição acima do solo, em m

### `eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)`
Equação de Penman-Monteith padronizada pelo boletim FAO-56 (Allen et al.,
1998), referenciada a uma cultura hipotética (grama, 0,12 m, albedo 0,23):

```math
ETo = \frac{0{,}408 \, \Delta (R_n - G) + \gamma \cdot \frac{900}{T_{med}+273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma (1 + 0{,}34 \, u_2)} \quad [\text{mm/dia}]
```

**Onde:**
- $ETo$ (retorno) — evapotranspiração de referência, em mm/dia
- $R_n$ (`Rn`) — saldo de radiação, em MJ/m² dia
- $G$ — fluxo de calor no solo, em MJ/m² dia ($G = 0$ para escala diária)
- $T_{med}$ (`Tmed`) — temperatura média diária do ar, em °C
- $u_2$ (`u2`) — velocidade do vento a 2 m de altura, em m/s
- $e_s$ (`es`) — pressão de saturação de vapor, em kPa
- $e_a$ (`ea`) — pressão parcial (atual) de vapor d'água, em kPa
- $\Delta$ (`Delta`) — declive da curva de pressão de saturação de vapor, em kPa/°C
- $\gamma$ (`gamma`) — constante psicrométrica, em kPa/°C

---

## 7. Grau-Dias

Módulo `agrometeorologiapy.grau_dias`.

### Regra do grau-dia diário (GDi)
Usada tanto por `data_maturacao_fisiologica` quanto por `data_semeadura`:

- $GD_i = T_{med} - T_b$, se $T_b < T_{min}$
- $GD_i = \dfrac{(T_{max} - T_b)^2}{2(T_{max} - T_{min})}$, se $T_b \ge T_{min}$

**Onde:**
- $GD_i$ — grau-dia do período (dia, decêndio ou mês), em °C·dia
- $T_{med}$, $T_{max}$, $T_{min}$ — temperaturas média, máxima e mínima do
  período, em °C
- $T_b$ — temperatura base da cultura (abaixo da qual não há
  desenvolvimento), em °C

### `data_maturacao_fisiologica(df, Tb, CT, dia_semeadura, mes_semeadura, intervalo='d', ano=2023)`
A partir da data de semeadura, acumula $GD_i \times n_{período}$
período a período (diário, decendial ou mensal) até que o acumulado
atinja a constante térmica do ciclo:

```math
\sum GD_i \cdot n_{período} \ge CT
```

Retorna a data em que isso ocorre — a maturação fisiológica. No último
período, conta-se apenas os dias necessários para completar `CT`
($\lceil (CT - GDA_{anterior}) / GD_i \rceil$), de modo que a última linha
do DataFrame é a própria data de maturação. O ano avança ao cruzar
dezembro → janeiro.
O resultado tem uma linha por período, com as colunas `Meses` (mês
abreviado), `data`, `Tmed`, `GDi` (grau-dia diário), `GDA_mes` ($GD_i \times$
dias contados no período) e `GDA_ciclo` (acumulado no ciclo).

**Onde:**
- `df` — série climática (colunas `dia`, `mes`, `Tmed`, `Tmax`, `Tmin`), em ordem cronológica
- `Tb` — temperatura base da cultura, em °C
- `CT` — constante térmica total do ciclo, em °C·dia
- `dia_semeadura`, `mes_semeadura` — data de semeadura, inteiros
- `intervalo` — `'d'` (diário), `'dec'` (decendial, 10 dias) ou `'M'` (mensal)
- $n_{período}$ — número de dias do período (1 no diário, 10 no decendial,
  dias do mês no mensal)
- `ano` — ano de referência, inteiro (padrão 2023, não bissexto)

### `data_semeadura(df, Tb, CT, dia_maturacao, mes_maturacao, intervalo='d', ano=2023)`
O mesmo acúmulo de graus-dia, mas percorrendo o calendário **de trás para
frente** a partir de uma data de maturação conhecida (ex.: colheita-alvo),
até acumular `CT` — retornando a data de semeadura necessária.

**Onde:** mesmos parâmetros de `data_maturacao_fisiologica`, trocando
`dia_semeadura`/`mes_semeadura` (data de partida, conhecida) por
`dia_maturacao`/`mes_maturacao` (data de referência a partir da qual se
retrocede).

---

## 8. Balanço Hídrico

Módulo `agrometeorologiapy.balanco_hidrico`. Ambas as funções implementam o
método de contabilidade sequencial de Thornthwaite & Mather (1955).

### `balanco_hidrico_climatologico(df, CAD=100.0, ciclico=True)`
Para cada período (mês) $i$:

```math
P - ETP
```

- Se $P - ETP < 0$ (déficit): acumula o negativo e recalcula o
  armazenamento por via exponencial —

  $`\displaystyle \text{NEG.ACUM}_i = \text{NEG.ACUM}_{i-1} + (P - ETP)`$

  $`\displaystyle ARM_i = CAD \cdot e^{\,\text{NEG.ACUM}_i / CAD}`$

- Se $P - ETP \ge 0$ (reposição): o solo recebe água até no máximo `CAD` —

  $`\displaystyle ARM_i = \min(ARM_{i-1} + (P - ETP),\; CAD)`$

  e, se $ARM_i < CAD$, o NEG.ACUM é recalculado por inversão:

  $`\displaystyle \text{NEG.ACUM}_i = CAD \cdot \ln(ARM_i / CAD)`$

A partir do armazenamento, derivam-se:

```math
ALT_i = ARM_i - ARM_{i-1}
```

- $ETR_i = P_i + |ALT_i|$, se $P_i - ETP_i < 0$
- $ETR_i = ETP_i$, caso contrário

```math
DEF_i = ETP_i - ETR_i
```

- $EXC_i = (P_i - ETP_i) - ALT_i$, se $P_i - ETP_i > 0$ e $ARM_i = CAD$
- $EXC_i = 0$, caso contrário

**Onde:**
- `df` — colunas `Meses`, `P (mm/mês)` (precipitação) e `ETP (mm/mês)`
  (evapotranspiração potencial)
- `CAD` — capacidade de água disponível no solo, em mm (padrão 100,0)
- $P_i$, $ETP_i$ — precipitação e evapotranspiração potencial do período
  $i$, em mm/mês
- $ARM_i$ (coluna `ARM (mm/mês)`) — armazenamento de água no solo ao fim do
  período $i$, em mm; $i-1$ é o período anterior. Com `ciclico=True`
  (padrão), o ARM anterior a janeiro é o de dezembro no equilíbrio do ciclo
  anual (os 12 meses são repetidos até o ARM de dezembro convergir); com
  `ciclico=False`, parte de $ARM = CAD$ (solo cheio)
- $\text{NEG.ACUM}_i$ (coluna `NEG.ACUM (mm)`) — negativo acumulado de
  $P-ETP$, variável auxiliar em mm, usada no modelo exponencial de secagem
- $ALT_i$ (coluna `ALT (mm/mês)`) — variação do armazenamento entre
  períodos, em mm
- $ETR_i$ (coluna `ETR (mm/mês)`) — evapotranspiração real, em mm/mês
- $DEF_i$ (coluna `DEF (mm/mês)`) — deficiência hídrica, em mm/mês
- $EXC_i$ (coluna `EXC (mm/mês)`) — excedente hídrico, em mm/mês

### `balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0, ciclico=True, tol=0.01, max_iter=100)`
Mesmas equações de `balanco_hidrico_climatologico`, vetorizadas com numpy
para muitos locais de uma vez (ex.: todos os pixels de um raster com a
normal mensal).

**Onde:**
- `P`, `ETP` — arrays de formato `(12, ...)`, em mm/mês, com janeiro a
  dezembro no primeiro eixo
- `CAD` — escalar ou array com o formato das dimensões espaciais, em mm
  (> 0; use NaN para locais sem dado)
- `ciclico` — armazenamento inicial de equilíbrio do ciclo anual (padrão)
  ou solo cheio antes de janeiro
- `tol`, `max_iter` — tolerância (mm) e número máximo de ciclos anuais na
  busca do equilíbrio
- retorno — `dict` com arrays `(12, ...)`: `P-ETP`, `ARM`, `NEG.ACUM`,
  `ALT`, `ETR`, `DEF`, `EXC`

### `balanco_hidrico_cultura(df)`
Mesmo método, aplicado a uma cultura específica com `Chuva`, `ETc` e `CAD`
pré-calculados pelo usuário para cada período (mesmas fórmulas de `ARM`,
`ALT`, `ETR`, `DEF` e `EXC` da função anterior, trocando $P \to$ `Chuva` e
$ETP \to$ `ETc`), mais o Índice de Satisfação das Necessidades de Água:

```math
ISNA = \frac{ETR}{ETc}
```

**Onde:**
- `df` — colunas `Chuva` (precipitação, mm/período), `ETc` (evapotranspiração
  da cultura, mm/período) e `CAD` (capacidade de água disponível no
  período, mm), em ordem cronológica
- `ETc` — evapotranspiração da cultura, mm/período; pré-calculada pelo
  usuário como $K_c \times ETo$ (coeficiente de cultura × evapotranspiração
  de referência), já considerando a fase fenológica
- `CAD` — capacidade de água disponível no período, mm; pré-calculada pelo
  usuário como $z \times DTA$ (profundidade radicular × disponibilidade
  total de água), já considerando o avanço da profundidade das raízes
- $ISNA$ (retorno, coluna `ISNA`) — Índice de Satisfação das Necessidades
  de Água, adimensional (0–1: quanto mais próximo de 1, menor o estresse
  hídrico da cultura)

---

## 9. Classificação Climática

Módulo `agrometeorologiapy.classificacao_climatica`. Cada método tem duas
versões: uma para **um local** (listas de 12 meses, jan–dez) e outra
`_grade` para **muitos locais** (ex.: todos os pixels de um raster), com os
12 meses no primeiro eixo, formato `(12, ...)`. Locais com NaN em qualquer
mês recebem código 0 e classe `''`.

**Estações do ano** (Thornthwaite e Camargo), astronômicas e ponderadas pela
fração do mês em cada uma (hemisfério sul; no norte, verão ↔ inverno e
outono ↔ primavera):

| Estação | Meses |
|---|---|
| Verão | ⅓ dez + jan + fev + ⅔ mar |
| Outono | ⅓ mar + abr + mai + ⅔ jun |
| Inverno | ⅓ jun + jul + ago + ⅔ set |
| Primavera | ⅓ set + out + nov + ⅔ dez |

> A Tabela 1 de Aparecido et al. (2016) traz ⅓ jun no outono e ⅔ jun no
> inverno; aqui junho segue o calendário astronômico (o inverno começa em 21/06).

### `classificacao_koppen(T, P, lat)` · `classificacao_koppen_grade(T, P, lat)`
Köppen-Geiger pela metodologia de Alvares et al. (2013), com a sazonalidade
f/s/w dos grupos C e D de Kottek et al. (2006). Variáveis: $T_{ano}$ (média
anual), $T_{frio}$ e $T_{quente}$ (mês mais frio e mais quente), $P_{ano}$ e
$P_{seco}$ (total anual e mês mais seco); verão = out–mar no hemisfério sul.

- **A** (tropical): $T_{frio} \ge 18$ °C — Af se $P_{seco} \ge 60$ mm; Am se
  $P_{seco} \ge 100 - P_{ano}/25$; senão As (seca no verão) ou Aw (no inverno)
- **B** (seco, prevalece sobre os demais): $P_{ano} < 10 P_{lim}$, com
  $P_{lim} = 2T_{ano} + 14$ (ou $2T_{ano}$ se ≥ 70% da chuva cai no inverno;
  $2T_{ano} + 28$ se ≥ 70% cai no verão) — BW se $P_{ano} < 5 P_{lim}$, senão
  BS; h se $T_{ano} \ge 18$ °C, senão k
- **C**: $-3 < T_{frio} < 18$ °C e $T_{quente} > 10$ °C; **D**: $T_{frio} \le -3$ °C;
  **E**: $T_{quente} \le 10$ °C (ET se $T_{quente} > 0$, senão EF)
- C/D — s: seca no verão ($P_{s,seco} < P_{w,seco}$, $P_{w,úmido} > 3P_{s,seco}$ e
  $P_{s,seco} < 40$ mm); w: seca no inverno ($P_{w,seco} < P_{s,seco}$ e
  $P_{s,úmido} > 10P_{w,seco}$); f: nem s nem w. a: $T_{quente} \ge 22$ °C;
  b: ≥ 4 meses acima de 10 °C; c: 1–3 meses; d: $T_{frio} < -38$ °C

**Onde:**
- `T` — temperatura média mensal, em °C
- `P` — precipitação mensal, em mm/mês
- `lat` — latitude, em graus (negativa no hemisfério sul)
- retorno — `id` (1–31) e `classe` (ex.: `'Cfa'`)

### `classificacao_thornthwaite(P, ETP, lat, CAD=100)` · `classificacao_thornthwaite_grade(...)`
Thornthwaite (1948), a partir do balanço hídrico climatológico
(`balanco_hidrico_climatologico_grade`):

```math
I_h = 100\,\frac{EXC}{ETP} \qquad I_a = 100\,\frac{DEF}{ETP} \qquad I_m = I_h - 0{,}6\,I_a
```

```math
ETP_{verão}\,(\%) = 100\,\frac{ETP_{verão}}{ETP_{anual}}
```

Classe = umidade ($I_m$) + subtipo + eficiência térmica ($ETP_{anual}$) +
concentração estival ($ETP_{verão}$, %), ex.: `B1rA'a'`.

| $I_m$ | Umidade | | $ETP_{anual}$ (mm) | Térmica |
|---|---|---|---|---|
| ≥ 100 | A (superúmido) | | ≥ 1140 | A' (megatérmico) |
| 80–100 / 60–80 / 40–60 / 20–40 | B4 / B3 / B2 / B1 (úmido) | | 997–1140 / 855–997 / 712–855 / 570–712 | B'4 / B'3 / B'2 / B'1 (mesotérmico) |
| 0–20 | C2 (subúmido) | | 427–570 / 285–427 | C'2 / C'1 (microtérmico) |
| −20–0 | C1 (subúmido seco) | | 142–285 | D' (tundra) |
| −40– −20 | D (semiárido) | | < 142 | E' (gelo perpétuo) |
| < −40 | E (árido) | | | |

- **Subtipo, climas úmidos** ($I_m \ge 0$), pela deficiência: r ($I_a < 16{,}7$);
  s/w ($16{,}7 \le I_a < 33{,}3$); s2/w2 ($I_a \ge 33{,}3$) — s se a DEF do verão
  for maior que a do inverno, senão w
- **Subtipo, climas secos** ($I_m < 0$), pelo excedente: d ($I_h < 10$); s/w
  ($10 \le I_h < 20$); s2/w2 ($I_h \ge 20$) — w se o EXC do verão for maior que
  o do inverno, senão s
- **Concentração estival** ($ETP_{verão}$, %): a' < 48; b'4 < 51,9; b'3 < 56,3;
  b'2 < 61,6; b'1 < 68; c'2 < 76,3; c'1 < 88; d' ≥ 88

> Tabelas de Aparecido et al. (2016), corrigidas pelo original de Thornthwaite
> (1948) em três pontos em que o artigo traz erro de digitação: limite
> B'3/B'2 = 855 mm (o artigo traz 885); secos s2/w2 com $I_h \ge 20$ (o artigo
> traz 33,3); secos com excedente maior no verão = w (o artigo traz s).

**Onde:**
- `P`, `ETP` — precipitação e evapotranspiração potencial mensais, em mm/mês
- `lat` — latitude, em graus (define verão e inverno)
- `CAD` — capacidade de água disponível no solo, em mm (padrão 100)
- retorno — `classe`, códigos (`umidade`, `subtipo`, `termica`,
  `concentracao`), `Ih`, `Ia`, `Im`, `ETP_anual`, `DEF_anual`, `EXC_anual`,
  `ETP_verao_pct` (e `descricao`, na versão de um local)

### `classificacao_camargo(T, P, ETP, lat, CAD=100, tabela_termica='coerente')` · `classificacao_camargo_grade(...)`
Camargo (1991) modificada por Maluf (2000), Tabelas 6–8 de Aparecido et al.
(2016). Classe = térmica + `-` + hídrica + letra da estação seca, ex.: `ST-UMi`.

| $T_{ano}$ (°C) | $T_{frio}$ (°C) | Térmica |
|---|---|---|
| ≤ 3 | | GL (glacial) |
| 3–7 | | FR (frio) |
| 7–12 | | CO (frio moderado) |
| 12–18 | | TE (temperado) |
| 18–22 | ≤ 13 | STE (subtemperado) |
| 18–22 | 13–20 | ST (subtropical) |
| 18–22 | > 20 | TR (tropical) |
| 22–25 | | TR (tropical) |
| > 25 | | EQ (equatorial) |

| DEF anual (mm) | EXC anual (mm) | Hídrica |
|---|---|---|
| > 800 | 0 | DE (desértico) |
| 150–800 | 0 | AR (árido) |
| > 150 | 0–200 | SE (seco) |
| > 150 | > 200 | MO (monçônico) |
| 0–150 | 0–200 | SB (subúmido) |
| 0–150 (> 0) | > 200 | UM (úmido) |
| 0 | 200–1000 | PU (superúmido) |
| 0 | > 1000 | SU (extremamente úmido) |

Nas classes SE, MO, SB e UM, acrescenta-se a letra da estação com maior
deficiência: v (verão), o (outono), i (inverno) ou p (primavera).

> **Leitura da Tabela 6.** O artigo traz "22 < $T_{ano}$ ≤ 25 **ou**
> $T_{frio}$ > 20 → TR", o que, lido ao pé da letra, impede a classe EQ onde o
> mês mais frio passa de 20 °C (praticamente toda a região equatorial). O
> padrão `tabela_termica='coerente'` usa $T_{frio}$ só para subdividir a faixa
> 18–22 °C (tabela acima); `'literal'` aplica a regra $T_{frio} > 20 \rightarrow TR$ a todos os
> casos, como no script do repositório `climas_brasil`.
>
> **Casos fora da Tabela 7.** O artigo não cobre 0 < DEF ≤ 150 mm com EXC = 0,
> nem DEF = 0 com EXC ≤ 200 mm; aqui eles entram em SB. DEF e EXC abaixo de
> 0,05 mm contam como zero.

**Onde:**
- `T` — temperatura média mensal, em °C
- `P`, `ETP` — precipitação e evapotranspiração potencial mensais, em mm/mês
- `lat` — latitude, em graus
- `CAD` — capacidade de água disponível no solo, em mm (padrão 100)
- `tabela_termica` — `'coerente'` (padrão) ou `'literal'`
- retorno — `classe`, códigos (`termica`, `hidrica`, `estacao_seca`),
  `T_anual`, `T_mes_frio`, `DEF_anual`, `EXC_anual` (e `descricao`, na versão
  de um local)

### `classificacao_holdridge(T, P, lat, ETP=None, limiar_correcao=24)` · `classificacao_holdridge_grade(...)`
Zonas de vida de Holdridge (38 zonas; numeração e nomes de Jungkunst et al.,
2021, a partir de Leemans, 1990). Biotemperatura:

```math
t^*_m = t_m - \frac{3\,|\phi|}{100}\,(t_m - 24)^2 \;\;(\text{se } t_m > 24\ °C), \qquad BT = \overline{\min(\max(t^*_m, 0), 30)}
```

```math
ETP_{anual} = 58{,}93 \cdot BT \qquad R = \frac{ETP_{anual}}{P_{anual}}
```

A zona sai da combinação da faixa de biotemperatura (polar < 1,5; subpolar
< 3; boreal < 6; temperado frio < 12; temperado quente < 18; subtropical
< 24; tropical ≥ 24 °C) com a província de umidade, em faixas de $R$ que
dobram a cada classe (0,125; 0,25; 0,5; 1; 2; 4; 8; 16; 32).

**Onde:**
- `T` — temperatura média mensal, em °C
- `P` — precipitação mensal, em mm/mês
- `lat` — latitude ($\phi$), em graus
- `ETP` — ETP mensal, em mm/mês (opcional; se `None`, usa $58{,}93 \cdot BT$)
- `limiar_correcao` — temperatura mensal acima da qual se corrige pela
  latitude, em °C (padrão 24; `None` corrige todos os meses)
- retorno — `id` (1–38), `classe` (nome da zona em português), `classe_en`
  (nome original em inglês, como na fonte),
  `biotemperatura` (°C), `P_anual`, `ETP_anual` (mm) e `razao_ETP`

**Referências:** Alvares, C. A. et al. (2013) *Meteorol. Z.* 22:711–728 ·
Kottek, M. et al. (2006) *Meteorol. Z.* 15:259–263 · Thornthwaite, C. W. (1948)
*Geogr. Rev.* 38:55–94 · Camargo, A. P. (1991) · Maluf, J. R. T. (2000)
*Rev. Bras. Agrometeorol.* 8:141–150 · Aparecido, L. E. O. et al. (2016)
*Ciênc. Agrotec.* 40:405–417 · Holdridge, L. R. (1967) *Life zone ecology* ·
Jungkunst, H. F. et al. (2021) *J. Plant Nutr. Soil Sci.* 184:5–11

---

## 10. Produtividade Potencial e Atingível

Módulo `agrometeorologiapy.produtividade`. Modelo da Zona Agroecológica da
FAO (Doorenbos & Kassam, 1979), com as correções de temperatura de Barbieri &
Tuon (1992), em Pereira et al. (2002). As constantes e os valores tabelados
têm os valores usuais como padrão, mas todos são parâmetros da função.

### n/N
Quando não há heliógrafo, a razão de insolação vem de `relacao_n_N` (seção 2),
a equação de Angström-Prescott invertida. `produtividade_potencial` a usa
quando o df traz `Qg` em vez de `nN`.

### `cTn(T, rota, T_limiar=16.5)` · `cTc(T, rota, T_limiar=16.5)`
Correções de temperatura da fotossíntese para céu nublado ($cT_n$) e claro
($cT_c$), limitadas a zero.

**C3** (soja, trigo, feijão, girassol):

```math
cT_n = -0{,}0425 + 0{,}035\,T + 0{,}00325\,T^2 - 0{,}0000925\,T^3
```

```math
cT_c = 0{,}583 + 0{,}014\,T + 0{,}0013\,T^2 - 0{,}000037\,T^3
```

**C4** (milho, sorgo, cana):

| Condição | $cT_n$ | $cT_c$ |
|---|---|---|
| $T \ge 16{,}5$ °C | $-1{,}064 + 0{,}173 T - 0{,}0029 T^2$ | $-4{,}16 + 0{,}4325 T - 0{,}00725 T^2$ |
| $T < 16{,}5$ °C | $-4{,}16 + 0{,}4325 T - 0{,}00725 T^2$ | $-9{,}32 + 0{,}865 T - 0{,}0145 T^2$ |

**Onde:**
- `T` — temperatura média do ar no período, em °C
- `rota` — `'C3'` ou `'C4'` (obrigatório)
- `T_limiar` — troca de polinômio das C4, em °C (padrão 16,5)

### `PPBp(Qo, nN, cTc, cTn, a_c=107.2, b_c=8.604, a_n=31.7, b_n=5.234)`
Produtividade potencial bruta padrão da cultura padrão (IAF = 5):

```math
PPB_c = (107{,}2 + 8{,}604\,Q_o)\, cT_c\, \frac{n}{N} \qquad PPB_n = (31{,}7 + 5{,}234\,Q_o)\, cT_n \left(1 - \frac{n}{N}\right)
```

```math
PPBp = PPB_c + PPB_n
```

**Onde:**
- `Qo` — irradiância extraterrestre, em MJ m⁻² d⁻¹
- `nN` — razão de insolação
- $PPBp$ (retorno) — em kg MS ha⁻¹ d⁻¹
- os coeficientes são os de Doorenbos & Kassam (0,36 e 0,219, com $Q_o$ em
  cal cm⁻² d⁻¹) convertidos para MJ m⁻² d⁻¹

### `CIAF(IAF, IAF_padrao=5.0)` · `CR(T, T_limiar=20.0, CR_frio=0.6, CR_quente=0.5)`

```math
CIAF = 0{,}0093 + 0{,}185\,IAF - 0{,}0175\,IAF^2 \;\;(IAF < 5), \qquad CIAF = 0{,}5 \;\;(IAF \ge 5)
```

```math
CR = 0{,}6 \;\;(T < 20\ °C), \qquad CR = 0{,}5 \;\;(T \ge 20\ °C)
```

### `produtividade_potencial(df, lat, rota, IAF=5.0, Cc=0.35, U=13.0, intervalo='dec', semeadura=None, ciclo=None, ...)`
Acumula a produtividade potencial período a período ao longo do ciclo:

```math
PPf = \frac{\sum_{p} PPBp_p \cdot ND_p \cdot CIAF \cdot CR_p \cdot C_c}{1 - 0{,}01\,U}
```

**Onde:**
- `df` — uma linha por período: `Tmed` (°C); radiação como `nN`, como `n` e
  `N`, ou como `Qg` (MJ m⁻² d⁻¹); `Qo` (opcional, calculada a partir de `lat`
  e das datas); datas em `data`, DatetimeIndex ou `dia` e `mes` (+ `ano`);
  `ND` (opcional)
- `intervalo` — escala de tempo: `'d'` (diário), `'dec'` (decendial do
  calendário: 1–10, 11–20 e 21–fim do mês, com 8 a 11 dias) ou `'M'` (mensal)
- $ND_p$ — dias do período com a cultura no campo. Com `semeadura` e `ciclo`
  (dias, contando o da semeadura), os períodos das pontas ficam parciais e os
  de fora do ciclo são descartados
- `IAF`, `Cc`, `U` — IAF máximo, índice de colheita e umidade do produto, em %
- retorno — DataFrame com `ND`, `Qo`, `nN`, `cTn`, `cTc`, `PPBc`, `PPBn`,
  `PPBp`, `CR`, `PPf_periodo` e `PPf_acum` (kg ha⁻¹); `attrs['PPf']`

### `produtividade_atingivel(PPf, ETR, ETc, ky, fase=None, metodo='produtorio')`
Penalização pelo déficit hídrico:

```math
1 - \frac{PA}{PPf} = k_y \left(1 - \frac{ETR}{ETc}\right)
```

- **Produtório** (Battisti et al., 2013), fase a fase em sequência:

  $`\displaystyle PA = PPf \cdot \prod_i \left[1 - k_{y,i}\left(1 - \frac{ETR_i}{ETc_i}\right)\right]`$

- **Etapa única**, com o $k_y$ do ciclo:

  $`\displaystyle PA = PPf \cdot \left[1 - k_{y,ciclo}\left(1 - \frac{\sum ETR}{\sum ETc}\right)\right]`$

```math
EC = \frac{PA}{PPf}
```

**Onde:**
- `ETR`, `ETc` — em mm, por fase ou por período (ex.: saída de
  `balanco_hidrico_cultura`); com `fase`, são somados por fase
- `ky` — um por fase (sequência ou dict `{fase: ky}`) ou o do ciclo
- cada fator é limitado a zero; fases com $ETc = 0$ não penalizam
- retorno — uma linha por fase com `ETc`, `ETR`, `ETR/ETc`, `ky`, `fator`,
  `PA_inicial`, `PA_final`, `perda` e `EC`; `attrs['PA']` e `attrs['EC']`

**Referências:** Doorenbos, J.; Kassam, A. H. (1979) *Yield response to water*,
FAO Irrigation and Drainage Paper 33 · Barbieri, V.; Tuon, R. L. (1992) ESALQ/USP ·
Pereira, A. R.; Angelocci, L. R.; Sentelhas, P. C. (2002) *Agrometeorologia* ·
Battisti, R. et al. (2013) *Ciência Rural* 43:390–396
