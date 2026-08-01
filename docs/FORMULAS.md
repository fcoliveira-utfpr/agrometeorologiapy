# Fórmulas de Agrometeorologia — Referência das Funções

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

---

## 1. Funções auxiliares de trigonometria

Módulo `agrometeorologiapy._trig`. Atalhos usados internamente por quase
todas as demais funções, para evitar conversões manuais grau↔radiano.

### `sind(x)`
$$ \text{sind}(x) = \sin\left(x \cdot \frac{\pi}{180}\right) $$

### `cosd(x)`
$$ \text{cosd}(x) = \cos\left(x \cdot \frac{\pi}{180}\right) $$

### `tand(x)`
$$ \text{tand}(x) = \tan\left(x \cdot \frac{\pi}{180}\right) $$

### `acosd(x)`
$$ \text{acosd}(x) = \arccos(x) \cdot \frac{180}{\pi} $$

---

## 2. Radiação Solar

Módulo `agrometeorologiapy.radiacao`.

### `nda(dia, mes, ano=2023)`
Número do Dia do Ano (NDA): posição ordinal do dia dentro do ano
(1 de janeiro = 1, 31 de dezembro = 365 ou 366). Calculado via
`datetime.timetuple().tm_yday`, sem fórmula fechada.

### `declinacao_solar(NDA)`
Aproximação senoidal da declinação solar δ, decorrente da inclinação do
eixo terrestre (23,45°):
$$ \delta = 23{,}45 \cdot \sin\left(\frac{360 (NDA - 80)}{365}\right) \quad [\text{graus}] $$

### `angulo_horario(hora, minuto=0)`
$$ h = (\text{hora} + \tfrac{\text{minuto}}{60} - 12) \times 15 \quad [\text{graus}] $$
(15°/hora, nulo ao meio-dia solar, negativo pela manhã, positivo à tarde)

### `angulo_zenital(lat, declinacao)`
Ângulo zenital do Sol ao meio-dia local (h = 0°):
$$ \cos Z = \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta) $$
$$ Z = \arccos(\cos Z) $$
(o valor de $\cos Z$ é limitado ao intervalo $[-1, 1]$ antes do arco-cosseno,
para evitar erro de domínio por arredondamento)

### `azimute_solar(lat, declinacao, Z)`
Direção horizontal do Sol em relação à linha Norte-Sul:
$$ \cos \alpha = \frac{\sin(\varphi)\cos(Z) - \sin(\delta)}{\cos(\varphi)\sin(Z)} $$
$$ \alpha = \arccos(\cos \alpha) $$

### `comprimento_sombra(d, Z)`
Comprimento da sombra projetada por um objeto de altura `d`:
$$ S = d \cdot \tan(Z) $$

### `fotoperiodo(Hn)`
Duração do dia (dobro do ângulo horário do nascer, convertido para horas):
$$ N = \frac{2 \, H_n}{15} \quad [\text{horas}] $$

### `angulo_horario_nascer(lat, declinacao)`
Obtida impondo-se $Z = 90°$ (cos Z = 0) na equação do ângulo zenital:
$$ H_n = \arccos(-\tan(\varphi)\tan(\delta)) $$

### `fator_correcao_distancia(NDA)`
Correção $(d/D)^2$ da excentricidade da órbita terrestre:
$$ \left(\frac{d}{D}\right)^2 = 1 + 0{,}033 \cdot \cos\left(\frac{360 \cdot NDA}{365}\right) $$

### `irradiancia_extraterrestre(lat, declinacao, Hn, dD2)`
Irradiância solar global no topo da atmosfera (Qo):
$$ Q_o = 37{,}6 \cdot \left(\frac{d}{D}\right)^2 \cdot \left[ H_{n,\text{rad}} \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)\sin(H_n) \right] \quad [\text{MJ/m}^2\text{dia}] $$
onde $H_{n,\text{rad}}$ é $H_n$ convertido para radianos apenas nesse termo.

### `insolacao(N, Tmax, Tmin, lat, k=0.19)`
Estimativa do número de horas de brilho solar (n) a partir da amplitude
térmica diária:
$$ n = \frac{N}{0{,}52} \cdot \left[ k \sqrt{T_{max} - T_{min}} - 0{,}29 \cos(\varphi) \right] \quad [\text{horas}] $$
`k = 0,16` para regiões interioranas, `k = 0,19` (padrão) para regiões
costeiras.

### `Qg_angstrom(insolacao, N, Qo, lat, b=0.52)`
Equação de Angström-Prescott (variante Glover-McCulloch):
$$ a = 0{,}29 \cos(\varphi) $$
$$ Q_g = Q_o \cdot \left( a + b \cdot \frac{n}{N} \right) \quad [\text{MJ/m}^2] $$

### `Qg_hargreaves(Tmax, Tmin, Qo, k=0.16)`
Equação de Hargreaves-Samani, sem depender de dados de insolação:
$$ Q_g = k \sqrt{T_{max} - T_{min}} \cdot Q_o \quad [\text{MJ/m}^2] $$

---

## 3. Temperatura

Módulo `agrometeorologiapy.temperatura`.

### `temp_media_extremos(Tmax, Tmin)`
$$ T_{med} = \frac{T_{max} + T_{min}}{2} $$

### `temp_media_estacao_automatica(temperaturas)`
Média aritmética simples de $n$ observações no período:
$$ T_{med} = \frac{1}{n}\sum_{i=1}^{n} T_i $$

---

## 4. Umidade do Ar

Módulo `agrometeorologiapy.umidade`.

### `es_tetens(T_ar)`
Equação de Tetens — pressão de saturação de vapor d'água:
$$ e_s = 0{,}6108 \cdot 10^{\frac{7{,}5 \, T}{237{,}3 + T}} \quad [\text{kPa}] $$

### `ea_umidade(es, UR)`
Pressão parcial (atual) de vapor:
$$ e_a = e_s \cdot \frac{UR}{100} $$

### `deficit_saturacao(es, ea)`
$$ \Delta e = e_s - e_a $$

### `patm_altitude(A)`
Pressão atmosférica local em função da altitude:
$$ P_{atm} = 101{,}3 \cdot \left( \frac{293 - 0{,}0065 A}{293} \right)^{5{,}26} \quad [\text{kPa}] $$

### `umidade_absoluta(ea, T_ar_C)`
$$ UA = \frac{2168 \cdot e_a}{T + 273{,}15} \quad [\text{g de H}_2\text{O / m}^3] $$

### `umidade_saturacao(es, T_ar_C)`
$$ US = \frac{2168 \cdot e_s}{T + 273{,}15} \quad [\text{g de H}_2\text{O / m}^3] $$

### `umidade_relativa(ea, es)`
$$ UR = 100 \cdot \frac{e_a}{e_s} \quad [\%] $$

### `ponto_orvalho(ea)`
Inversão algébrica da equação de Tetens:
$$ x = \log_{10}\left(\frac{e_a}{0{,}6108}\right) $$
$$ T_o = \frac{237{,}3 \, x}{7{,}5 - x} \quad [°C] $$

### `constante_psicrometrica(Patm)`
$$ \gamma = 0{,}665 \times 10^{-3} \cdot P_{atm} \quad [\text{kPa}/°C] $$

---

## 5. Balanço de Energia

Módulo `agrometeorologiapy.energia`.

### `boc_saldo(Qg, r=0.25)`
Balanço de ondas curtas — radiação global incidente menos a fração refletida
(albedo `r`; 0,25 é a média para gramado):
$$ BOC = Q_g (1 - r) $$

### `bol_saldo(Tmax, Tmin, ea, Qg, Qg_cs)`
Balanço de ondas longas, equação de Stefan-Boltzmann corrigida pela FAO-56.
`Tmax`/`Tmin` são recebidos em °C (consistente com o resto do pacote) e
convertidos internamente para Kelvin, porque o termo de Stefan-Boltzmann
exige temperatura absoluta:
$$ T_{max,K} = T_{max} + 273{,}15, \qquad T_{min,K} = T_{min} + 273{,}15 $$
$$ BOL = -\left[ 4{,}903 \times 10^{-9} \cdot \frac{T_{max,K}^4 + T_{min,K}^4}{2} \right] \cdot \left[ 0{,}34 - 0{,}14\sqrt{e_a} \right] \cdot \left[ 1{,}35 \frac{Q_g}{Q_{g,cs}} - 0{,}35 \right] $$

> 🐛 **Nota:** o notebook original de origem aplicava esse termo diretamente
> sobre Tmax/Tmin em °C (sem converter para Kelvin), o que subestimava o BOL
> em várias ordens de grandeza. Foi corrigido nesta biblioteca — veja o
> [CHANGELOG](../CHANGELOG.md).

### `saldo_radiacao(BOC, BOL)`
$$ R_n = BOC + BOL \quad [\text{MJ/m}^2\text{dia}] $$

---

## 6. Evapotranspiração

Módulo `agrometeorologiapy.evapotranspiracao`.

### `thornthwaite_mensal(df, col_T='T_media_C', lat=None)`
Método de Thornthwaite (1948), mensal. Para cada um dos 12 meses:
$$ i = \left(\frac{T}{5}\right)^{1{,}514} \quad \text{se } T > 0, \text{ senão } 0 $$
$$ I = \sum_{m=1}^{12} i_m \qquad \text{(índice de calor anual)} $$
$$ a = 6{,}75 \times 10^{-7} I^3 - 7{,}71 \times 10^{-5} I^2 + 1{,}792 \times 10^{-2} I + 0{,}49239 $$
$$ ETP_{nc} = 16 \left( \frac{10 T}{I} \right)^{a} \quad \text{(não corrigida, mm/mês de 30 dias e 12h de sol)} $$

Correção pelo fotoperíodo real do mês (via declinação solar do dia juliano
médio de cada mês) e pelo número real de dias do mês:
$$ K = \frac{N}{12} \cdot \frac{\text{dias\_mes}}{30} $$
$$ ETP = ETP_{nc} \cdot K \quad [\text{mm/mês}] $$

### `camargo_maluf_mensal(df, col_T='T_media_C', lat=None)`
Método de Camargo (1971), coeficiente F modificado por Maluf (Camargo et
al., 1999). Função autocontida: recalcula internamente a radiação
extraterrestre mensal ($Q_o$, a partir da declinação solar do dia juliano
médio de cada mês) e converte para equivalente de evaporação:
$$ Q_{o,mm} = 0{,}408 \cdot Q_{o} $$

O coeficiente F depende da temperatura média **anual** ($\bar{T}_{anual}$):

- $F = 0{,}0100$ se $\bar{T}_{anual} < 23°C$
- $F = 0{,}0105$ se $23°C \le \bar{T}_{anual} < 24°C$
- $F = 0{,}0110$ se $\bar{T}_{anual} \ge 24°C$

$$ ETP = F \cdot Q_{o,mm} \cdot T \cdot \text{dias\_mes} \quad [\text{mm/mês}] $$

### `etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)`
$$ Q_{o,mm} = 0{,}408 \cdot Q_o $$
$$ ETP = 0{,}0023 \cdot Q_{o,mm} \cdot \sqrt{T_{max} - T_{min}} \cdot (T_{med} + 17{,}8) \quad [\text{mm/dia}] $$

### `declive_pressao_vapor(T_ar)`
Derivada da equação de Tetens em relação à temperatura:
$$ \Delta = \frac{4098 \cdot e_s(T)}{(T + 237{,}3)^2} \quad [\text{kPa}/°C] $$

### `etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26)`
$$ ETP = \alpha \cdot \frac{\Delta}{\Delta + \gamma} \cdot \frac{R_n - G}{\lambda} \quad [\text{mm/dia}] $$
onde $\lambda = 2{,}45$ MJ/kg (calor latente de vaporização) e
$\alpha = 1{,}26$ (coeficiente padrão de Priestley-Taylor).

### `eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)`
Equação de Penman-Monteith padronizada pelo boletim FAO-56 (Allen et al.,
1998), referenciada a uma cultura hipotética (grama, 0,12 m, albedo 0,23):
$$ ETo = \frac{0{,}408 \, \Delta (R_n - G) + \gamma \cdot \frac{900}{T_{med}+273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma (1 + 0{,}34 \, u_2)} \quad [\text{mm/dia}] $$

---

## 7. Grau-Dias

Módulo `agrometeorologiapy.grau_dias`.

### Regra do grau-dia diário (GDi)
Usada tanto por `data_maturacao_fisiologica` quanto por `data_semeadura`,
onde $T_b$ é a temperatura base da cultura:

- $GD_i = T_{med} - T_b$, se $T_b < T_{min}$
- $GD_i = \dfrac{(T_{max} - T_b)^2}{2(T_{max} - T_{min})}$, se $T_b \ge T_{min}$

### `data_maturacao_fisiologica(df, Tb, CT, dia_semeadura, mes_semeadura, intervalo='d', ano=2023)`
A partir da data de semeadura, acumula $GD_i \times n_{período}$
período a período (diário, decendial ou mensal) até que o acumulado
atinja a constante térmica do ciclo:
$$ \sum GD_i \cdot n_{período} \ge CT $$
Retorna a data em que isso ocorre — a maturação fisiológica.

### `data_semeadura(df, Tb, CT, dia_maturacao, mes_maturacao, intervalo='d', ano=2023)`
O mesmo acúmulo de graus-dia, mas percorrendo o calendário **de trás para
frente** a partir de uma data de maturação conhecida (ex.: colheita-alvo),
até acumular `CT` — retornando a data de semeadura necessária.

---

## 8. Balanço Hídrico

Módulo `agrometeorologiapy.balanco_hidrico`. Ambas as funções implementam o
método de contabilidade sequencial de Thornthwaite & Mather (1955).

### `balanco_hidrico_climatologico(df, CAD=100.0)`
Para cada período (mês), com $CAD$ = capacidade de água disponível no solo:
$$ P - ETP $$

- Se $P - ETP < 0$ (déficit): acumula o negativo e recalcula o
  armazenamento por via exponencial —
  $$ \text{NEG.ACUM}_i = \text{NEG.ACUM}_{i-1} + (P - ETP) $$
  $$ ARM_i = CAD \cdot e^{\,\text{NEG.ACUM}_i / CAD} $$
- Se $P - ETP \ge 0$ (reposição): o solo recebe água até no máximo `CAD` —
  $$ ARM_i = \min(ARM_{i-1} + (P - ETP),\; CAD) $$
  e, se $ARM_i < CAD$, o NEG.ACUM é recalculado por inversão:
  $$ \text{NEG.ACUM}_i = CAD \cdot \ln(ARM_i / CAD) $$

A partir do armazenamento, derivam-se:
$$ ALT_i = ARM_i - ARM_{i-1} $$

- $ETR_i = P_i + |ALT_i|$, se $P_i - ETP_i < 0$
- $ETR_i = ETP_i$, caso contrário

$$ DEF_i = ETP_i - ETR_i $$

- $EXC_i = (P_i - ETP_i) - ALT_i$, se $P_i - ETP_i > 0$ e $ARM_i = CAD$
- $EXC_i = 0$, caso contrário

### `balanco_hidrico_cultura(df)`
Mesmo método, aplicado a uma cultura específica com `Chuva`, `ETc` (=
$K_c \times ETo$, já considerando a fase fenológica) e `CAD` (=
$z \times DTA$, já considerando o avanço da profundidade radicular)
pré-calculados pelo usuário para cada período. Além de `ARM`, `ALT`, `ETR`,
`DEF` e `EXC` (mesmas fórmulas acima, trocando P→Chuva e ETP→ETc), calcula
também o Índice de Satisfação das Necessidades de Água:
$$ ISNA = \frac{ETR}{ETc} $$
