# Agrometeorology Formulas — Function Reference

*[🇧🇷 Versão em português](FORMULAS.md)*

Mathematical documentation for every public function in `agrometeorologiapy`:
the original formula, what each variable means, and the expected unit. For
installation and usage examples, see the [README](../README.md); for a
runnable walkthrough, see
[`examples/tutorial_colab.ipynb`](../examples/tutorial_colab.ipynb).

> ℹ️ **Note on naming:** the package's function and parameter names stay in
> Portuguese (e.g. `declinacao_solar`, `Tmax`, `umidade_relativa`) — that is
> the actual public API, matching the source textbook. This document
> translates the surrounding explanation into English, not the identifiers
> you call in code.

**Primary source:** Pereira, Angelocci & Sentelhas (2002) — *Agrometeorologia:
fundamentos e aplicações práticas* ("Agrometeorology: fundamentals and
practical applications"), ESALQ/USP, complemented by Allen et al. (1998,
FAO-56) for the evapotranspiration methods.

> ⚠️ **Angle convention:** every trigonometric function in the package
> (`sind`, `cosd`, `tand`, `acosd`, and by extension any function that
> computes angles) works in **degrees**, not radians.

## Contents

- [1. Trigonometry helper functions](#1-trigonometry-helper-functions)
- [2. Solar Radiation](#2-solar-radiation)
- [3. Temperature](#3-temperature)
- [4. Air Humidity](#4-air-humidity)
- [5. Energy Balance](#5-energy-balance)
- [6. Evapotranspiration](#6-evapotranspiration)
- [7. Growing Degree-Days](#7-growing-degree-days)
- [8. Water Balance](#8-water-balance)
- [9. Climate Classification](#9-climate-classification)
- [10. Potential and Attainable Yield](#10-potential-and-attainable-yield)

---

## 1. Trigonometry helper functions

Module `agrometeorologiapy._trig`. Shortcuts used internally by almost every
other function, to avoid manual degree↔radian conversions.

### `sind(x)`

```math
\text{sind}(x) = \sin\left(x \cdot \frac{\pi}{180}\right)
```

**Where:** `x` — angle, in degrees.

### `cosd(x)`

```math
\text{cosd}(x) = \cos\left(x \cdot \frac{\pi}{180}\right)
```

**Where:** `x` — angle, in degrees.

### `tand(x)`

```math
\text{tand}(x) = \tan\left(x \cdot \frac{\pi}{180}\right)
```

**Where:** `x` — angle, in degrees.

### `acosd(x)`

```math
\text{acosd}(x) = \arccos(x) \cdot \frac{180}{\pi}
```

**Where:** `x` — cosine of an angle, dimensionless (between -1 and 1).
Returns the corresponding angle in degrees.

---

## 2. Solar Radiation

Module `agrometeorologiapy.radiacao`.

### `nda(dia, mes, ano=2023)`
Day of the Year number (NDA, from the Portuguese *Número do Dia do Ano*):
the ordinal position of the day within the year (January 1st = 1, December
31st = 365 or 366). Computed via `datetime.timetuple().tm_yday`, no closed
formula.

**Where:**
- `dia` — day of the month, integer (1–31)
- `mes` — month, integer (1–12)
- `ano` — year, integer (default 2023, non-leap)
- `NDA` (return) — day-of-year number, integer (1–366)

### `declinacao_solar(NDA)`
Sinusoidal approximation of solar declination δ, resulting from the tilt of
Earth's rotation axis (23.45°):

```math
\delta = 23{.}45 \cdot \sin\left(\frac{360 (NDA - 80)}{365}\right) \quad [\text{degrees}]
```

**Where:**
- $\delta$ (return) — solar declination for the day, in degrees
- $NDA$ — day-of-year number (1–366)

### `angulo_horario(hora, minuto=0)`

```math
h = (\text{hora} + \tfrac{\text{minuto}}{60} - 12) \times 15 \quad [\text{degrees}]
```

(15°/hour, zero at local solar noon, negative in the morning, positive in
the afternoon)

**Where:**
- $h$ (return) — hour angle, in degrees
- `hora` — hour of the day, integer (0–23)
- `minuto` — minute of the hour, integer (0–59)

### `angulo_zenital(lat, declinacao)`
Solar zenith angle at local solar noon (h = 0°):

```math
\cos Z = \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)
```

```math
Z = \arccos(\cos Z)
```

($\cos Z$ is clamped to $[-1, 1]$ before the arc-cosine, to avoid a domain
error from rounding)

**Where:**
- $Z$ (return) — zenith angle, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees (negative in the
  southern hemisphere)
- $\delta$ (`declinacao`) — solar declination for the day, in degrees

### `azimute_solar(lat, declinacao, Z)`
Horizontal direction of the Sun relative to the North-South line:

```math
\cos \alpha = \frac{\sin(\varphi)\cos(Z) - \sin(\delta)}{\cos(\varphi)\sin(Z)}
```

```math
\alpha = \arccos(\cos \alpha)
```

**Where:**
- $\alpha$ (return) — solar azimuth, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declinacao`) — solar declination for the day, in degrees
- $Z$ — zenith angle, in degrees (previous equation)

### `comprimento_sombra(d, Z)`
Length of the shadow cast by an object of height `d`:

```math
S = d \cdot \tan(Z)
```

**Where:**
- $S$ (return) — shadow length, in the same unit as `d`
- `d` — height of the object, in meters
- $Z$ — zenith angle at the instant considered, in degrees

### `fotoperiodo(Hn)`
Day length (twice the sunrise hour angle, converted to hours):

```math
N = \frac{2 \, H_n}{15} \quad [\text{hours}]
```

**Where:**
- $N$ (return) — photoperiod (day length), in hours
- $H_n$ (`Hn`) — sunrise hour angle, in degrees

### `angulo_horario_nascer(lat, declinacao)`
Obtained by setting $Z = 90°$ (cos Z = 0) in the general zenith-angle
equation:

```math
H_n = \arccos(-\tan(\varphi)\tan(\delta))
```

**Where:**
- $H_n$ (return) — sunrise hour angle, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declinacao`) — solar declination for the day, in degrees

### `fator_correcao_distancia(NDA)`
Correction $(d/D)^2$ for the eccentricity of Earth's orbit:

```math
\left(\frac{d}{D}\right)^2 = 1 + 0{.}033 \cdot \cos\left(\frac{360 \cdot NDA}{365}\right)
```

**Where:**
- $(d/D)^2$ (return) — Earth-Sun distance correction factor, dimensionless
- $NDA$ — day-of-year number (1–366)

### `irradiancia_extraterrestre(lat, declinacao, Hn, dD2)`
Extraterrestrial global solar irradiance at the top of the atmosphere (Qo):

```math
Q_o = 37{.}6 \cdot \left(\frac{d}{D}\right)^2 \cdot \left[ H_{n,\text{rad}} \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)\sin(H_n) \right] \quad [\text{MJ/m}^2\text{day}]
```

**Where:**
- $Q_o$ (return) — daily extraterrestrial global solar irradiance, in
  MJ/m² day
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declinacao`) — solar declination for the day, in degrees
- $H_n$ (`Hn`) — sunrise hour angle, in degrees; $H_{n,\text{rad}}$ is the
  same value converted to radians, used only in that term
- $(d/D)^2$ (`dD2`) — Earth-Sun distance correction factor, dimensionless

### `insolacao(N, Tmax, Tmin, lat, k=0.19)`
Estimate of the number of hours of sunshine (n) from the daily thermal
amplitude:

```math
n = \frac{N}{0{.}52} \cdot \left[ k \sqrt{T_{max} - T_{min}} - 0{.}29 \cos(\varphi) \right] \quad [\text{hours}]
```

**Where:**
- $n$ (return) — estimated number of hours of sunshine, in hours
- $N$ — photoperiod for the day, in hours
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $k$ — empirical coefficient, dimensionless (0.16 for inland regions, 0.19
  — the function's default — for coastal regions)

### `Qg_angstrom(insolacao, N, Qo, lat, b=0.52)`
Ångström-Prescott equation (Glover-McCulloch variant):

```math
a = 0{.}29 \cos(\varphi)
```

```math
Q_g = Q_o \cdot \left( a + b \cdot \frac{n}{N} \right) \quad [\text{MJ/m}^2]
```

**Where:**
- $Q_g$ (return) — global solar irradiance, in MJ/m²
- $n$ (`insolacao`) — measured hours of sunshine for the day, in hours
- $N$ — photoperiod for the day, in hours
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m²
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $a$ — regression coefficient derived from latitude, dimensionless
- $b$ — empirical regression coefficient, dimensionless (default 0.52)

### `relacao_n_N(Qg, Qo, lat, b=0.52, a=None, limitar=True)`
The same Ångström-Prescott equation as `Qg_angstrom`, solved for the sunshine
ratio — for when there is $Q_g$ (pyranometer or gridded data such as
BR-DWGD) but no sunshine recorder:

```math
\frac{n}{N} = \frac{Q_g / Q_o - a}{b}, \qquad a = 0{.}29 \cos(\varphi)
```

**Where:**
- $n/N$ (return) — sunshine ratio, dimensionless; with `limitar=True`
  (default), clipped to 0–1
- $Q_g$ (`Qg`) — global solar irradiance, in MJ/m² day
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m² day
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $b$ — empirical regression coefficient, dimensionless (default 0.52)
- $a$ — intercept, dimensionless (default `None`: $0{.}29 \cos\varphi$)

### `Qg_hargreaves(Tmax, Tmin, Qo, k=0.16)`
Hargreaves-Samani equation, without requiring sunshine data:

```math
Q_g = k \sqrt{T_{max} - T_{min}} \cdot Q_o \quad [\text{MJ/m}^2]
```

**Where:**
- $Q_g$ (return) — global solar irradiance, in MJ/m²
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m²
- $k$ — empirical coefficient, dimensionless (0.16 — the function's default
  — for inland regions, 0.19 for coastal regions)

---

## 3. Temperature

Module `agrometeorologiapy.temperatura`.

### `temp_media_extremos(Tmax, Tmin)`

```math
T_{med} = \frac{T_{max} + T_{min}}{2}
```

**Where:**
- $T_{med}$ (return) — estimated daily mean temperature, in °C
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C

### `temp_media_estacao_automatica(temperaturas)`
Simple arithmetic mean of $n$ observations over the period:

```math
T_{med} = \frac{1}{n}\sum_{i=1}^{n} T_i
```

**Where:**
- $T_{med}$ (return) — mean temperature for the period, in °C
- $T_i$ (`temperaturas`) — each temperature observation, in °C
- $n$ — number of observations in the list

---

## 4. Air Humidity

Module `agrometeorologiapy.umidade`.

### `es_tetens(T_ar)`
Tetens equation — saturation vapor pressure:

```math
e_s = 0{.}6108 \cdot 10^{\frac{7{.}5 \, T}{237{.}3 + T}} \quad [\text{kPa}]
```

**Where:**
- $e_s$ (return) — saturation vapor pressure, in kPa
- $T$ (`T_ar`) — air temperature, in °C

### `ea_umidade(es, UR)`
Actual (partial) vapor pressure:

```math
e_a = e_s \cdot \frac{UR}{100}
```

**Where:**
- $e_a$ (return) — actual (partial) vapor pressure, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $UR$ — relative humidity of the air, in % (e.g. 65)

### `deficit_saturacao(es, ea)`

```math
\Delta e = e_s - e_a
```

**Where:**
- $\Delta e$ (return) — vapor saturation deficit, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa

### `patm_altitude(A)`
Local atmospheric pressure as a function of altitude:

```math
P_{atm} = 101{.}3 \cdot \left( \frac{293 - 0{.}0065 A}{293} \right)^{5{.}26} \quad [\text{kPa}]
```

**Where:**
- $P_{atm}$ (return) — local atmospheric pressure, in kPa
- $A$ — altitude of the location, in meters

### `umidade_absoluta(ea, T_ar_C)`

```math
UA = \frac{2168 \cdot e_a}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3]
```

**Where:**
- $UA$ (return) — absolute humidity of the air, in g of H₂O per m³ of air
- $e_a$ (`ea`) — partial vapor pressure, in kPa
- $T$ (`T_ar_C`) — air temperature, in °C (converted internally to Kelvin,
  $T + 273{.}15$)

### `umidade_saturacao(es, T_ar_C)`

```math
US = \frac{2168 \cdot e_s}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3]
```

**Where:**
- $US$ (return) — saturation humidity of the air, in g of H₂O per m³ of air
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $T$ (`T_ar_C`) — air temperature, in °C

### `umidade_relativa(ea, es)`

```math
UR = 100 \cdot \frac{e_a}{e_s} \quad [\%]
```

**Where:**
- $UR$ (return) — relative humidity of the air, in %
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa

### `ponto_orvalho(ea)`
Algebraic inversion of the Tetens equation:

```math
x = \log_{10}\left(\frac{e_a}{0{.}6108}\right)
```

```math
T_o = \frac{237{.}3 \, x}{7{.}5 - x} \quad [°C]
```

**Where:**
- $T_o$ (return) — dew point temperature, in °C
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $x$ — auxiliary variable (logarithm of the ratio $e_a/0{.}6108$),
  dimensionless

### `constante_psicrometrica(Patm)`

```math
\gamma = 0{.}665 \times 10^{-3} \cdot P_{atm} \quad [\text{kPa}/°C]
```

**Where:**
- $\gamma$ (return) — psychrometric constant, in kPa/°C
- $P_{atm}$ (`Patm`) — local atmospheric pressure, in kPa

---

## 5. Energy Balance

Module `agrometeorologiapy.energia`.

### `boc_saldo(Qg, r=0.25)`
Net shortwave radiation — incoming global radiation minus the reflected
fraction (albedo `r`; 0.25 is the average for grass):

```math
BOC = Q_g (1 - r)
```

**Where:**
- $BOC$ (return) — net shortwave radiation, in the same unit as $Q_g$
- $Q_g$ (`Qg`) — global solar irradiance incident on the surface, in
  MJ/m² day (or W/m²)
- $r$ — surface reflection coefficient (albedo), dimensionless (0–1;
  default 0.25, typical for grass)

### `bol_saldo(Tmax, Tmin, ea, Qg, Qg_cs)`
Net longwave radiation, Stefan-Boltzmann equation corrected per FAO-56.
`Tmax`/`Tmin` are received in °C (consistent with the rest of the package)
and converted internally to Kelvin, because the Stefan-Boltzmann term
requires absolute temperature:

```math
T_{max,K} = T_{max} + 273{.}15, \qquad T_{min,K} = T_{min} + 273{.}15
```

```math
BOL = -\left[ 4{.}903 \times 10^{-9} \cdot \frac{T_{max,K}^4 + T_{min,K}^4}{2} \right] \cdot \left[ 0{.}34 - 0{.}14\sqrt{e_a} \right] \cdot \left[ 1{.}35 \frac{Q_g}{Q_{g,cs}} - 0{.}35 \right]
```

**Where:**
- $BOL$ (return) — net longwave radiation, in MJ/m² (negative: net energy
  loss through terrestrial emission)
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C; $T_{max,K}$, $T_{min,K}$ are the same in Kelvin
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $Q_g$ (`Qg`) — measured/estimated global solar radiation for the day, in
  MJ/m²
- $Q_{g,cs}$ (`Qg_cs`) — clear-sky solar radiation, in MJ/m²
- $4{.}903 \times 10^{-9}$ — Stefan-Boltzmann constant, in MJ K⁻⁴ m⁻² day⁻¹

> 🐛 **Note:** the original source notebook applied this term directly to
> Tmax/Tmin in °C (without converting to Kelvin), which underestimated BOL
> by several orders of magnitude. This was fixed in this library — see the
> [CHANGELOG](../CHANGELOG.md).

### `saldo_radiacao(BOC, BOL)`

```math
R_n = BOC + BOL \quad [\text{MJ/m}^2\text{day}]
```

**Where:**
- $R_n$ (return) — net radiation, in MJ/m² day (or W/m²)
- $BOC$ — net shortwave radiation, same unit as $R_n$
- $BOL$ — net longwave radiation, same unit as $R_n$

---

## 6. Evapotranspiration

Module `agrometeorologiapy.evapotranspiracao`.

### `thornthwaite_mensal(df, col_T='T_media_C', lat=None)`
Thornthwaite (1948) method, monthly. For each of the 12 months:

```math
i = \left(\frac{T}{5}\right)^{1{.}514} \quad \text{if } T > 0, \text{ else } 0
```

```math
I = \sum_{m=1}^{12} i_m \qquad \text{(annual heat index)}
```

```math
a = 6{.}75 \times 10^{-7} I^3 - 7{.}71 \times 10^{-5} I^2 + 1{.}792 \times 10^{-2} I + 0{.}49239
```

```math
ETP_{nc} = 16 \left( \frac{10 T}{I} \right)^{a} \quad \text{(uncorrected, mm/month of 30 days and 12h of sunlight)}
```

For $T \ge 26{.}5$ °C the original method uses a table instead of the equation above:

```math
ETP_{nc} = -415{.}85 + 32{.}24\,T - 0{.}43\,T^2
```

Corrected by the month's real photoperiod (via the solar declination for
the month's mean Julian day) and by the actual number of days in the
month:

```math
K = \frac{N}{12} \cdot \frac{\text{days}_{\text{month}}}{30}
```

```math
ETP = ETP_{nc} \cdot K \quad [\text{mm/month}]
```

**Where:**
- $ETP$ (return, column `ETP_mm_mes`) — corrected potential
  evapotranspiration, in mm/month
- $T$ (`col_T` column of `df`) — mean temperature for each month, in °C
- $i$ — monthly heat index, dimensionless; $I$ — sum of the 12 monthly
  indices (annual heat index), dimensionless
- $a$ — empirical exponent, a function of $I$, dimensionless
- $ETP_{nc}$ — uncorrected potential evapotranspiration, in mm/month
- $N$ — photoperiod for the month (computed internally from `lat` and the
  solar declination for the month's mean Julian day), in hours
- `dias_mes` — actual number of days in the month (28–31)

### `camargo_maluf_mensal(df, col_T='T_media_C', lat=None)`
Camargo (1971) method, with the F coefficient modified by Maluf (Camargo
et al., 1999). Self-contained function: it recomputes monthly
extraterrestrial radiation internally ($Q_o$, from the solar declination
for the month's mean Julian day) and converts it to evaporation
equivalent:

```math
Q_{o,mm} = 0{.}408 \cdot Q_{o}
```

The F coefficient depends on the **annual** mean temperature
($\bar{T}_{annual}$):

- $F = 0{.}0100$ if $\bar{T}_{annual} < 23°C$
- $F = 0{.}0105$ if $23°C \le \bar{T}_{annual} < 24°C$
- $F = 0{.}0110$ if $\bar{T}_{annual} \ge 24°C$

```math
ETP = F \cdot Q_{o,mm} \cdot T \cdot \text{days}_{\text{month}} \quad [\text{mm/month}]
```

**Where:**
- $ETP$ (return, column `ETP_mm_mes`) — monthly potential
  evapotranspiration, in mm/month
- $Q_o$ (column `Qo_MJ_m2dia`) — monthly extraterrestrial radiation, in
  MJ/m² day; $Q_{o,mm}$ (column `Qo_mm_dia`) — the same, in evaporation
  equivalent, mm/day
- $F$ — empirical coefficient, dimensionless, chosen from the annual mean
  temperature
- $\bar{T}_{annual}$ — annual mean temperature (mean of the 12 input
  monthly temperatures), in °C
- $T$ (`col_T` column of `df`) — mean temperature for each month, in °C
- `dias_mes` — actual number of days in the month (28–31)

### `etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)`

```math
Q_{o,mm} = 0{.}408 \cdot Q_o
```

```math
ETP = 0{.}0023 \cdot Q_{o,mm} \cdot \sqrt{T_{max} - T_{min}} \cdot (T_{med} + 17{.}8) \quad [\text{mm/day}]
```

**Where:**
- $ETP$ (return) — potential evapotranspiration, in mm/day
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m² day;
  $Q_{o,mm}$ — the same, in evaporation equivalent, mm/day
- $T_{max}$, $T_{min}$, $T_{med}$ (`Tmax`, `Tmin`, `Tmed`) — maximum,
  minimum and mean air temperatures, in °C

### `declive_pressao_vapor(T_ar)`
Derivative of the Tetens equation with respect to temperature:

```math
\Delta = \frac{4098 \cdot e_s(T)}{(T + 237{.}3)^2} \quad [\text{kPa}/°C]
```

**Where:**
- $\Delta$ (return) — slope of the saturation vapor pressure curve, in
  kPa/°C
- $T$ (`T_ar`) — air temperature, in °C (typically the daily mean
  temperature)
- $e_s(T)$ — saturation vapor pressure at temperature $T$ (`es_tetens`), in
  kPa

### `etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26)`

```math
ETP = \alpha \cdot \frac{\Delta}{\Delta + \gamma} \cdot \frac{R_n - G}{\lambda} \quad [\text{mm/day}]
```

**Where:**
- $ETP$ (return) — potential evapotranspiration, in mm/day
- $R_n$ (`Rn`) — net radiation, in MJ/m² day
- $G$ — soil heat flux, in MJ/m² day (usually $G = 0$ on a daily scale)
- $\Delta$ (`Delta`) — slope of the saturation vapor pressure curve, in
  kPa/°C
- $\gamma$ (`gamma`) — psychrometric constant, in kPa/°C
- $\lambda = 2{.}45$ — latent heat of vaporization, in MJ/kg (fixed constant
  in the code)
- $\alpha$ (`alfa`) — Priestley-Taylor coefficient, dimensionless (default
  1.26)

### `vento_2m(uz, z)`
Converts wind speed measured at any height to 2 m, using the logarithmic
wind profile (FAO-56, eq. 47):

```math
u_2 = u_z \cdot \frac{4{.}87}{\ln(67{.}8 \, z - 5{.}42)} \quad [\text{m/s}]
```

**Where:**
- $u_2$ (return) — wind speed at 2 m height, in m/s
- $u_z$ (`uz`) — wind speed measured at height $z$, in m/s
- $z$ (`z`) — measurement height above the ground, in m

### `eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)`
Penman-Monteith equation as standardized by the FAO-56 bulletin (Allen et
al., 1998), referenced to a hypothetical reference crop (grass, 0.12 m,
albedo 0.23):

```math
ETo = \frac{0{.}408 \, \Delta (R_n - G) + \gamma \cdot \frac{900}{T_{med}+273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma (1 + 0{.}34 \, u_2)} \quad [\text{mm/day}]
```

**Where:**
- $ETo$ (return) — reference evapotranspiration, in mm/day
- $R_n$ (`Rn`) — net radiation, in MJ/m² day
- $G$ — soil heat flux, in MJ/m² day ($G = 0$ on a daily scale)
- $T_{med}$ (`Tmed`) — daily mean air temperature, in °C
- $u_2$ (`u2`) — wind speed at 2 m height, in m/s
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $\Delta$ (`Delta`) — slope of the saturation vapor pressure curve, in
  kPa/°C
- $\gamma$ (`gamma`) — psychrometric constant, in kPa/°C

---

## 7. Growing Degree-Days

Module `agrometeorologiapy.grau_dias`.

### Daily growing-degree-day rule (GDi)
Used by both `data_maturacao_fisiologica` and `data_semeadura`:

- $GD_i = T_{med} - T_b$, if $T_b < T_{min}$
- $GD_i = \dfrac{(T_{max} - T_b)^2}{2(T_{max} - T_{min})}$, if $T_b \ge T_{min}$

**Where:**
- $GD_i$ — growing degree-days for the period (day, dekad or month), in
  °C·day
- $T_{med}$, $T_{max}$, $T_{min}$ — mean, maximum and minimum temperatures
  for the period, in °C
- $T_b$ — crop base temperature (below which there is no development), in
  °C

### `data_maturacao_fisiologica(df, Tb, CT, dia_semeadura, mes_semeadura, intervalo='d', ano=2023)`
Starting from the sowing date, accumulates $GD_i \times n_{period}$
period by period (daily, dekadal or monthly) until the accumulated sum
reaches the cycle's thermal constant:

```math
\sum GD_i \cdot n_{period} \ge CT
```

Returns the date on which that happens — physiological maturity. In the
last period only the days needed to complete `CT` are counted
($\lceil (CT - GDA_{previous}) / GD_i \rceil$), so the last row of the
DataFrame is the maturity date itself. The year rolls over when the cycle
crosses December → January.
The result has one row per period with the columns `Meses` (abbreviated
month), `data`, `Tmed`, `GDi` (daily degree-days), `GDA_mes` ($GD_i \times$
days counted in the period) and `GDA_ciclo` (cycle accumulation).

**Where:**
- `df` — climate series (columns `dia`, `mes`, `Tmed`, `Tmax`, `Tmin`), in
  chronological order
- `Tb` — crop base temperature, in °C
- `CT` — total thermal constant for the cycle, in °C·day
- `dia_semeadura`, `mes_semeadura` — sowing date, integers
- `intervalo` — `'d'` (daily), `'dec'` (dekadal, 10 days) or `'M'` (monthly)
- $n_{period}$ — number of days in the period (1 for daily, 10 for
  dekadal, days in the month for monthly)
- `ano` — reference year, integer (default 2023, non-leap)

### `data_semeadura(df, Tb, CT, dia_maturacao, mes_maturacao, intervalo='d', ano=2023)`
The same growing-degree-day accumulation, but walking the calendar
**backwards** from a known maturity date (e.g. a target harvest date),
until `CT` is accumulated — returning the required sowing date.

**Where:** same parameters as `data_maturacao_fisiologica`, swapping
`dia_semeadura`/`mes_semeadura` (known starting date) for
`dia_maturacao`/`mes_maturacao` (reference date from which the calendar is
walked backwards).

---

## 8. Water Balance

Module `agrometeorologiapy.balanco_hidrico`. Both functions implement the
sequential accounting method of Thornthwaite & Mather (1955).

### `balanco_hidrico_climatologico(df, CAD=100.0, ciclico=True)`
For each period (month) $i$:

```math
P - ETP
```

- If $P - ETP < 0$ (deficit): accumulate the negative value and recompute
  soil water storage through the exponential model —

  $`\displaystyle \text{NEG.ACUM}_i = \text{NEG.ACUM}_{i-1} + (P - ETP)`$

  $`\displaystyle ARM_i = CAD \cdot e^{\,\text{NEG.ACUM}_i / CAD}`$

- If $P - ETP \ge 0$ (replenishment): the soil receives water up to a
  maximum of `CAD` —

  $`\displaystyle ARM_i = \min(ARM_{i-1} + (P - ETP),\; CAD)`$

  and, if $ARM_i < CAD$, NEG.ACUM is recomputed by inversion:

  $`\displaystyle \text{NEG.ACUM}_i = CAD \cdot \ln(ARM_i / CAD)`$

From soil water storage, the following are derived:

```math
ALT_i = ARM_i - ARM_{i-1}
```

- $ETR_i = P_i + |ALT_i|$, if $P_i - ETP_i < 0$
- $ETR_i = ETP_i$, otherwise

```math
DEF_i = ETP_i - ETR_i
```

- $EXC_i = (P_i - ETP_i) - ALT_i$, if $P_i - ETP_i > 0$ and $ARM_i = CAD$
- $EXC_i = 0$, otherwise

**Where:**
- `df` — columns `Meses` (months), `P (mm/mês)` (precipitation) and
  `ETP (mm/mês)` (potential evapotranspiration)
- `CAD` — soil water holding capacity, in mm (default 100.0)
- $P_i$, $ETP_i$ — precipitation and potential evapotranspiration for
  period $i$, in mm/month
- $ARM_i$ (column `ARM (mm/mês)`) — soil water storage at the end of
  period $i$, in mm; $i-1$ is the previous period. With `ciclico=True`
  (default), the storage before January is December's storage at the
  annual-cycle equilibrium (the 12 months are repeated until December's ARM
  converges); with `ciclico=False`, it starts from $ARM = CAD$ (full soil)
- $\text{NEG.ACUM}_i$ (column `NEG.ACUM (mm)`) — accumulated negative
  $P-ETP$, an auxiliary variable in mm, used in the exponential drying
  model
- $ALT_i$ (column `ALT (mm/mês)`) — change in storage between periods, in
  mm
- $ETR_i$ (column `ETR (mm/mês)`) — actual evapotranspiration, in mm/month
- $DEF_i$ (column `DEF (mm/mês)`) — water deficiency, in mm/month
- $EXC_i$ (column `EXC (mm/mês)`) — water surplus, in mm/month

### `balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0, ciclico=True, tol=0.01, max_iter=100)`
Same equations as `balanco_hidrico_climatologico`, vectorized with numpy
for many locations at once (e.g. every pixel of a raster with the monthly
climate normal).

**Where:**
- `P`, `ETP` — arrays of shape `(12, ...)`, in mm/month, January to
  December along the first axis
- `CAD` — scalar or array with the shape of the spatial dimensions, in mm
  (> 0; use NaN for locations without data)
- `ciclico` — equilibrium initial storage of the annual cycle (default) or
  full soil before January
- `tol`, `max_iter` — tolerance (mm) and maximum number of annual cycles
  in the search for equilibrium
- return — `dict` of `(12, ...)` arrays: `P-ETP`, `ARM`, `NEG.ACUM`,
  `ALT`, `ETR`, `DEF`, `EXC`

### `balanco_hidrico_cultura(df)`
The same method, applied to a specific crop with `Chuva` (rainfall),
`ETc` and `CAD` pre-computed by the user for each period (same formulas
for `ARM`, `ALT`, `ETR`, `DEF` and `EXC` as the function above, swapping
$P \to$ `Chuva` and $ETP \to$ `ETc`), plus the Water Requirement
Satisfaction Index:

```math
ISNA = \frac{ETR}{ETc}
```

**Where:**
- `df` — columns `Chuva` (rainfall, mm/period), `ETc` (crop
  evapotranspiration, mm/period) and `CAD` (available water capacity for
  the period, mm), in chronological order
- `ETc` — crop evapotranspiration, mm/period; pre-computed by the user as
  $K_c \times ETo$ (crop coefficient × reference evapotranspiration),
  already accounting for the phenological stage
- `CAD` — available water capacity for the period, mm; pre-computed by the
  user as $z \times DTA$ (root depth × total available water), already
  accounting for root depth advancement
- $ISNA$ (return, column `ISNA`) — Water Requirement Satisfaction Index
  (*Índice de Satisfação das Necessidades de Água*), dimensionless (0–1:
  the closer to 1, the lower the crop's water stress)

---

## 9. Climate Classification

Module `agrometeorologiapy.classificacao_climatica`. Each method has two
versions: one for **a single site** (lists of 12 months, Jan–Dec) and a
`_grade` one for **many sites** (e.g. every pixel of a raster), with the
12 months on the first axis, shape `(12, ...)`. Sites with NaN in any month
get code 0 and class `''`.

**Seasons** (Thornthwaite and Camargo): astronomical, weighted by the
fraction of each month in the season (southern hemisphere; in the north,
summer ↔ winter and autumn ↔ spring):

| Season | Months |
|---|---|
| Summer | ⅓ Dec + Jan + Feb + ⅔ Mar |
| Autumn | ⅓ Mar + Apr + May + ⅔ Jun |
| Winter | ⅓ Jun + Jul + Aug + ⅔ Sep |
| Spring | ⅓ Sep + Oct + Nov + ⅔ Dec |

> Table 1 of Aparecido et al. (2016) has ⅓ Jun in autumn and ⅔ Jun in
> winter; here June follows the astronomical calendar (winter starts on 21 June).

### `classificacao_koppen(T, P, lat)` · `classificacao_koppen_grade(T, P, lat)`
Köppen-Geiger following Alvares et al. (2013), with the f/s/w seasonality of
groups C and D from Kottek et al. (2006). Variables: $T_{ann}$ (annual mean),
$T_{cold}$ and $T_{hot}$ (coldest and warmest month), $P_{ann}$ and $P_{dry}$
(annual total and driest month); summer = Oct–Mar in the southern hemisphere.

- **A** (tropical): $T_{cold} \ge 18$ °C — Af if $P_{dry} \ge 60$ mm; Am if
  $P_{dry} \ge 100 - P_{ann}/25$; otherwise As (dry summer) or Aw (dry winter)
- **B** (arid, overrides all others): $P_{ann} < 10 P_{th}$, with
  $P_{th} = 2T_{ann} + 14$ (or $2T_{ann}$ if ≥ 70% of rain falls in winter;
  $2T_{ann} + 28$ if ≥ 70% falls in summer) — BW if $P_{ann} < 5 P_{th}$,
  otherwise BS; h if $T_{ann} \ge 18$ °C, otherwise k
- **C**: $-3 < T_{cold} < 18$ °C and $T_{hot} > 10$ °C; **D**: $T_{cold} \le -3$ °C;
  **E**: $T_{hot} \le 10$ °C (ET if $T_{hot} > 0$, otherwise EF)
- C/D — s: dry summer ($P_{s,dry} < P_{w,dry}$, $P_{w,wet} > 3P_{s,dry}$ and
  $P_{s,dry} < 40$ mm); w: dry winter ($P_{w,dry} < P_{s,dry}$ and
  $P_{s,wet} > 10P_{w,dry}$); f: neither s nor w. a: $T_{hot} \ge 22$ °C;
  b: ≥ 4 months above 10 °C; c: 1–3 months; d: $T_{cold} < -38$ °C

**Where:**
- `T` — monthly mean temperature, in °C
- `P` — monthly precipitation, in mm/month
- `lat` — latitude, in degrees (negative in the southern hemisphere)
- return — `id` (1–31) and `classe` (e.g. `'Cfa'`)

### `classificacao_thornthwaite(P, ETP, lat, CAD=100)` · `classificacao_thornthwaite_grade(...)`
Thornthwaite (1948), from the climatological water balance
(`balanco_hidrico_climatologico_grade`):

```math
I_h = 100\,\frac{EXC}{ETP} \qquad I_a = 100\,\frac{DEF}{ETP} \qquad I_m = I_h - 0{.}6\,I_a
```

```math
ETP_{summer}\,(\%) = 100\,\frac{ETP_{summer}}{ETP_{annual}}
```

Class = moisture ($I_m$) + subtype + thermal efficiency ($ETP_{annual}$) +
summer concentration ($ETP_{summer}$, %), e.g. `B1rA'a'`.

| $I_m$ | Moisture | | $ETP_{annual}$ (mm) | Thermal |
|---|---|---|---|---|
| ≥ 100 | A (perhumid) | | ≥ 1140 | A' (megathermal) |
| 80–100 / 60–80 / 40–60 / 20–40 | B4 / B3 / B2 / B1 (humid) | | 997–1140 / 855–997 / 712–855 / 570–712 | B'4 / B'3 / B'2 / B'1 (mesothermal) |
| 0–20 | C2 (moist subhumid) | | 427–570 / 285–427 | C'2 / C'1 (microthermal) |
| −20–0 | C1 (dry subhumid) | | 142–285 | D' (tundra) |
| −40– −20 | D (semiarid) | | < 142 | E' (perpetual frost) |
| < −40 | E (arid) | | | |

- **Subtype, humid climates** ($I_m \ge 0$), by deficit: r ($I_a < 16{.}7$);
  s/w ($16{.}7 \le I_a < 33{.}3$); s2/w2 ($I_a \ge 33{.}3$) — s if the summer
  DEF exceeds the winter DEF, otherwise w
- **Subtype, dry climates** ($I_m < 0$), by surplus: d ($I_h < 10$); s/w
  ($10 \le I_h < 20$); s2/w2 ($I_h \ge 20$) — w if the summer EXC exceeds the
  winter EXC, otherwise s
- **Summer concentration** ($ETP_{summer}$, %): a' < 48; b'4 < 51.9;
  b'3 < 56.3; b'2 < 61.6; b'1 < 68; c'2 < 76.3; c'1 < 88; d' ≥ 88

> Tables from Aparecido et al. (2016), corrected against the original
> Thornthwaite (1948) at three points where the paper has typos: B'3/B'2
> limit = 855 mm (paper: 885); dry s2/w2 with $I_h \ge 20$ (paper: 33.3); dry
> climates with larger summer surplus = w (paper: s).

**Where:**
- `P`, `ETP` — monthly precipitation and potential evapotranspiration, in mm/month
- `lat` — latitude, in degrees (sets summer and winter)
- `CAD` — soil available water capacity, in mm (default 100)
- return — `classe`, codes (`umidade`, `subtipo`, `termica`,
  `concentracao`), `Ih`, `Ia`, `Im`, `ETP_anual`, `DEF_anual`, `EXC_anual`,
  `ETP_verao_pct` (plus `descricao` in the single-site version)

### `classificacao_camargo(T, P, ETP, lat, CAD=100, tabela_termica='coerente')` · `classificacao_camargo_grade(...)`
Camargo (1991) as modified by Maluf (2000), Tables 6–8 of Aparecido et al.
(2016). Class = thermal + `-` + water + dry-season letter, e.g. `ST-UMi`.

| $T_{ann}$ (°C) | $T_{cold}$ (°C) | Thermal |
|---|---|---|
| ≤ 3 | | GL (glacial) |
| 3–7 | | FR (frigid) |
| 7–12 | | CO (cold) |
| 12–18 | | TE (temperate) |
| 18–22 | ≤ 13 | STE (subtemperate) |
| 18–22 | 13–20 | ST (subtropical) |
| 18–22 | > 20 | TR (tropical) |
| 22–25 | | TR (tropical) |
| > 25 | | EQ (equatorial) |

| Annual DEF (mm) | Annual EXC (mm) | Water |
|---|---|---|
| > 800 | 0 | DE (desert) |
| 150–800 | 0 | AR (arid) |
| > 150 | 0–200 | SE (dry) |
| > 150 | > 200 | MO (monsoonal) |
| 0–150 | 0–200 | SB (subhumid) |
| 0–150 (> 0) | > 200 | UM (humid) |
| 0 | 200–1000 | PU (very humid) |
| 0 | > 1000 | SU (extremely humid) |

For classes SE, MO, SB and UM, the letter of the season with the largest
deficit is appended: v (summer), o (autumn), i (winter) or p (spring).

> **Reading Table 6.** The paper has "22 < $T_{ann}$ ≤ 25 **or**
> $T_{cold}$ > 20 → TR", which, read literally, rules out EQ wherever the
> coldest month exceeds 20 °C (nearly the whole equatorial region). The
> default `tabela_termica='coerente'` uses $T_{cold}$ only to split the
> 18–22 °C range (table above); `'literal'` applies the rule $T_{cold} > 20 \rightarrow TR$ in
> every case, as in the `climas_brasil` repository script.
>
> **Cases outside Table 7.** The paper does not cover 0 < DEF ≤ 150 mm with
> EXC = 0, nor DEF = 0 with EXC ≤ 200 mm; here they fall into SB. DEF and EXC
> below 0.05 mm count as zero.

**Where:**
- `T` — monthly mean temperature, in °C
- `P`, `ETP` — monthly precipitation and potential evapotranspiration, in mm/month
- `lat` — latitude, in degrees
- `CAD` — soil available water capacity, in mm (default 100)
- `tabela_termica` — `'coerente'` (default) or `'literal'`
- return — `classe`, codes (`termica`, `hidrica`, `estacao_seca`),
  `T_anual`, `T_mes_frio`, `DEF_anual`, `EXC_anual` (plus `descricao` in the
  single-site version)

### `classificacao_holdridge(T, P, lat, ETP=None, limiar_correcao=24)` · `classificacao_holdridge_grade(...)`
Holdridge life zones (38 zones; numbering and names from Jungkunst et al.,
2021, based on Leemans, 1990). Biotemperature:

```math
t^*_m = t_m - \frac{3\,|\phi|}{100}\,(t_m - 24)^2 \;\;(\text{if } t_m > 24\ °C), \qquad BT = \overline{\min(\max(t^*_m, 0), 30)}
```

```math
ETP_{annual} = 58{.}93 \cdot BT \qquad R = \frac{ETP_{annual}}{P_{annual}}
```

The zone comes from combining the biotemperature belt (polar < 1.5;
subpolar < 3; boreal < 6; cool temperate < 12; warm temperate < 18;
subtropical < 24; tropical ≥ 24 °C) with the humidity province, in ranges of
$R$ that double at each class (0.125; 0.25; 0.5; 1; 2; 4; 8; 16; 32).

**Where:**
- `T` — monthly mean temperature, in °C
- `P` — monthly precipitation, in mm/month
- `lat` — latitude ($\phi$), in degrees
- `ETP` — monthly PET, in mm/month (optional; if `None`, uses $58{.}93 \cdot BT$)
- `limiar_correcao` — monthly temperature above which the latitude
  correction applies, in °C (default 24; `None` corrects every month)
- return — `id` (1–38), `classe` (zone name in Portuguese), `classe_en`
  (original English name, as in the source),
  `biotemperatura` (°C), `P_anual`, `ETP_anual` (mm) and `razao_ETP`

**References:** Alvares, C. A. et al. (2013) *Meteorol. Z.* 22:711–728 ·
Kottek, M. et al. (2006) *Meteorol. Z.* 15:259–263 · Thornthwaite, C. W. (1948)
*Geogr. Rev.* 38:55–94 · Camargo, A. P. (1991) · Maluf, J. R. T. (2000)
*Rev. Bras. Agrometeorol.* 8:141–150 · Aparecido, L. E. O. et al. (2016)
*Ciênc. Agrotec.* 40:405–417 · Holdridge, L. R. (1967) *Life zone ecology* ·
Jungkunst, H. F. et al. (2021) *J. Plant Nutr. Soil Sci.* 184:5–11

---

## 10. Potential and Attainable Yield

Module `agrometeorologiapy.produtividade`. FAO Agro-Ecological Zone model
(Doorenbos & Kassam, 1979), with the temperature corrections of Barbieri &
Tuon (1992), in Pereira et al. (2002). Constants and tabulated values default
to the usual values, but all of them are function parameters.

### n/N
Without a sunshine recorder, the sunshine ratio comes from `relacao_n_N`
(section 2), the inverted Ångström-Prescott equation. `produtividade_potencial`
uses it when the df has `Qg` instead of `nN`.

### `cTn(T, rota, T_limiar=16.5)` · `cTc(T, rota, T_limiar=16.5)`
Temperature corrections of photosynthesis for overcast ($cT_n$) and clear
($cT_c$) sky, clipped at zero.

**C3** (soybean, wheat, bean, sunflower):

```math
cT_n = -0.0425 + 0.035\,T + 0.00325\,T^2 - 0.0000925\,T^3
```

```math
cT_c = 0.583 + 0.014\,T + 0.0013\,T^2 - 0.000037\,T^3
```

**C4** (maize, sorghum, sugarcane):

| Condition | $cT_n$ | $cT_c$ |
|---|---|---|
| $T \ge 16.5$ °C | $-1.064 + 0.173 T - 0.0029 T^2$ | $-4.16 + 0.4325 T - 0.00725 T^2$ |
| $T < 16.5$ °C | $-4.16 + 0.4325 T - 0.00725 T^2$ | $-9.32 + 0.865 T - 0.0145 T^2$ |

**Where:**
- `T` — mean air temperature of the period, °C
- `rota` — photosynthetic pathway, `'C3'` or `'C4'` (required)
- `T_limiar` — polynomial switch for C4 crops, °C (default 16.5)

### `PPBp(Qo, nN, cTc, cTn, a_c=107.2, b_c=8.604, a_n=31.7, b_n=5.234)`
Standard gross potential productivity of the standard crop (LAI = 5):

```math
PPB_c = (107.2 + 8.604\,Q_o)\, cT_c\, \frac{n}{N} \qquad PPB_n = (31.7 + 5.234\,Q_o)\, cT_n \left(1 - \frac{n}{N}\right)
```

```math
PPBp = PPB_c + PPB_n
```

**Where:**
- `Qo` — extraterrestrial irradiance, MJ m⁻² d⁻¹
- `nN` — sunshine ratio
- $PPBp$ (return) — kg DM ha⁻¹ d⁻¹
- the coefficients are Doorenbos & Kassam's (0.36 and 0.219, with $Q_o$ in
  cal cm⁻² d⁻¹) converted to MJ m⁻² d⁻¹

### `CIAF(IAF, IAF_padrao=5.0)` · `CR(T, T_limiar=20.0, CR_frio=0.6, CR_quente=0.5)`
Leaf area index correction and maintenance respiration correction:

```math
CIAF = 0.0093 + 0.185\,LAI - 0.0175\,LAI^2 \;\;(LAI < 5), \qquad CIAF = 0.5 \;\;(LAI \ge 5)
```

```math
CR = 0.6 \;\;(T < 20\ °C), \qquad CR = 0.5 \;\;(T \ge 20\ °C)
```

### `produtividade_potencial(df, lat, rota, IAF=5.0, Cc=0.35, U=13.0, intervalo='dec', semeadura=None, ciclo=None, ...)`
Accumulates potential yield period by period over the crop cycle:

```math
PPf = \frac{\sum_{p} PPBp_p \cdot ND_p \cdot CIAF \cdot CR_p \cdot C_c}{1 - 0.01\,U}
```

**Where:**
- `df` — one row per period: `Tmed` (°C); radiation as `nN`, as `n` and `N`,
  or as `Qg` (MJ m⁻² d⁻¹); `Qo` (optional, computed from `lat` and the dates);
  dates in `data`, a DatetimeIndex or `dia` and `mes` (+ `ano`); `ND` (optional)
- `intervalo` — time step: `'d'` (daily), `'dec'` (calendar 10-day period:
  1–10, 11–20 and 21–end of month, 8 to 11 days) or `'M'` (monthly)
- $ND_p$ — days of the period with the crop in the field. With `semeadura`
  (sowing date) and `ciclo` (days, counting the sowing day), the first and
  last periods are partial and periods outside the cycle are dropped
- `IAF`, `Cc`, `U` — maximum LAI, harvest index and product moisture, %
- return — DataFrame with `ND`, `Qo`, `nN`, `cTn`, `cTc`, `PPBc`, `PPBn`,
  `PPBp`, `CR`, `PPf_periodo` and `PPf_acum` (kg ha⁻¹); `attrs['PPf']`

### `produtividade_atingivel(PPf, ETR, ETc, ky, fase=None, metodo='produtorio')`
Water-deficit penalty:

```math
1 - \frac{PA}{PPf} = k_y \left(1 - \frac{ETR}{ETc}\right)
```

- **Product** (`'produtorio'`, Battisti et al., 2013), phase by phase in sequence:

  $`\displaystyle PA = PPf \cdot \prod_i \left[1 - k_{y,i}\left(1 - \frac{ETR_i}{ETc_i}\right)\right]`$

- **Single step** (`'etapa_unica'`), with the whole-cycle $k_y$:

  $`\displaystyle PA = PPf \cdot \left[1 - k_{y,cycle}\left(1 - \frac{\sum ETR}{\sum ETc}\right)\right]`$

```math
EC = \frac{PA}{PPf}
```

**Where:**
- `ETR`, `ETc` — mm, per phase or per period (e.g. output of
  `balanco_hidrico_cultura`); with `fase`, they are summed per phase
- `ky` — one per phase (sequence or dict `{phase: ky}`) or the whole-cycle value
- each factor is clipped at zero; phases with $ETc = 0$ are not penalized
- return — one row per phase with `ETc`, `ETR`, `ETR/ETc`, `ky`, `fator`,
  `PA_inicial`, `PA_final`, `perda` and `EC`; `attrs['PA']` and `attrs['EC']`

**References:** Doorenbos, J.; Kassam, A. H. (1979) *Yield response to water*,
FAO Irrigation and Drainage Paper 33 · Barbieri, V.; Tuon, R. L. (1992) ESALQ/USP ·
Pereira, A. R.; Angelocci, L. R.; Sentelhas, P. C. (2002) *Agrometeorologia* ·
Battisti, R. et al. (2013) *Ciência Rural* 43:390–396
