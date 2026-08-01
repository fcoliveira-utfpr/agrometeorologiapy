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

---

## 1. Trigonometry helper functions

Module `agrometeorologiapy._trig`. Shortcuts used internally by almost every
other function, to avoid manual degree↔radian conversions.

### `sind(x)`
$$ \text{sind}(x) = \sin\left(x \cdot \frac{\pi}{180}\right) $$
**Where:** `x` — angle, in degrees.

### `cosd(x)`
$$ \text{cosd}(x) = \cos\left(x \cdot \frac{\pi}{180}\right) $$
**Where:** `x` — angle, in degrees.

### `tand(x)`
$$ \text{tand}(x) = \tan\left(x \cdot \frac{\pi}{180}\right) $$
**Where:** `x` — angle, in degrees.

### `acosd(x)`
$$ \text{acosd}(x) = \arccos(x) \cdot \frac{180}{\pi} $$
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
$$ \delta = 23{.}45 \cdot \sin\left(\frac{360 (NDA - 80)}{365}\right) \quad [\text{degrees}] $$
**Where:**
- $\delta$ (return) — solar declination for the day, in degrees
- $NDA$ — day-of-year number (1–366)

### `angulo_horario(hora, minuto=0)`
$$ h = (\text{hora} + \tfrac{\text{minuto}}{60} - 12) \times 15 \quad [\text{degrees}] $$
(15°/hour, zero at local solar noon, negative in the morning, positive in
the afternoon)

**Where:**
- $h$ (return) — hour angle, in degrees
- `hora` — hour of the day, integer (0–23)
- `minuto` — minute of the hour, integer (0–59)

### `angulo_zenital(lat, declinacao)`
Solar zenith angle at local solar noon (h = 0°):
$$ \cos Z = \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta) $$
$$ Z = \arccos(\cos Z) $$
($\cos Z$ is clamped to $[-1, 1]$ before the arc-cosine, to avoid a domain
error from rounding)

**Where:**
- $Z$ (return) — zenith angle, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees (negative in the
  southern hemisphere)
- $\delta$ (`declinacao`) — solar declination for the day, in degrees

### `azimute_solar(lat, declinacao, Z)`
Horizontal direction of the Sun relative to the North-South line:
$$ \cos \alpha = \frac{\sin(\varphi)\cos(Z) - \sin(\delta)}{\cos(\varphi)\sin(Z)} $$
$$ \alpha = \arccos(\cos \alpha) $$
**Where:**
- $\alpha$ (return) — solar azimuth, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declinacao`) — solar declination for the day, in degrees
- $Z$ — zenith angle, in degrees (previous equation)

### `comprimento_sombra(d, Z)`
Length of the shadow cast by an object of height `d`:
$$ S = d \cdot \tan(Z) $$
**Where:**
- $S$ (return) — shadow length, in the same unit as `d`
- `d` — height of the object, in meters
- $Z$ — zenith angle at the instant considered, in degrees

### `fotoperiodo(Hn)`
Day length (twice the sunrise hour angle, converted to hours):
$$ N = \frac{2 \, H_n}{15} \quad [\text{hours}] $$
**Where:**
- $N$ (return) — photoperiod (day length), in hours
- $H_n$ (`Hn`) — sunrise hour angle, in degrees

### `angulo_horario_nascer(lat, declinacao)`
Obtained by setting $Z = 90°$ (cos Z = 0) in the general zenith-angle
equation:
$$ H_n = \arccos(-\tan(\varphi)\tan(\delta)) $$
**Where:**
- $H_n$ (return) — sunrise hour angle, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declinacao`) — solar declination for the day, in degrees

### `fator_correcao_distancia(NDA)`
Correction $(d/D)^2$ for the eccentricity of Earth's orbit:
$$ \left(\frac{d}{D}\right)^2 = 1 + 0{.}033 \cdot \cos\left(\frac{360 \cdot NDA}{365}\right) $$
**Where:**
- $(d/D)^2$ (return) — Earth-Sun distance correction factor, dimensionless
- $NDA$ — day-of-year number (1–366)

### `irradiancia_extraterrestre(lat, declinacao, Hn, dD2)`
Extraterrestrial global solar irradiance at the top of the atmosphere (Qo):
$$ Q_o = 37{.}6 \cdot \left(\frac{d}{D}\right)^2 \cdot \left[ H_{n,\text{rad}} \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)\sin(H_n) \right] \quad [\text{MJ/m}^2\text{day}] $$
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
$$ n = \frac{N}{0{.}52} \cdot \left[ k \sqrt{T_{max} - T_{min}} - 0{.}29 \cos(\varphi) \right] \quad [\text{hours}] $$
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
$$ a = 0{.}29 \cos(\varphi) $$
$$ Q_g = Q_o \cdot \left( a + b \cdot \frac{n}{N} \right) \quad [\text{MJ/m}^2] $$
**Where:**
- $Q_g$ (return) — global solar irradiance, in MJ/m²
- $n$ (`insolacao`) — measured hours of sunshine for the day, in hours
- $N$ — photoperiod for the day, in hours
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m²
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $a$ — regression coefficient derived from latitude, dimensionless
- $b$ — empirical regression coefficient, dimensionless (default 0.52)

### `Qg_hargreaves(Tmax, Tmin, Qo, k=0.16)`
Hargreaves-Samani equation, without requiring sunshine data:
$$ Q_g = k \sqrt{T_{max} - T_{min}} \cdot Q_o \quad [\text{MJ/m}^2] $$
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
$$ T_{med} = \frac{T_{max} + T_{min}}{2} $$
**Where:**
- $T_{med}$ (return) — estimated daily mean temperature, in °C
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C

### `temp_media_estacao_automatica(temperaturas)`
Simple arithmetic mean of $n$ observations over the period:
$$ T_{med} = \frac{1}{n}\sum_{i=1}^{n} T_i $$
**Where:**
- $T_{med}$ (return) — mean temperature for the period, in °C
- $T_i$ (`temperaturas`) — each temperature observation, in °C
- $n$ — number of observations in the list

---

## 4. Air Humidity

Module `agrometeorologiapy.umidade`.

### `es_tetens(T_ar)`
Tetens equation — saturation vapor pressure:
$$ e_s = 0{.}6108 \cdot 10^{\frac{7{.}5 \, T}{237{.}3 + T}} \quad [\text{kPa}] $$
**Where:**
- $e_s$ (return) — saturation vapor pressure, in kPa
- $T$ (`T_ar`) — air temperature, in °C

### `ea_umidade(es, UR)`
Actual (partial) vapor pressure:
$$ e_a = e_s \cdot \frac{UR}{100} $$
**Where:**
- $e_a$ (return) — actual (partial) vapor pressure, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $UR$ — relative humidity of the air, in % (e.g. 65)

### `deficit_saturacao(es, ea)`
$$ \Delta e = e_s - e_a $$
**Where:**
- $\Delta e$ (return) — vapor saturation deficit, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa

### `patm_altitude(A)`
Local atmospheric pressure as a function of altitude:
$$ P_{atm} = 101{.}3 \cdot \left( \frac{293 - 0{.}0065 A}{293} \right)^{5{.}26} \quad [\text{kPa}] $$
**Where:**
- $P_{atm}$ (return) — local atmospheric pressure, in kPa
- $A$ — altitude of the location, in meters

### `umidade_absoluta(ea, T_ar_C)`
$$ UA = \frac{2168 \cdot e_a}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3] $$
**Where:**
- $UA$ (return) — absolute humidity of the air, in g of H₂O per m³ of air
- $e_a$ (`ea`) — partial vapor pressure, in kPa
- $T$ (`T_ar_C`) — air temperature, in °C (converted internally to Kelvin,
  $T + 273{.}15$)

### `umidade_saturacao(es, T_ar_C)`
$$ US = \frac{2168 \cdot e_s}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3] $$
**Where:**
- $US$ (return) — saturation humidity of the air, in g of H₂O per m³ of air
- $e_s$ (`es`) — saturation vapor pressure, in kPa
- $T$ (`T_ar_C`) — air temperature, in °C

### `umidade_relativa(ea, es)`
$$ UR = 100 \cdot \frac{e_a}{e_s} \quad [\%] $$
**Where:**
- $UR$ (return) — relative humidity of the air, in %
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $e_s$ (`es`) — saturation vapor pressure, in kPa

### `ponto_orvalho(ea)`
Algebraic inversion of the Tetens equation:
$$ x = \log_{10}\left(\frac{e_a}{0{.}6108}\right) $$
$$ T_o = \frac{237{.}3 \, x}{7{.}5 - x} \quad [°C] $$
**Where:**
- $T_o$ (return) — dew point temperature, in °C
- $e_a$ (`ea`) — actual (partial) vapor pressure, in kPa
- $x$ — auxiliary variable (logarithm of the ratio $e_a/0{.}6108$),
  dimensionless

### `constante_psicrometrica(Patm)`
$$ \gamma = 0{.}665 \times 10^{-3} \cdot P_{atm} \quad [\text{kPa}/°C] $$
**Where:**
- $\gamma$ (return) — psychrometric constant, in kPa/°C
- $P_{atm}$ (`Patm`) — local atmospheric pressure, in kPa

---

## 5. Energy Balance

Module `agrometeorologiapy.energia`.

### `boc_saldo(Qg, r=0.25)`
Net shortwave radiation — incoming global radiation minus the reflected
fraction (albedo `r`; 0.25 is the average for grass):
$$ BOC = Q_g (1 - r) $$
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
$$ T_{max,K} = T_{max} + 273{.}15, \qquad T_{min,K} = T_{min} + 273{.}15 $$
$$ BOL = -\left[ 4{.}903 \times 10^{-9} \cdot \frac{T_{max,K}^4 + T_{min,K}^4}{2} \right] \cdot \left[ 0{.}34 - 0{.}14\sqrt{e_a} \right] \cdot \left[ 1{.}35 \frac{Q_g}{Q_{g,cs}} - 0{.}35 \right] $$
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
$$ R_n = BOC + BOL \quad [\text{MJ/m}^2\text{day}] $$
**Where:**
- $R_n$ (return) — net radiation, in MJ/m² day (or W/m²)
- $BOC$ — net shortwave radiation, same unit as $R_n$
- $BOL$ — net longwave radiation, same unit as $R_n$

---

## 6. Evapotranspiration

Module `agrometeorologiapy.evapotranspiracao`.

### `thornthwaite_mensal(df, col_T='T_media_C', lat=None)`
Thornthwaite (1948) method, monthly. For each of the 12 months:
$$ i = \left(\frac{T}{5}\right)^{1{.}514} \quad \text{if } T > 0, \text{ else } 0 $$
$$ I = \sum_{m=1}^{12} i_m \qquad \text{(annual heat index)} $$
$$ a = 6{.}75 \times 10^{-7} I^3 - 7{.}71 \times 10^{-5} I^2 + 1{.}792 \times 10^{-2} I + 0{.}49239 $$
$$ ETP_{nc} = 16 \left( \frac{10 T}{I} \right)^{a} \quad \text{(uncorrected, mm/month of 30 days and 12h of sunlight)} $$

Corrected by the month's real photoperiod (via the solar declination for
the month's mean Julian day) and by the actual number of days in the
month:
$$ K = \frac{N}{12} \cdot \frac{\text{dias\_mes}}{30} $$
$$ ETP = ETP_{nc} \cdot K \quad [\text{mm/month}] $$
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
$$ Q_{o,mm} = 0{.}408 \cdot Q_{o} $$

The F coefficient depends on the **annual** mean temperature
($\bar{T}_{annual}$):

- $F = 0{.}0100$ if $\bar{T}_{annual} < 23°C$
- $F = 0{.}0105$ if $23°C \le \bar{T}_{annual} < 24°C$
- $F = 0{.}0110$ if $\bar{T}_{annual} \ge 24°C$

$$ ETP = F \cdot Q_{o,mm} \cdot T \cdot \text{dias\_mes} \quad [\text{mm/month}] $$
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
$$ Q_{o,mm} = 0{.}408 \cdot Q_o $$
$$ ETP = 0{.}0023 \cdot Q_{o,mm} \cdot \sqrt{T_{max} - T_{min}} \cdot (T_{med} + 17{.}8) \quad [\text{mm/day}] $$
**Where:**
- $ETP$ (return) — potential evapotranspiration, in mm/day
- $Q_o$ (`Qo`) — extraterrestrial solar irradiance, in MJ/m² day;
  $Q_{o,mm}$ — the same, in evaporation equivalent, mm/day
- $T_{max}$, $T_{min}$, $T_{med}$ (`Tmax`, `Tmin`, `Tmed`) — maximum,
  minimum and mean air temperatures, in °C

### `declive_pressao_vapor(T_ar)`
Derivative of the Tetens equation with respect to temperature:
$$ \Delta = \frac{4098 \cdot e_s(T)}{(T + 237{.}3)^2} \quad [\text{kPa}/°C] $$
**Where:**
- $\Delta$ (return) — slope of the saturation vapor pressure curve, in
  kPa/°C
- $T$ (`T_ar`) — air temperature, in °C (typically the daily mean
  temperature)
- $e_s(T)$ — saturation vapor pressure at temperature $T$ (`es_tetens`), in
  kPa

### `etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26)`
$$ ETP = \alpha \cdot \frac{\Delta}{\Delta + \gamma} \cdot \frac{R_n - G}{\lambda} \quad [\text{mm/day}] $$
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

### `eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)`
Penman-Monteith equation as standardized by the FAO-56 bulletin (Allen et
al., 1998), referenced to a hypothetical reference crop (grass, 0.12 m,
albedo 0.23):
$$ ETo = \frac{0{.}408 \, \Delta (R_n - G) + \gamma \cdot \frac{900}{T_{med}+273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma (1 + 0{.}34 \, u_2)} \quad [\text{mm/day}] $$
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
$$ \sum GD_i \cdot n_{period} \ge CT $$
Returns the date on which that happens — physiological maturity.

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

### `balanco_hidrico_climatologico(df, CAD=100.0)`
For each period (month) $i$:
$$ P - ETP $$

- If $P - ETP < 0$ (deficit): accumulate the negative value and recompute
  soil water storage through the exponential model —
  $$ \text{NEG.ACUM}_i = \text{NEG.ACUM}_{i-1} + (P - ETP) $$
  $$ ARM_i = CAD \cdot e^{\,\text{NEG.ACUM}_i / CAD} $$
- If $P - ETP \ge 0$ (replenishment): the soil receives water up to a
  maximum of `CAD` —
  $$ ARM_i = \min(ARM_{i-1} + (P - ETP),\; CAD) $$
  and, if $ARM_i < CAD$, NEG.ACUM is recomputed by inversion:
  $$ \text{NEG.ACUM}_i = CAD \cdot \ln(ARM_i / CAD) $$

From soil water storage, the following are derived:
$$ ALT_i = ARM_i - ARM_{i-1} $$

- $ETR_i = P_i + |ALT_i|$, if $P_i - ETP_i < 0$
- $ETR_i = ETP_i$, otherwise

$$ DEF_i = ETP_i - ETR_i $$

- $EXC_i = (P_i - ETP_i) - ALT_i$, if $P_i - ETP_i > 0$ and $ARM_i = CAD$
- $EXC_i = 0$, otherwise

**Where:**
- `df` — columns `Meses` (months), `P (mm/mês)` (precipitation) and
  `ETP (mm/mês)` (potential evapotranspiration)
- `CAD` — soil water holding capacity, in mm (default 100.0)
- $P_i$, $ETP_i$ — precipitation and potential evapotranspiration for
  period $i$, in mm/month
- $ARM_i$ (column `ARM (mm/mês)`) — soil water storage at the end of
  period $i$, in mm; $i-1$ is the previous period (the first period starts
  from $ARM = CAD$, full soil)
- $\text{NEG.ACUM}_i$ (column `NEG.ACUM (mm)`) — accumulated negative
  $P-ETP$, an auxiliary variable in mm, used in the exponential drying
  model
- $ALT_i$ (column `ALT (mm/mês)`) — change in storage between periods, in
  mm
- $ETR_i$ (column `ETR (mm/mês)`) — actual evapotranspiration, in mm/month
- $DEF_i$ (column `DEF (mm/mês)`) — water deficiency, in mm/month
- $EXC_i$ (column `EXC (mm/mês)`) — water surplus, in mm/month

### `balanco_hidrico_cultura(df)`
The same method, applied to a specific crop with `Chuva` (rainfall),
`ETc` and `CAD` pre-computed by the user for each period (same formulas
for `ARM`, `ALT`, `ETR`, `DEF` and `EXC` as the function above, swapping
$P \to$ `Chuva` and $ETP \to$ `ETc`), plus the Water Requirement
Satisfaction Index:
$$ ISNA = \frac{ETR}{ETc} $$
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
