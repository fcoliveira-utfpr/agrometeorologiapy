# Agrometeorology Formulas — Function Reference

*[🇧🇷 Versão em português](FORMULAS.md)*

Mathematical documentation for every public function in `agrometeorologiapy`:
the original formula, what each variable means, and the expected unit. For
installation and usage examples, see the [README](../README.md); for a
runnable walkthrough, see
[`examples/tutorial_colab_en.ipynb`](../examples/tutorial_colab_en.ipynb).

> ℹ️ **Note on naming:** this document describes the English API,
> `agrometeorologiapy.en` (`from agrometeorologiapy import en as agm`), the
> core implementation of the package. Every function is also available with
> a Portuguese name at the top level (`import agrometeorologiapy as amp`),
> shown next to each heading as "pt: ...". Terms and symbols follow FAO-56
> (Allen et al., 1998) wherever they exist. See [API_PT_EN.md](API_PT_EN.md)
> for the full name mapping.

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
- [7. Growing Degree Days](#7-growing-degree-days)
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

Module `agrometeorologiapy.en.radiation` (Portuguese API:
`agrometeorologiapy.radiacao`; names in parentheses). Terms and symbols
follow FAO-56 (Allen et al., 1998, ch. 3).

### `day_of_year(day, month, year=2023)` · pt: `nda`
Ordinal position of the day within the year (J; January 1st = 1, December
31st = 365 or 366). Computed via `datetime.timetuple().tm_yday`, no closed
formula.

**Where:**
- `day` — day of the month, integer (1–31)
- `month` — month, integer (1–12)
- `year` — year, integer (default 2023, non-leap)
- $J$ (return) — day of the year, integer (1–366)

### `solar_declination(doy)` · pt: `declinacao_solar`
Sinusoidal approximation of solar declination δ, resulting from the tilt of
Earth's rotation axis (23.45°):

```math
\delta = 23{.}45 \cdot \sin\left(\frac{360 (J - 80)}{365}\right) \quad [\text{degrees}]
```

**Where:**
- $\delta$ (return) — solar declination for the day, in degrees
- $J$ (`doy`) — day of the year (1–366)

### `hour_angle(hour, minute=0)` · pt: `angulo_horario`

```math
\omega = (\text{hour} + \tfrac{\text{minute}}{60} - 12) \times 15 \quad [\text{degrees}]
```

(15°/hour, zero at local solar noon, negative in the morning, positive in
the afternoon)

**Where:**
- $\omega$ (return) — hour angle, in degrees
- `hour` — hour of the day, integer (0–23)
- `minute` — minute of the hour, integer (0–59)

### `solar_zenith_angle(lat, declination)` · pt: `angulo_zenital`
Solar zenith angle at local solar noon ($\omega = 0°$):

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
- $\delta$ (`declination`) — solar declination for the day, in degrees

### `solar_azimuth(lat, declination, zenith)` · pt: `azimute_solar`
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
- $\delta$ (`declination`) — solar declination for the day, in degrees
- $Z$ (`zenith`) — zenith angle, in degrees (previous equation)

### `shadow_length(height, zenith)` · pt: `comprimento_sombra`
Length of the shadow cast by an object of height $d$:

```math
S = d \cdot \tan(Z)
```

**Where:**
- $S$ (return) — shadow length, in the same unit as `height`
- $d$ (`height`) — height of the object, e.g. in meters
- $Z$ (`zenith`) — zenith angle at the instant considered, in degrees

### `daylight_hours(ws)` · pt: `fotoperiodo`
Daylight hours, the maximum possible duration of sunshine (FAO-56 eq. 34,
with $\omega_s$ in degrees):

```math
N = \frac{2 \, \omega_s}{15} \quad [\text{hours}]
```

**Where:**
- $N$ (return) — daylight hours, in hours
- $\omega_s$ (`ws`) — sunset hour angle, in degrees

### `sunset_hour_angle(lat, declination)` · pt: `angulo_horario_nascer`
Obtained by setting $Z = 90°$ (cos Z = 0) in the general zenith-angle
equation (FAO-56 eq. 25). By symmetry, it equals the sunrise hour angle in
absolute value:

```math
\omega_s = \arccos(-\tan(\varphi)\tan(\delta))
```

**Where:**
- $\omega_s$ (return) — sunset hour angle, in degrees
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declination`) — solar declination for the day, in degrees

### `inverse_relative_distance(doy)` · pt: `fator_correcao_distancia`
Inverse relative distance Earth-Sun (FAO-56 eq. 23), the correction
$(d/D)^2$ for the eccentricity of Earth's orbit:

```math
d_r = 1 + 0{.}033 \cdot \cos\left(\frac{360 \cdot J}{365}\right)
```

**Where:**
- $d_r$ (return) — inverse relative distance Earth-Sun, dimensionless
- $J$ (`doy`) — day of the year (1–366)

### `extraterrestrial_radiation(lat, declination, ws, dr)` · pt: `irradiancia_extraterrestre`
Daily extraterrestrial radiation on a horizontal surface at the top of the
atmosphere (FAO-56 eq. 21, with $37{.}6 \approx 24 \cdot 60 \cdot 0{.}0820/\pi$):

```math
R_a = 37{.}6 \cdot d_r \cdot \left[ \omega_{s,\text{rad}} \sin(\varphi)\sin(\delta) + \cos(\varphi)\cos(\delta)\sin(\omega_s) \right] \quad [\text{MJ m}^{-2}\text{day}^{-1}]
```

**Where:**
- $R_a$ (return) — extraterrestrial radiation, in MJ m⁻² day⁻¹
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $\delta$ (`declination`) — solar declination for the day, in degrees
- $\omega_s$ (`ws`) — sunset hour angle, in degrees; $\omega_{s,\text{rad}}$
  is the same value in radians, used only in that term
- $d_r$ (`dr`) — inverse relative distance Earth-Sun, dimensionless

### `sunshine_duration(N, Tmax, Tmin, lat, k=0.19)` · pt: `insolacao`
Estimate of the actual duration of sunshine (n) from the daily temperature
range:

```math
n = \frac{N}{0{.}52} \cdot \left[ k \sqrt{T_{max} - T_{min}} - 0{.}29 \cos(\varphi) \right] \quad [\text{hours}]
```

**Where:**
- $n$ (return) — estimated actual duration of sunshine, in hours
- $N$ — daylight hours, in hours
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $k$ — empirical coefficient, dimensionless (0.16 for inland regions, 0.19
  — the function's default — for coastal regions)

### `solar_radiation_angstrom(n, N, Ra, lat, b=0.52)` · pt: `Qg_angstrom`
Angstrom formula (FAO-56 eq. 35), in the Glover-McCulloch variant:

```math
a = 0{.}29 \cos(\varphi)
```

```math
R_s = R_a \cdot \left( a + b \cdot \frac{n}{N} \right) \quad [\text{MJ m}^{-2}\text{day}^{-1}]
```

**Where:**
- $R_s$ (return) — solar radiation, in MJ m⁻² day⁻¹
- $n$ — actual duration of sunshine for the day, in hours
- $N$ — daylight hours, in hours
- $R_a$ (`Ra`) — extraterrestrial radiation, in MJ m⁻² day⁻¹
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $a$ — regression coefficient derived from latitude, dimensionless
- $b$ — empirical regression coefficient, dimensionless (default 0.52)

### `relative_sunshine_duration(Rs, Ra, lat, b=0.52, a=None, clip=True)` · pt: `relacao_n_N`
The same Angstrom formula as `solar_radiation_angstrom`, solved for the
relative sunshine duration — for when there is $R_s$ (pyranometer or
gridded data such as BR-DWGD) but no sunshine recorder:

```math
\frac{n}{N} = \frac{R_s / R_a - a}{b}, \qquad a = 0{.}29 \cos(\varphi)
```

**Where:**
- $n/N$ (return) — relative sunshine duration, dimensionless; with
  `clip=True` (default), clipped to 0–1
- $R_s$ (`Rs`) — solar radiation, in MJ m⁻² day⁻¹
- $R_a$ (`Ra`) — extraterrestrial radiation, in MJ m⁻² day⁻¹
- $\varphi$ (`lat`) — latitude of the location, in degrees
- $b$ — empirical regression coefficient, dimensionless (default 0.52)
- $a$ — intercept, dimensionless (default `None`: $0{.}29 \cos\varphi$)

### `solar_radiation_hargreaves(Tmax, Tmin, Ra, k=0.16)` · pt: `Qg_hargreaves`
Hargreaves' radiation formula (FAO-56 eq. 50), without requiring sunshine
data:

```math
R_s = k_{Rs} \sqrt{T_{max} - T_{min}} \cdot R_a \quad [\text{MJ m}^{-2}\text{day}^{-1}]
```

**Where:**
- $R_s$ (return) — solar radiation, in MJ m⁻² day⁻¹
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C
- $R_a$ (`Ra`) — extraterrestrial radiation, in MJ m⁻² day⁻¹
- $k_{Rs}$ (`k`) — adjustment coefficient, dimensionless (0.16 — the
  function's default — for inland regions, 0.19 for coastal regions)

---


## 3. Temperature

Module `agrometeorologiapy.en.temperature` (Portuguese API:
`agrometeorologiapy.temperatura`).

### `mean_temperature(Tmax, Tmin)` · pt: `temp_media_extremos`
Mean daily air temperature (FAO-56 eq. 9):

```math
T_{mean} = \frac{T_{max} + T_{min}}{2}
```

**Where:**
- $T_{mean}$ (return) — estimated daily mean temperature, in °C
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C

### `mean_temperature_observations(temperatures)` · pt: `temp_media_estacao_automatica`
Simple arithmetic mean of $n$ observations over the period:

```math
T_{mean} = \frac{1}{n}\sum_{i=1}^{n} T_i
```

**Where:**
- $T_{mean}$ (return) — mean temperature for the period, in °C
- $T_i$ (`temperatures`) — each temperature observation, in °C
- $n$ — number of observations in the list

---

## 4. Air Humidity

Module `agrometeorologiapy.en.humidity` (Portuguese API:
`agrometeorologiapy.umidade`).

### `saturation_vapour_pressure(T)` · pt: `es_tetens`
Tetens equation (FAO-56 eq. 11):

```math
e_s = 0{.}6108 \cdot 10^{\frac{7{.}5 \, T}{237{.}3 + T}} \quad [\text{kPa}]
```

**Where:**
- $e_s$ (return) — saturation vapour pressure, in kPa
- $T$ (`T`) — air temperature, in °C

### `actual_vapour_pressure(es, RH)` · pt: `ea_umidade`

```math
e_a = e_s \cdot \frac{RH}{100}
```

**Where:**
- $e_a$ (return) — actual vapour pressure, in kPa
- $e_s$ (`es`) — saturation vapour pressure, in kPa
- $RH$ (`RH`) — relative humidity, in % (e.g. 65)

### `vapour_pressure_deficit(es, ea)` · pt: `deficit_saturacao`

```math
\Delta e = e_s - e_a
```

**Where:**
- $\Delta e$ (return) — vapour pressure deficit, in kPa
- $e_s$ (`es`) — saturation vapour pressure, in kPa
- $e_a$ (`ea`) — actual vapour pressure, in kPa

### `atmospheric_pressure(z)` · pt: `patm_altitude`
Atmospheric pressure as a function of elevation (FAO-56 eq. 7):

```math
P_{atm} = 101{.}3 \cdot \left( \frac{293 - 0{.}0065 z}{293} \right)^{5{.}26} \quad [\text{kPa}]
```

**Where:**
- $P_{atm}$ (return) — atmospheric pressure, in kPa
- $z$ (`z`) — elevation above sea level, in m

### `absolute_humidity(ea, T)` · pt: `umidade_absoluta`

```math
AH = \frac{2168 \cdot e_a}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3]
```

**Where:**
- $AH$ (return) — absolute humidity, in g of H₂O per m³ of air
- $e_a$ (`ea`) — actual vapour pressure, in kPa
- $T$ (`T`) — air temperature, in °C (converted internally to kelvin,
  $T + 273{.}15$)

### `saturation_humidity(es, T)` · pt: `umidade_saturacao`

```math
SH = \frac{2168 \cdot e_s}{T + 273{.}15} \quad [\text{g of H}_2\text{O / m}^3]
```

**Where:**
- $SH$ (return) — saturation humidity, in g of H₂O per m³ of air
- $e_s$ (`es`) — saturation vapour pressure, in kPa
- $T$ (`T`) — air temperature, in °C

### `relative_humidity(ea, es)` · pt: `umidade_relativa`
FAO-56 eq. 10:

```math
RH = 100 \cdot \frac{e_a}{e_s} \quad [\%]
```

**Where:**
- $RH$ (return) — relative humidity, in %
- $e_a$ (`ea`) — actual vapour pressure, in kPa
- $e_s$ (`es`) — saturation vapour pressure, in kPa

### `dewpoint_temperature(ea)` · pt: `ponto_orvalho`
Algebraic inversion of the Tetens equation:

```math
x = \log_{10}\left(\frac{e_a}{0{.}6108}\right)
```

```math
T_{dew} = \frac{237{.}3 \, x}{7{.}5 - x} \quad [°C]
```

**Where:**
- $T_{dew}$ (return) — dewpoint temperature, in °C
- $e_a$ (`ea`) — actual vapour pressure, in kPa
- $x$ — auxiliary variable (logarithm of the ratio $e_a/0{.}6108$),
  dimensionless

### `psychrometric_constant(Patm)` · pt: `constante_psicrometrica`
FAO-56 eq. 8:

```math
\gamma = 0{.}665 \times 10^{-3} \cdot P_{atm} \quad [\text{kPa}/°C]
```

**Where:**
- $\gamma$ (return) — psychrometric constant, in kPa/°C
- $P_{atm}$ (`Patm`) — atmospheric pressure, in kPa

---

## 5. Energy Balance

Module `agrometeorologiapy.en.energy` (Portuguese API:
`agrometeorologiapy.energia`).

> ℹ️ **Sign convention:** the package returns the net longwave radiation as
> a negative value (net loss) and computes $R_n = R_{ns} + R_{nl}$. FAO-56
> writes $R_{nl}$ as a positive outgoing flux and $R_n = R_{ns} - R_{nl}$;
> the net radiation is the same.

### `net_shortwave_radiation(Rs, albedo=0.25)` · pt: `boc_saldo`
Incoming solar radiation minus the reflected fraction (FAO-56 eq. 38):

```math
R_{ns} = R_s (1 - \alpha)
```

**Where:**
- $R_{ns}$ (return) — net shortwave radiation, in the same unit as $R_s$
- $R_s$ (`Rs`) — solar radiation incident on the surface, in
  MJ m⁻² day⁻¹ (or W/m²)
- $\alpha$ (`albedo`) — surface albedo, dimensionless (0–1; default 0.25,
  typical for grass; 0.23 for the FAO-56 reference crop)

### `net_longwave_radiation(Tmax, Tmin, ea, Rs, Rso)` · pt: `bol_saldo`
Stefan-Boltzmann law corrected for humidity and cloudiness (FAO-56 eq. 39).
`Tmax`/`Tmin` are received in °C (consistent with the rest of the package)
and converted internally to kelvin, because the Stefan-Boltzmann term
requires absolute temperature:

```math
T_{max,K} = T_{max} + 273{.}15, \qquad T_{min,K} = T_{min} + 273{.}15
```

```math
R_{nl} = -\left[ 4{.}903 \times 10^{-9} \cdot \frac{T_{max,K}^4 + T_{min,K}^4}{2} \right] \cdot \left[ 0{.}34 - 0{.}14\sqrt{e_a} \right] \cdot \left[ 1{.}35 \frac{R_s}{R_{so}} - 0{.}35 \right]
```

**Where:**
- $R_{nl}$ (return) — net longwave radiation, in MJ m⁻² day⁻¹ (negative:
  net energy loss through terrestrial emission)
- $T_{max}$, $T_{min}$ (`Tmax`, `Tmin`) — daily maximum and minimum
  temperatures, in °C; $T_{max,K}$, $T_{min,K}$ are the same in kelvin
- $e_a$ (`ea`) — actual vapour pressure, in kPa
- $R_s$ (`Rs`) — measured/estimated solar radiation for the day, in
  MJ m⁻² day⁻¹
- $R_{so}$ (`Rso`) — clear-sky solar radiation, in MJ m⁻² day⁻¹
- $4{.}903 \times 10^{-9}$ — Stefan-Boltzmann constant, in MJ K⁻⁴ m⁻² day⁻¹

> 🐛 **Note:** the original source notebook applied this term directly to
> Tmax/Tmin in °C (without converting to kelvin), which underestimated it
> by several orders of magnitude. This was fixed in this library — see the
> [CHANGELOG](../CHANGELOG.md).

### `net_radiation(Rns, Rnl)` · pt: `saldo_radiacao`
FAO-56 eq. 40, with the sign convention above:

```math
R_n = R_{ns} + R_{nl} \quad [\text{MJ m}^{-2}\text{day}^{-1}]
```

**Where:**
- $R_n$ (return) — net radiation, in MJ m⁻² day⁻¹ (or W/m²)
- $R_{ns}$ (`Rns`) — net shortwave radiation, same unit as $R_n$
- $R_{nl}$ (`Rnl`) — net longwave radiation (negative), same unit as $R_n$

---

## 6. Evapotranspiration

Module `agrometeorologiapy.en.evapotranspiration` (Portuguese API:
`agrometeorologiapy.evapotranspiracao`).

### `pet_thornthwaite(df, col_T='T_mean_C', lat=None)` · pt: `thornthwaite_mensal`
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
PET_{u} = 16 \left( \frac{10 T}{I} \right)^{a} \quad \text{(unadjusted, mm/month of 30 days and 12h of sunlight)}
```

For $T \ge 26{.}5$ °C the original method uses a table instead of the equation above:

```math
PET_{u} = -415{.}85 + 32{.}24\,T - 0{.}43\,T^2
```

Adjusted by the month's daylight hours (via the solar declination for the
month's mean day of the year) and by the actual number of days in the
month:

```math
K = \frac{N}{12} \cdot \frac{\text{days}_{\text{month}}}{30}
```

```math
PET = PET_{u} \cdot K \quad [\text{mm/month}]
```

**Where:**
- $PET$ (return, column `PET_mm_month`) — adjusted potential
  evapotranspiration, in mm/month
- $T$ (`col_T` column of `df`) — mean temperature for each month, in °C
- $i$ — monthly heat index, dimensionless; $I$ — sum of the 12 monthly
  indices (annual heat index), dimensionless
- $a$ — empirical exponent, a function of $I$, dimensionless
- $PET_{u}$ — unadjusted potential evapotranspiration, in mm/month
- $N$ — daylight hours for the month (computed internally from `lat` and
  the solar declination for the month's mean day of the year), in hours
- $\text{days}_{\text{month}}$ — actual number of days in the month (28–31)

### `pet_camargo(df, col_T='T_mean_C', lat=None)` · pt: `camargo_maluf_mensal`
Camargo (1971) method, with the F coefficient modified by Maluf (Camargo
et al., 1999). Self-contained function: it recomputes the monthly
extraterrestrial radiation internally ($R_a$, from the solar declination
for the month's mean day of the year) and converts it to evaporation
equivalent:

```math
R_{a,mm} = 0{.}408 \cdot R_{a}
```

The F coefficient depends on the **annual** mean temperature
($\bar{T}_{annual}$):

- $F = 0{.}0100$ if $\bar{T}_{annual} < 23°C$
- $F = 0{.}0105$ if $23°C \le \bar{T}_{annual} < 24°C$
- $F = 0{.}0110$ if $\bar{T}_{annual} \ge 24°C$

```math
PET = F \cdot R_{a,mm} \cdot T \cdot \text{days}_{\text{month}} \quad [\text{mm/month}]
```

**Where:**
- $PET$ (return, column `PET_mm_month`) — monthly potential
  evapotranspiration, in mm/month
- $R_a$ (column `Ra_MJ_m2day`) — monthly extraterrestrial radiation, in
  MJ m⁻² day⁻¹; $R_{a,mm}$ (column `Ra_mm_day`) — the same, in evaporation
  equivalent, mm/day
- $F$ — empirical coefficient, dimensionless, chosen from the annual mean
  temperature
- $\bar{T}_{annual}$ — annual mean temperature (mean of the 12 input
  monthly temperatures), in °C
- $T$ (`col_T` column of `df`) — mean temperature for each month, in °C
- $\text{days}_{\text{month}}$ — actual number of days in the month (28–31)

### `pet_hargreaves_samani(Ra, Tmax, Tmin, Tmean)` · pt: `etp_hargreaves_samani`
Hargreaves & Samani (1985) (FAO-56 eq. 52):

```math
R_{a,mm} = 0{.}408 \cdot R_a
```

```math
PET = 0{.}0023 \cdot R_{a,mm} \cdot \sqrt{T_{max} - T_{min}} \cdot (T_{mean} + 17{.}8) \quad [\text{mm/day}]
```

**Where:**
- $PET$ (return) — potential evapotranspiration, in mm/day
- $R_a$ (`Ra`) — extraterrestrial radiation, in MJ m⁻² day⁻¹; $R_{a,mm}$ —
  the same, in evaporation equivalent, mm/day
- $T_{max}$, $T_{min}$, $T_{mean}$ (`Tmax`, `Tmin`, `Tmean`) — maximum,
  minimum and mean air temperatures, in °C

### `slope_vapour_pressure_curve(T)` · pt: `declive_pressao_vapor`
Derivative of the Tetens equation with respect to temperature (FAO-56
eq. 13):

```math
\Delta = \frac{4098 \cdot e_s(T)}{(T + 237{.}3)^2} \quad [\text{kPa}/°C]
```

**Where:**
- $\Delta$ (return) — slope of the saturation vapour pressure curve, in
  kPa/°C
- $T$ (`T`) — air temperature, in °C (typically the daily mean
  temperature)
- $e_s(T)$ — saturation vapour pressure at temperature $T$
  (`saturation_vapour_pressure`), in kPa

### `pet_priestley_taylor(Rn, G, Delta, gamma, alpha=1.26)` · pt: `etp_priestley_taylor`

```math
PET = \alpha \cdot \frac{\Delta}{\Delta + \gamma} \cdot \frac{R_n - G}{\lambda} \quad [\text{mm/day}]
```

**Where:**
- $PET$ (return) — potential evapotranspiration, in mm/day
- $R_n$ (`Rn`) — net radiation, in MJ m⁻² day⁻¹
- $G$ — soil heat flux, in MJ m⁻² day⁻¹ (usually $G = 0$ on a daily scale)
- $\Delta$ (`Delta`) — slope of the saturation vapour pressure curve, in
  kPa/°C
- $\gamma$ (`gamma`) — psychrometric constant, in kPa/°C
- $\lambda = 2{.}45$ — latent heat of vaporization, in MJ/kg (fixed constant
  in the code)
- $\alpha$ (`alpha`) — Priestley-Taylor coefficient, dimensionless (default
  1.26)

### `wind_speed_2m(uz, z)` · pt: `vento_2m`
Converts wind speed measured at any height to 2 m, using the logarithmic
wind profile (FAO-56 eq. 47):

```math
u_2 = u_z \cdot \frac{4{.}87}{\ln(67{.}8 \, z - 5{.}42)} \quad [\text{m/s}]
```

**Where:**
- $u_2$ (return) — wind speed at 2 m height, in m/s
- $u_z$ (`uz`) — wind speed measured at height $z$, in m/s
- $z$ (`z`) — measurement height above the ground, in m

### `eto_penman_monteith(Rn, G, Tmean, u2, es, ea, Delta, gamma)` · pt: `eto_penman_monteith_fao56`
FAO Penman-Monteith equation (FAO-56 eq. 6; Allen et al., 1998), for the
hypothetical reference crop (grass, 0.12 m, albedo 0.23):

```math
ET_o = \frac{0{.}408 \, \Delta (R_n - G) + \gamma \cdot \frac{900}{T_{mean}+273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma (1 + 0{.}34 \, u_2)} \quad [\text{mm/day}]
```

**Where:**
- $ET_o$ (return) — reference evapotranspiration, in mm/day
- $R_n$ (`Rn`) — net radiation, in MJ m⁻² day⁻¹
- $G$ — soil heat flux, in MJ m⁻² day⁻¹ ($G = 0$ on a daily scale)
- $T_{mean}$ (`Tmean`) — daily mean air temperature, in °C
- $u_2$ (`u2`) — wind speed at 2 m height, in m/s
- $e_s$ (`es`) — saturation vapour pressure, in kPa
- $e_a$ (`ea`) — actual vapour pressure, in kPa
- $\Delta$ (`Delta`) — slope of the saturation vapour pressure curve, in
  kPa/°C
- $\gamma$ (`gamma`) — psychrometric constant, in kPa/°C

---

## 7. Growing Degree Days

Module `agrometeorologiapy.en.degree_days` (Portuguese API:
`agrometeorologiapy.grau_dias`).

### Daily growing degree days (GDD)
Used by both `maturity_date` and `sowing_date`:

- $GDD = T_{mean} - T_b$, if $T_b < T_{min}$
- $GDD = \dfrac{(T_{max} - T_b)^2}{2(T_{max} - T_{min})}$, if $T_b \ge T_{min}$

**Where:**
- $GDD$ — daily growing degree days for the period (day, dekad or month),
  in °C·day
- $T_{mean}$, $T_{max}$, $T_{min}$ — mean, maximum and minimum temperatures
  for the period, in °C
- $T_b$ — crop base temperature (below which there is no development), in
  °C

### `maturity_date(df, Tb, TT, sowing_day, sowing_month, interval='d', year=2023, verbose=True)` · pt: `data_maturacao_fisiologica`
Starting from the sowing date, accumulates $GDD \times n_{period}$ period
by period (daily, dekadal or monthly) until the sum reaches the thermal
time requirement of the cycle:

```math
\sum GDD \cdot n_{period} \ge TT
```

Returns the date on which that happens — physiological maturity. In the
last period only the days needed to complete `TT` are counted
($\lceil (TT - GDD_{previous}) / GDD \rceil$), so the last row of the
DataFrame is the maturity date itself. The year rolls over when the cycle
crosses December → January. The result has one row per period with the
columns `month` (abbreviated), `date`, `Tmean`, `GDD` (daily degree days),
`GDD_period` ($GDD \times$ days counted in the period) and `GDD_cycle`
(cycle accumulation). With `verbose=True` (default), the sowing and
maturity dates are printed.

**Where:**
- `df` — climate series (columns `day`, `month`, `Tmean`, `Tmax`, `Tmin`),
  in chronological order
- `Tb` — crop base temperature, in °C
- `TT` — thermal time requirement of the cycle, in °C·day
- `sowing_day`, `sowing_month` — sowing date, integers
- `interval` — `'d'` (daily), `'dekad'` (dekadal, 10 days) or `'M'`
  (monthly)
- $n_{period}$ — number of days in the period (1 for daily, 10 for
  dekadal, days in the month for monthly)
- `year` — reference year, integer (default 2023, non-leap)

### `sowing_date(df, Tb, TT, maturity_day, maturity_month, interval='d', year=2023, verbose=True)` · pt: `data_semeadura`
The same growing degree day accumulation, but walking the calendar
**backwards** from a known maturity date (e.g. a target harvest date),
until `TT` is accumulated — returning the required sowing date.

**Where:** same parameters as `maturity_date`, swapping
`sowing_day`/`sowing_month` (known starting date) for
`maturity_day`/`maturity_month` (reference date from which the calendar is
walked backwards).

---

## 8. Water Balance

Module `agrometeorologiapy.en.water_balance` (Portuguese API:
`agrometeorologiapy.balanco_hidrico`). The functions implement the
sequential accounting method of Thornthwaite & Mather (1955).

### `climatic_water_balance(df, AWC=100.0, cyclic=True)` · pt: `balanco_hidrico_climatologico`
For each period (month) $i$:

```math
P - PET
```

- If $P - PET < 0$ (deficit): accumulate the potential water loss and
  recompute soil water storage through the exponential model —

  $`\displaystyle APWL_i = APWL_{i-1} + (P - PET)`$

  $`\displaystyle STO_i = AWC \cdot e^{\,APWL_i / AWC}`$

- If $P - PET \ge 0$ (recharge): the soil receives water up to a maximum
  of `AWC` —

  $`\displaystyle STO_i = \min(STO_{i-1} + (P - PET),\; AWC)`$

  and, if $STO_i < AWC$, APWL is recomputed by inversion:

  $`\displaystyle APWL_i = AWC \cdot \ln(STO_i / AWC)`$

From soil water storage, the following are derived:

```math
\Delta STO_i = STO_i - STO_{i-1}
```

- $ETa_i = P_i + |\Delta STO_i|$, if $P_i - PET_i < 0$
- $ETa_i = PET_i$, otherwise

```math
DEF_i = PET_i - ETa_i
```

- $SUR_i = (P_i - PET_i) - \Delta STO_i$, if $P_i - PET_i > 0$ and $STO_i = AWC$
- $SUR_i = 0$, otherwise

**Where:**
- `df` — columns `P (mm/month)` (precipitation) and `PET (mm/month)`
  (potential evapotranspiration); a `month` column, if present, is kept
- `AWC` — available water capacity of the soil, in mm (default 100.0)
- $P_i$, $PET_i$ — precipitation and potential evapotranspiration for
  period $i$, in mm/month
- $STO_i$ (column `STO (mm/month)`) — soil water storage at the end of
  period $i$, in mm; $i-1$ is the previous period. With `cyclic=True`
  (default), the storage before January is December's storage at the
  annual-cycle equilibrium (the 12 months are repeated until December's STO
  converges); with `cyclic=False`, it starts from $STO = AWC$ (full soil)
- $APWL_i$ (column `APWL (mm)`) — accumulated potential water loss, an
  auxiliary variable in mm, used in the exponential drying model
- $\Delta STO_i$ (column `dSTO (mm/month)`) — change in storage between
  periods, in mm
- $ETa_i$ (column `ETa (mm/month)`) — actual evapotranspiration, in mm/month
- $DEF_i$ (column `DEF (mm/month)`) — water deficit, in mm/month
- $SUR_i$ (column `SUR (mm/month)`) — water surplus, in mm/month

### `climatic_water_balance_grid(P, PET, AWC=100.0, cyclic=True, tol=0.01, max_iter=100)` · pt: `balanco_hidrico_climatologico_grade`
Same equations as `climatic_water_balance`, vectorized with numpy for many
locations at once (e.g. every pixel of a raster with the monthly climate
normal).

**Where:**
- `P`, `PET` — arrays of shape `(12, ...)`, in mm/month, January to
  December along the first axis
- `AWC` — scalar or array with the shape of the spatial dimensions, in mm
  (> 0; use NaN for locations without data)
- `cyclic` — equilibrium initial storage of the annual cycle (default) or
  full soil before January
- `tol`, `max_iter` — tolerance (mm) and maximum number of annual cycles
  in the search for equilibrium
- return — `dict` of `(12, ...)` arrays: `P-PET`, `STO`, `APWL`, `dSTO`,
  `ETa`, `DEF`, `SUR`

### `crop_water_balance(df)` · pt: `balanco_hidrico_cultura`
The same method, applied to a specific crop with `P`, `ETc` and `AWC`
pre-computed by the user for each period (same formulas for `STO`, `dSTO`,
`ETa`, `DEF` and `SUR` as above, with $PET \to$ `ETc`), plus the water
requirement satisfaction index (Frère & Popov, 1979):

```math
WRSI = \frac{ETa}{ETc}
```

**Where:**
- `df` — columns `P` (precipitation, mm/period), `ETc` (crop
  evapotranspiration, mm/period) and `AWC` (available water capacity for
  the period, mm), in chronological order
- `ETc` — crop evapotranspiration, mm/period; pre-computed by the user as
  $K_c \times ET_o$ (crop coefficient × reference evapotranspiration),
  already accounting for the phenological stage
- `AWC` — available water capacity for the period, mm; pre-computed by the
  user as root depth × available water of the soil, already accounting for
  root growth
- $WRSI$ (return, column `WRSI`) — water requirement satisfaction index,
  dimensionless (0–1: the closer to 1, the lower the crop's water stress)

---

## 9. Climate Classification

Module `agrometeorologiapy.en.climate_classification` (Portuguese API:
`agrometeorologiapy.classificacao_climatica`). Each method has two
versions: one for **a single site** (lists of 12 months, Jan–Dec) and a
`_grid` one for **many sites** (e.g. every pixel of a raster), with the 12
months on the first axis, shape `(12, ...)`. Sites with NaN in any month get
code 0 and class `''`. The class codes (`Cfa`, `B1rA'a'`, `ST-UMi`...) are
the same in both APIs; names and descriptions are in English here and in
Portuguese in the Portuguese API.

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

### `koppen_classification(T, P, lat)` · `koppen_classification_grid(T, P, lat)` · pt: `classificacao_koppen`
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
- return — `id` (1–31) and `class` (e.g. `'Cfa'`)

### `thornthwaite_classification(P, PET, lat, AWC=100)` · `thornthwaite_classification_grid(...)` · pt: `classificacao_thornthwaite`
Thornthwaite (1948), from the climatic water balance
(`climatic_water_balance_grid`):

```math
I_h = 100\,\frac{SUR}{PET} \qquad I_a = 100\,\frac{DEF}{PET} \qquad I_m = I_h - 0{.}6\,I_a
```

```math
PET_{summer}\,(\%) = 100\,\frac{PET_{summer}}{PET_{annual}}
```

Class = moisture ($I_m$) + subtype + thermal efficiency ($PET_{annual}$) +
summer concentration ($PET_{summer}$, %), e.g. `B1rA'a'`.

| $I_m$ | Moisture | | $PET_{annual}$ (mm) | Thermal |
|---|---|---|---|---|
| ≥ 100 | A (perhumid) | | ≥ 1140 | A' (megathermal) |
| 80–100 / 60–80 / 40–60 / 20–40 | B4 / B3 / B2 / B1 (humid) | | 997–1140 / 855–997 / 712–855 / 570–712 | B'4 / B'3 / B'2 / B'1 (mesothermal) |
| 0–20 | C2 (moist subhumid) | | 427–570 / 285–427 | C'2 / C'1 (microthermal) |
| −20–0 | C1 (dry subhumid) | | 142–285 | D' (tundra) |
| −40– −20 | D (semiarid) | | < 142 | E' (frost) |
| < −40 | E (arid) | | | |

- **Subtype, humid climates** ($I_m \ge 0$), by deficit: r ($I_a < 16{.}7$);
  s/w ($16{.}7 \le I_a < 33{.}3$); s2/w2 ($I_a \ge 33{.}3$) — s if the summer
  DEF exceeds the winter DEF, otherwise w
- **Subtype, dry climates** ($I_m < 0$), by surplus: d ($I_h < 10$); s/w
  ($10 \le I_h < 20$); s2/w2 ($I_h \ge 20$) — w if the summer SUR exceeds the
  winter SUR, otherwise s
- **Summer concentration** ($PET_{summer}$, %): a' < 48; b'4 < 51.9;
  b'3 < 56.3; b'2 < 61.6; b'1 < 68; c'2 < 76.3; c'1 < 88; d' ≥ 88

> Tables from Aparecido et al. (2016), corrected against the original
> Thornthwaite (1948) at three points where the paper has typos: B'3/B'2
> limit = 855 mm (paper: 885); dry s2/w2 with $I_h \ge 20$ (paper: 33.3); dry
> climates with larger summer surplus = w (paper: s).

**Where:**
- `P`, `PET` — monthly precipitation and potential evapotranspiration, in mm/month
- `lat` — latitude, in degrees (sets summer and winter)
- `AWC` — available water capacity of the soil, in mm (default 100)
- return — `class`, codes (`moisture`, `subtype`, `thermal`,
  `concentration`), `Ih`, `Ia`, `Im`, `PET_annual`, `DEF_annual`,
  `SUR_annual`, `PET_summer_pct` (plus `description` in the single-site
  version, e.g. `'Humid, mesothermal'`)

### `camargo_classification(T, P, PET, lat, AWC=100, thermal_table='consistent')` · `camargo_classification_grid(...)` · pt: `classificacao_camargo`
Camargo (1991) as modified by Maluf (2000), Tables 6–8 of Aparecido et al.
(2016). Class = thermal + `-` + moisture + dry-season letter, e.g. `ST-UMi`.

| $T_{ann}$ (°C) | $T_{cold}$ (°C) | Thermal |
|---|---|---|
| ≤ 3 | | GL (glacial) |
| 3–7 | | FR (cold) |
| 7–12 | | CO (moderately cold) |
| 12–18 | | TE (temperate) |
| 18–22 | ≤ 13 | STE (subtemperate) |
| 18–22 | 13–20 | ST (subtropical) |
| 18–22 | > 20 | TR (tropical) |
| 22–25 | | TR (tropical) |
| > 25 | | EQ (equatorial) |

| Annual DEF (mm) | Annual SUR (mm) | Moisture |
|---|---|---|
| > 800 | 0 | DE (desert) |
| 150–800 | 0 | AR (arid) |
| > 150 | 0–200 | SE (dry) |
| > 150 | > 200 | MO (monsoonal) |
| 0–150 | 0–200 | SB (subhumid) |
| 0–150 (> 0) | > 200 | UM (humid) |
| 0 | 200–1000 | PU (perhumid) |
| 0 | > 1000 | SU (extremely humid) |

For classes SE, MO, SB and UM, the letter of the season with the largest
deficit is appended — Portuguese initials: v (*verão*, summer), o (*outono*,
autumn), i (*inverno*, winter) or p (*primavera*, spring). The codes are
those of the original system.

> **Reading Table 6.** The paper has "22 < $T_{ann}$ ≤ 25 **or**
> $T_{cold}$ > 20 → TR", which, read literally, rules out EQ wherever the
> coldest month exceeds 20 °C (nearly the whole equatorial region). The
> default `thermal_table='consistent'` uses $T_{cold}$ only to split the
> 18–22 °C range (table above); `'literal'` applies the rule $T_{cold} > 20 \rightarrow TR$ in
> every case, as in the `climas_brasil` repository script.
>
> **Cases outside Table 7.** The paper does not cover 0 < DEF ≤ 150 mm with
> SUR = 0, nor DEF = 0 with SUR ≤ 200 mm; here they fall into SB. DEF and SUR
> below 0.05 mm count as zero.

**Where:**
- `T` — monthly mean temperature, in °C
- `P`, `PET` — monthly precipitation and potential evapotranspiration, in mm/month
- `lat` — latitude, in degrees
- `AWC` — available water capacity of the soil, in mm (default 100)
- `thermal_table` — `'consistent'` (default) or `'literal'`
- return — `class`, codes (`thermal`, `moisture`, `dry_season`),
  `T_annual`, `T_coldest`, `DEF_annual`, `SUR_annual` (plus `description` in
  the single-site version, e.g. `'Subtropical perhumid'`)

### `holdridge_classification(T, P, lat, PET=None, correction_threshold=24)` · `holdridge_classification_grid(...)` · pt: `classificacao_holdridge`
Holdridge life zones (38 zones; numbering and names from Jungkunst et al.,
2021, based on Leemans, 1990). Biotemperature:

```math
t^*_m = t_m - \frac{3\,|\phi|}{100}\,(t_m - 24)^2 \;\;(\text{if } t_m > 24\ °C), \qquad BT = \overline{\min(\max(t^*_m, 0), 30)}
```

```math
PET_{annual} = 58{.}93 \cdot BT \qquad R = \frac{PET_{annual}}{P_{annual}}
```

The zone comes from combining the biotemperature belt (polar < 1.5;
subpolar < 3; boreal < 6; cool temperate < 12; warm temperate < 18;
subtropical < 24; tropical ≥ 24 °C) with the humidity province, in ranges of
$R$ that double at each class (0.125; 0.25; 0.5; 1; 2; 4; 8; 16; 32).

**Where:**
- `T` — monthly mean temperature, in °C
- `P` — monthly precipitation, in mm/month
- `lat` — latitude ($\phi$), in degrees
- `PET` — monthly PET, in mm/month (optional; if `None`, uses $58{.}93 \cdot BT$)
- `correction_threshold` — monthly temperature above which the latitude
  correction applies, in °C (default 24; `None` corrects every month)
- return — `id` (1–38), `class` (zone name, as in the source),
  `biotemperature` (°C), `P_annual`, `PET_annual` (mm) and `PET_ratio`.
  The Portuguese API returns `classe` (name in Portuguese) and `classe_en`

**References:** Alvares, C. A. et al. (2013) *Meteorol. Z.* 22:711–728 ·
Kottek, M. et al. (2006) *Meteorol. Z.* 15:259–263 · Thornthwaite, C. W. (1948)
*Geogr. Rev.* 38:55–94 · Camargo, A. P. (1991) · Maluf, J. R. T. (2000)
*Rev. Bras. Agrometeorol.* 8:141–150 · Aparecido, L. E. O. et al. (2016)
*Ciênc. Agrotec.* 40:405–417 · Holdridge, L. R. (1967) *Life zone ecology* ·
Jungkunst, H. F. et al. (2021) *J. Plant Nutr. Soil Sci.* 184:5–11

---

## 10. Potential and Attainable Yield

Module `agrometeorologiapy.en.crop_yield` (Portuguese API:
`agrometeorologiapy.produtividade`). FAO Agro-Ecological Zone model
(Doorenbos & Kassam, 1979), with the temperature corrections of Barbieri &
Tuon (1992), in Pereira et al. (2002). Constants and tabulated values default
to the usual values, but all of them are function parameters. Names follow
the yield gap literature (van Ittersum et al., 2013): $Y_p$ is the potential
yield and $Y_w$ the attainable (water-limited) yield — not $Y_a$, which there
denotes the actual (observed) yield.

### n/N
Without a sunshine recorder, the relative sunshine duration comes from
`relative_sunshine_duration` (section 2), the inverted Angstrom formula.
`potential_yield` uses it when the df has `Rs` instead of `nN`.

### `temperature_correction_overcast(T, pathway, T_threshold=16.5)` · `temperature_correction_clear(...)` · pt: `cTn`, `cTc`
Temperature corrections of photosynthesis for overcast ($cT_o$) and clear
($cT_c$) sky, clipped at zero.

**C3** (soybean, wheat, bean, sunflower):

```math
cT_o = -0.0425 + 0.035\,T + 0.00325\,T^2 - 0.0000925\,T^3
```

```math
cT_c = 0.583 + 0.014\,T + 0.0013\,T^2 - 0.000037\,T^3
```

**C4** (maize, sorghum, sugarcane):

| Condition | $cT_o$ | $cT_c$ |
|---|---|---|
| $T \ge 16.5$ °C | $-1.064 + 0.173 T - 0.0029 T^2$ | $-4.16 + 0.4325 T - 0.00725 T^2$ |
| $T < 16.5$ °C | $-4.16 + 0.4325 T - 0.00725 T^2$ | $-9.32 + 0.865 T - 0.0145 T^2$ |

**Where:**
- `T` — mean air temperature of the period, °C
- `pathway` — photosynthetic pathway, `'C3'` or `'C4'` (required)
- `T_threshold` — polynomial switch for C4 crops, °C (default 16.5)

### `gross_potential_productivity(Ra, nN, cTc, cTo, a_c=107.2, b_c=8.604, a_o=31.7, b_o=5.234)` · pt: `PPBp`
Gross potential productivity of the standard crop (LAI = 5):

```math
GP_c = (107.2 + 8.604\,R_a)\, cT_c\, \frac{n}{N} \qquad GP_o = (31.7 + 5.234\,R_a)\, cT_o \left(1 - \frac{n}{N}\right)
```

```math
GP = GP_c + GP_o
```

**Where:**
- `Ra` — extraterrestrial radiation, MJ m⁻² day⁻¹
- `nN` — relative sunshine duration
- $GP$ (return) — kg DM ha⁻¹ day⁻¹
- the coefficients are Doorenbos & Kassam's (0.36 and 0.219, with $R_a$ in
  cal cm⁻² day⁻¹) converted to MJ m⁻² day⁻¹

### `lai_correction(LAI, LAI_standard=5.0)` · `respiration_correction(T, T_threshold=20.0, CR_cool=0.6, CR_warm=0.5)` · pt: `CIAF`, `CR`
Leaf area index correction and maintenance respiration correction:

```math
LAI_{corr} = 0.0093 + 0.185\,LAI - 0.0175\,LAI^2 \;\;(LAI < 5), \qquad LAI_{corr} = 0.5 \;\;(LAI \ge 5)
```

```math
CR = 0.6 \;\;(T < 20\ °C), \qquad CR = 0.5 \;\;(T \ge 20\ °C)
```

### `potential_yield(df, lat, pathway, LAI=5.0, HI=0.35, moisture=13.0, interval='dekad', sowing=None, cycle=None, ...)` · pt: `produtividade_potencial`
Accumulates the potential yield period by period over the crop cycle:

```math
Y_p = \frac{\sum_{p} GP_p \cdot ND_p \cdot LAI_{corr} \cdot CR_p \cdot HI}{1 - 0.01\,moisture}
```

**Where:**
- `df` — one row per period: `Tmean` (°C); radiation as `nN`, as `n` and
  `N`, or as `Rs` (MJ m⁻² day⁻¹); `Ra` (optional, computed from `lat` and the
  dates); dates in `date`, a DatetimeIndex or `day` and `month` (+ `year`);
  `ND` (optional)
- `interval` — time step: `'d'` (daily), `'dekad'` (calendar 10-day period:
  1–10, 11–20 and 21–end of month, 8 to 11 days) or `'M'` (monthly)
- $ND_p$ — days of the period with the crop in the field. With `sowing`
  (sowing date) and `cycle` (days, counting the sowing day), the first and
  last periods are partial and periods outside the cycle are dropped
- `LAI`, `HI`, `moisture` — maximum LAI, harvest index and moisture content
  of the product, %
- return — DataFrame with `ND`, `Ra`, `nN`, `cTo`, `cTc`, `GPc`, `GPo`,
  `GP`, `CR`, `Yp_period` and `Yp_cum` (kg ha⁻¹); `attrs['Yp']`

### `attainable_yield(Yp, ETa, ETc, ky, phase=None, method='product')` · pt: `produtividade_atingivel`
Water-deficit penalty:

```math
1 - \frac{Y_w}{Y_p} = k_y \left(1 - \frac{ETa}{ETc}\right)
```

- **Product** (`'product'`, Battisti et al., 2013), phase by phase in sequence:

  $`\displaystyle Y_w = Y_p \cdot \prod_i \left[1 - k_{y,i}\left(1 - \frac{ETa_i}{ETc_i}\right)\right]`$

- **Single step** (`'single'`), with the whole-cycle $k_y$:

  $`\displaystyle Y_w = Y_p \cdot \left[1 - k_{y,cycle}\left(1 - \frac{\sum ETa}{\sum ETc}\right)\right]`$

```math
CE = \frac{Y_w}{Y_p}
```

**Where:**
- `ETa`, `ETc` — mm, per phase or per period (e.g. output of
  `crop_water_balance`); with `phase`, they are summed per phase
- `ky` — one per phase (sequence or dict `{phase: ky}`) or the whole-cycle value
- each factor is clipped at zero; phases with $ETc = 0$ are not penalized
- return — one row per phase with `ETc`, `ETa`, `ETa/ETc`, `ky`, `factor`,
  `Yw_start`, `Yw_end`, `loss` and `CE` (climatic efficiency);
  `attrs['Yw']` and `attrs['CE']`

**References:** Doorenbos, J.; Kassam, A. H. (1979) *Yield response to water*,
FAO Irrigation and Drainage Paper 33 · Barbieri, V.; Tuon, R. L. (1992) ESALQ/USP ·
Pereira, A. R.; Angelocci, L. R.; Sentelhas, P. C. (2002) *Agrometeorologia* ·
Battisti, R. et al. (2013) *Ciência Rural* 43:390–396 · van Ittersum, M. K.
et al. (2013) *Field Crops Res.* 143:4–17
