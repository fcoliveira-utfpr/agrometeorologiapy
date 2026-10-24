# Nomes da API: português ↔ inglês

*Function names: Portuguese ↔ English*

O `agrometeorologiapy` tem **uma única implementação**, em inglês
(`agrometeorologiapy.en`), e duas APIs com o mesmo comportamento:

```python
import agrometeorologiapy as amp          # português
from agrometeorologiapy import en as agm  # English

amp.relacao_n_N(Qg=22.5, Qo=42.96, lat=-24.86)            # 0.501
agm.relative_sunshine_duration(Rs=22.5, Ra=42.96, lat=-24.86)  # 0.501
```

A API em português é permanente (não será descontinuada). Cada função em
português chama a função em inglês correspondente, traduzindo os parâmetros,
as colunas de DataFrame e as opções em texto.

## Convenções

- **Funções:** inglês, em `snake_case`, com os termos da FAO-56 (Allen et al.,
  1998) sempre que existem.
- **Parâmetros e colunas:** os símbolos da FAO-56 na grafia original (`Ra`,
  `Rs`, `Rso`, `Rn`, `G`, `N`, `n`, `Tmax`, `Tmin`, `Tmean`, `es`, `ea`, `u2`,
  `RH`). Na API em português continuam os símbolos de Pereira, Angelocci &
  Sentelhas (2002): `Qo`, `Qg`, `Tmed`, `UR` etc.
- **Ortografia:** a da FAO-56 ("vapour").
- `sind`, `cosd`, `tand` e `acosd` têm o mesmo nome nas duas APIs.

**Fontes da terminologia:** Allen et al. (1998), FAO-56, cap. 3 (radiação,
umidade, pressão, vento, ETo); Thornthwaite & Mather (1955) e Thornthwaite
(1948) (balanço hídrico e classificação); Frère & Popov (1979) (WRSI);
Doorenbos & Kassam (1979), FAO-33 (ky, produtividade), e van Ittersum et al.
(2013) (Yp, Yw); Battisti et al. (2013) e Battisti & Sentelhas (artigos em
inglês sobre o modelo da Zona Agroecológica no Brasil); Kassam (1977) (FAO
Agro-Ecological Zones).

**Situação:** ✅ implementado · ✅ proposto (aguardando aprovação)

---

## Radiação solar — `radiacao` → `en.radiation` ✅

| Português | Inglês | Parâmetros (pt → en) | Termo de referência |
|---|---|---|---|
| `nda(dia, mes, ano=2023)` | `day_of_year(day, month, year=2023)` | dia → day, mes → month, ano → year | day of the year, J (FAO-56) |
| `declinacao_solar(NDA)` | `solar_declination(doy)` | NDA → doy | solar declination, δ (eq. 24) |
| `angulo_horario(hora, minuto=0)` | `hour_angle(hour, minute=0)` | hora → hour, minuto → minute | hour angle, ω |
| `angulo_zenital(lat, declinacao)` | `solar_zenith_angle(lat, declination)` | declinacao → declination | solar zenith angle (ao meio-dia) |
| `azimute_solar(lat, declinacao, Z)` | `solar_azimuth(lat, declination, zenith)` | Z → zenith | solar azimuth |
| `comprimento_sombra(d, Z)` | `shadow_length(height, zenith)` | d → height, Z → zenith | — |
| `fotoperiodo(Hn)` | `daylight_hours(ws)` | Hn → ws | daylight hours, N (eq. 34) |
| `angulo_horario_nascer(lat, declinacao)` | `sunset_hour_angle(lat, declination)` | — | sunset hour angle, ωs (eq. 25)¹ |
| `fator_correcao_distancia(NDA)` | `inverse_relative_distance(doy)` | NDA → doy | inverse relative distance Earth-Sun, dr (eq. 23) |
| `irradiancia_extraterrestre(lat, declinacao, Hn, dD2)` | `extraterrestrial_radiation(lat, declination, ws, dr)` | Hn → ws, dD2 → dr | extraterrestrial radiation, Ra (eq. 21) |
| `insolacao(N, Tmax, Tmin, lat, k=0.19)` | `sunshine_duration(N, Tmax, Tmin, lat, k=0.19)` | — | actual duration of sunshine, n |
| `Qg_angstrom(insolacao, N, Qo, lat, b=0.52)` | `solar_radiation_angstrom(n, N, Ra, lat, b=0.52)` | insolacao → n, Qo → Ra | solar radiation, Rs; Angstrom formula (eq. 35) |
| `relacao_n_N(Qg, Qo, lat, b=0.52, a=None, limitar=True)` | `relative_sunshine_duration(Rs, Ra, lat, b=0.52, a=None, clip=True)` | Qg → Rs, Qo → Ra, limitar → clip | relative sunshine duration, n/N (eq. 35) |
| `Qg_hargreaves(Tmax, Tmin, Qo, k=0.16)` | `solar_radiation_hargreaves(Tmax, Tmin, Ra, k=0.16)` | Qo → Ra | Hargreaves' radiation formula, kRs (eq. 50) |

¹ Em português, "nascer" (Pereira et al., 2002); na FAO-56, "sunset". O valor é
o mesmo, por simetria.

## Temperatura — `temperatura` → `en.temperature` ✅

| Português | Inglês | Parâmetros | Termo |
|---|---|---|---|
| `temp_media_extremos(Tmax, Tmin)` | `mean_temperature(Tmax, Tmin)` | — | mean daily air temperature, Tmean (eq. 9) |
| `temp_media_estacao_automatica(temperaturas)` | `mean_temperature_observations(temperatures)` | temperaturas → temperatures | — |

## Umidade do ar — `umidade` → `en.humidity` ✅

| Português | Inglês | Parâmetros | Termo |
|---|---|---|---|
| `es_tetens(T_ar)` | `saturation_vapour_pressure(T)` | T_ar → T | saturation vapour pressure, e°(T) (eq. 11) |
| `ea_umidade(es, UR)` | `actual_vapour_pressure(es, RH)` | UR → RH | actual vapour pressure, ea |
| `deficit_saturacao(es, ea)` | `vapour_pressure_deficit(es, ea)` | — | vapour pressure deficit, es − ea |
| `patm_altitude(A)` | `atmospheric_pressure(z)` | A → z | atmospheric pressure, P (eq. 7) |
| `umidade_absoluta(ea, T_ar_C)` | `absolute_humidity(ea, T)` | T_ar_C → T | absolute humidity |
| `umidade_saturacao(es, T_ar_C)` | `saturation_humidity(es, T)` | T_ar_C → T | saturation (absolute) humidity |
| `umidade_relativa(ea, es)` | `relative_humidity(ea, es)` | — | relative humidity, RH (eq. 10) |
| `ponto_orvalho(ea)` | `dewpoint_temperature(ea)` | — | dewpoint temperature, Tdew |
| `constante_psicrometrica(Patm)` | `psychrometric_constant(Patm)` | — | psychrometric constant, γ (eq. 8) |

## Balanço de energia — `energia` → `en.energy` ✅

| Português | Inglês | Parâmetros | Termo |
|---|---|---|---|
| `boc_saldo(Qg, r=0.25)` | `net_shortwave_radiation(Rs, albedo=0.25)` | Qg → Rs, r → albedo | net shortwave radiation, Rns (eq. 38) |
| `bol_saldo(Tmax, Tmin, ea, Qg, Qg_cs)` | `net_longwave_radiation(Tmax, Tmin, ea, Rs, Rso)` | Qg → Rs, Qg_cs → Rso | net longwave radiation, Rnl (eq. 39)² |
| `saldo_radiacao(BOC, BOL)` | `net_radiation(Rns, Rnl)` | BOC → Rns, BOL → Rnl | net radiation, Rn (eq. 40)² |

² O sinal segue o pacote: o BOL é negativo (perda) e Rn = BOC + BOL. Na FAO-56,
Rnl é positivo e Rn = Rns − Rnl. A docstring em inglês deixa isso explícito.

## Evapotranspiração — `evapotranspiracao` → `en.evapotranspiration` ✅

| Português | Inglês | Parâmetros / colunas | Termo |
|---|---|---|---|
| `thornthwaite_mensal(df, col_T='T_media_C', lat=None)` | `pet_thornthwaite(df, col_T='T_mean_C', lat=None)` | coluna de saída ETP_mm_mes → PET_mm_month | potential evapotranspiration, PET (Thornthwaite, 1948) |
| `camargo_maluf_mensal(df, col_T='T_media_C', lat=None)` | `pet_camargo(df, col_T='T_mean_C', lat=None)` | Qo_MJ_m2dia → Ra_MJ_m2day, Qo_mm_dia → Ra_mm_day, ETP_mm_mes → PET_mm_month | — |
| `etp_hargreaves_samani(Qo, Tmax, Tmin, Tmed)` | `pet_hargreaves_samani(Ra, Tmax, Tmin, Tmean)` | Qo → Ra, Tmed → Tmean | Hargreaves equation (FAO-56 eq. 52) |
| `declive_pressao_vapor(T_ar)` | `slope_vapour_pressure_curve(T)` | T_ar → T | slope of saturation vapour pressure curve, Δ (eq. 13) |
| `etp_priestley_taylor(Rn, G, Delta, gamma, alfa=1.26)` | `pet_priestley_taylor(Rn, G, Delta, gamma, alpha=1.26)` | alfa → alpha | Priestley-Taylor (1972) |
| `vento_2m(uz, z)` | `wind_speed_2m(uz, z)` | — | wind speed at 2 m, u2 (eq. 47) |
| `eto_penman_monteith_fao56(Rn, G, Tmed, u2, es, ea, Delta, gamma)` | `eto_penman_monteith(Rn, G, Tmean, u2, es, ea, Delta, gamma)` | Tmed → Tmean | FAO Penman-Monteith, ETo (eq. 6) |

## Graus-dia — `grau_dias` → `en.degree_days` ✅

| Português | Inglês | Parâmetros |
|---|---|---|
| `data_maturacao_fisiologica(df, Tb, CT, dia_semeadura, mes_semeadura, intervalo='d', ano=2023)` | `maturity_date(df, Tb, TT, sowing_day, sowing_month, interval='d', year=2023)` | CT → TT (thermal time, °C·day) |
| `data_semeadura(df, Tb, CT, dia_maturacao, mes_maturacao, intervalo='d', ano=2023)` | `sowing_date(df, Tb, TT, maturity_day, maturity_month, interval='d', year=2023)` | — |

- `intervalo`: `'d'`, `'dec'`, `'M'` → `interval`: `'d'`, `'dekad'`, `'M'`
  ("dekad" é o termo da FAO/OMM para o decêndio).
- Colunas de entrada: dia, mes, Tmed, Tmax, Tmin → day, month, Tmean, Tmax, Tmin.
- Colunas de saída: Meses, data, Tmed, GDi, GDA_mes, GDA_ciclo → month, date,
  Tmean, GDD, GDD_period, GDD_cycle (growing degree days).
- As funções em inglês têm o parâmetro `verbose=True`, que imprime as datas em
  inglês. A API em português não tem esse parâmetro: imprime sempre as suas
  mensagens, em português, como antes.

## Balanço hídrico — `balanco_hidrico` → `en.water_balance` ✅

| Português | Inglês | Parâmetros |
|---|---|---|
| `balanco_hidrico_climatologico(df, CAD=100.0, ciclico=True)` | `climatic_water_balance(df, AWC=100.0, cyclic=True)` | CAD → AWC (available water capacity) |
| `balanco_hidrico_climatologico_grade(P, ETP, CAD=100.0, ciclico=True, tol=0.01, max_iter=100)` | `climatic_water_balance_grid(P, PET, AWC=100.0, cyclic=True, tol=0.01, max_iter=100)` | ETP → PET |
| `balanco_hidrico_cultura(df)` | `crop_water_balance(df)` | — |

Colunas (Thornthwaite & Mather, 1955; Frère & Popov, 1979):

| Português | Inglês | Termo |
|---|---|---|
| Meses, P (mm/mês), ETP (mm/mês) | month, P (mm/month), PET (mm/month) | precipitation, potential evapotranspiration |
| Chuva | P | precipitation |
| CAD | AWC | available water capacity |
| P-ETP | P-PET | — |
| ARM | STO | soil water storage |
| NEG.ACUM | APWL | accumulated potential water loss |
| ALT | dSTO | change in soil water storage |
| ETR | ETa | actual evapotranspiration |
| DEF | DEF | water deficit |
| EXC | SUR | water surplus |
| ISNA | WRSI | water requirement satisfaction index (ETa/ETc) |

## Classificação climática — `classificacao_climatica` → `en.climate_classification` ✅

| Português | Inglês | Parâmetros |
|---|---|---|
| `classificacao_koppen(T, P, lat)` (+ `_grade`) | `koppen_classification(T, P, lat)` (+ `_grid`) | — |
| `classificacao_thornthwaite(P, ETP, lat, CAD=100.0)` (+ `_grade`) | `thornthwaite_classification(P, PET, lat, AWC=100.0)` (+ `_grid`) | — |
| `classificacao_camargo(T, P, ETP, lat, CAD=100.0, tabela_termica='coerente')` (+ `_grade`) | `camargo_classification(T, P, PET, lat, AWC=100.0, thermal_table='consistent')` (+ `_grid`) | 'coerente' → 'consistent', 'literal' → 'literal' |
| `classificacao_holdridge(T, P, lat, ETP=None, limiar_correcao=24.0)` (+ `_grade`) | `holdridge_classification(T, P, lat, PET=None, correction_threshold=24.0)` (+ `_grid`) | — |

- Chaves do resultado: classe → class, descricao → description, umidade →
  moisture, subtipo → subtype, termica → thermal, concentracao →
  concentration, hidrica → moisture, estacao_seca → dry_season, T_anual →
  T_annual, T_mes_frio → T_coldest, ETP_anual → PET_annual, DEF_anual →
  DEF_annual, EXC_anual → SUR_annual, ETP_verao_pct → PET_summer_pct,
  biotemperatura → biotemperature, P_anual → P_annual, razao_ETP → PET_ratio.
- Os códigos das classes (Cfa, B2rB'4a', ST-PU...) não mudam. As descrições
  (`description`) ficam em inglês no núcleo e em português na API em
  português, com os termos de Thornthwaite (1948) (Perhumid, Humid, Moist
  subhumid, Dry subhumid, Semiarid, Arid; Megathermal, Mesothermal,
  Microthermal, Tundra, Frost) e, no Camargo, Glacial, Cold, Moderately cold,
  Temperate, Subtemperate, Subtropical, Tropical, Equatorial / Desert, Arid,
  Dry, Monsoonal, Subhumid, Humid, Perhumid, Extremely humid.
- No Holdridge, `class` é o nome original em inglês (Jungkunst et al., 2021).
  A API em português continua devolvendo `classe` (em português) e
  `classe_en`.
- As constantes públicas do módulo em português (`KOPPEN_CLASSES`,
  `TH_UMIDADE`, `HOLDRIDGE_ZONAS` etc.) continuam disponíveis.

## Produtividade — `produtividade` → `en.crop_yield` ✅

| Português | Inglês | Parâmetros |
|---|---|---|
| `cTn(T, rota, T_limiar=16.5)` | `temperature_correction_overcast(T, pathway, T_threshold=16.5)` | rota → pathway ('C3'/'C4') |
| `cTc(T, rota, T_limiar=16.5)` | `temperature_correction_clear(T, pathway, T_threshold=16.5)` | — |
| `PPBp(Qo, nN, cTc, cTn, a_c, b_c, a_n, b_n)` | `gross_potential_productivity(Ra, nN, cTc, cTo, a_c, b_c, a_o, b_o)` | Qo → Ra; n (nublado) → o (overcast) |
| `CIAF(IAF, IAF_padrao=5.0)` | `lai_correction(LAI, LAI_standard=5.0)` | IAF → LAI |
| `CR(T, T_limiar=20.0, CR_frio=0.6, CR_quente=0.5)` | `respiration_correction(T, T_threshold=20.0, CR_cool=0.6, CR_warm=0.5)` | — |
| `produtividade_potencial(df, lat, rota, IAF=5.0, Cc=0.35, U=13.0, intervalo='dec', semeadura=None, ciclo=None, ...)` | `potential_yield(df, lat, pathway, LAI=5.0, HI=0.35, moisture=13.0, interval='dekad', sowing=None, cycle=None, ...)` | Cc → HI (harvest index), U → moisture |
| `produtividade_atingivel(PPf, ETR, ETc, ky, fase=None, metodo='produtorio')` | `attainable_yield(Yp, ETa, ETc, ky, phase=None, method='product')` | 'produtorio' → 'product', 'etapa_unica' → 'single' |

Colunas e `attrs`:

| Português | Inglês | Termo |
|---|---|---|
| data, dia, mes, ano | date, day, month, year | — |
| Tmed, Qg, Qo, nN, n, N | Tmean, Rs, Ra, nN, n, N | FAO-56 |
| ND | ND | days of the period in the field |
| cTn, cTc | cTo, cTc | temperature corrections, overcast/clear |
| PPBn, PPBc, PPBp | GPo, GPc, GP | gross potential productivity (standard crop) |
| CR, CIAF | CR, LAI_corr | — |
| PPf_periodo, PPf_acum, PPf | Yp_period, Yp_cum, Yp | potential yield, Yp |
| ETR, ETR/ETc, fator | ETa, ETa/ETc, factor | — |
| PA_inicial, PA_final, PA, perda | Yw_start, Yw_end, Yw, loss | attainable (water-limited) yield, Yw³ |
| EC | CE | climatic efficiency (Yw/Yp) |

³ Na literatura de *yield gap* (van Ittersum et al., 2013), Ya é a produtividade
**real** (observada). A produtividade atingível pelo déficit hídrico é a
water-limited yield, Yw; por isso PA → Yw, e não Ya.
