"""Agrometeorology formulas in Python — English API.

This is the core of the package: every function is implemented here, with
FAO-56 terminology and symbols (Allen et al., 1998) wherever they exist.
The same functions are available with Portuguese names at the top level of
``agrometeorologiapy`` (see docs/API_PT_EN.md for the name mapping).

>>> from agrometeorologiapy import en as agm
>>> round(agm.relative_sunshine_duration(Rs=22.5, Ra=42.96, lat=-24.86), 3)
0.501
"""

from .._trig import acosd, cosd, sind, tand
from .radiation import (
    day_of_year,
    daylight_hours,
    extraterrestrial_radiation,
    hour_angle,
    inverse_relative_distance,
    relative_sunshine_duration,
    shadow_length,
    solar_azimuth,
    solar_declination,
    solar_radiation_angstrom,
    solar_radiation_hargreaves,
    solar_zenith_angle,
    sunset_hour_angle,
    sunshine_duration,
)
from .temperature import (
    mean_temperature,
    mean_temperature_observations,
)
from .humidity import (
    absolute_humidity,
    actual_vapour_pressure,
    atmospheric_pressure,
    dewpoint_temperature,
    psychrometric_constant,
    relative_humidity,
    saturation_humidity,
    saturation_vapour_pressure,
    vapour_pressure_deficit,
)
from .energy import (
    net_longwave_radiation,
    net_radiation,
    net_shortwave_radiation,
)
from .evapotranspiration import (
    eto_penman_monteith,
    pet_camargo,
    pet_hargreaves_samani,
    pet_priestley_taylor,
    pet_thornthwaite,
    slope_vapour_pressure_curve,
    wind_speed_2m,
)
from .degree_days import (
    maturity_date,
    sowing_date,
)
from .water_balance import (
    climatic_water_balance,
    climatic_water_balance_grid,
    crop_water_balance,
)
from .climate_classification import (
    camargo_classification,
    camargo_classification_grid,
    holdridge_classification,
    holdridge_classification_grid,
    koppen_classification,
    koppen_classification_grid,
    thornthwaite_classification,
    thornthwaite_classification_grid,
)
from .crop_yield import (
    attainable_yield,
    gross_potential_productivity,
    lai_correction,
    potential_yield,
    respiration_correction,
    temperature_correction_clear,
    temperature_correction_overcast,
)
from .climate_risk import (
    frost_probability,
    heat_stress_days,
    heat_stress_degree_days,
    leaf_wetness_duration,
    return_period,
    temperature_humidity_index,
    thi_category,
)
from .model_evaluation import (
    confidence_index,
    confidence_index_category,
    evaluate_model,
    kge,
    mae,
    mbe,
    nse,
    pearson_r,
    rmse,
    willmott_d,
)

__all__ = [
    "sind", "cosd", "tand", "acosd",
    # radiation
    "day_of_year",
    "solar_declination",
    "hour_angle",
    "solar_zenith_angle",
    "solar_azimuth",
    "shadow_length",
    "daylight_hours",
    "sunset_hour_angle",
    "inverse_relative_distance",
    "extraterrestrial_radiation",
    "sunshine_duration",
    "solar_radiation_angstrom",
    "relative_sunshine_duration",
    "solar_radiation_hargreaves",
    # temperature
    "mean_temperature",
    "mean_temperature_observations",
    # humidity
    "saturation_vapour_pressure",
    "actual_vapour_pressure",
    "vapour_pressure_deficit",
    "atmospheric_pressure",
    "absolute_humidity",
    "saturation_humidity",
    "relative_humidity",
    "dewpoint_temperature",
    "psychrometric_constant",
    # energy
    "net_radiation",
    "net_shortwave_radiation",
    "net_longwave_radiation",
    # evapotranspiration
    "pet_thornthwaite",
    "pet_camargo",
    "pet_hargreaves_samani",
    "slope_vapour_pressure_curve",
    "pet_priestley_taylor",
    "wind_speed_2m",
    "eto_penman_monteith",
    # degree_days
    "maturity_date",
    "sowing_date",
    # water_balance
    "climatic_water_balance",
    "climatic_water_balance_grid",
    "crop_water_balance",
    # climate_classification
    "koppen_classification",
    "koppen_classification_grid",
    "thornthwaite_classification",
    "thornthwaite_classification_grid",
    "camargo_classification",
    "camargo_classification_grid",
    "holdridge_classification",
    "holdridge_classification_grid",
    # crop_yield
    "temperature_correction_overcast",
    "temperature_correction_clear",
    "gross_potential_productivity",
    "lai_correction",
    "respiration_correction",
    "potential_yield",
    "attainable_yield",
    # climate_risk
    "frost_probability",
    "return_period",
    "temperature_humidity_index",
    "thi_category",
    "heat_stress_days",
    "heat_stress_degree_days",
    "leaf_wetness_duration",
    # model_evaluation
    "mbe",
    "mae",
    "rmse",
    "pearson_r",
    "willmott_d",
    "confidence_index",
    "confidence_index_category",
    "nse",
    "kge",
    "evaluate_model",
]
