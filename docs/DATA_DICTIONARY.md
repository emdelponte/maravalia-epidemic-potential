# Data Dictionary

This document provides complete descriptions, formats, and units for all variables across datasets in `data/metadata/`, `data/derived/`, and numerical summary tables in `outputs/tables/`.

---

## 1. Derived Data Artifacts (`data/derived/`)

### `phase2_episodes_rh90.rds`
Extracted continuous wetness episodes (RH ≥ 90%, $W \ge 6\text{ h}$) across 65 weather grid cells in Northeastern Brazil (2001–2025).

| Column Name | Type | Units / Format | Description |
|---|---|---|---|
| `cell_id` | String | `lat_XX.X_lon_YY.Y` | Weather grid cell identifier |
| `lat` | Double | Decimal degrees (°S) | Grid cell center latitude |
| `lon` | Double | Decimal degrees (°W) | Grid cell center longitude |
| `start_year` | Integer | YYYY (2001–2025) | Calendar year of episode initiation |
| `start_month` | Integer | 1–12 | Calendar month of episode initiation |
| `start_day` | Integer | 1–31 | Day of month of episode initiation |
| `start_hour` | Integer | 0–23 (UTC) | Hour of episode start |
| `end_hour` | Integer | 0–23 (UTC) | Hour of episode termination |
| `start_idx` | Integer | Count | Row index in continuous hourly weather record |
| `end_idx` | Integer | Count | Row index of episode termination |
| `wet_h` | Integer | Hours (h) | Total duration of continuous wetness ($W = \text{end} - \text{start} + 1$) |
| `T_wet` | Double | Degrees Celsius (°C) | Mean air temperature during wetness episode |
| `I_pot` | Double | Dimensionless [0, 1] | Infection potential calculated via Model 3 ($I_3(T_{\text{wet}}, W)$) |

### `phase2_latent_cohorts.rds`
Hourly tracked physiological development cohorts for all infection episodes with $I_{\text{pot}} > 0$. Contains all columns from `phase2_episodes_rh90.rds`, plus:

| Column Name | Type | Units / Format | Description |
|---|---|---|---|
| `init_idx` | Integer | Index | Hourly index where latent incubation initiates ($\text{end\_idx} + 1$) |
| `spore_idx` | Integer | Index | Hourly index where cumulative development $D(t) \ge 1.0$ (NA if incomplete) |
| `completed` | Logical | TRUE / FALSE | Whether sporulation completed within the weather record |
| `lp_hours` | Double | Hours (h) | Latent period duration in hours |
| `lp_days` | Double | Days | Latent period duration in days ($\text{lp\_hours} / 24$) |
| `spore_year` | Integer | YYYY | Calendar year of sporulation completion |
| `spore_month` | Integer | 1–12 | Calendar month of sporulation completion |
| `spore_day` | Integer | 1–31 | Day of month of sporulation completion |

### `phase2_annual_site_metrics.csv`
Annual aggregate epidemic potential metrics per weather grid cell.

| Column Name | Type | Units / Format | Description |
|---|---|---|---|
| `cell_id` | String | Identifier | Grid cell identifier |
| `lat` | Double | Decimal degrees (°S) | Grid cell latitude |
| `lon` | Double | Decimal degrees (°W) | Grid cell longitude |
| `start_year` | Integer | YYYY (2001–2025) | Calendar year |
| `episodes_count` | Integer | Episodes / year | Number of wetness episodes initiated in year |
| `high_inf_episodes`| Integer | Episodes / year | Number of episodes with high infection potential ($I_{\text{pot}} \ge 0.5$) |
| `EPI_inf` | Double | Cumulative index | Annual sum of infection potential |
| `mean_lp_days` | Double | Days | Mean latent period for cohorts initiated in year |
| `completed_count` | Integer | Cohorts / year | Number of cohorts completing development in year |
| `EPI_spore` | Double | Cumulative index | Annual sum of completed sporulation potential |

### `bootstrap_gamma_samples.rds`
R list object caching parametric bootstrap null replicates under $H_0$ ($B = 2,000$).
- `valid_delta`: Vector of length 2,000 containing simulated $\Delta\text{RSS} = \text{RSS}_{H0} - \text{RSS}_{H1}$.
- `valid_gamma`: Vector of length 2,000 containing fitted dynamic breadth parameter $\gamma$ under $H_0$ data.
- `p_value`: Empirical bootstrap $P$-value.

---

## 2. Output Summary Tables (`outputs/tables/`)

### `site_epidemic_rankings.csv`
25-year climatological ranking of 65 Brazilian weather cells sorted by annual infection pressure.

| Column Name | Type | Units | Description |
|---|---|---|---|
| `cell_id` | String | - | Grid cell code |
| `lat` | Double | °S | Center latitude |
| `lon` | Double | °W | Center longitude |
| `Mean_Annual_EPI_inf` | Double | EPI units / year | 25-year mean annual infection potential (divided by $N = 25$) |
| `Mean_Annual_EPI_spore` | Double | EPI units / year | 25-year mean annual sporulation potential |
| `Mean_Annual_Episodes` | Double | Episodes / year | 25-year mean annual wetness episodes |
| `Mean_Annual_High_Inf_Episodes` | Double | Episodes / year | 25-year mean annual episodes with $I \ge 0.5$ |
| `Mean_Annual_Completed_Cohorts` | Double | Cohorts / year | 25-year mean annual completed sporulation cohorts |
| `Mean_LP_days` | Double | Days | 25-year mean latent period duration |

### `monthly_climatology_summary.csv`
Monthly climatology across all 65 grid cells in Northeastern Brazil (2001–2025).

| Column Name | Type | Units | Description |
|---|---|---|---|
| `start_month` | Integer | 1–12 | Calendar month (1 = January, ..., 12 = December) |
| `Mean_Monthly_EPI_inf` | Double | EPI units / month | Multi-cell average monthly infection potential |
| `SD_Monthly_EPI_inf` | Double | EPI units / month | Standard deviation across 65 grid cells |
| `Mean_Monthly_EPI_spore` | Double | EPI units / month | Multi-cell average monthly sporulation potential |
| `SD_Monthly_EPI_spore` | Double | EPI units / month | Standard deviation of sporulation potential |
| `Mean_Latent_Period_days`| Double | Days | Mean latent period duration for cohorts initiated in month |

### `australian_model_validation_summary.csv`
Validation metrics and historical field impact classifications across 13 Queensland field release sites (Tomley & Evans 2004).

| Column Name | Type | Units | Description |
|---|---|---|---|
| `site_id` | String | - | Station identifier |
| `site_name` | String | - | Geographic locality name |
| `region` | String | - | Queensland subregion |
| `environment_type` | String | - | Vegetation and ecological environment |
| `Mean_Annual_EPI_inf` | Double | EPI units / year | Predicted 25-year mean annual infection potential |
| `Mean_Annual_EPI_spore` | Double | EPI units / year | Predicted 25-year mean annual sporulation potential |
| `Mean_Annual_Episodes` | Double | Episodes / year | Mean annual wetness episodes ($W \ge 6\text{ h}$) |
| `Mean_Annual_High_Inf_Episodes` | Double | Episodes / year | Mean annual high-conduciveness episodes |
| `Mean_LP_days` | Double | Days | Mean latent development duration |
| `lat` | Double | °S | Latitude decimal degrees |
| `lon` | Double | °E | Longitude decimal degrees |
| `observed_field_impact` | String | - | Historical field release observations (Tomley & Evans 2004) |

### `infection_model_metrics.csv`
Goodness-of-fit and cross-validation comparison across primary candidate models.

| Column Name | Type | Units | Description |
|---|---|---|---|
| `Model` | String | - | Model description |
| `Parameters` | Integer | Count | Number of estimated parameters ($p$) |
| `RMSE` | Double | Dimensionless | Root mean squared error on relative infection response |
| `R2` | Double | Dimensionless | Coefficient of determination ($R^2$) |
| `AICc` | Double | Dimensionless | Second-order Akaike Information Criterion |
| `deltaAICc` | Double | Dimensionless | $\Delta\text{AICc} = \text{AICc} - \min(\text{AICc})$ |
| `LOWO_RMSE` | Double | Dimensionless | Leave-one-wetness-level-out cross-validation RMSE |
| `LOTO_RMSE` | Double | Dimensionless | Leave-one-temperature-level-out cross-validation RMSE |

### `infection_model_parameters.csv`
Nonlinear least-squares parameter estimates for candidate infection models.

| Column Name | Type | Units | Description |
|---|---|---|---|
| `Model` | String | - | Model name |
| `Parameter` | String | - | Parameter symbol ($T_{\text{opt}}, \sigma_{\bar{W}}, \gamma, \lambda, k, a, b, c_0, c_1$) |
| `Estimate` | Double | Model-specific | Point estimate |
| `Std_Error` | Double | Model-specific | Asymptotic standard error |
| `t_value` | Double | - | $t$-statistic |
| `p_value` | Double | - | Two-tailed $P$-value |

### `bootstrap_gamma_results.csv`
Summary statistics of parametric bootstrap test for dynamic breadth parameter $\gamma$.

| Column Name | Type | Description |
|---|---|---|
| `Metric` | String | Statistic name (Observed $\Delta\text{RSS}$, 95th Percentile $H_0$, Bootstrap $P$-value, etc.) |
| `Value` | Double | Numerical value |

### `rh_sensitivity_comparison.csv`
Epidemic conduciveness under varying relative humidity thresholds (85%, 90%, 95%).

| Column Name | Type | Description |
|---|---|---|
| `RH_Threshold` | String | Relative humidity threshold string ("85%", "90%", "95%") |
| `Total_Episodes` | Integer | Total continuous wetness episodes recorded |
| `Episodes_ge_6h` | Integer | Episodes lasting $\ge 6$ hours |
| `Episodes_ge_8h` | Integer | Episodes lasting $\ge 8$ hours |
| `Episodes_ge_12h`| Integer | Episodes lasting $\ge 12$ hours |
| `Mean_Wetness_Duration_h` | Double | Average episode duration (hours) |
| `Mean_Wetness_Temp_C` | Double | Average incubation temperature during wetness (°C) |
| `Total_Infection_Potential` | Double | Sum of $I_{\text{pot}}$ across all episodes |
| `Mean_Annual_Infection_Potential_per_Site` | Double | Mean annual $EPI_{\text{inf}}$ per site |

### `digitization_sensitivity_params.csv`
Parameter confidence intervals generated by stochastic perturbations ($\pm 0.05\text{ pustules cm}^{-2}$).

| Column Name | Type | Description |
|---|---|---|
| `Submodel` | String | Parameter identifier |
| `Baseline_Estimate` | Double | Point estimate from unperturbed dataset |
| `P05_Sensitivity` | Double | 5th percentile under perturbation |
| `P50_Sensitivity` | Double | Median under perturbation |
| `P95_Sensitivity` | Double | 95th percentile under perturbation |
