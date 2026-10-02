# Pipeline v2 Run Report

**Date**: 2026-09-28  
**Working Directory**: `<project root>` (portable relative path)  
**Execution Environment**: R version 4.4.1 (macOS aarch64)

---

## 1. Sanity-Check Summary (Expected vs. Obtained)

| Check / Metric | Target / Expected | Obtained (v2 Default) | Status | Difference / Notes |
|---|---|---|---|---|
| **Model 3 $T_{\text{opt}}$** | 21.26 | 21.2597 (21.26) | **PASS** | Identical to `tables/` (<0.01% diff) |
| **Model 3 $\bar{\sigma}$ ($\sigma_{W_{\text{bar}}}$)** | 2.280 | 2.2798 (2.280) | **PASS** | Identical to `tables/` (<0.01% diff) |
| **Model 3 $\gamma$** | 0.2984 | 0.2984 | **PASS** | Identical to `tables/` (<0.01% diff) |
| **Model 3 $\lambda$** | 8.367 | 8.3666 (8.367) | **PASS** | Identical to `tables/` (<0.01% diff) |
| **Model 2 $c_0$ bound** | Not at 20; not on bound | 26.2017 | **PASS** | Bounds [0.01, 200]; interior solution |
| **Bootstrap Observed $\Delta\text{RSS}$** | ~0.260 (not 0.258) | 0.2601 | **PASS** | Exact: 0.260125 |
| **Bootstrap Negative $\Delta\text{RSS}$** | None ($\ge 0$) | Min = 0.0000 | **PASS** | 0 negative values across 2,000 replicates |
| **Bootstrap Failed Fits** | 0 | 0 / 2,000 | **PASS** | 100% convergence rate |
| **Bootstrap $P$-value** | Low ($<0.01$) | 0.001999 | **PASS** | $(3 + 1) / (2000 + 1) = 4 / 2001$ |
| **Bootstrap 95th Percentile** | - | 0.06788 | **PASS** | Observed 0.2601 exceeds 95th percentile |
| **RH 90% Mean Episode Duration** | ~9 h (not 31 h) | 9.41 h | **PASS** | Previous was 31.2 h with $P>0$ clause |
| **RH 90% Mean Annual $\text{EPI}_{\text{inf}}$ (BR)** | ~63 | 62.83 | **PASS** | Exact site mean: 62.83 per site/year |
| **Brazil Peak Month** | March–June | April (13.02) | **PASS** | March (11.50), May (11.91), June (9.43) |
| **Brazil Sep–Dec $\text{EPI}_{\text{inf}}$** | Near zero | Sep: 0.00, Oct: 0.015, Nov: 0.013, Dec: 0.22 | **PASS** | Off-season average <0.07/month |
| **Brazil May–July Share** | - | 42.68% | **PASS** | 26.88 / 62.98 annual total |
| **Brazil March–May Share** | - | 57.85% | **PASS** | 36.43 / 62.98 annual total |
| **Brazil March–June Share** | - | 72.82% | **PASS** | 45.86 / 62.98 annual total |
| **Australia: McLeod River** | ~161 | 161.00 | **PASS** | Python reference: 161.03 (diff = 0.02%) |
| **Australia: Inkerman** | ~6 | 5.74 | **PASS** | Python reference: 5.74 (diff = 0.00%) |
| **Australia: Delta Downs** | ~0.7 | 0.82 | **FLAGGED (>2%)** | Python reference: 0.72 (diff = 13.6%). Cause: R averages over 22 non-zero episode years ($0.818 \times 22 / 25 = 0.720$) |
| **Australia: Welch $t$** | - | $t = 1.78$, $\text{df} = 4.6$, $P = 0.140$ | **PASS** | Script 12 output |
| **Australia: Mann–Whitney $W$** | - | $W = 32$, $P = 0.093$ | **PASS** | $n_1 = 5, n_2 = 8$ |

---

## 2. Model 2 & Bootstrap Results

### 2.1 Infection Model Parameters (`tables_v2/infection_model_parameters.csv`)

#### Model 2: Wetness-Expanded Analytis $\times$ Weibull
- Formula: $I_{\text{obs}} \sim \left[a \left(\frac{T - 10}{25}\right)^b \left(1 - \frac{T - 10}{25}\right)\right]^{c_0 + c_1 \cdot \text{dew\_h}} \left(1 - \exp\left(-\left(\frac{\text{dew\_h}}{\lambda}\right)^k\right)\right)$
- Bounds: $a \in [0.01, 100], b \in [0.01, 20], c_0 \in [0.01, 200], c_1 \in [-1.0, 5.0], \lambda \in [0.1, 50], k \in [0.1, 20]$

| Parameter | Estimate | Std. Error | $t$-value | $P$-value | Boundary Status |
|---|---|---|---|---|---|
| $a$ | 3.394190 | 0.123720 | 27.434 | $1.43 \times 10^{-13}$ | Interior |
| $b$ | 0.792513 | 0.043275 | 18.313 | $3.54 \times 10^{-11}$ | Interior |
| $c_0$ | 26.201733 | 7.603607 | 3.446 | 0.00394 | **Interior (was bounded at 20 in v1)** |
| $c_1$ | -1.000000 | 0.316932 | -3.155 | 0.00702 | On lower bound (-1.0) |
| $\lambda$ | 8.785261 | 1.047417 | 8.388 | $7.86 \times 10^{-7}$ | Interior |
| $k$ | 8.406592 | 9.418183 | 0.893 | 0.38716 | Interior |

#### Model 3: Centered Gaussian Dynamic Breadth $\times$ Weibull (Primary Model)
- Formula: $I_{\text{obs}} \sim \exp\left(-0.5 \left(\frac{T - T_{\text{opt}}}{\max(1.0, \sigma_{\bar{W}} + \gamma (\text{dew\_h} - 12.5))}\right)^2\right) \left(1 - \exp\left(-\left(\frac{\text{dew\_h}}{\lambda}\right)^k\right)\right)$

| Parameter | Estimate | Std. Error | $t$-value | $P$-value | Boundary Status |
|---|---|---|---|---|---|
| $T_{\text{opt}}$ | 21.259724 | 0.206830 | 102.788 | $8.79 \times 10^{-23}$ | Interior |
| $\sigma_{\bar{W}}$ | 2.279779 | 0.253794 | 8.983 | $2.01 \times 10^{-7}$ | Interior |
| $\gamma$ | 0.298371 | 0.063040 | 4.733 | $2.67 \times 10^{-4}$ | Interior |
| $\lambda$ | 8.366565 | 0.555938 | 15.049 | $1.85 \times 10^{-10}$ | Interior |
| $k$ | 8.928839 | 9.138353 | 0.977 | 0.34403 | Interior |

### 2.2 Parametric Bootstrap Results (`tables_v2/bootstrap_gamma_results.csv`)

| Metric | Value |
|---|---|
| Replicates ($B$) | 2,000 |
| Failed fits | 0 / 2,000 (0.0%) |
| Observed RSS under $H_0$ (fixed breadth, $\gamma = 0$) | 0.4495 |
| Observed RSS under $H_1$ (dynamic breadth, $\gamma > 0$) | 0.1894 |
| Observed $\Delta\text{RSS}$ ($\text{RSS}_{H_0} - \text{RSS}_{H_1}$) | **0.260125** |
| Minimum $\Delta\text{RSS}$ in replicates | **0.000000** |
| Replicates with negative $\Delta\text{RSS}$ | **0** (0.0%) |
| Mean $H_0$ $\Delta\text{RSS}$ | 0.011640 |
| 95th Percentile $H_0$ $\Delta\text{RSS}$ | **0.067885** |
| Bootstrap $P$-value: $(x + 1) / (B + 1)$ | **0.001999** |
| $\gamma$ under $H_0$: Median | 0.000000 |
| $\gamma$ under $H_0$: 2.5% quantile | 0.000000 |
| $\gamma$ under $H_0$: 97.5% quantile | 0.301319 |

---

## 3. Brazil Pipeline Outputs

### 3.1 Monthly Climatology Summary (`tables_v2/monthly_climatology_summary.csv`)

| Month (`start_month`) | Mean Monthly $\text{EPI}_{\text{inf}}$ | SD Monthly $\text{EPI}_{\text{inf}}$ | Mean Monthly $\text{EPI}_{\text{spore}}$ | SD Monthly $\text{EPI}_{\text{spore}}$ | Mean Latent Period (days) |
|---|---|---|---|---|---|
| 1 (Jan) | 3.3824 | 2.4333 | 1.5936 | 1.2615 | 9.59 |
| 2 (Feb) | 6.6626 | 4.3737 | 5.3682 | 3.6830 | 9.63 |
| 3 (Mar) | 11.5039 | 6.2379 | 10.2490 | 5.8100 | 9.69 |
| 4 (Apr) | **13.0159** | 6.3705 | 12.6869 | 6.3502 | 9.73 |
| 5 (May) | 11.9110 | 6.6333 | **12.7275** | 6.6620 | 9.83 |
| 6 (Jun) | 9.4277 | 5.9391 | 10.3771 | 6.1943 | 10.06 |
| 7 (Jul) | 5.5375 | 4.0954 | 7.3883 | 5.1969 | 10.09 |
| 8 (Aug) | 1.2849 | 1.0169 | 2.2715 | 1.8154 | 9.97 |
| 9 (Sep) | 0.0000 (NA) | NA | 0.3263 | 0.2813 | 9.76 |
| 10 (Oct) | 0.0153 | 0.0174 | 0.0200 | 0.0224 | 9.41 |
| 11 (Nov) | 0.0130 | 0.0154 | 0.0145 | 0.0176 | 9.33 |
| 12 (Dec) | 0.2238 | 0.2189 | 0.1418 | 0.1401 | 9.41 |

- **Annual Total Mean $\text{EPI}_{\text{inf}}$**: 62.98
- **May–July Cumulative $\text{EPI}_{\text{inf}}$**: 26.88 (42.68% of annual total)
- **March–May Cumulative $\text{EPI}_{\text{inf}}$**: 36.43 (57.85% of annual total)
- **March–June Cumulative $\text{EPI}_{\text{inf}}$**: 45.86 (72.82% of annual total)

### 3.2 10 Highest-EPI Grid Cells (`tables_v2/site_epidemic_rankings.csv`)

| Rank | Cell ID | Latitude (°S) | Longitude (°W) | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual $\text{EPI}_{\text{spore}}$ | Mean Annual Episodes | Mean LP (days) |
|---|---|---|---|---|---|---|---|
| 1 | `lat_-3.8_lon_-40.8` | -3.8 | -40.8 | 117.84 | 117.84 | 188.52 | 9.54 |
| 2 | `lat_-3.9_lon_-40.4` | -3.9 | -40.4 | 117.84 | 117.84 | 188.52 | 9.54 |
| 3 | `lat_-3.9_lon_-40.5` | -3.9 | -40.5 | 117.84 | 117.84 | 188.52 | 9.54 |
| 4 | `lat_-4.2_lon_-40.4` | -4.2 | -40.4 | 117.84 | 117.84 | 188.52 | 9.54 |
| 5 | `lat_-4.3_lon_-39.3` | -4.3 | -39.3 | 111.37 | 111.37 | 186.84 | 9.60 |
| 6 | `lat_-4.7_lon_-39.5` | -4.7 | -39.5 | 111.37 | 111.37 | 186.84 | 9.60 |
| 7 | `lat_-3.3_lon_-40.5` | -3.3 | -40.5 | 108.44 | 108.44 | 197.44 | 9.48 |
| 8 | `lat_-3.4_lon_-40.7` | -3.4 | -40.7 | 108.44 | 108.44 | 197.44 | 9.48 |
| 9 | `lat_-3.6_lon_-40.6` | -3.6 | -40.6 | 108.44 | 108.44 | 197.44 | 9.48 |
| 10 | `lat_-4.0_lon_-40.1` | -4.0 | -40.1 | 105.40 | 105.40 | 196.20 | 9.48 |

### 3.3 10 Lowest-EPI Grid Cells (`tables_v2/site_epidemic_rankings.csv`)

| Rank | Cell ID | Latitude (°S) | Longitude (°W) | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual $\text{EPI}_{\text{spore}}$ | Mean Annual Episodes | Mean LP (days) |
|---|---|---|---|---|---|---|---|
| 56 | `lat_-3.1_lon_-40.2` | -3.1 | -40.2 | 25.09 | 25.09 | 150.20 | 9.53 |
| 57 | `lat_-3.2_lon_-39.9` | -3.2 | -39.9 | 25.09 | 25.09 | 150.20 | 9.53 |
| 58 | `lat_-3.2_lon_-40.1` | -3.2 | -40.1 | 25.09 | 25.09 | 150.20 | 9.53 |
| 59 | `lat_-3.0_lon_-41.1` | -3.0 | -41.1 | 23.04 | 23.04 | 154.80 | 9.38 |
| 60 | `lat_-3.0_lon_-41.2` | -3.0 | -41.2 | 23.04 | 23.04 | 154.80 | 9.38 |
| 61 | `lat_-4.2_lon_-38.2` | -4.2 | -38.2 | 2.04 | 2.04 | 32.92 | 9.77 |
| 62 | `lat_-3.6_lon_-39.0` | -3.6 | -39.0 | 0.94 | 0.94 | 23.72 | 9.88 |
| 63 | `lat_-3.7_lon_-38.8` | -3.7 | -38.8 | 0.94 | 0.94 | 23.72 | 9.88 |
| 64 | `lat_-4.6_lon_-37.8` | -4.6 | -37.8 | 0.02 | 0.02 | 2.58 | 9.71 |
| 65 | `lat_-4.7_lon_-37.8` | -4.7 | -37.8 | 0.02 | 0.02 | 2.58 | 9.71 |

---

## 4. Australia Field Validation Outputs (`tables_v2/australian_model_validation_summary.csv`)

### 4.1 Site Rankings and Observed Field Impacts

| Site ID | Site Name | Region | Environment Type | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual $\text{EPI}_{\text{spore}}$ | Mean Annual Episodes | Mean LP (days) | Observed Field Impact | Impact Category |
|---|---|---|---|---|---|---|---|---|---|
| `aus_site_08` | McLeod River | Cooktown Hinterland | Riparian / Watercourse | 161.00 | 160.69 | 332.12 | 13.93 | Heavy rust (47 pustules cm^-2); total defoliation | High impact |
| `aus_site_13` | Laura Lea | Cape York Peninsula | Moist South-East Airflow | 80.09 | 80.02 | 202.84 | 10.81 | Enhanced rust intensity; 20% mortality by 1999 | High impact |
| `aus_site_02` | Lakefield NP | Cape York Peninsula | Moist Coastal Airflow | 69.30 | 69.27 | 190.92 | 10.34 | High defoliation & 20-75% weed mortality | High impact |
| `aus_site_01` | Chillagoe NP | Tropical North QLD | Inland / Karst | 57.53 | 57.46 | 202.60 | 12.06 | Spread to 80 km by Aug 1996 | Moderate/low impact |
| `aus_site_11` | Charters Towers | Inland Central QLD | Inland Research Hub (TWRC) | 42.50 | 42.37 | 176.92 | 13.58 | Key monitoring hub; heavy rust in wet years | High impact |
| `aus_site_03` | Wrotham Park | Gulf Region | Gulf Savanna | 41.44 | 41.39 | 154.32 | 10.68 | Rapid spread to 50 km in 8 months | Moderate/low impact |
| `aus_site_10` | Rockhampton | Central QLD Coast | Subtropical Coast | 35.83 | 35.82 | 157.32 | 14.64 | Rapid spread over 100 km in 14 months | Moderate/low impact |
| `aus_site_07` | Georgetown | Etheridge Inland | Elevated Dry Inland | 23.03 | 22.88 | 86.84 | 12.23 | Moderate rust; shorter wet season duration | Moderate/low impact |
| `aus_site_04` | Rutland Plains | Gulf Coast | Gulf Lowlands | 15.42 | 15.42 | 100.88 | 9.39 | Satellite spread >20 km detected | Moderate/low impact |
| `aus_site_09` | Strathmore | Gulf / Croydon | Gulf Savanna | 11.68 | 11.63 | 52.04 | 9.65 | Spread to 45 km by Aug 1996 | Moderate/low impact |
| `aus_site_12` | Hughenden | Dry Inland QLD | Semi-Arid (<600mm rain) | 8.05 | 8.00 | 36.52 | 14.24 | Inadequate control in most years due to dry air | Moderate/low impact |
| `aus_site_05` | Inkerman | Burdekin Coast | Sand Ridge / Clay Floodplain | 5.74 | 5.74 | 21.92 | 10.95 | Severe rust damage; up to 75% mortality | High impact |
| `aus_site_06` | Delta Downs | Gulf of Carpentaria | Coastal Gulf | 0.82 | 0.82 | 4.05 | 10.31 | Spread to 70 km in 18 months | Moderate/low impact |

### 4.2 Statistical Tests

- **High-Impact Sites ($n = 5$) Mean $\text{EPI}_{\text{inf}}$**: 71.73
- **Other Sites ($n = 8$) Mean $\text{EPI}_{\text{inf}}$**: 24.22
- **Difference**: 47.51 EPI units
- **Welch Two-Sample $t$-test**:
  - $t = 1.78$
  - $\text{df} = 4.6$
  - $P = 0.140$
- **Mann–Whitney $U$ test (Wilcoxon rank sum)**:
  - $W = 32$
  - $P = 0.093$

---

## 5. Comparison: Default Run vs. Sensitivity Runs (a) and (b)

- **Default (v2)**: `PRECIP_THRESH = NA`, `MAX_EPISODE_H = NA` (outputs in `tables_v2/`)
- **Sensitivity (a)**: `PRECIP_THRESH = 1`, `MAX_EPISODE_H = NA` (outputs in `tables_v2_precip1/`)
- **Sensitivity (b)**: `PRECIP_THRESH = NA`, `MAX_EPISODE_H = 24` (outputs in `tables_v2_maxh24/`)

| Metric | Default v2 (`tables_v2`) | Sensitivity (a) $P \ge 1$ (`tables_v2_precip1`) | Sensitivity (b) $W \le 24$ (`tables_v2_maxh24`) |
|---|---|---|---|
| **Brazil Mean Annual $\text{EPI}_{\text{inf}}$** | 62.83 | 51.13 | 62.85 |
| **Brazil Mean Annual $\text{EPI}_{\text{spore}}$** | 62.83 | 51.12 | 62.85 |
| **Brazil Peak Month** | Month 4 / April (13.02) | Month 5 / May (10.13) | Month 4 / April (13.02) |
| **Australia All Sites Mean $\text{EPI}_{\text{inf}}$** | 42.49 | 41.38 | 42.49 |
| **Australia High-Impact Sites Mean** | 71.73 | 60.71 | 71.73 |
| **Australia Other Sites Mean** | 24.22 | 29.29 | 24.22 |
| **Australia Difference (High - Other)** | 47.51 | 31.42 | 47.51 |
| **Welch $t$-statistic** | 1.78 | 2.08 | 1.78 |
| **Welch $\text{df}$** | 4.6 | 5.5 | 4.6 |
| **Welch $P$-value** | 0.140 | 0.087 | 0.140 |
| **Mann–Whitney $W$** | 32 | 34 | 32 |
| **Mann–Whitney $P$-value** | 0.093 | 0.045 | 0.093 |

*Note: In `12_australian_validation_pipeline_v2.R`, wet episodes are identified via `rle(is_wet)` rather than `wet_runs(is_wet)`, so `MAX_EPISODE_H = 24` operates on Brazil (where episodes are segmented) while Australian site metrics remain unchanged from default.*

---

## 6. Code Fix Log

| Date / Step | File | Line | Error / Issue Encountered | Action Taken |
|---|---|---|---|---|
| Main pipeline execution | All v2 scripts (`02`, `04`, `07`, `08`, `09`, `10`, `12`) | N/A | None | **No script failures occurred.** All required CRAN packages were present and all scripts executed to completion with return code 0 on the first attempt without bug fixes. |
| Sensitivity run (a) | `R/maravalia_functions_v2.R` | 9 | Parameter modification per protocol | Changed `PRECIP_THRESH <- NA` to `PRECIP_THRESH <- 1`. Re-ran scripts 07→12. Copied `tables_v2/` to `tables_v2_precip1/`. |
| Sensitivity run (b) | `R/maravalia_functions_v2.R` | 9, 15 | Parameter modification per protocol | Changed `PRECIP_THRESH <- NA` and `MAX_EPISODE_H <- 24`. Re-ran scripts 07→12. Copied `tables_v2/` to `tables_v2_maxh24/`. |
| Default restoration | `R/maravalia_functions_v2.R` | 15 | Protocol completion | Restored `MAX_EPISODE_H <- NA` (with `PRECIP_THRESH <- NA`). Re-ran scripts 07→12 so that `tables_v2/` and `data_v2/` remain in the default state. |

---

## 7. Figure Files in `figures_v2/`

All figures were saved at 300 DPI publication quality to `figures_v2/`:

| Filename | Dimensions / Resolution | File Size | Description |
|---|---|---|---|
| `fig1_infection_submodel_main.png` | $10 \times 8.5$ in, 300 dpi | 523,880 bytes | 4-panel (2x2) main infection submodel figure (A: raw pustules, B: $T$ slices, C: $W$ slices, D: contour surface) |
| `fig5_combined_wetness_episodes.png` | $10 \times 4.2$ in, 300 dpi | 293,334 bytes | Combined wetness episode diagnostic (A: RH sensitivity column chart, B: scatter of $I_{\text{pot}}$ vs. duration) |
| `fig6_combined_brazil_climatology.png` | $8 \times 7.5$ in, 300 dpi | 184,787 bytes | Composite Brazil climatology (A: monthly infection and sporulation pressure, B: latent period duration boxplot) |
| `fig6a_latent_period_seasonal_variation.png` | $7 \times 4.5$ in, 300 dpi | 97,458 bytes | Standalone boxplot of completed cohort latent period (days) across initiation months |
| `fig7a_monthly_epidemic_pressure_climatology.png` | $7.5 \times 4.2$ in, 300 dpi | 97,386 bytes | Standalone monthly climatology of $\text{EPI}_{\text{inf}}$ bars and $\text{EPI}_{\text{spore}}$ line |
| `fig8_combined_spatial_maps.png` | $10 \times 9.5$ in, 300 dpi | 616,165 bytes | Multi-panel spatial composite: A ($\text{EPI}_{\text{inf}}$ map), B ($\text{EPI}_{\text{spore}}$ map), C (infestation risk overlay) |
| `fig8a_spatial_map_infection_pressure.png` | $6.5 \times 5.5$ in, 300 dpi | 292,702 bytes | Spatial map of Northeastern Brazil showing $\text{EPI}_{\text{inf}}$ with state boundaries and postal labels |
| `fig8b_spatial_map_sporulation_pressure.png` | $6.5 \times 5.5$ in, 300 dpi | 299,306 bytes | Spatial map of Northeastern Brazil showing $\text{EPI}_{\text{spore}}$ with state boundaries and postal labels |
| `fig8c_host_infestation_risk_overlay.png` | $7 \times 4.5$ in, 300 dpi | 148,277 bytes | Scatterplot and linear trend of nearest grid cell $\text{EPI}_{\text{inf}}$ vs. $\log_{10}(\text{infested area})$ |
| `fig9_combined_australian_validation.png` | $8.5 \times 14$ in, 300 dpi | 363,477 bytes | 3x1 vertical stack: A (monthly climatology), B (site ranking), C (Queensland spatial map) |
| `fig9a_australian_seasonal_climatology.png` | $7.5 \times 4.2$ in, 300 dpi | 99,468 bytes | Monthly climatology across Australian field sites ($\text{EPI}_{\text{inf}}$ and $\text{EPI}_{\text{spore}}$) |
| `fig9b_australian_site_risk_ranking.png` | $8 \times 5$ in, 300 dpi | 97,727 bytes | Horizontal bar chart of Queensland release sites ranked by $\text{EPI}_{\text{inf}}$, colored by field impact |
| `fig9c_australian_site_map.png` | $7 \times 6$ in, 300 dpi | 161,733 bytes | Geographic map of Queensland field sites with repelled labels and $\text{EPI}_{\text{inf}}$ color fill |
| `figS1_all_infection_models_comparison.png` | $8 \times 6.5$ in, 300 dpi | 268,429 bytes | Multi-model comparison across wetness slices (Models 1, 2, 3, and Duthie) |
| `figS2_infection_residuals_diagnostics.png` | $11 \times 3.8$ in, 300 dpi | 104,512 bytes | Residual diagnostic plots across fitted values for Models 1, 2, and 3 |
| `figS3_bootstrap_gamma_distribution.png` | $7 \times 4.5$ in, 300 dpi | 59,334 bytes | Histogram of bootstrap $\Delta\text{RSS}$ replicates under $H_0$ with 95th percentile and observed statistic |

---

## 8. Round 2 Execution & Post-Fix Validation

This section documents the results following the **Round 2 fixes** applied on 2026-09-28:
1. **Model 2 Boundary Expansion** (`R/02_infection_submodel_v2.R`): $c_1$ lower bound extended from $-1.0$ to $-10.0$ with supplementary grid starts.
2. **Annual Climatology Normalization** (`R/09_epidemic_pressure_indices_v2.R` and `R/12_australian_validation_pipeline_v2.R`): Annual site summaries now divide by $N_{\text{YEARS}} = 25$ explicitly rather than averaging only over non-zero epidemic years.
3. **Monthly Climatology Zero-Completion** (`R/09_epidemic_pressure_indices_v2.R`): Monthly cell combinations with 0 infection events are now filled with zeros across all 12 months, resolving undefined (NA) standard deviations for dry off-season months (e.g. September).
4. **Shared Episode Splitting in Australia** (`R/12_australian_validation_pipeline_v2.R`): Australian episodes are now extracted using the shared `wet_runs(is_wet)` function, allowing `MAX_EPISODE_H = 24` to partition multi-day wetness episodes.

---

### 8.1 Model 2 Parameters & Candidate Model Information Criteria

#### Model 2 Parameter Estimates (Bound $[-10, 5]$ on $c_1$)
- Formula: $I_{\text{obs}} \sim \left[a \left(\frac{T - 10}{25}\right)^b \left(1 - \frac{T - 10}{25}\right)\right]^{c_0 + c_1 \cdot \text{dew\_h}} \left(1 - \exp\left(-\left(\frac{\text{dew\_h}}{\lambda}\right)^k\right)\right)$

| Parameter | Estimate | Std. Error | $t$-value | $P$-value | Bounds | Boundary Status |
|---|---|---|---|---|---|---|
| $a$ | 3.395972 | 0.119029 | 28.531 | $8.34 \times 10^{-14}$ | $[0.01, 100]$ | **Interior** |
| $b$ | 0.790922 | 0.041705 | 18.964 | $2.21 \times 10^{-11}$ | $[0.01, 20]$ | **Interior** |
| $c_0$ | 29.186176 | 8.618326 | 3.387 | 0.00443 | $[0.01, 200]$ | **Interior** |
| $c_1$ | -1.124675 | 0.359248 | -3.131 | 0.00737 | $[-10.0, 5.0]$ | **Interior (moved off lower bound $-1.0$)** |
| $\lambda$ | 8.796228 | 1.062667 | 8.278 | $9.19 \times 10^{-7}$ | $[0.1, 50]$ | **Interior** |
| $k$ | 8.423668 | 9.476344 | 0.889 | 0.38907 | $[0.1, 20]$ | **Interior** |

*Note: With the expanded lower bound on $c_1$, all parameters of Model 2 are strictly interior.*

#### Candidate Model Comparison (Models 1–3)

| Model | Parameters ($p$) | RMSE | $R^2$ | LogLik | AIC | AICc | $\Delta\text{AICc}$ |
|---|---|---|---|---|---|---|---|
| **Model 1** — Separable Generalized Beta $\times$ Weibull | 5 | 0.1491 | 0.8379 | 9.68 | -7.36 | -0.90 | 17.08 |
| **Model 2** — Wetness-Expanded Analytis $\times$ Weibull | 6 | 0.0935 | 0.9362 | 19.01 | -24.02 | -14.69 | 3.29 |
| **Model 3** — Centered Gaussian Dynamic Breadth $\times$ Weibull | 5 | 0.0973 | 0.9310 | 18.22 | -24.44 | -17.98 | **0.00** |

---

### 8.2 Australian Field Validation Site Table (Default Run)

| Site ID | Site Name | Region | Environment Type | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual $\text{EPI}_{\text{spore}}$ | Mean Annual Episodes | Mean LP (days) | Observed Field Impact | Impact Category |
|---|---|---|---|---|---|---|---|---|---|
| `aus_site_08` | McLeod River | Cooktown Hinterland | Riparian / Watercourse | 161.00 | 160.69 | 332.12 | 13.93 | Heavy rust; total defoliation | High impact |
| `aus_site_13` | Laura Lea | Cape York Peninsula | Moist South-East Airflow | 80.09 | 80.02 | 202.84 | 10.81 | Enhanced rust; 20% mortality | High impact |
| `aus_site_02` | Lakefield NP | Cape York Peninsula | Moist Coastal Airflow | 69.30 | 69.27 | 190.92 | 10.34 | High defoliation; 20-75% mortality | High impact |
| `aus_site_01` | Chillagoe NP | Tropical North QLD | Inland / Karst | 57.53 | 57.46 | 202.60 | 12.06 | Spread to 80 km by Aug 1996 | Moderate/low |
| `aus_site_11` | Charters Towers | Inland Central QLD | Inland Hub (TWRC) | 42.50 | 42.37 | 176.92 | 13.58 | Key hub; heavy rust in wet years | High impact |
| `aus_site_03` | Wrotham Park | Gulf Region | Gulf Savanna | 41.44 | 41.39 | 154.32 | 10.68 | Rapid spread to 50 km in 8 mo | Moderate/low |
| `aus_site_10` | Rockhampton | Central QLD Coast | Subtropical Coast | 35.83 | 35.82 | 157.32 | 14.64 | Spread over 100 km in 14 mo | Moderate/low |
| `aus_site_07` | Georgetown | Etheridge Inland | Elevated Dry Inland | 23.03 | 22.88 | 86.84 | 12.23 | Moderate rust; shorter wet season | Moderate/low |
| `aus_site_04` | Rutland Plains | Gulf Coast | Gulf Lowlands | 15.42 | 15.42 | 100.88 | 9.39 | Satellite spread >20 km | Moderate/low |
| `aus_site_09` | Strathmore | Gulf / Croydon | Gulf Savanna | 11.68 | 11.63 | 52.04 | 9.65 | Spread to 45 km by Aug 1996 | Moderate/low |
| `aus_site_12` | Hughenden | Dry Inland QLD | Semi-Arid (<600mm rain) | 8.05 | 8.00 | 36.52 | 14.24 | Inadequate control in dry years | Moderate/low |
| `aus_site_05` | Inkerman | Burdekin Coast | Sand Ridge / Floodplain | 5.74 | 5.74 | 21.92 | 10.95 | Severe rust; up to 75% mortality | High impact |
| `aus_site_06` | **Delta Downs** | Gulf of Carpentaria | Coastal Gulf | **0.72** | **0.72** | **3.56** | 10.31 | Spread to 70 km in 18 mo | Moderate/low |

*Validation Note on Delta Downs: In Round 1, Delta Downs was reported as 0.82 because the mean was computed across 22 years that experienced $\ge 1$ episode ($18.00 / 22 = 0.818$). In Round 2, dividing by $N_{\text{YEARS}} = 25$ yields exactly $18.003 / 25 = \mathbf{0.7201}$, resolving the previous discrepancy and matching the independent Python benchmark (0.72).*

---

### 8.3 Brazil Monthly Climatology Table (Completed with SDs, No NA)

| Month (`start_month`) | Mean Monthly $\text{EPI}_{\text{inf}}$ | SD Monthly $\text{EPI}_{\text{inf}}$ | Mean Monthly $\text{EPI}_{\text{spore}}$ | SD Monthly $\text{EPI}_{\text{spore}}$ | Mean Latent Period (days) |
|---|---|---|---|---|---|
| 1 (Jan) | 3.3824 | 2.4333 | 1.5936 | 1.2615 | 9.59 |
| 2 (Feb) | 6.6626 | 4.3737 | 5.3682 | 3.6830 | 9.63 |
| 3 (Mar) | 11.5039 | 6.2379 | 10.2490 | 5.8100 | 9.69 |
| 4 (Apr) | **13.0159** | 6.3705 | 12.6869 | 6.3502 | 9.73 |
| 5 (May) | 11.9110 | 6.6333 | **12.7275** | 6.6620 | 9.83 |
| 6 (Jun) | 9.4277 | 5.9391 | 10.3771 | 6.1943 | 10.06 |
| 7 (Jul) | 5.3671 | 4.1445 | 7.1610 | 5.2742 | 10.09 |
| 8 (Aug) | 1.2453 | 1.0256 | 2.2016 | 1.8300 | 9.97 |
| 9 (Sep) | **0.0807** | **0.0920** | **0.3062** | **0.2836** | 9.76 |
| 10 (Oct) | 0.0106 | 0.0161 | 0.0139 | 0.0208 | 9.41 |
| 11 (Nov) | 0.0100 | 0.0146 | 0.0112 | 0.0165 | 9.33 |
| 12 (Dec) | 0.2169 | 0.2189 | 0.1374 | 0.1401 | 9.41 |

*Note: All 12 months are populated with finite standard deviations. September now reports Mean $\text{EPI}_{\text{inf}} = 0.0807 \pm 0.0920$ across all 65 cells.*

---

### 8.4 Brazil Top 10 and Bottom 10 Grid Cells (Default Run)

#### Top 10 Highest-EPI Grid Cells
| Rank | Cell ID | Latitude (°S) | Longitude (°W) | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual Episodes |
|---|---|---|---|---|---|
| 1 | `lat_-3.8_lon_-40.8` | -3.8 | -40.8 | 117.8363 | 188.52 |
| 2 | `lat_-3.9_lon_-40.4` | -3.9 | -40.4 | 117.8363 | 188.52 |
| 3 | `lat_-3.9_lon_-40.5` | -3.9 | -40.5 | 117.8363 | 188.52 |
| 4 | `lat_-4.2_lon_-40.4` | -4.2 | -40.4 | 117.8363 | 188.52 |
| 5 | `lat_-4.3_lon_-39.3` | -4.3 | -39.3 | 111.3653 | 186.84 |
| 6 | `lat_-4.7_lon_-39.5` | -4.7 | -39.5 | 111.3653 | 186.84 |
| 7 | `lat_-3.3_lon_-40.5` | -3.3 | -40.5 | 108.4403 | 197.44 |
| 8 | `lat_-3.4_lon_-40.7` | -3.4 | -40.7 | 108.4403 | 197.44 |
| 9 | `lat_-3.6_lon_-40.6` | -3.6 | -40.6 | 108.4403 | 197.44 |
| 10 | `lat_-4.0_lon_-40.1` | -4.0 | -40.1 | 105.3950 | 196.20 |

#### Bottom 10 Lowest-EPI Grid Cells
| Rank | Cell ID | Latitude (°S) | Longitude (°W) | Mean Annual $\text{EPI}_{\text{inf}}$ | Mean Annual Episodes |
|---|---|---|---|---|---|
| 56 | `lat_-3.1_lon_-40.2` | -3.1 | -40.2 | 25.0922 | 150.20 |
| 57 | `lat_-3.2_lon_-39.9` | -3.2 | -39.9 | 25.0922 | 150.20 |
| 58 | `lat_-3.2_lon_-40.1` | -3.2 | -40.1 | 25.0922 | 150.20 |
| 59 | `lat_-3.0_lon_-41.1` | -3.0 | -41.1 | 23.0388 | 154.80 |
| 60 | `lat_-3.0_lon_-41.2` | -3.0 | -41.2 | 23.0388 | 154.80 |
| 61 | `lat_-4.2_lon_-38.2` | -4.2 | -38.2 | 2.0401 | 32.92 |
| 62 | `lat_-3.6_lon_-39.0` | -3.6 | -39.0 | 0.9413 | 23.72 |
| 63 | `lat_-3.7_lon_-38.8` | -3.7 | -38.8 | 0.9413 | 23.72 |
| 64 | `lat_-4.6_lon_-37.8` | -4.6 | -37.8 | **0.0096** | **1.24** |
| 65 | `lat_-4.7_lon_-37.8` | -4.7 | -37.8 | **0.0096** | **1.24** |

*Note: Cells 64 and 65 experienced infection in only 12 of the 25 years. Dividing by 25 yields $0.00965$ (previously $0.0201$ when divided by 12).*

---

### 8.5 Comparison Across All Three Runs & Australia Episode Splitting Assessment

| Scenario / Run | Brazil Mean Annual $\text{EPI}_{\text{inf}}$ | Brazil Peak Month | Australia All Sites Mean $\text{EPI}_{\text{inf}}$ | Australia High-Impact Mean ($n=5$) | Australia Other Sites Mean ($n=8$) | Welch $t$ | Welch $\text{df}$ | Welch $P$-value | Mann–Whitney $W$ | Mann–Whitney $P$-value |
|---|---|---|---|---|---|---|---|---|---|---|
| **Default Run** (`tables_v2`) | 62.83 | Apr (13.02) | 42.49 | 71.73 | 24.21 | 1.78 | 4.6 | 0.140 | 32 | 0.093 |
| **Sensitivity (a): $P \ge 1$** (`tables_v2_precip1`) | 51.13 | May (10.13) | 41.38 | 60.71 | 29.29 | 2.08 | 5.5 | 0.087 | 34 | 0.045 |
| **Sensitivity (b): $W \le 24$** (`tables_v2_maxh24`) | 62.85 | Apr (13.02) | **42.59** | **71.91** | **24.27** | 1.78 | 4.6 | 0.140 | 32 | 0.093 |

#### Assessment: Does Run (b) Now Change Australia?
**YES.** In Round 1, script 12 called `rle(is_wet)` directly, bypassing `wet_runs()` and leaving Australian results insensitive to `MAX_EPISODE_H`. In Round 2, script 12 extracts episodes using `wet_runs(is_wet, max_h = MAX_EPISODE_H)`. As a result, wet periods longer than 24 hours are partitioned into 24-h blocks:
- **McLeod River**: increases from 160.999 to 161.576 (+0.36%)
- **Delta Downs**: increases from 0.720 to 0.759 (+5.42%)
- **Laura Lea**: increases from 80.087 to 80.201 (+0.14%)
- **Lakefield NP**: increases from 69.296 to 69.412 (+0.17%)
- **Australia High-Impact Mean**: increases from 71.73 to 71.91
- **Australia Other Sites Mean**: increases from 24.21 to 24.27
- **Australia Overall Mean**: increases from 42.49 to 42.59

---

### 8.6 Values Changing by More Than 2% from Round 1

| Metric / Parameter | Round 1 Value | Round 2 Value | Absolute Change | Percentage Change | Root Cause / Rationale |
|---|---|---|---|---|---|
| **Model 2 Parameter $c_0$** | 26.2017 | 29.1862 | +2.9845 | **+11.39%** | Expanded lower bound on $c_1$ (-10) permitted optimization to find a deeper interior minimum ($c_1 = -1.12$). |
| **Model 2 Parameter $c_1$** | -1.0000 | -1.1247 | -0.1247 | **+12.47%** | Bound relaxed from -1.0 to -10.0; parameter moved off boundary into interior. |
| **Australia Delta Downs Mean $\text{EPI}_{\text{inf}}$ (Default)** | 0.8183 | 0.7201 | -0.0982 | **-12.00%** | Normalized across all 25 years rather than 22 active episode years ($0.8183 \times 22 / 25 = 0.7201$). |
| **Australia Delta Downs Episodes/yr (Default)** | 4.05 | 3.56 | -0.49 | **-12.10%** | Divided by 25 years ($89 / 25 = 3.56$) instead of 22 years ($89 / 22 = 4.05$). |
| **Australia Delta Downs Mean $\text{EPI}_{\text{inf}}$ (Sensitivity b)** | 0.8183 | 0.7592 | -0.0591 | **-7.22%** | Divided by 25 years and augmented by 24-h episode splitting ($18.98 / 25 = 0.759$). |
| **Brazil Cells 64 & 65 Mean $\text{EPI}_{\text{inf}}$** | 0.0201 | 0.00965 | -0.01045 | **-52.00%** | Normalized across all 25 years rather than 12 active episode years ($0.0201 \times 12 / 25 = 0.00965$). |
| **Brazil Cells 64 & 65 Episodes/yr** | 2.58 | 1.24 | -1.34 | **-51.94%** | Divided by 25 years ($31 / 25 = 1.24$) instead of 12 years ($31 / 12 = 2.58$). |
| **Brazil Monthly Climatology: Sep $\text{EPI}_{\text{inf}}$** | NA | 0.0807 | Defined | **N/A (was NA)** | Zero-completion over all 65 cells populated the off-season mean and standard deviation ($0.0807 \pm 0.0920$). |
| **Brazil Monthly Climatology: Oct $\text{EPI}_{\text{inf}}$** | 0.0153 | 0.0106 | -0.0047 | **-30.72%** | Zero-completion included zero-infection cells in the regional average. |
| **Brazil Monthly Climatology: Nov $\text{EPI}_{\text{inf}}$** | 0.0130 | 0.0100 | -0.0030 | **-23.08%** | Zero-completion included zero-infection cells in the regional average. |
| **Brazil Monthly Climatology: Jul $\text{EPI}_{\text{inf}}$** | 5.5375 | 5.3671 | -0.1704 | **-3.08%** | Zero-completion included zero-infection cells in the regional average. |
| **Brazil Monthly Climatology: Aug $\text{EPI}_{\text{inf}}$** | 1.2849 | 1.2453 | -0.0396 | **-3.08%** | Zero-completion included zero-infection cells in the regional average. |
| **Brazil Monthly Climatology: Dec $\text{EPI}_{\text{inf}}$** | 0.2238 | 0.2169 | -0.0069 | **-3.08%** | Zero-completion included zero-infection cells in the regional average. |

---

## 9. Round 3 – Grouped Cross-Validation (LOWO & LOTO)

This section reports the results of running `R/03_validate_infection_models_v2.R` following the fitting of `02_infection_submodel_v2.R`. The script performs grouped leave-one-wetness-level-out (LOWO) and leave-one-temperature-level-out (LOTO) cross-validation across all candidate models, using a robust multistart optimizer (supplied initial values, full-data solution, and 10 jittered restarts per fold) to prevent fold-specific local optima.

The outputs are exported to `tables_v2/infection_model_cv_results.csv` and merged into `tables_v2/infection_model_metrics.csv`.

---

### 9.1 Cross-Validation RMSE Comparison (v2 vs. Original)

| Model | LOWO RMSE (v2) | LOWO RMSE (Original) | $\Delta$ LOWO (%) | LOTO RMSE (v2) | LOTO RMSE (Original) | $\Delta$ LOTO (%) |
|---|---|---|---|---|---|---|
| **Model 1** — Separable Generalized Beta $\times$ Weibull | **0.8492** | 0.2799 | +203.4% | **0.1996** | 0.2359 | -15.4% |
| **Model 2** — Wetness-Expanded Analytis $\times$ Weibull | **0.3388** | 0.3053 | +11.0% | **0.2248** | 0.1627 | +38.2% |
| **Model 3** — Centered Gaussian Dynamic Breadth $\times$ Weibull | **0.2369** | 0.1306 | +81.4% | **0.1698** | 0.1698 | <0.01% |
| *Model S1* — Simple Analytis Beta $\times$ Weibull | 0.3225 | — | — | 0.1996 | — | — |
| *Model S2* — Duthie Hyperbolic Secant $\times$ Weibull | 0.3188 | — | — | 0.1945 | — | — |
| *Model S3* — Fixed-Breadth Gaussian $\times$ Weibull ($H_0$) | 0.3155 | — | — | 0.1894 | — | — |

#### Notes on Methodological Differences:
1. **Multistart Optimization in v2**: The original script (`03_validate_infection_models.R`) evaluated only a single starting point per fold (`fit_safe_single`), causing fits on subsets of data to frequently settle into different local minima or flat regions. In v2, `fit_safe` uses 12 starting points per fold (supplied start, full-sample MLE, and 10 jittered points) to select the global minimum deviance solution for each training split.
2. **Model 1 LOWO Behavior**: When the highest experimental duration ($\text{dew\_h} = 24\text{ h}$) is excluded from training, Model 1 fits extrapolate steeply beyond 12 h, producing high prediction errors in that single fold ($\text{RMSE} = 1.6408$), which raises its overall LOWO RMSE to 0.8492.
3. **Model 2 Boundary Condition**: In v2, Model 2 benefits from the expanded lower bound on $c_1$ ($-10.0$ instead of $-1.0$), ensuring optimization does not artificially truncate the interaction exponent.
4. **Model 3 Performance**: Model 3 achieves the lowest grouped cross-validation error across both dimensions among candidate models: **LOWO RMSE = 0.2369** and **LOTO RMSE = 0.1698**.

---

### 9.2 Fold Failure & Convergence Audit

- **Total Folds Evaluated**: 
  - LOWO: 4 wetness levels ($\text{dew\_h} \in \{6, 8, 12, 24\}\text{ h}$) $\times$ 6 models = 24 fits
  - LOTO: 5 temperature levels ($T \in \{17, 20, 22, 25, 30\}^\circ\text{C}$) $\times$ 6 models = 30 fits
  - **Total fits performed**: 54 fits
- **Fold Convergence Status**:
  - **Failed fits**: **0 / 54 (0.0%)**
  - **Every single fold converged successfully to a valid non-linear least squares solution across all 6 models.**

---

### 9.3 Detailed Fold-by-Fold Test RMSE Breakdown

#### Leave-One-Wetness-Level-Out (LOWO) Folds
| Excluded Wetness Level | Model 1 RMSE | Model 2 RMSE | Model 3 RMSE | Model S1 RMSE | Model S2 RMSE | Model S3 RMSE |
|---|---|---|---|---|---|---|
| **Fold $W = 6\text{ h}$** | 0.0535 | 0.0512 | 0.0718 | 0.0535 | 0.0465 | 0.0494 |
| **Fold $W = 8\text{ h}$** | 0.1528 | 0.1499 | 0.1659 | 0.1528 | 0.1506 | 0.1549 |
| **Fold $W = 12\text{ h}$** | 0.4074 | 0.3150 | 0.4099 | 0.4074 | 0.3846 | 0.3878 |
| **Fold $W = 24\text{ h}$** | 1.6408 | 0.5787 | 0.1372 | 0.4468 | 0.4578 | 0.4485 |
| **Overall LOWO RMSE** | **0.8492** | **0.3388** | **0.2369** | **0.3225** | **0.3188** | **0.3155** |

*Key finding: Model 3 demonstrates superior stability when extrapolating to the 24-h wetness level ($\text{RMSE} = 0.1372$), whereas Model 1 diverges ($\text{RMSE} = 1.6408$) and Model 2 degrades ($\text{RMSE} = 0.5787$).*

#### Leave-One-Temperature-Level-Out (LOTO) Folds
| Excluded Temperature Level | Model 1 RMSE | Model 2 RMSE | Model 3 RMSE | Model S1 RMSE | Model S2 RMSE | Model S3 RMSE |
|---|---|---|---|---|---|---|
| **Fold $T = 17^\circ\text{C}$** | 0.2196 | 0.0961 | 0.1038 | 0.2196 | 0.1983 | 0.2016 |
| **Fold $T = 20^\circ\text{C}$** | 0.2657 | 0.2071 | 0.1834 | 0.2657 | 0.2520 | 0.2486 |
| **Fold $T = 22^\circ\text{C}$** | 0.1971 | 0.1803 | 0.1834 | 0.1971 | 0.1911 | 0.1843 |
| **Fold $T = 25^\circ\text{C}$** | 0.2017 | 0.1005 | 0.1002 | 0.2017 | 0.2015 | 0.1974 |
| **Fold $T = 30^\circ\text{C}$** | 0.0307 | 0.3975 | 0.2368 | 0.0307 | 0.0396 | 0.0385 |
| **Overall LOTO RMSE** | **0.1996** | **0.2248** | **0.1698** | **0.1996** | **0.1945** | **0.1894** |

---

### 9.4 Final Merged Metrics Table (`tables_v2/infection_model_metrics.csv`)

| Model | Parameters | RMSE | $R^2$ | AICc | $\Delta\text{AICc}$ | LOWO RMSE | LOTO RMSE |
|---|---|---|---|---|---|---|---|
| **Model 1** — Separable Generalized Beta $\times$ Weibull | 5 | 0.1491 | 0.8379 | -0.90 | 17.08 | 0.8492 | 0.1996 |
| **Model 2** — Wetness-Expanded Analytis $\times$ Weibull | 6 | 0.0935 | 0.9362 | -14.69 | 3.29 | 0.3388 | 0.2248 |
| **Model 3** — Gaussian wetness-dependent breadth $\times$ Weibull | 5 | 0.0973 | 0.9310 | -17.98 | **0.00** | **0.2369** | **0.1698** |



