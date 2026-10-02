# Manuscript Figures Revision v3 Report (`figures_v3/`)
**Date:** 2026-09-28  
**Scope:** Removal of all sporulation pressure ($EPI_{\text{spore}}$) lines, panels, legends, and labels from manuscript figures. Creation of publication-ready 300 DPI PNGs and vector/raster PDFs in `figures_v3/`. Independent verification of table calculations against `tables_v2/` via `review_2026-09-28/v3_check/`.

---

## 1. Summary of Figures in `figures_v3/`

| File Name (PNG & PDF) | Manuscript Fig. | Panels Description | Source Script | New vs Copied |
|---|---|---|---|---|
| `Fig1_infection_submodel` | Figure 1 | (A) Raw pustule density vs temp by dew duration; (B) Model 3 fitted infection potential vs temp; (C) Model 3 fitted response vs wetness duration; (D) Infection response surface contour. | `R/02_infection_submodel_v2.R` | Copied from `figures_v2/fig1_infection_submodel_main.png` |
| `Fig2_latent_period` | Figure 2 | (A) Observed latent period vs temp; (B) Development rate regression $r(T)$; (C) Fitted back-transformed latent curves; (D) Residuals vs fitted latent period. | `R/05_latent_submodel.R` | Copied from `figures/fig3_combined_latent_submodel.png` |
| `Fig3_seasonal_climatology` | Figure 3 | (A) Mean monthly $EPI_{\text{inf}}$ across 65 grid cells with $\pm\text{SD}$ error bars (no sporulation line, no legend); (B) Latent period duration (days) boxplots by month of infection initiation. | `R/09_epidemic_pressure_indices_v3.R` | **New** (regenerated without $EPI_{\text{spore}}$) |
| `Fig4_spatial_maps` | Figure 4 | (A) Spatial map of mean annual $EPI_{\text{inf}}$ across Northeastern Brazil (viridis scale 0–120) with 73 host localities and visible state labels (CE, RN, PB, PI); (B) $\log_{10}(\text{infested area, ha})$ vs $EPI_{\text{inf}}$ linear trend with dual correlation annotations. | `R/10_seasonal_spatial_analysis_v3.R` | **New** (regenerated without $EPI_{\text{spore}}$ map) |
| `Fig5_australia` | Figure 5 | (A) Queensland monthly mean $EPI_{\text{inf}}$ bars only; (B) Australian site ranking bar chart colored by "High impact" vs "Other sites"; (C) Queensland spatial map with release sites sized/colored by $EPI_{\text{inf}}$. | `R/12_australian_validation_pipeline_v3.R` | **New** (regenerated without $EPI_{\text{spore}}$, updated legend) |
| `FigS1_bootstrap_gamma` | Figure S1 | Parametric bootstrap null distribution of $\Delta\text{RSS}$ ($B=2,000$), 95th percentile cut-off, and observed $\Delta\text{RSS}$ annotated with $P = 0.0020$. | `R/04_bootstrap_gamma_v2.R` | Copied from `figures_v2/figS3_bootstrap_gamma_distribution.png` |
| `FigS2_model_comparison` | Figure S2 | Multi-panel comparison of Model 1, Model 2, Model 3, and Duthie against observed pustule infection data across wetness duration slices ($6, 8, 12, 24\text{ h}$). | `R/02_infection_submodel_v2.R` | Copied from `figures_v2/figS1_all_infection_models_comparison.png` |
| `FigS3_residuals` | Figure S3 | Infection submodel residual diagnostics for (A) Model 1, (B) Model 2, and (C) Model 3 vs fitted response. | `R/02_infection_submodel_v2.R` | Copied from `figures_v2/figS2_infection_residuals_diagnostics.png` |
| `FigS4_digitization_sensitivity` | Figure S4 | Parameter and curve sensitivity under $\pm 5\%$ digitization noise ($B=500$) for (A) infection submodel and (B) latent period submodel. | `R/06_digitization_sensitivity.R` | Copied from `figures/fig4_combined_digitization_sensitivity.png` |

---

## 2. Correlation Values Printed on Fig 4B

The correlation between recorded rubber vine patch size ($\log_{10}$-transformed infested area in ha) and nearest weather grid cell annual infection potential ($EPI_{\text{inf}}$) was evaluated at two spatial resolutions:

- **Locality Level ($n = 73$ municipal records):**
  - **Pearson $r$:** $-0.487$ ($P < 0.001$, exact $P = 1.40 \times 10^{-5}$)
  - **Spearman $\rho$:** $-0.466$ ($P < 0.001$, exact $P = 3.96 \times 10^{-5}$)
- **Grid Cell Level ($n = 28$ unique weather cells, mean $\log_{10}$ area per cell):**
  - **Pearson $r$:** $-0.615$ ($P < 0.001$, exact $P = 0.0005$)
  - **Spearman $\rho$:** $-0.622$ ($P < 0.001$, exact $P = 0.0004$)

---

## 3. Table-Identity Verification (`v3_check/` vs `tables_v2/`)

To confirm that the `_v3.R` pipeline did not alter any underlying calculations, model parameters, wetness rules, or summary statistics, all data outputs were redirected to `review_2026-09-28/v3_check/` and evaluated against their `tables_v2/` counterparts using `all.equal()` in R:

```r
site_epidemic_rankings.csv:              TRUE (Identical)
monthly_climatology_summary.csv:         TRUE (Identical)
australian_model_validation_summary.csv: TRUE (Identical)
phase2_annual_site_metrics.csv:          TRUE (Identical to data_v2/)
```

**Result:** Zero numerical or structural discrepancies detected across all exported tables.

---

## 4. Code Fixes & Implementation Details

1. **`R/09_epidemic_pressure_indices_v3.R`**:
   - Redirected all table exports (`phase2_annual_site_metrics.csv`, `site_epidemic_rankings.csv`, `monthly_climatology_summary.csv`) to `review_2026-09-28/v3_check/`.
   - In Panel A, removed the sporulation line/points and legend entirely. Formatted y-axis label using `expression("Mean monthly " * EPI[inf])`.
   - Exported `Fig3_seasonal_climatology.png` (300 DPI) and `Fig3_seasonal_climatology.pdf` to `figures_v3/`.

2. **`R/10_seasonal_spatial_analysis_v3.R`**:
   - Removed the $EPI_{\text{spore}}$ spatial map panel.
   - Updated Panel A color scale to sequential viridis with fixed limits `c(0, 120)` and restricted state abbreviation labels to only those within the map view extent (`CE`, `RN`, `PB`, `PI`).
   - In Panel B, rendered the linear regression of $\log_{10}(\text{infested area})$ vs $EPI_{\text{inf}}$ with annotations for both locality ($n=73$) and grid-cell ($n=28$) Pearson $r$ and Spearman $\rho$.
   - Utilized `cairo_pdf` in `ggsave()` to ensure that the Greek symbol $\rho$ (U+03C1) renders cleanly without font conversion warnings.
   - Exported `Fig4_spatial_maps.png` (300 DPI) and `Fig4_spatial_maps.pdf` to `figures_v3/`.

3. **`R/12_australian_validation_pipeline_v3.R`**:
   - Redirected `australian_model_validation_summary.csv` to `review_2026-09-28/v3_check/`.
   - In Panel A, removed the sporulation pressure line and point markers, retaining only the monthly mean $EPI_{\text{inf}}$ bar chart.
   - In Panel B, recoded impact categories from `"Moderate/low impact"` to `"Other sites"` to ensure the legend displays `"High impact"` vs `"Other sites"`.
   - In Panel C, maintained the Queensland release site spatial map with $EPI_{\text{inf}}$ scaling, with no references to sporulation.
   - Exported `Fig5_australia.png` (300 DPI) and `Fig5_australia.pdf` to `figures_v3/`.

---

## 5. Verification Checks

- **Grep Inspection:** Running `grep -n -i "spore"` on `R/09_epidemic_pressure_indices_v3.R`, `R/10_seasonal_spatial_analysis_v3.R`, and `R/12_australian_validation_pipeline_v3.R` confirmed that `spore` appears only in internal calculation variables and script header comments. Not a single ggplot layer, legend, scale, or label contains `spore` or `sporulation`.
- **Visual Inspection:** Every PNG in `figures_v3/` was opened and verified. No panel, axis, or legend mentions sporulation or $EPI_{\text{spore}}$.
- **FigS1 Verification:** `figures_v3/FigS1_bootstrap_gamma.png` was inspected; its annotation correctly displays `Observed ΔRSS = 0.2601 (P = 0.0020)`.
- **DPI Verification:** All 9 PNG files in `figures_v3/` have resolution verified at exactly 300.000 DPI.
