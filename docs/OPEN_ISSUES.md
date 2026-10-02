# Open Scientific & Technical Issues

This document catalogs open questions, pending empirical validations, and editorial items identified during the review and reproducibility audit. In accordance with the project guidelines, **no scientific logic was altered** during repository packaging; these items are documented here for the authors' review and future revisions.

---

## 1. INMET Ground Station Humidity Validation for Coastal/Low-EPI Cells

- **Status**: Pending empirical validation.
- **Description**: Several coastal lowland municipalities (e.g., Aracati, Caucaia, Cascavel in Ceará) exhibit surprisingly low predicted annual infection potential ($EPI_{\text{inf}} < 3.0$), in contrast to interior upland stations ($EPI_{\text{inf}} > 100$). A key hypothesis is that NASA POWER 0.5° $\times$ 0.5° satellite-model reanalysis relative humidity (`RH2M`) may underestimate nocturnal boundary-layer humidity or sea-breeze moisture in coastal littoral zones due to spatial averaging over ocean/land grid cells.
- **Recommended Action**: Validate NASA POWER hourly relative humidity and dew point depression against ground-based automatic weather station records from INMET (Instituto Nacional de Meteorologia) for Fortaleza, Aracati, and nearby coastal stations (2001–2025).

---

## 2. Magnitude of Perturbation Noise in Digitization Sensitivity Analysis

- **Status**: Methodological consideration for sensitivity testing.
- **Description**: In `R/06_digitization_sensitivity.R`, stochastic perturbation of digitized infection data was implemented by adding uniform additive noise $\delta \sim \text{Uniform}(-0.05, +0.05)\text{ pustules cm}^{-2}$. Because observed pustule densities range up to $14.35\text{ pustules cm}^{-2}$, an additive perturbation of $\pm 0.05$ represents $<0.4\%$ noise at high pustule counts, although it represents a larger fraction at the lowest non-zero densities ($0.05\text{ pustules cm}^{-2}$).
- **Recommended Action**: In future sensitivity runs, evaluate proportional (multiplicative) perturbation noise (e.g., $\pm 3\%$ to $\pm 5\%$ relative noise across all data points) to confirm that Model 3 parameter estimates ($\gamma$, $\sigma_{\bar{W}}$, $T_{\text{opt}}$) remain robust across the entire dynamic range.

---

## 3. Supplementary Figures S2–S4 Missing Captions

- **Status**: Editorial / Manuscript Supplement.
- **Description**: While Supplementary Figure S1 ($H_0$ bootstrap distribution of $\Delta\text{RSS}$) is fully documented in the manuscript text, detailed figure captions for Supplementary Figures S2 (Multi-Model Comparison across Wetness Slices), S3 (Infection Submodel Residual Diagnostics), and S4 (Digitization Sensitivity Confidence Envelopes) are currently missing from the supplement document.
- **Recommended Action**: Add formal captions defining panels A–D, axes, and model legends to the Supplementary Information document before final submission.

---

## 4. Duplicate Grid Coordinates & Spatial Non-Independence in Table S3

- **Status**: Spatial aggregation property.
- **Description**: Table S3 (`data/metadata/tableS3_v2.csv` / `outputs/tables/site_epidemic_rankings.csv`) reports predicted annual epidemic pressure indices across 73 municipal occurrence records from Bonilla et al. (2023). However, multiple municipalities share identical coordinates because they fall within the same 0.5° $\times$ 0.5° NASA POWER grid cell (for example, Acaraú, Bela Cruz, and Cruz all map to `lat_-2.9_lon_-40.1`, yielding identical $EPI_{\text{inf}} = 25.09$).
- **Impact**: In regression analyses of infested patch area versus climatic suitability, treating all 73 municipal records as independent observations introduces spatial pseudoreplication. In `R/10_seasonal_spatial_analysis.R`, this was addressed by reporting both locality-level ($n = 73$, Pearson $r = -0.487$) and cell-level aggregated ($n = 28$ unique weather cells, Pearson $r = -0.615$) statistics.
- **Recommended Action**: Ensure the manuscript text explicitly notes that 73 municipal records map to 28 unique meteorological grid cells, and highlight the cell-aggregated correlation ($r = -0.615, P = 0.0005$) as the primary spatial test.

---

## 5. Discrepancy in Digitization Noise Labeling in `FIGURES_v3.md`

- **Status**: Documentation error in revision notes.
- **Description**: In [`docs/FIGURES_v3.md`](FIGURES_v3.md) (Table 1, row `FigS4_digitization_sensitivity`), the panel description states: *"Parameter and curve sensitivity under $\pm 5\%$ digitization noise ($B=500$)"*. However, the underlying analysis script (`R/06_digitization_sensitivity.R`, line 46) implements additive noise of $\pm 0.05\text{ pustules cm}^{-2}$ (`noise <- runif(nrow(inf_data), min = -0.05, max = 0.05)`).
- **Recommended Action**: Correct the caption and documentation to clarify that the analysis used additive noise of $\pm 0.05\text{ pustules cm}^{-2}$, or re-run the sensitivity analysis under explicit $\pm 5\%$ multiplicative noise.

---

## 6. Author Names, Initials, and ORCIDs Confirmation

- **Status**: Editorial / Metadata confirmation.
- **Description**: In `CITATION.cff`, authors are recorded as "Del Ponte, Emerson M." (ORCID: 0000-0003-4398-409X), "de Costa, J. S.", and "Barreto, R. W.". Author initials and ORCIDs must be confirmed by the co-authors prior to final publication and Zenodo DOI minting. The DOI field is currently left as a placeholder (`10.5281/zenodo.placeholder`).
- **Recommended Action**: Co-authors should confirm their exact preferred name formatting, initials, and ORCID identifiers for the final Zenodo release.
