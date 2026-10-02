# Epidemic Potential of *Maravalia cryptostegiae*, a Classical Biological Control Agent of Rubber Vine, in Northeastern Brazil

[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](LICENSE-data)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23111213-blue.svg)](https://doi.org/10.5281/zenodo.23111213)
[![Data](https://img.shields.io/badge/Data-10.5281%2Fzenodo.23111223-blue.svg)](https://doi.org/10.5281/zenodo.23111223)
[![Version](https://img.shields.io/badge/version-2.0.1-green.svg)](https://github.com/emdelponte/maravalia-epidemic-potential/releases/tag/v2.0.1)
[![R 4.4+](https://img.shields.io/badge/R-4.4%2B-blue.svg)](https://www.r-project.org/)

This repository contains the complete, reproducible computational workflow, nonlinear epidemiological models, and evaluation pipelines for the manuscript:

> **Del Ponte, E. M., de Costa, J. H., and Barreto, R. W. (2026).**  
> *Epidemic Potential of Maravalia cryptostegiae, a Classical Biological Control Agent of Rubber Vine, in Northeastern Brazil.*

The study reconstructs and validates temperature- and wetness-dependent infection and latent development submodels for the rust fungus *Maravalia cryptostegiae* (Evans & Fleureau 1993), couples them to 25-year (2001–2025) continuous hourly meteorological reanalysis from NASA POWER across 65 grid cells in Northeastern Brazil covering 73 municipal infestations of rubber vine (*Cryptostegia madagascariensis*), and validates cross-continental accuracy against 10 years of historical field release and establishment monitoring across 13 field sites in Queensland, Australia (Tomley & Evans 2004).

---

## Quick Start (One-Command Reproducibility)

Clone this repository and regenerate every table, model, and publication figure from a clean R terminal:

```bash
# 1. Clone repository
git clone https://github.com/emdelponte/maravalia-epidemic-potential.git
cd maravalia-epidemic-potential

# 2. Check and verify meteorological data integrity (or retrieve from Zenodo/NASA POWER)
Rscript scripts/get_data.R

# 3. Execute master pipeline (regenerates all models, tables, and figures)
Rscript run_all.R

# 4. Verify numerical identity against reference publication benchmarks
Rscript scripts/verify_outputs.R
```

Alternatively, use the provided `Makefile`:

```bash
make all      # Runs get_data -> analysis -> verify
make verify   # Compares outputs against reference outputs
make quick    # Fast test run (reduced bootstrap B = 100)
```

**Expected Runtime**:
- Full reproducible run from scratch (including $B = 2,000$ bootstrap iterations and sensitivity pipelines): **~4.5 minutes**.
- Quick test mode (`Rscript run_all.R --quick`): **~1.5 minutes**.

### Meteorological Data Policy & Reference Baseline

Archived NASA POWER data are the reference. The `PRECIP_THRESH = 1` sensitivity run depends on `PRECTOTCORR`, which NASA revises, so it may not be reproduced exactly from a fresh download. Main results use `RH2M` and `T2M` only, which matched exactly.

---

## Key Findings & Benchmark Metrics

The table below summarizes the key numerical benchmarks verified by `scripts/verify_outputs.R`:

| Epidemiological Benchmark | Metric / Target | Obtained Value | Verification Status |
|---|---|---|---|
| **Model 3 Optimum Temp ($T_{\text{opt}}$)** | $21.26\text{ °C}$ | **$21.2597\text{ °C}$** | PASS ($<0.001\%$) |
| **Model 3 Thermal Breadth ($\sigma_{\bar{W}}$)** | $2.280\text{ °C}$ | **$2.2798\text{ °C}$** | PASS ($<0.01\%$) |
| **Model 3 Dynamic Breadth ($\gamma$)** | $0.2984\text{ °C h}^{-1}$ | **$0.2984\text{ °C h}^{-1}$** | PASS ($<0.001\%$) |
| **Model 3 Weibull Scale ($\lambda$)** | $8.367\text{ h}$ | **$8.3666\text{ h}$** | PASS ($<0.01\%$) |
| **Model 3 Weibull Shape ($k$)** | $8.929$ | **$8.9288$** | PASS ($<0.01\%$) |
| **Bootstrap Dynamic Breadth Test ($\Delta\text{RSS}$)** | $0.2601$ | **$0.2601$** | PASS |
| **Bootstrap Empirical $P$-value** | $0.0020$ | **$0.0020$** ($4 / 2001$) | PASS ($P < 0.01$) |
| **Brazil Mean Annual $EPI_{\text{inf}}$** | $62.83$ | **$62.83$** per site/year | PASS |
| **Brazil Peak Conduciveness Window** | March–June | **$72.82\%$** of annual $EPI$ | PASS |
| **Australia: McLeod River (High Impact)** | $161.00$ | **$160.999$** | PASS |
| **Australia: Laura Lea (High Impact)** | $80.12$ | **$80.087$** | PASS |
| **Australia: Inkerman (Other Sites)** | $5.74$ | **$5.745$** | PASS |
| **Australia: Delta Downs (Other Sites)** | $0.72$ | **$0.720$** | PASS |
| **Australia Impact Discrimination (Welch $t$)** | $P = 0.140$ | **$t = 1.78, \text{df} = 4.6, P = 0.140$** | PASS |
| **Australia Impact Discrimination (Mann-Whitney)** | $P = 0.093$ | **$W = 32, P = 0.093$** | PASS |

---

## Manuscript Figure & Table Mapping

Every figure and table in the manuscript and supplement maps to an explicit script and output file:

| Manuscript Item | Description | Generating Script | Output Artifact |
|---|---|---|---|
| **Figure 1** | Infection submodel response curves and surface | `R/02_infection_submodel.R` | `outputs/figures/Fig1_infection_submodel.{png,pdf}` |
| **Figure 2** | Latent development rate and temperature curves | `R/05_latent_submodel.R` | `outputs/figures/Fig2_latent_period.{png,pdf}` |
| **Figure 3** | Monthly $EPI_{\text{inf}}$ climatology and latent periods | `R/09_epidemic_pressure_indices.R` | `outputs/figures/Fig3_seasonal_climatology.{png,pdf}` |
| **Figure 4** | Spatial risk map and weed area correlation | `R/10_seasonal_spatial_analysis.R` | `outputs/figures/Fig4_spatial_maps.{png,pdf}` |
| **Figure 5** | Australian 10-year validation and site rankings | `R/12_australian_validation_pipeline.R` | `outputs/figures/Fig5_australia.{png,pdf}` |
| **Figure S1** | Parametric bootstrap null distribution of $\Delta\text{RSS}$ | `R/04_bootstrap_gamma.R` | `outputs/figures/FigS1_bootstrap_gamma.{png,pdf}` |
| **Figure S2** | Comparison of candidate infection models across wetness | `R/02_infection_submodel.R` | `outputs/figures/FigS2_model_comparison.{png,pdf}` |
| **Figure S3** | Infection submodel residual diagnostics | `R/02_infection_submodel.R` | `outputs/figures/FigS3_residuals.{png,pdf}` |
| **Figure S4** | Digitization sensitivity perturbation envelopes | `R/06_digitization_sensitivity.R` | `outputs/figures/FigS4_digitization_sensitivity.{png,pdf}` |
| **Table 1** | Infection submodel comparison and information criteria | `R/02_infection_submodel.R`, `R/03_validate_infection_models.R` | `outputs/tables/infection_model_metrics.csv` |
| **Table 2** | Parameter estimates for candidate infection models | `R/02_infection_submodel.R` | `outputs/tables/infection_model_parameters.csv` |
| **Table 3** | Latent development model metrics and parameters | `R/05_latent_submodel.R` | `outputs/tables/latent_model_metrics.csv` |
| **Table 4** | 25-year climatological ranking of 65 Brazilian sites | `R/09_epidemic_pressure_indices.R` | `outputs/tables/site_epidemic_rankings.csv` |
| **Table 5** | Australian field release site validation summary | `R/12_australian_validation_pipeline.R` | `outputs/tables/australian_model_validation_summary.csv` |
| **Table S1** | Grouped cross-validation (LOWO and LOTO) results | `R/03_validate_infection_models.R` | `outputs/tables/infection_model_cv_results.csv` |
| **Table S2** | Relative humidity threshold sensitivity comparison | `R/07_wetness_episodes.R` | `outputs/tables/rh_sensitivity_comparison.csv` |
| **Table S3** | Predicted epidemic potential for 73 municipal records | `R/09_epidemic_pressure_indices.R` | `data/metadata/tableS3_v2.csv` |

---

## Repository Structure

```
maravalia-epidemic-potential/
├── README.md                      # Project documentation, quick-start, and mapping
├── LICENSE                        # MIT License for software code
├── LICENSE-data                   # Creative Commons Attribution 4.0 for data & tables
├── CITATION.cff                   # Citation metadata (Del Ponte et al. 2026)
├── Makefile                       # Workflow targets (all, data, analysis, figures, verify)
├── run_all.R                      # Master end-to-end reproducible pipeline driver
├── renv.lock                      # Exact R package lockfile for reproduction
├── sessionInfo.txt                # Full R environment session specifications
├── requirements.txt               # Python dependencies for independent cross-checks
├── .gitignore                     # Excludes raw NASA POWER RDS files and OS artifacts
├── config/
│   └── settings.R                 # Global constants, thresholds (RH=90%, W_CAP=24), and seeds
├── R/
│   ├── maravalia_functions.R      # Shared core functions (wetness, Model 3 infection, latent rate)
│   ├── 00_download_nasapower_brazil.R      # Automated NASA POWER retrieval for Brazil (newly written)
│   ├── 01_qc_and_data_prep.R               # Structural QC and digitized experimental dataset prep
│   ├── 02_infection_submodel.R             # Model fitting, AICc, Fig 1, Fig S2, Fig S3
│   ├── 03_validate_infection_models.R      # LOWO and LOTO cross-validation
│   ├── 04_bootstrap_gamma.R                # Parametric bootstrap test for dynamic breadth, Fig S1
│   ├── 05_latent_submodel.R                # Latent rate model fitting, Fig 2
│   ├── 06_digitization_sensitivity.R       # Stochastic perturbation sensitivity analysis, Fig S4
│   ├── 07_wetness_episodes.R               # Continuous wetness episode extraction across Brazil
│   ├── 08_latent_tracking.R                # Latent development cohort tracking across Brazil
│   ├── 09_epidemic_pressure_indices.R      # EPI calculation, rankings table, Fig 3
│   ├── 10_seasonal_spatial_analysis.R      # Spatial risk mapping, host patch correlation, Fig 4
│   ├── 11_download_australian_nasapower.R  # NASA POWER retrieval for 13 Australian sites
│   └── 12_australian_validation_pipeline.R # Australian 10-year validation, t-tests, Fig 5
├── data/
│   ├── raw/
│   │   ├── MANIFEST.csv           # File registry with SHA-256 checksums and record counts
│   │   ├── nasapower_maravalia/   # 65 Brazilian hourly RDS weather files (gitignored)
│   │   ├── nasapower_australia/   # 13 Australian hourly RDS weather files (gitignored)
│   │   └── nasapower_subset_example/ # 2-site minimal example for rapid local testing
│   ├── derived/                   # Derived caches (episodes, cohorts, bootstrap samples)
│   └── metadata/                  # Digitized experimental inputs and geographic coordinates
├── outputs/
│   ├── tables/                    # Regenerated publication CSV tables
│   ├── figures/                   # Regenerated 300 DPI PNG and vector PDF publication figures
│   ├── models/                    # Fitted RDS model objects
│   └── logs/                      # Execution runtime logs
├── reference_outputs/             # Frozen benchmark tables and figures for verify_outputs.R
├── scripts/
│   ├── get_data.R                 # Checks local data existence and verifies SHA-256 integrity
│   └── verify_outputs.R           # Automated audit suite asserting numerical table identity
├── python/                        # Independent Python cross-checks (scenarios, host localities)
└── docs/
    ├── DATA_DICTIONARY.md         # Full column definitions and units
    ├── METHODS_SUMMARY.md         # Mathematical equations, thresholds, and v2 corrections
    ├── MECHANICAL_CHANGES.md      # Detailed audit log of all mechanical script modifications
    ├── OPEN_ISSUES.md             # Documented open scientific questions and editorial items
    ├── FILE_INVENTORY.csv         # Full 410-file inventory and classification log
    ├── CHANGES_v2.md              # Original audit log from review (2026-09-28)
    ├── RUN_REPORT_v2.md           # Benchmark execution report from review
    └── FIGURES_v3.md              # Specification report for final manuscript figures
```

---

## Data Policy & Zenodo Archive

- **Raw Hourly Weather Data**: The complete 25-year (2001–2025) continuous hourly weather record consists of 65 grid cells for Northeastern Brazil (~106 MB) and 13 monitoring stations for Queensland, Australia (~22 MB). Because these files exceed standard Git repository best-practice sizes, they are excluded via `.gitignore`.
- **Integrity Validation**: All 78 raw files are cataloged in [`data/raw/MANIFEST.csv`](data/raw/MANIFEST.csv) with their exact file sizes, record counts, date ranges, and SHA-256 cryptographic checksums.
- **Archived Copy**: An immutable archive of the complete raw dataset is deposited on Zenodo with DOI: [`10.5281/zenodo.23111223`](https://doi.org/10.5281/zenodo.23111223).
- **Note on ERA5 Exclusion**: Early preliminary analyses explored ERA5-Land reanalysis (`era5_maravalia/`, `download_era5_maravalia.R`). Due to spatial interpolation artifacts and lower relative humidity fidelity in coastal semi-arid zones, ERA5 was excluded from the final research and superseded entirely by 0.5° NASA POWER hourly meteorological reanalysis.

---

## Citation & Authorship

To cite this computational repository or its findings:

```bibtex
@article{delponte2026maravalia,
  author    = {Del Ponte, Emerson M. and de Costa, J. H. and Barreto, Robert W.},
  title     = {Epidemic Potential of \textit{Maravalia cryptostegiae}, a Classical Biological Control Agent of Rubber Vine, in Northeastern Brazil},
  year      = {2026},
  version   = {2.0.1},
  doi       = {10.5281/zenodo.23111213},
  url       = {https://github.com/emdelponte/maravalia-epidemic-potential}
}
```

---

## Licenses

- **Code**: [MIT License](LICENSE) applies to all R, Python, and Makefile code.
- **Data & Tables**: [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-data) applies to all derived data, summaries, and documentation.
