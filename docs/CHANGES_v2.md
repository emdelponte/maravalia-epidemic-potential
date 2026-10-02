# v2 pipeline — changes from the review (2026-09-28)

    # Set working directory to project root:
    setwd(".")
    source("R/run_phase2_v2.R")

Outputs go to `data_v2/`, `tables_v2/`, `figures_v2/`, `models_v2/`. Original outputs are untouched.
All settings live in `R/maravalia_functions_v2.R`.

| Script | Change |
|---|---|
| maravalia_functions_v2.R (new) | One wetness rule, infection and latent functions shared by Brazil and Australia |
| 07_wetness_episodes_v2 | Wet hour = RH >= 90 (precipitation clause off; `P > 0` was true in 50–74% of hours). sigma(W) capped at W = 24 h. Optional splitting of long episodes (`MAX_EPISODE_H`) |
| 08_latent_tracking_v2 | Latent rate held at r(27 °C) above 27 °C (`LAT_T_PLATEAU`) |
| 09, 10 _v2 | Paths only |
| 12_australian_validation_pipeline_v2 | Same functions as Brazil; reports Welch t/df/P and Mann–Whitney; group label corrected ("other sites", not "dry inland") |
| 04_bootstrap_gamma_v2 | H1 refits start from the H0 solution (gamma = 0) + extra starts, so Delta RSS >= 0; H0 bounds match H1 (nested); P = (x+1)/(B+1); new cache in data_v2/ |
| 02_infection_submodel_v2 | Model 2 c0 upper bound 20 -> 200 (was on the bound); AIC = -2logL + 2(p+1), consistent with the LogLik column |

Sensitivity runs: change `RH_THRESH`, `PRECIP_THRESH` (e.g. 1), `MAX_EPISODE_H` (e.g. 24) or `LAT_T_PLATEAU` in the functions file and rerun from 07.

Expected (from the Python check, `epi_scenarios.py`): Brazil mean annual EPI_inf about 63 (range 0–118 across cells), peak March–June;
Australia peak Jan–Apr; McLeod River about 161, Inkerman about 6, Delta Downs about 0.7.

Not changed (manuscript/decision items): EPI_spore definition; Model 1 `a` unidentified; kappa (k) unidentified; the σ >= 1 floor;
digitization noise size; spatial non-independence of municipalities sharing a grid cell.

## Round 2 fixes (after RUN_REPORT_v2)
- 09 and 12: annual means were averaged only over years that had at least one infection, which inflated dry sites (Delta Downs 0.82 instead of 0.72). They now divide by `N_YEARS` (25). Monthly climatologies are completed with zeros for cell/site × month combinations that had no infection (the Sep SD was NA).
- 12: episodes now come from the shared `wet_runs()`, so `MAX_EPISODE_H` also applies to Australia.
- 02: Model 2 c1 had moved onto its lower bound (−1); the bound is now −10 and there are extra starting values.
