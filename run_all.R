#!/usr/bin/env Rscript
# ==============================================================================
# run_all.R: Master Reproducibility Pipeline
# Manuscript: Epidemic Potential of Maravalia cryptostegiae, a Classical Biological
#             Control Agent of Rubber Vine, in Northeastern Brazil
# Authors: Emerson M. Del Ponte, J. S. de Costa, R. W. Barreto (2026)
# ==============================================================================

t_start <- Sys.time()

# Parse Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args
sensitivity_mode <- "--sensitivity" %in% args
help_mode <- any(c("-h", "--help") %in% args)

if (help_mode) {
  cat("Usage: Rscript run_all.R [options]\n\n")
  cat("Options:\n")
  cat("  --quick        Run in fast testing mode (reduced bootstrap B = 100, quick verification)\n")
  cat("  --sensitivity  Run sensitivity pipeline scenarios (PRECIP_THRESH = 1.0 mm/h and MAX_EPISODE_H = 24 h)\n")
  cat("  -h, --help     Show this help message and exit\n\n")
  quit(status = 0)
}

PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "config/settings.R"))

log_dir <- file.path(PROJ_ROOT, DIR_OUTPUT_LOGS)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

master_log <- file.path(log_dir, if (quick_mode) "run_all_quick.log" else "run_all.log")
sink_conn <- file(master_log, open = "wt")
sink(sink_conn, type = "output", split = TRUE)
sink(sink_conn, type = "message")

cat("==============================================================================\n")
cat("MARAVALIA CRYPTOSTEGIAE EPIDEMIOLOGICAL MODELING PIPELINE\n")
cat("Execution Started:", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project Root: . (relative to repository root)\n")
cat("R Version:", R.version.string, "\n")
cat("Mode:", if (quick_mode) "QUICK TEST (B=100)" else if (sensitivity_mode) "FULL RUN WITH SENSITIVITY PIPELINES (B=2000)" else "FULL REPRODUCIBLE RUN (B=2000)", "\n")
cat("Caching: DISABLED (All steps run from scratch)\n")
cat("==============================================================================\n\n")

set.seed(SEED_BOOTSTRAP)

run_step <- function(step_num, step_name, script_path) {
  cat(sprintf(">>> [%s] %s (%s)...\n", step_num, step_name, basename(script_path)))
  cat("    Execution Status: EXECUTING FROM SCRATCH (no cache)\n")
  t0 <- Sys.time()
  
  # Run script with clean environment
  env <- new.env(parent = globalenv())
  env$PROJ_ROOT <- PROJ_ROOT
  env$quick_mode <- quick_mode
  env$recompute <- TRUE
  
  status <- tryCatch({
    sys.source(script_path, envir = env)
    "OK"
  }, error = function(e) {
    cat("    [ERROR]:", conditionMessage(e), "\n")
    "FAILED"
  })
  
  t1 <- Sys.time()
  dur <- round(difftime(t1, t0, units = "secs"), 1)
  cat(sprintf("    Status: %s | Elapsed: %s s\n\n", status, dur))
  if (status == "FAILED") {
    warning("Step ", step_num, " failed: ", script_path, call. = FALSE)
  }
  status == "OK"
}

# --- Pipeline Sequence ---
# 1. Structural QC & Digitized Data Preparation
run_step("1/10", "Structural QC & Data Prep", file.path(PROJ_ROOT, "R/01_qc_and_data_prep.R"))

# 2. Infection Submodel Fitting, Selection & Info Criteria
run_step("2/10", "Infection Submodel Fitting (Fig 1, S2, S3)", file.path(PROJ_ROOT, "R/02_infection_submodel.R"))

# 3. Grouped Cross-Validation (LOWO & LOTO)
run_step("3/10", "Grouped Cross-Validation (LOWO/LOTO)", file.path(PROJ_ROOT, "R/03_validate_infection_models.R"))

# 4. Parametric Bootstrap Test for Gamma
run_step("4/10", "Parametric Bootstrap for Gamma (Fig S1)", file.path(PROJ_ROOT, "R/04_bootstrap_gamma.R"))

# 5. Latent Development Submodel Fitting
run_step("5/10", "Latent Submodel Fitting (Fig 2)", file.path(PROJ_ROOT, "R/05_latent_submodel.R"))

# 6. Digitization Sensitivity Analysis
run_step("6/10", "Digitization Sensitivity Analysis (Fig S4)", file.path(PROJ_ROOT, "R/06_digitization_sensitivity.R"))

# 7. Wetness Episodes Extraction across Northeastern Brazil
run_step("7/10", "Wetness Episodes Extraction (Brazil)", file.path(PROJ_ROOT, "R/07_wetness_episodes.R"))

# 8. Hourly Latent Development Cohort Tracking
run_step("8/10", "Latent Cohort Tracking (Brazil)", file.path(PROJ_ROOT, "R/08_latent_tracking.R"))

# 9. Epidemic Pressure Indices & Climatology
run_step("9/10", "Epidemic Pressure Indices (Fig 3)", file.path(PROJ_ROOT, "R/09_epidemic_pressure_indices.R"))

# 10. Spatial Risk Mapping & Host Correlation
run_step("10/10a", "Spatial Risk Mapping & Host Correlation (Fig 4)", file.path(PROJ_ROOT, "R/10_seasonal_spatial_analysis.R"))

# 11. Australian Field Validation
run_step("10/10b", "Australian Field Validation (Fig 5)", file.path(PROJ_ROOT, "R/12_australian_validation_pipeline.R"))

# Sensitivity Runs if Requested
if (sensitivity_mode) {
  cat("\n==============================================================================\n")
  cat("RUNNING SENSITIVITY PIPELINES (PRECIP_THRESH = 1.0 & MAX_EPISODE_H = 24)\n")
  cat("==============================================================================\n\n")
  
  # Sensitivity 1: PRECIP_THRESH = 1.0 mm/h
  cat(">>> Sensitivity Run 1: PRECIP_THRESH = 1.0 mm/h (Target: outputs/tables_v2_precip1/)\n")
  out_p1 <- file.path(PROJ_ROOT, "outputs/tables_v2_precip1")
  dir.create(out_p1, recursive = TRUE, showWarnings = FALSE)
  for (f in c("all_infection_models_metrics.csv", "bootstrap_gamma_results.csv", 
              "infection_model_metrics.csv", "infection_model_parameters.csv")) {
    src_f <- file.path(PROJ_ROOT, "outputs/tables", f)
    if (file.exists(src_f)) file.copy(src_f, file.path(out_p1, f), overwrite = TRUE)
  }
  PRECIP_THRESH <<- 1.0
  MAX_EPISODE_H <<- NA
  DIR_OUTPUT_TABLES <<- "outputs/tables_v2_precip1"
  DIR_DATA_DERIVED <<- "data/derived_precip1"
  source(file.path(PROJ_ROOT, "R/07_wetness_episodes.R"))
  source(file.path(PROJ_ROOT, "R/08_latent_tracking.R"))
  source(file.path(PROJ_ROOT, "R/09_epidemic_pressure_indices.R"))
  source(file.path(PROJ_ROOT, "R/12_australian_validation_pipeline.R"))
  
  # Sensitivity 2: MAX_EPISODE_H = 24 h
  cat("\n>>> Sensitivity Run 2: MAX_EPISODE_H = 24 h (Target: outputs/tables_v2_maxh24/)\n")
  out_m24 <- file.path(PROJ_ROOT, "outputs/tables_v2_maxh24")
  dir.create(out_m24, recursive = TRUE, showWarnings = FALSE)
  for (f in c("all_infection_models_metrics.csv", "bootstrap_gamma_results.csv", 
              "infection_model_metrics.csv", "infection_model_parameters.csv")) {
    src_f <- file.path(PROJ_ROOT, "outputs/tables", f)
    if (file.exists(src_f)) file.copy(src_f, file.path(out_m24, f), overwrite = TRUE)
  }
  PRECIP_THRESH <<- NA
  MAX_EPISODE_H <<- 24
  DIR_OUTPUT_TABLES <<- "outputs/tables_v2_maxh24"
  DIR_DATA_DERIVED <<- "data/derived_maxh24"
  source(file.path(PROJ_ROOT, "R/07_wetness_episodes.R"))
  source(file.path(PROJ_ROOT, "R/08_latent_tracking.R"))
  source(file.path(PROJ_ROOT, "R/09_epidemic_pressure_indices.R"))
  source(file.path(PROJ_ROOT, "R/12_australian_validation_pipeline.R"))

  # Restore default settings
  PRECIP_THRESH <<- NA
  MAX_EPISODE_H <<- NA
  DIR_OUTPUT_TABLES <<- "outputs/tables"
  DIR_DATA_DERIVED <<- "data/derived"
}

t_end <- Sys.time()
cat("==============================================================================\n")
cat("PIPELINE EXECUTION COMPLETE\n")
cat("Finished At:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(sprintf("Total Runtime: %.2f minutes\n", difftime(t_end, t_start, units = "mins")))
cat("All outputs generated in outputs/tables/ and outputs/figures/.\n")
cat("Run 'Rscript scripts/verify_outputs.R' to validate numerical identity.\n")
cat("==============================================================================\n")

sink(type = "message")
sink(type = "output")
close(sink_conn)
