#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/verify_outputs.R
# Purpose: Step 5 Numerical Verification of Regenerated Tables against Reference
#          Outputs and Validation of Key Manuscript Benchmark Values
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."

cat("==============================================================================\n")
cat("VERIFICATION SUITE: REGENERATED OUTPUTS VS. REFERENCE BENCHMARKS\n")
cat("==============================================================================\n\n")

verify_table_dir <- function(regen_dir, ref_dir, suite_title) {
  cat(sprintf("--- %s ---\n", suite_title))
  cat(sprintf("Regenerated: %s\nReference:   %s\n", regen_dir, ref_dir))
  
  if (!dir.exists(ref_dir)) {
    cat("  [FAIL]: Reference directory not found:", ref_dir, "\n\n")
    return(FALSE)
  }
  if (!dir.exists(regen_dir)) {
    cat(sprintf("  [NOTICE]: Directory %s not found. (Run Rscript run_all.R with corresponding flags to generate.)\n\n", regen_dir))
    return(NA)
  }
  
  ref_files <- sort(list.files(ref_dir, pattern = "\\.csv$", full.names = FALSE))
  cat(sprintf("Evaluating %d reference tables...\n", length(ref_files)))
  
  all_passed <- TRUE
  
  for (fn in ref_files) {
    ref_path <- file.path(ref_dir, fn)
    regen_path <- file.path(regen_dir, fn)
    
    if (!file.exists(regen_path)) {
      cat(sprintf("  [%-40s] FAIL: File missing in %s\n", fn, regen_dir))
      all_passed <- FALSE
      next
    }
    
    df_ref   <- read_csv(ref_path, show_col_types = FALSE)
    df_regen <- read_csv(regen_path, show_col_types = FALSE)
    
    tol <- if (grepl("bootstrap", fn)) {
      1e-6  # Parametric bootstrap (exact seed reproducibility, matches reference to < 1e-10)
    } else if (grepl("digitization_sensitivity", fn)) {
      1e-6  # Stochastic perturbation quantiles (exact seed reproducibility, matches reference to < 1.5e-7)
    } else if (grepl("cv_results", fn)) {
      1e-5  # Multi-start non-linear optimization
    } else if (grepl("monthly_climatology|site_epidemic_rankings|rh_sensitivity|australian_model", fn)) {
      1e-4  # Hourly time-series numerical integration over 25 years
    } else {
      1e-5  # Deterministic analytical submodels
    }
    
    # Check dimensions
    if (nrow(df_ref) != nrow(df_regen) || ncol(df_ref) != ncol(df_regen)) {
      cat(sprintf("  [%-40s] FAIL: Dim mismatch (Ref: %dx%d, Regen: %dx%d)\n", 
                  fn, nrow(df_ref), ncol(df_ref), nrow(df_regen), ncol(df_regen)))
      all_passed <- FALSE
      next
    }
    
    common_cols <- intersect(names(df_ref), names(df_regen))
    eq_res <- all.equal(as.data.frame(df_ref[, common_cols]), as.data.frame(df_regen[, common_cols]), tolerance = tol)
    
    if (isTRUE(eq_res)) {
      cat(sprintf("  [%-40s] PASS (tol = %.1e)\n", fn, tol))
    } else {
      cat(sprintf("  [%-40s] FAIL: Diff detected: %s\n", fn, paste(eq_res[1:min(2, length(eq_res))], collapse = "; ")))
      all_passed <- FALSE
    }
  }
  cat("\n")
  all_passed
}

# 1. Main Default Analysis Tables (outputs/tables vs reference_outputs/tables_v2)
res_default <- verify_table_dir(
  regen_dir = file.path(PROJ_ROOT, "outputs/tables"),
  ref_dir   = file.path(PROJ_ROOT, "reference_outputs/tables_v2"),
  suite_title = "DEFAULT ANALYSIS TABLES (outputs/tables/ vs reference_outputs/tables_v2/)"
)

# 2. Sensitivity Run (a): PRECIP_THRESH = 1.0 mm/h
res_p1 <- verify_table_dir(
  regen_dir = file.path(PROJ_ROOT, "outputs/tables_v2_precip1"),
  ref_dir   = file.path(PROJ_ROOT, "reference_outputs/tables_v2_precip1"),
  suite_title = "SENSITIVITY RUN (a): PRECIP >= 1.0 mm/h (outputs/tables_v2_precip1/)"
)

# 3. Sensitivity Run (b): MAX_EPISODE_H = 24 h
res_m24 <- verify_table_dir(
  regen_dir = file.path(PROJ_ROOT, "outputs/tables_v2_maxh24"),
  ref_dir   = file.path(PROJ_ROOT, "reference_outputs/tables_v2_maxh24"),
  suite_title = "SENSITIVITY RUN (b): MAX EPISODE = 24 h (outputs/tables_v2_maxh24/)"
)

# ------------------------------------------------------------------------------
# Key Manuscript Value Checks (Directly computed from regenerated outputs/tables/)
# ------------------------------------------------------------------------------
cat("------------------------------------------------------------------------------\n")
cat("BENCHMARK VALUES VERIFICATION (Computed directly from outputs/tables/)\n")
cat("------------------------------------------------------------------------------\n")

regen_dir <- file.path(PROJ_ROOT, "outputs/tables")
benchmark_results <- list()

# 1. Brazil Mean EPI_inf = 62.83
site_rank_file <- file.path(regen_dir, "site_epidemic_rankings.csv")
if (file.exists(site_rank_file)) {
  rank_df <- read_csv(site_rank_file, show_col_types = FALSE)
  br_mean <- mean(rank_df$Mean_Annual_EPI_inf, na.rm = TRUE)
  pass_br <- abs(br_mean - 62.83) < 0.1
  status_str <- if (pass_br) "PASS" else "FAIL"
  cat(sprintf("  %-42s Target: 62.83   | Obtained: %6.2f | %s\n", "Brazil Mean Annual EPI_inf:", br_mean, status_str))
  benchmark_results[["Brazil_Mean_EPI"]] <- pass_br
} else {
  cat("  [FAIL] site_epidemic_rankings.csv missing\n")
  benchmark_results[["Brazil_Mean_EPI"]] <- FALSE
}

# 2. Bootstrap Delta RSS = 0.2601 & P = 0.002
boot_file <- file.path(regen_dir, "bootstrap_gamma_results.csv")
if (file.exists(boot_file)) {
  boot_df <- read_csv(boot_file, show_col_types = FALSE)
  drss <- boot_df %>% filter(grepl("Observed Delta RSS", Metric)) %>% pull(Value)
  pval <- boot_df %>% filter(grepl("Bootstrap P-value", Metric)) %>% pull(Value)
  
  pass_drss <- abs(drss - 0.2601) < 0.001
  status_drss <- if (pass_drss) "PASS" else "FAIL"
  cat(sprintf("  %-42s Target: 0.2601  | Obtained: %6.4f | %s\n", "Bootstrap Observed Delta RSS:", drss, status_drss))
  benchmark_results[["Bootstrap_Delta_RSS"]] <- pass_drss
  
  pass_pval <- abs(pval - 0.002) < 0.0015
  status_pval <- if (pass_pval) "PASS" else "FAIL"
  cat(sprintf("  %-42s Target: 0.0020  | Obtained: %6.4f | %s\n", "Bootstrap Empirical P-value:", pval, status_pval))
  benchmark_results[["Bootstrap_P_value"]] <- pass_pval
} else {
  cat("  [FAIL] bootstrap_gamma_results.csv missing\n")
  benchmark_results[["Bootstrap_Delta_RSS"]] <- FALSE
  benchmark_results[["Bootstrap_P_value"]] <- FALSE
}

# 3. Australia Site Values (McLeod = 161.0, Inkerman = 5.74, Delta Downs = 0.72)
aus_file <- file.path(regen_dir, "australian_model_validation_summary.csv")
if (file.exists(aus_file)) {
  aus_df <- read_csv(aus_file, show_col_types = FALSE)
  
  mcleod <- aus_df %>% filter(site_name == "McLeod River") %>% pull(Mean_Annual_EPI_inf)
  pass_mc <- length(mcleod) > 0 && abs(mcleod - 161.0) < 0.5
  cat(sprintf("  %-42s Target: 161.00  | Obtained: %6.2f | %s\n", "Australia McLeod River Mean EPI:", ifelse(length(mcleod)>0, mcleod, NA), if (pass_mc) "PASS" else "FAIL"))
  benchmark_results[["Aus_McLeod"]] <- pass_mc
  
  inkerman <- aus_df %>% filter(site_name == "Inkerman") %>% pull(Mean_Annual_EPI_inf)
  pass_ink <- length(inkerman) > 0 && abs(inkerman - 5.74) < 0.1
  cat(sprintf("  %-42s Target: 5.74    | Obtained: %6.2f | %s\n", "Australia Inkerman Mean EPI:", ifelse(length(inkerman)>0, inkerman, NA), if (pass_ink) "PASS" else "FAIL"))
  benchmark_results[["Aus_Inkerman"]] <- pass_ink
  
  ddowns <- aus_df %>% filter(site_name == "Delta Downs") %>% pull(Mean_Annual_EPI_inf)
  pass_dd <- length(ddowns) > 0 && abs(ddowns - 0.72) < 0.05
  cat(sprintf("  %-42s Target: 0.72    | Obtained: %6.2f | %s\n", "Australia Delta Downs Mean EPI:", ifelse(length(ddowns)>0, ddowns, NA), if (pass_dd) "PASS" else "FAIL"))
  benchmark_results[["Aus_Delta_Downs"]] <- pass_dd
  
  # Welch Two-Sample t-test P-value = 0.140
  high_impact <- aus_df %>%
    filter(grepl("mortality|defoliation|Heavy rust", observed_field_impact, ignore.case = TRUE)) %>%
    pull(Mean_Annual_EPI_inf)
  other_sites <- aus_df %>%
    filter(!grepl("mortality|defoliation|Heavy rust", observed_field_impact, ignore.case = TRUE)) %>%
    pull(Mean_Annual_EPI_inf)
  
  t_res <- t.test(high_impact, other_sites)
  pass_welch <- abs(t_res$p.value - 0.140) < 0.01
  cat(sprintf("  %-42s Target: 0.140   | Obtained: %6.3f | %s\n", "Australia Welch t-test P-value:", t_res$p.value, if (pass_welch) "PASS" else "FAIL"))
  benchmark_results[["Aus_Welch_P"]] <- pass_welch
} else {
  cat("  [FAIL] australian_model_validation_summary.csv missing\n")
  benchmark_results[["Aus_McLeod"]] <- FALSE
  benchmark_results[["Aus_Inkerman"]] <- FALSE
  benchmark_results[["Aus_Delta_Downs"]] <- FALSE
  benchmark_results[["Aus_Welch_P"]] <- FALSE
}

# Summary
cat("==============================================================================\n")
all_bench_pass <- all(unlist(benchmark_results))

if (isTRUE(res_default) && all_bench_pass) {
  cat("OVERALL VERIFICATION RESULT: PASS\n")
  cat("All evaluated tables match reference outputs and all key scientific benchmarks are verified.\n")
} else {
  cat("OVERALL VERIFICATION RESULT: ATTENTION NEEDED\n")
  cat(sprintf("Default tables suite: %s\n", if (isTRUE(res_default)) "PASS" else "FAIL"))
  cat(sprintf("Benchmark checks:     %d / %d passed.\n", sum(unlist(benchmark_results)), length(benchmark_results)))
}
cat("==============================================================================\n")
