# ==============================================================================
# Script: R/03_validate_infection_models.R
# Purpose: Step 3 Grouped Cross-Validation (LOWO and LOTO) for Candidate Infection Models
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

table_dir <- file.path(PROJ_ROOT, DIR_OUTPUT_TABLES)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("STEP 3: GROUPED CROSS-VALIDATION (LOWO & LOTO)\n")
cat("==============================================================================\n\n")

data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_infection_data.csv")
if (!file.exists(data_path)) data_path <- file.path(PROJ_ROOT, "data/digitized_infection_data.csv")

inf_data <- read_csv(data_path, show_col_types = FALSE)
N <- nrow(inf_data)

# Model Formulas & Starting Values
m1_formula <- I_obs ~ a * ((T - 10)^d) * ((35 - T)^e) * (1 - exp(-(dew_h / lambda)^k))
m1_start   <- list(a = 0.0001, d = 2.7, e = 3.3, lambda = 12.0, k = 3.7)
m1_lower   <- c(a = 1e-9, d = 0.01, e = 0.01, lambda = 0.1, k = 0.1)
m1_upper   <- c(a = 1, d = 20, e = 20, lambda = 50, k = 20)

m2_formula <- I_obs ~ (a * ((T - 10)/25)^b * (1 - ((T - 10)/25)))^(c0 + c1 * dew_h) * (1 - exp(-(dew_h / lambda)^k))
m2_start   <- list(a = 3.5, b = 0.8, c0 = 10.0, c1 = -0.3, lambda = 11.0, k = 4.0)
m2_lower   <- c(a = 0.01, b = 0.01, c0 = 0.01, c1 = -10, lambda = 0.1, k = 0.1)
m2_upper   <- c(a = 100, b = 20, c0 = 200, c1 = 5.0, lambda = 50, k = 20)

m3_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / pmax(1.0, sigma_Wbar + gamma * (dew_h - W_BAR)))^2) * (1 - exp(-(dew_h / lambda)^k))
m3_start   <- list(Topt = 21.2, sigma_Wbar = 2.4, gamma = 0.29, lambda = 9.6, k = 4.3)
m3_lower   <- c(Topt = 15, sigma_Wbar = 0.1, gamma = 0, lambda = 1.0, k = 1.0)
m3_upper   <- c(Topt = 35, sigma_Wbar = 10, gamma = 5, lambda = 50, k = 20)

s1_formula <- I_obs ~ (a * ((T - 10)/25)^b * (1 - ((T - 10)/25)))^c * (1 - exp(-(dew_h / lambda)^k))
s1_start   <- list(a = 3.0, b = 1.5, c = 1.0, lambda = 10.0, k = 4.0)
s1_lower   <- c(a = 0.01, b = 0.01, c = 0.01, lambda = 0.1, k = 0.1)
s1_upper   <- c(a = 100, b = 20, c = 20, lambda = 50, k = 20)

s2_formula <- I_obs ~ (1 - exp(-(B * dew_h)^D)) / cosh((T - F) * G / 2)
s2_start   <- list(B = 0.085, D = 3.88, F = 21.24, G = 0.470)
s2_lower   <- c(B = 0.001, D = 0.1, F = 15, G = 0.01)
s2_upper   <- c(B = 1, D = 20, F = 35, G = 5)

s3_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / sigma_Wbar)^2) * (1 - exp(-(dew_h / lambda)^k))
s3_start   <- list(Topt = 21.2, sigma_Wbar = 4.8, lambda = 11.8, k = 3.8)
s3_lower   <- c(Topt = 15, sigma_Wbar = 0.1, lambda = 0.1, k = 0.1)
s3_upper   <- c(Topt = 35, sigma_Wbar = 10, lambda = 50, k = 20)

fit_safe_single <- function(formula, data, start, lower, upper) {
  try(nlsLM(formula, data = data, start = start, lower = lower, upper = upper,
            control = nls.lm.control(maxiter = 500)), silent = TRUE)
}

fit_safe <- function(formula, data, start, lower, upper) {
  full <- fit_safe_single(formula, inf_data, start, lower, upper)
  starts <- list(start)
  if (!inherits(full, "try-error")) starts[[2]] <- as.list(coef(full))
  set.seed(SEED_OPTIM_CV)
  for (j in 1:10) {
    starts[[length(starts) + 1]] <- mapply(
      function(v, lo, hi) min(hi, max(lo, v * exp(runif(1, -0.5, 0.5)))),
      start, lower[names(start)], upper[names(start)], SIMPLIFY = FALSE
    )
  }
  best <- NULL
  for (st in starts) {
    f <- fit_safe_single(formula, data, st, lower, upper)
    if (!inherits(f, "try-error") && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  if (is.null(best)) structure("fail", class = "try-error") else best
}

# --- 1. Leave-One-Wetness-Level-Out (LOWO) CV ---
cat("[1/3] Running Leave-One-Wetness-Level-Out (LOWO) Cross-Validation...\n")
wetness_levels <- unique(inf_data$dew_h) # c(6, 8, 12, 24)

calc_lowo_rmse <- function(formula, start, lower, upper) {
  sq_errs <- c()
  for (w in wetness_levels) {
    train_d <- inf_data %>% filter(dew_h != w)
    test_d  <- inf_data %>% filter(dew_h == w)
    fit <- fit_safe(formula, train_d, start, lower, upper)
    if (!inherits(fit, "try-error")) {
      preds <- predict(fit, newdata = test_d)
      sq_errs <- c(sq_errs, (test_d$I_obs - preds)^2)
    }
  }
  sqrt(mean(sq_errs))
}

lowo_m1 <- calc_lowo_rmse(m1_formula, m1_start, m1_lower, m1_upper)
lowo_m2 <- calc_lowo_rmse(m2_formula, m2_start, m2_lower, m2_upper)
lowo_m3 <- calc_lowo_rmse(m3_formula, m3_start, m3_lower, m3_upper)
lowo_s1 <- calc_lowo_rmse(s1_formula, s1_start, s1_lower, s1_upper)
lowo_s2 <- calc_lowo_rmse(s2_formula, s2_start, s2_lower, s2_upper)
lowo_s3 <- calc_lowo_rmse(s3_formula, s3_start, s3_lower, s3_upper)

cat("  - LOWO RMSE Model 1 (Beta x Weibull):", round(lowo_m1, 4), "\n")
cat("  - LOWO RMSE Model 2 (Expanded Analytis):", round(lowo_m2, 4), "\n")
cat("  - LOWO RMSE Model 3 (Centered Dynamic Breadth):", round(lowo_m3, 4), "\n\n")

# --- 2. Leave-One-Temperature-Level-Out (LOTO) CV ---
cat("[2/3] Running Leave-One-Temperature-Level-Out (LOTO) Cross-Validation...\n")
temp_levels <- unique(inf_data$T) # c(17, 20, 22, 25, 30)

calc_loto_rmse <- function(formula, start, lower, upper) {
  sq_errs <- c()
  for (temp in temp_levels) {
    train_d <- inf_data %>% filter(T != temp)
    test_d  <- inf_data %>% filter(T == temp)
    fit <- fit_safe(formula, train_d, start, lower, upper)
    if (!inherits(fit, "try-error")) {
      preds <- predict(fit, newdata = test_d)
      sq_errs <- c(sq_errs, (test_d$I_obs - preds)^2)
    }
  }
  sqrt(mean(sq_errs))
}

loto_m1 <- calc_loto_rmse(m1_formula, m1_start, m1_lower, m1_upper)
loto_m2 <- calc_loto_rmse(m2_formula, m2_start, m2_lower, m2_upper)
loto_m3 <- calc_loto_rmse(m3_formula, m3_start, m3_lower, m3_upper)
loto_s1 <- calc_loto_rmse(s1_formula, s1_start, s1_lower, s1_upper)
loto_s2 <- calc_loto_rmse(s2_formula, s2_start, s2_lower, s2_upper)
loto_s3 <- calc_loto_rmse(s3_formula, s3_start, s3_lower, s3_upper)

cat("  - LOTO RMSE Model 1 (Beta x Weibull):", round(loto_m1, 4), "\n")
cat("  - LOTO RMSE Model 2 (Expanded Analytis):", round(loto_m2, 4), "\n")
cat("  - LOTO RMSE Model 3 (Centered Dynamic Breadth):", round(loto_m3, 4), "\n\n")

# --- 3. Save Combined Metrics Table ---
cat("[3/3] Exporting Grouped Cross-Validation Summary Table...\n")

cv_table <- tibble(
  Model = c(
    "Model 1 — Separable Generalized Beta × Weibull",
    "Model 2 — Wetness-expanded Analytis × Weibull",
    "Model 3 — Gaussian wetness-dependent breadth × Weibull"
  ),
  LOWO_RMSE = c(lowo_m1, lowo_m2, lowo_m3),
  LOTO_RMSE = c(loto_m1, loto_m2, loto_m3)
)

write_csv(cv_table, file.path(table_dir, "infection_model_cv_results.csv"))

main_metrics_file <- file.path(table_dir, "infection_model_metrics.csv")
if (file.exists(main_metrics_file)) {
  main_metrics <- read_csv(main_metrics_file, show_col_types = FALSE)
  merged_metrics <- left_join(main_metrics, cv_table, by = "Model") %>%
    select(Model, Parameters, RMSE, R2, AICc, deltaAICc, LOWO_RMSE, LOTO_RMSE)
  write_csv(merged_metrics, main_metrics_file)
  cat("  -> Updated outputs/tables/infection_model_metrics.csv with LOWO & LOTO RMSE\n\n")
  print(merged_metrics)
}

cat("\nStep 3 completed successfully.\n")
cat("==============================================================================\n")
