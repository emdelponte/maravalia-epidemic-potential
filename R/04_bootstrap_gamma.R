# ==============================================================================
# Script: R/04_bootstrap_gamma.R
# Purpose: Step 4 Parametric Bootstrap Likelihood Ratio / RSS Test for Dynamic
#          Breadth Parameter gamma in Model 3 vs Fixed Breadth (Null H0)
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
  library(ggplot2)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS),
            file.path(PROJ_ROOT, DIR_DATA_DERIVED))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 4: PARAMETRIC BOOTSTRAP TEST FOR GAMMA (DYNAMIC BREADTH)\n")
cat("==============================================================================\n\n")

data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_infection_data.csv")
if (!file.exists(data_path)) data_path <- file.path(PROJ_ROOT, "data/digitized_infection_data.csv")

inf_data <- read_csv(data_path, show_col_types = FALSE)
N <- nrow(inf_data)

# 1. Fit Null H0 (Fixed breadth, gamma = 0)
m_h0_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / sigma_Wbar)^2) * (1 - exp(-(dew_h / lambda)^k))
st_h0 <- list(Topt = 21.2, sigma_Wbar = 4.8, lambda = 11.8, k = 3.8)

fit_h0 <- try(nlsLM(
  m_h0_formula,
  data = inf_data,
  start = st_h0,
  lower = c(Topt = 15, sigma_Wbar = 0.1, lambda = 1.0, k = 1.0),
  upper = c(Topt = 35, sigma_Wbar = 10, lambda = 50, k = 20),
  control = nls.lm.control(maxiter = 500)
), silent = TRUE)

# 2. Fit Alternative H1 (Dynamic breadth, gamma >= 0)
m_h1_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / pmax(1.0, sigma_Wbar + gamma * (dew_h - W_BAR)))^2) * (1 - exp(-(dew_h / lambda)^k))
st_h1 <- list(Topt = 21.2, sigma_Wbar = 2.4, gamma = 0.29, lambda = 9.6, k = 4.3)

fit_h1_best <- function(df, h0) {
  starts <- list(st_h1)
  if (!inherits(h0, "try-error")) {
    c0 <- coef(h0)
    for (g in c(0, 0.1, 0.3)) {
      starts[[length(starts) + 1]] <- list(
        Topt = c0[["Topt"]],
        sigma_Wbar = c0[["sigma_Wbar"]],
        gamma = g,
        lambda = max(1, c0[["lambda"]]),
        k = max(1, c0[["k"]])
      )
    }
  }
  for (j in 1:4) {
    starts[[length(starts) + 1]] <- lapply(st_h1, function(v) v * exp(runif(1, -0.4, 0.4)))
  }
  best <- NULL
  for (st in starts) {
    f <- try(nlsLM(
      m_h1_formula,
      data = df,
      start = st,
      lower = c(Topt = 15, sigma_Wbar = 0.1, gamma = 0, lambda = 1.0, k = 1.0),
      upper = c(Topt = 35, sigma_Wbar = 10, gamma = 5, lambda = 50, k = 20),
      control = nls.lm.control(maxiter = 300)
    ), silent = TRUE)
    if (!inherits(f, "try-error") && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  if (is.null(best)) structure("fail", class = "try-error") else best
}

set.seed(SEED_H1_FIT)
fit_h1 <- fit_h1_best(inf_data, fit_h0)

rss_h0_obs <- deviance(fit_h0)
rss_h1_obs <- deviance(fit_h1)
delta_rss_obs <- rss_h0_obs - rss_h1_obs

cat(sprintf("Observed RSS under H0 (Fixed breadth): %.4f\n", rss_h0_obs))
cat(sprintf("Observed RSS under H1 (Dynamic breadth): %.4f\n", rss_h1_obs))
cat(sprintf("Observed Delta RSS (RSS_H0 - RSS_H1): %.4f\n\n", delta_rss_obs))

# Simulated data generation under H0
pred_h0 <- predict(fit_h0)
res_h0  <- residuals(fit_h0)
sigma_res_h0 <- sqrt(sum(res_h0^2) / (N - 4))

cache_file <- file.path(PROJ_ROOT, DIR_DATA_DERIVED, "bootstrap_gamma_samples.rds")

# Check for command line flag --quick or --recompute
args <- commandArgs(trailingOnly = TRUE)
if (!exists("quick_mode")) quick_mode <- "--quick" %in% args
if (!exists("recompute")) recompute <- "--recompute" %in% args

B <- if (quick_mode) 100 else BOOTSTRAP_B

if (file.exists(cache_file) && !recompute && !quick_mode) {
  cat("Loading cached bootstrap replicates from:", cache_file, "...\n")
  cached_boot <- readRDS(cache_file)
  valid_delta <- cached_boot$valid_delta
  valid_gamma <- cached_boot$valid_gamma
  p_value     <- cached_boot$p_value
} else {
  set.seed(SEED_BOOTSTRAP)
  boot_delta_rss <- numeric(B)
  boot_gamma     <- numeric(B)
  fails <- 0

  cat(sprintf("Running Parametric Bootstrap under H0 (B = %d replicates)...\n", B))
  t0 <- Sys.time()

  for (b in 1:B) {
    y_sim <- pred_h0 + rnorm(N, mean = 0, sd = sigma_res_h0)
    y_sim <- pmin(1.0, pmax(0.0, y_sim))
    df_sim <- inf_data %>% mutate(I_obs = y_sim)

    fit_h0_b <- try(nlsLM(
      m_h0_formula,
      data = df_sim,
      start = st_h0,
      lower = c(Topt = 15, sigma_Wbar = 0.1, lambda = 1.0, k = 1.0),
      upper = c(Topt = 35, sigma_Wbar = 10, lambda = 50, k = 20),
      control = nls.lm.control(maxiter = 300)
    ), silent = TRUE)

    fit_h1_b <- fit_h1_best(df_sim, fit_h0_b)

    if (!inherits(fit_h0_b, "try-error") && !inherits(fit_h1_b, "try-error")) {
      boot_delta_rss[b] <- deviance(fit_h0_b) - deviance(fit_h1_b)
      boot_gamma[b]     <- coef(fit_h1_b)["gamma"]
    } else {
      fails <- fails + 1
      boot_delta_rss[b] <- NA
      boot_gamma[b]     <- NA
    }

    if (b %% 500 == 0 || (quick_mode && b %% 25 == 0)) {
      cat(sprintf("   Replicate %d / %d completed...\n", b, B))
    }
  }
  t1 <- Sys.time()
  cat(sprintf("  -> Completed in %.2f seconds (Failed fits: %d / %d).\n", difftime(t1, t0, units="secs"), fails, B))

  valid_delta <- boot_delta_rss[!is.na(boot_delta_rss)]
  p_value <- (sum(valid_delta >= delta_rss_obs) + 1) / (length(valid_delta) + 1)
  valid_gamma <- boot_gamma[!is.na(boot_gamma)]

  if (!quick_mode) {
    saveRDS(list(valid_delta = valid_delta, valid_gamma = valid_gamma, p_value = p_value), cache_file)
  }
}

gamma_quantiles <- quantile(valid_gamma, probs = c(0.025, 0.50, 0.975))
q95 <- quantile(valid_delta, 0.95)

boot_summary <- tibble(
  Metric = c("Observed Delta RSS", "Mean H0 Delta RSS", "95th Percentile H0 Delta RSS", "Bootstrap P-value", "Gamma H0 Median", "Gamma H0 2.5%", "Gamma H0 97.5%"),
  Value = c(delta_rss_obs, mean(valid_delta), q95, p_value, gamma_quantiles["50%"], gamma_quantiles["2.5%"], gamma_quantiles["97.5%"])
)

target_tbl_file <- if (quick_mode) file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "bootstrap_gamma_results_quick.csv") else file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "bootstrap_gamma_results.csv")
write_csv(boot_summary, target_tbl_file)

# Save Supplementary Figure S1
df_boot <- tibble(Delta_RSS = valid_delta)

fig_s1 <- ggplot(df_boot, aes(x = Delta_RSS)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = q95, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = delta_rss_obs, color = "red", linetype = "solid", linewidth = 1.1) +
  annotate(
    "text", x = q95 + 0.015, y = 650,
    label = sprintf("95th Percentile\n(%.4f)", q95),
    color = "darkblue", fontface = "bold", hjust = 0, size = 3.6
  ) +
  annotate(
    "text", x = delta_rss_obs - 0.015, y = 350,
    label = sprintf("Observed ΔRSS\n= %.4f\n(P = %.4f)", delta_rss_obs, p_value),
    color = "red", fontface = "bold", hjust = 1, size = 3.6
  ) +
  labs(
    x = expression("Improvement statistic"~Delta*RSS~(RSS[H0] - RSS[H1])),
    y = "Frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS1_bootstrap_gamma.png"), fig_s1, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS1_bootstrap_gamma.pdf"), fig_s1, width = 7, height = 4.5)

cat("  -> Saved FigS1_bootstrap_gamma.png and .pdf\n")
cat(sprintf("Observed Delta RSS: %.4f | 95th Percentile: %.4f | Bootstrap P-value: %.4f\n\n",
            delta_rss_obs, q95, p_value))

cat("Step 4 completed successfully.\n")
cat("==============================================================================\n")
