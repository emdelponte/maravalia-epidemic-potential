# ==============================================================================
# Script: R/02_infection_submodel.R
# Purpose: Step 2 Infection Submodel Comparison, Fitting, and Evaluation
#          Fits candidate non-linear infection submodels across temperature
#          and dew period regimes.
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(minpack.lm)
  library(ggplot2)
  library(patchwork)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS),
            file.path(PROJ_ROOT, DIR_OUTPUT_MODELS))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 2: INFECTION SUBMODEL FITTING & EVALUATION\n")
cat("==============================================================================\n\n")

# Load digitized infection dataset
data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_infection_data.csv")
if (!file.exists(data_path)) data_path <- file.path(PROJ_ROOT, "data/digitized_infection_data.csv")

inf_data <- read_csv(data_path, show_col_types = FALSE)
N <- nrow(inf_data)
cat("Loaded digitized infection dataset with N =", N, "observations.\n\n")

# Multi-start optimization helper
fit_multistart <- function(formula, data, start_grid, lower, upper) {
  best_fit <- NULL
  best_rss <- Inf
  
  for (i in seq_len(nrow(start_grid))) {
    st <- as.list(start_grid[i, ])
    fit <- try(nlsLM(
      formula,
      data = data,
      start = st,
      lower = lower,
      upper = upper,
      control = nls.lm.control(maxiter = 500)
    ), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      current_rss <- deviance(fit)
      if (current_rss < best_rss) {
        best_rss <- current_rss
        best_fit <- fit
      }
    }
  }
  best_fit
}

# ------------------------------------------------------------------------------
# 2.1 Model 1: Separable Generalized Beta x Weibull
# ------------------------------------------------------------------------------
cat("[1/6] Fitting Model 1: Separable Generalized Beta x Weibull...\n")
m1_formula <- I_obs ~ a * ((T - 10)^d) * ((35 - T)^e) * (1 - exp(-(dew_h / lambda)^k))
m1_lower <- c(a = 1e-9, d = 0.01, e = 0.01, lambda = 0.1, k = 0.1)
m1_upper <- c(a = 1, d = 20, e = 20, lambda = 50, k = 20)
m1_grid <- expand.grid(
  a = c(1e-5, 1e-4, 1e-3),
  d = c(1.5, 2.7, 4.0),
  e = c(1.5, 3.3, 4.0),
  lambda = c(8, 12, 16),
  k = c(2, 3.7, 5)
)
fit_m1 <- fit_multistart(m1_formula, inf_data, m1_grid, m1_lower, m1_upper)

# ------------------------------------------------------------------------------
# 2.2 Model 2: Wetness-expanded Analytis x Weibull
# ------------------------------------------------------------------------------
cat("[2/6] Fitting Model 2: Wetness-expanded Analytis x Weibull...\n")
m2_formula <- I_obs ~ (a * ((T - 10)/25)^b * (1 - ((T - 10)/25)))^(c0 + c1 * dew_h) * (1 - exp(-(dew_h / lambda)^k))
m2_lower <- c(a = 0.01, b = 0.01, c0 = 0.01, c1 = -10.0, lambda = 0.1, k = 0.1)
m2_upper <- c(a = 100, b = 20, c0 = 200, c1 = 5.0, lambda = 50, k = 20)
m2_grid <- expand.grid(
  a = c(1.5, 3.5, 5.0),
  b = c(0.5, 0.8, 1.2),
  c0 = c(1.0, 5.0, 10.0, 25, 50),
  c1 = c(-2, -1, -0.4, -0.2, -0.05),
  lambda = c(8, 11, 15),
  k = c(2, 4, 6)
)
fit_m2 <- fit_multistart(m2_formula, inf_data, m2_grid, m2_lower, m2_upper)

# ------------------------------------------------------------------------------
# 2.3 Model 3: Centered Gaussian Dynamic Breadth x Weibull (Primary Model)
# ------------------------------------------------------------------------------
cat("[3/6] Fitting Model 3: Centered Gaussian Dynamic Breadth x Weibull (Primary)...\n")
m3_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / pmax(1.0, sigma_Wbar + gamma * (dew_h - W_BAR)))^2) * (1 - exp(-(dew_h / lambda)^k))
m3_lower <- c(Topt = 15, sigma_Wbar = 0.1, gamma = 0, lambda = 1.0, k = 1.0)
m3_upper <- c(Topt = 35, sigma_Wbar = 10, gamma = 5, lambda = 50, k = 20)
m3_grid <- expand.grid(
  Topt = c(20.0, 21.2, 22.5),
  sigma_Wbar = c(1.5, 2.4, 3.2),
  gamma = c(0.15, 0.29, 0.45),
  lambda = c(5.0, 9.5, 15.0),
  k = c(1.5, 3.5, 5.0)
)
fit_m3 <- fit_multistart(m3_formula, inf_data, m3_grid, m3_lower, m3_upper)

# ------------------------------------------------------------------------------
# 2.4 Supplementary Models S1, S2, S3
# ------------------------------------------------------------------------------
cat("[4/6] Fitting Supplementary Candidate Models (S1, S2, S3)...\n")
s1_formula <- I_obs ~ (a * ((T - 10)/25)^b * (1 - ((T - 10)/25)))^c * (1 - exp(-(dew_h / lambda)^k))
s1_lower <- c(a = 0.01, b = 0.01, c = 0.01, lambda = 0.1, k = 0.1)
s1_upper <- c(a = 100, b = 20, c = 20, lambda = 50, k = 20)
s1_grid <- expand.grid(a = c(2, 3, 5), b = c(0.5, 1.5), c = c(0.5, 1.0, 2.0), lambda = c(8, 10), k = c(2, 4))
fit_s1 <- fit_multistart(s1_formula, inf_data, s1_grid, s1_lower, s1_upper)

s2_formula <- I_obs ~ (1 - exp(-(B * dew_h)^D)) / cosh((T - F) * G / 2)
s2_lower <- c(B = 0.001, D = 0.1, F = 15, G = 0.01)
s2_upper <- c(B = 1, D = 20, F = 35, G = 5)
s2_grid <- expand.grid(B = c(0.05, 0.085, 0.12), D = c(2.0, 3.88, 5.0), F = c(18, 21.24, 24), G = c(0.2, 0.47, 0.8))
fit_s2 <- fit_multistart(s2_formula, inf_data, s2_grid, s2_lower, s2_upper)

s3_formula <- I_obs ~ exp(-0.5 * ((T - Topt) / sigma_Wbar)^2) * (1 - exp(-(dew_h / lambda)^k))
s3_lower <- c(Topt = 15, sigma_Wbar = 0.1, lambda = 0.1, k = 0.1)
s3_upper <- c(Topt = 35, sigma_Wbar = 10, lambda = 50, k = 20)
s3_grid <- expand.grid(Topt = c(19, 21.2, 23), sigma_Wbar = c(2.0, 4.8, 6.0), lambda = c(8, 11.8, 15), k = c(2, 3.8, 5))
fit_s3 <- fit_multistart(s3_formula, inf_data, s3_grid, s3_lower, s3_upper)

# Save models
saveRDS(fit_m1, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_genbeta.rds"))
saveRDS(fit_m2, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_analytis_exp.rds"))
saveRDS(fit_m3, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_primary.rds"))
saveRDS(fit_s1, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_analytis.rds"))
saveRDS(fit_s2, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_duthie.rds"))
saveRDS(fit_s3, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_fixed_breadth.rds"))

# ------------------------------------------------------------------------------
# 2.5 Information Criteria & Performance Metrics
# ------------------------------------------------------------------------------
cat("[5/6] Computing Performance Metrics and Parameter Tables...\n")
calc_metrics <- function(fit) {
  res <- residuals(fit)
  p_val <- length(coef(fit))
  rss_val <- deviance(fit)
  rmse_val <- sqrt(rss_val / N)
  r2_val <- 1 - (rss_val / sum((inf_data$I_obs - mean(inf_data$I_obs))^2))
  loglik <- as.numeric(logLik(fit))
  aic_val <- -2 * loglik + 2 * (p_val + 1)
  aicc_val <- aic_val + (2 * (p_val + 1) * (p_val + 2)) / (N - p_val - 2)
  list(RMSE = rmse_val, R2 = r2_val, LogLik = loglik, AIC = aic_val, AICc = aicc_val, p = p_val)
}

ic_m1 <- calc_metrics(fit_m1)
ic_m2 <- calc_metrics(fit_m2)
ic_m3 <- calc_metrics(fit_m3)
ic_s1 <- calc_metrics(fit_s1)
ic_s2 <- calc_metrics(fit_s2)
ic_s3 <- calc_metrics(fit_s3)

main_metrics_table <- tibble(
  Model = c(
    "Model 1 — Separable Generalized Beta × Weibull",
    "Model 2 — Wetness-expanded Analytis × Weibull",
    "Model 3 — Gaussian wetness-dependent breadth × Weibull"
  ),
  Parameters = c(ic_m1$p, ic_m2$p, ic_m3$p),
  RMSE = c(ic_m1$RMSE, ic_m2$RMSE, ic_m3$RMSE),
  R2 = c(ic_m1$R2, ic_m2$R2, ic_m3$R2),
  AICc = c(ic_m1$AICc, ic_m2$AICc, ic_m3$AICc)
) %>%
  mutate(deltaAICc = AICc - min(AICc))

write_csv(main_metrics_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "infection_model_metrics.csv"))

all_metrics_table <- tibble(
  Model = c(
    "Model 3 — Gaussian wetness-dependent breadth × Weibull (Primary)",
    "Model 2 — Wetness-expanded Analytis × Weibull (General interaction)",
    "Model 1 — Separable Generalized Beta × Weibull (Reference)",
    "Model S1 — Simple Analytis Beta × Weibull (Supplementary)",
    "Model S2 — Duthie Hyperbolic Secant × Weibull (Supplementary)",
    "Model S3 — Fixed-breadth Gaussian × Weibull (Null H0)"
  ),
  Parameters = c(ic_m3$p, ic_m2$p, ic_m1$p, ic_s1$p, ic_s2$p, ic_s3$p),
  RMSE = c(ic_m3$RMSE, ic_m2$RMSE, ic_m1$RMSE, ic_s1$RMSE, ic_s2$RMSE, ic_s3$RMSE),
  R2 = c(ic_m3$R2, ic_m2$R2, ic_m1$R2, ic_s1$R2, ic_s2$R2, ic_s3$R2),
  LogLik = c(ic_m3$LogLik, ic_m2$LogLik, ic_m1$LogLik, ic_s1$LogLik, ic_s2$LogLik, ic_s3$LogLik),
  AICc = c(ic_m3$AICc, ic_m2$AICc, ic_m1$AICc, ic_s1$AICc, ic_s2$AICc, ic_s3$AICc)
) %>%
  mutate(deltaAICc = AICc - min(AICc))

write_csv(all_metrics_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "all_infection_models_metrics.csv"))

extract_par_df <- function(fit, mod_name) {
  sm <- summary(fit)$coefficients
  tibble(
    Model = mod_name,
    Parameter = rownames(sm),
    Estimate = sm[, 1],
    Std_Error = sm[, 2],
    t_value = sm[, 3],
    p_value = sm[, 4]
  )
}

param_table <- bind_rows(
  extract_par_df(fit_m1, "Model 1: Generalized Beta x Weibull"),
  extract_par_df(fit_m2, "Model 2: Wetness-expanded Analytis x Weibull"),
  extract_par_df(fit_m3, "Model 3: Centered Gaussian Dynamic Breadth x Weibull")
)
write_csv(param_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "infection_model_parameters.csv"))

# ------------------------------------------------------------------------------
# 2.6 Export Main Figure 1 & Supplementary Figures S2 & S3
# ------------------------------------------------------------------------------
cat("[6/6] Generating Main Figure 1 & Supplementary Figures S2 & S3...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )

# --- Figure 1: Main Infection Submodel (2 x 2 Panel Grid) ---
fig1a <- ggplot(inf_data, aes(x = T, y = pustules_cm2, color = factor(dew_h), group = factor(dew_h))) +
  geom_line(linewidth = 0.8, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_color_viridis_d(name = "Dew (h)") +
  labs(
    x = "Temperature (°C)",
    y = expression("Pustule density (pustules cm"^-2*")")
  ) +
  theme_clean +
  theme(legend.position = "right")

grid_temp <- expand.grid(
  T = seq(12, 34, length.out = 100),
  dew_h = c(6, 8, 12, 24)
)
grid_temp$I_fit <- predict(fit_m3, newdata = grid_temp)

fig1b <- ggplot() +
  geom_line(data = grid_temp, aes(x = T, y = I_fit, color = factor(dew_h)), linewidth = 0.9) +
  geom_point(data = inf_data, aes(x = T, y = I_obs, color = factor(dew_h)), size = 2.2) +
  scale_color_viridis_d(name = "Dew (h)") +
  labs(
    x = "Temperature (°C)",
    y = "Relative infection I(T,W)"
  ) +
  theme_clean +
  theme(legend.position = "right")

grid_wet <- expand.grid(
  dew_h = seq(4, 26, length.out = 100),
  T = c(17, 20, 22, 25, 30)
)
grid_wet$I_fit <- predict(fit_m3, newdata = grid_wet)

fig1c <- ggplot() +
  geom_line(data = grid_wet, aes(x = dew_h, y = I_fit, color = factor(T)), linewidth = 0.9) +
  geom_point(data = inf_data, aes(x = dew_h, y = I_obs, color = factor(T)), size = 2.2) +
  scale_color_brewer(palette = "Set1", name = "Temp (°C)") +
  labs(
    x = "Wetness duration (h)",
    y = "Relative infection I(T,W)"
  ) +
  theme_clean +
  theme(legend.position = "right")

grid_surf <- expand.grid(
  T = seq(12, 34, length.out = 100),
  dew_h = seq(4, 26, length.out = 100)
)
grid_surf$I_fit <- predict(fit_m3, newdata = grid_surf)

fig1d <- ggplot() +
  geom_contour_filled(data = grid_surf, aes(x = T, y = dew_h, z = I_fit), alpha = 0.90) +
  geom_rect(aes(xmin = 17, xmax = 30, ymin = 6, ymax = 24), fill = NA, color = "white", linetype = "dashed", linewidth = 0.8) +
  geom_point(data = inf_data, aes(x = T, y = dew_h), color = "red", size = 2.2, shape = 4, stroke = 1.3) +
  scale_fill_viridis_d(name = "Fitted I(T,W)", option = "viridis") +
  labs(
    x = "Temperature (°C)",
    y = "Wetness duration (h)"
  ) +
  theme_clean +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 5)
  )

fig1_main <- (fig1a + fig1b) / (fig1c + fig1d) + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig1_infection_submodel.png"), fig1_main, width = 10, height = 8.5, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig1_infection_submodel.pdf"), fig1_main, width = 10, height = 8.5)
cat("  -> Saved Fig1_infection_submodel.png and .pdf\n")

# --- Supplementary Figure S2: Multi-Model Comparison ---
grid_comp <- expand.grid(
  T = seq(12, 34, length.out = 100),
  dew_h = c(6, 8, 12, 24)
)
grid_comp <- grid_comp %>%
  mutate(
    Model1 = predict(fit_m1, newdata = grid_comp),
    Model2 = predict(fit_m2, newdata = grid_comp),
    Model3 = predict(fit_m3, newdata = grid_comp),
    Duthie = predict(fit_s2, newdata = grid_comp)
  ) %>%
  pivot_longer(cols = c(Model1, Model2, Model3, Duthie), names_to = "Model", values_to = "I_fit") %>%
  mutate(Model = factor(Model, levels = c("Model1", "Model2", "Model3", "Duthie"),
                        labels = c("Model 1 (Beta)", "Model 2 (Expanded)", "Model 3 (Dynamic)", "Duthie")))

fig_s2 <- ggplot() +
  geom_line(data = grid_comp, aes(x = T, y = I_fit, color = Model, linetype = Model), linewidth = 0.8) +
  geom_point(data = inf_data, aes(x = T, y = I_obs), size = 1.8, color = "black") +
  facet_wrap(~ paste0("Wetness = ", dew_h, " h"), ncol = 2) +
  scale_color_brewer(palette = "Set1", name = "Submodel") +
  scale_linetype_discrete(name = "Submodel") +
  labs(
    x = "Temperature (°C)",
    y = "Relative infection I(T,W)"
  ) +
  theme_clean +
  theme(legend.position = "bottom")

fig_s2 <- fig_s2 + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS2_model_comparison.png"), fig_s2, width = 8, height = 6.5, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS2_model_comparison.pdf"), fig_s2, width = 8, height = 6.5)
cat("  -> Saved FigS2_model_comparison.png and .pdf\n")

# --- Supplementary Figure S3: Residual Diagnostics ---
plot_res <- inf_data %>%
  mutate(
    Fitted_M1 = predict(fit_m1),
    Res_M1    = I_obs - Fitted_M1,
    Fitted_M2 = predict(fit_m2),
    Res_M2    = I_obs - Fitted_M2,
    Fitted_M3 = predict(fit_m3),
    Res_M3    = I_obs - Fitted_M3
  )

p_res1 <- ggplot(plot_res, aes(x = Fitted_M1, y = Res_M1, color = factor(dew_h))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 2.2) +
  scale_color_viridis_d(name = "Dew (h)") +
  labs(x = "Fitted response (Model 1)", y = "Residual") +
  theme_clean

p_res2 <- ggplot(plot_res, aes(x = Fitted_M2, y = Res_M2, color = factor(dew_h))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 2.2) +
  scale_color_viridis_d(name = "Dew (h)") +
  labs(x = "Fitted response (Model 2)", y = "Residual") +
  theme_clean

p_res3 <- ggplot(plot_res, aes(x = Fitted_M3, y = Res_M3, color = factor(dew_h))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 2.2) +
  scale_color_viridis_d(name = "Dew (h)") +
  labs(x = "Fitted response (Model 3)", y = "Residual") +
  theme_clean

fig_s3 <- (p_res1 + p_res2 + p_res3) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

fig_s3 <- fig_s3 + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS3_residuals.png"), fig_s3, width = 11, height = 3.8, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS3_residuals.pdf"), fig_s3, width = 11, height = 3.8)
cat("  -> Saved FigS3_residuals.png and .pdf\n\n")

cat("Step 2 completed successfully.\n")
cat("==============================================================================\n")
