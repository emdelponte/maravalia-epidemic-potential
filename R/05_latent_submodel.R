# ==============================================================================
# Script: R/05_latent_submodel.R
# Purpose: Step 5 Latent Development Submodel Fitting, Alternative Model
#          Comparison, Diagnostics, and Manuscript Figure 2 Export
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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
cat("STEP 5: LATENT DEVELOPMENT SUBMODEL FITTING & EVALUATION\n")
cat("==============================================================================\n\n")

data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_latent_period.csv")
if (!file.exists(data_path)) data_path <- file.path(PROJ_ROOT, "data/digitized_latent_period.csv")

latent_data <- read_csv(data_path, show_col_types = FALSE)
N_lat <- nrow(latent_data)

cat("Loaded digitized latent period dataset with N =", N_lat, "observations.\n")
cat("Experimental temperature domain:", min(latent_data$T), "to", max(latent_data$T), "°C.\n\n")

# 1. Primary Model: Linear Development Rate r(T) = a + b*T
cat("[1/3] Fitting Primary Linear Development Rate Model r(T) = a + b*T...\n")

fit_lin_rate <- lm(r ~ T, data = latent_data)
sum_lin <- summary(fit_lin_rate)

coef_lin <- coef(fit_lin_rate)
a_hat <- coef_lin["(Intercept)"]
b_hat <- coef_lin["T"]

r2_rate_lin <- sum_lin$r.squared
rmse_rate_lin <- sqrt(mean(residuals(fit_lin_rate)^2))

latent_data <- latent_data %>%
  mutate(
    r_pred_lin = predict(fit_lin_rate),
    LP_pred_lin = 1 / r_pred_lin,
    res_LP_lin = LP_days - LP_pred_lin
  )

rmse_LP_lin <- sqrt(mean(latent_data$res_LP_lin^2))
r2_LP_lin <- 1 - (sum(latent_data$res_LP_lin^2) / sum((latent_data$LP_days - mean(latent_data$LP_days))^2))

saveRDS(fit_lin_rate, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_latent_linear.rds"))

# 2. Alternative Model: Quadratic Development Rate r(T) = a + b*T + c*T^2
cat("[2/3] Fitting Alternative Quadratic Development Rate Model...\n")

fit_quad_rate <- lm(r ~ T + I(T^2), data = latent_data)
sum_quad <- summary(fit_quad_rate)

r2_rate_quad <- sum_quad$r.squared
rmse_rate_quad <- sqrt(mean(residuals(fit_quad_rate)^2))

latent_data <- latent_data %>%
  mutate(
    r_pred_quad = predict(fit_quad_rate),
    LP_pred_quad = 1 / r_pred_quad,
    res_LP_quad = LP_days - LP_pred_quad
  )

rmse_LP_quad <- sqrt(mean(latent_data$res_LP_quad^2))
r2_LP_quad <- 1 - (sum(latent_data$res_LP_quad^2) / sum((latent_data$LP_days - mean(latent_data$LP_days))^2))

saveRDS(fit_quad_rate, file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_latent_nonlinear.rds"))

N_lat <- nrow(latent_data)
p_lin <- 2
p_quad <- 3

aic_lin <- N_lat * log(sum(residuals(fit_lin_rate)^2) / N_lat) + 2 * p_lin
aicc_lin <- aic_lin + (2 * p_lin * (p_lin + 1)) / (N_lat - p_lin - 1)

aic_quad <- N_lat * log(sum(residuals(fit_quad_rate)^2) / N_lat) + 2 * p_quad
aicc_quad <- aic_quad + (2 * p_quad * (p_quad + 1)) / (N_lat - p_quad - 1)

# Save Metrics Table
latent_metrics_table <- tibble(
  Model = c("Primary Linear Rate [r(T) = a + bT]", "Alternative Quadratic Rate [r(T) = a + bT + cT^2]"),
  Num_Params = c(p_lin, p_quad),
  R2_rate = c(r2_rate_lin, r2_rate_quad),
  RMSE_rate = c(rmse_rate_lin, rmse_rate_quad),
  RMSE_latent_days = c(rmse_LP_lin, rmse_LP_quad),
  AIC = c(aic_lin, aic_quad),
  AICc = c(aicc_lin, aicc_quad)
)

write_csv(latent_metrics_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "latent_model_metrics.csv"))

# ------------------------------------------------------------------------------
# 3. Export Manuscript Figure 2 (4-panel composite)
# ------------------------------------------------------------------------------
cat("[3/3] Exporting Manuscript Figure 2 (Fig2_latent_period)...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )

grid_lat <- tibble(T = seq(18, 27, length.out = 100)) %>%
  mutate(
    r_fit_lin = predict(fit_lin_rate, newdata = .),
    LP_fit_lin = 1 / r_fit_lin,
    r_fit_quad = predict(fit_quad_rate, newdata = .),
    LP_fit_quad = 1 / r_fit_quad
  )

# Panel A: Observed LP vs Temp
fig2a <- ggplot(latent_data, aes(x = T, y = LP_days)) +
  geom_point(size = 2.8, color = "darkblue", alpha = 0.8) +
  labs(x = "Temperature (°C)", y = "Latent period (days)") +
  theme_clean

# Panel B: Development Rate r(T)
fig2b <- ggplot(latent_data, aes(x = T, y = r)) +
  geom_point(size = 2.8, color = "firebrick", alpha = 0.8) +
  geom_abline(intercept = a_hat, slope = b_hat, color = "firebrick", linewidth = 0.9) +
  annotate(
    "text", x = 18.2, y = 0.12,
    label = sprintf("r(T) == %.4f + %.5f*T~~(R^2 == %.3f)", a_hat, b_hat, r2_rate_lin),
    parse = TRUE, hjust = 0, size = 3.2, color = "firebrick"
  ) +
  labs(x = "Temperature (°C)", y = expression("Development rate r (day"^-1*")")) +
  theme_clean

# Panel C: Fitted Back-Transformed LP Curves
fig2c <- ggplot() +
  geom_point(data = latent_data, aes(x = T, y = LP_days), size = 2.8, alpha = 0.8) +
  geom_line(data = grid_lat, aes(x = T, y = LP_fit_lin, color = "Linear rate [1/r(T)]", linetype = "Linear rate [1/r(T)]"), linewidth = 1) +
  geom_line(data = grid_lat, aes(x = T, y = LP_fit_quad, color = "Quadratic rate [1/r(T)]", linetype = "Quadratic rate [1/r(T)]"), linewidth = 0.9) +
  scale_color_manual(name = "Rate model [back-transformed]:", values = c("Linear rate [1/r(T)]" = "blue", "Quadratic rate [1/r(T)]" = "red")) +
  scale_linetype_manual(name = "Rate model [back-transformed]:", values = c("Linear rate [1/r(T)]" = "solid", "Quadratic rate [1/r(T)]" = "dashed")) +
  labs(x = "Temperature (°C)", y = "Latent period (days)") +
  theme_clean

# Panel D: Residuals vs Fitted LP
fig2d <- ggplot(latent_data, aes(x = LP_pred_lin, y = res_LP_lin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 2.8, color = "purple", alpha = 0.8) +
  labs(x = "Fitted latent period (days)", y = "Residual (days)") +
  theme_clean

fig2_combined <- (fig2a + fig2b) / (fig2c + fig2d) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

fig2_combined <- fig2_combined + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig2_latent_period.png"), fig2_combined, width = 9, height = 7.5, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig2_latent_period.pdf"), fig2_combined, width = 9, height = 7.5)
cat("  -> Saved Fig2_latent_period.png and .pdf\n\n")

cat("Step 5 completed successfully.\n")
cat("==============================================================================\n")
