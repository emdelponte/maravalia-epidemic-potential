# ==============================================================================
# Script: R/06_digitization_sensitivity.R
# Purpose: Step 6 Digitization Sensitivity Analysis via Stochastic Perturbations
#          and Generation of Sensitivity Bands for Model 3 and Latent Model
#          Exports Manuscript Figure S4 to outputs/figures/
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
  library(ggplot2)
  library(patchwork)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 6: DIGITIZATION SENSITIVITY ANALYSIS\n")
cat("==============================================================================\n\n")

set.seed(SEED_SENSITIVITY)

inf_data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_infection_data.csv")
if (!file.exists(inf_data_path)) inf_data_path <- file.path(PROJ_ROOT, "data/digitized_infection_data.csv")
inf_data <- read_csv(inf_data_path, show_col_types = FALSE)

latent_data_path <- file.path(PROJ_ROOT, "data/metadata/digitized_latent_period.csv")
if (!file.exists(latent_data_path)) latent_data_path <- file.path(PROJ_ROOT, "data/digitized_latent_period.csv")
latent_data <- read_csv(latent_data_path, show_col_types = FALSE)

max_pustules <- max(inf_data$pustules_cm2)

inf_model_path <- file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_primary.rds")
if (!file.exists(inf_model_path)) inf_model_path <- file.path(PROJ_ROOT, "models_v2/model_infection_primary.rds")
fit_primary <- readRDS(inf_model_path)
coef_prim <- coef(fit_primary)

lat_model_path <- file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_latent_linear.rds")
if (!file.exists(lat_model_path)) lat_model_path <- file.path(PROJ_ROOT, "models/model_latent_linear.rds")
fit_lin_rate <- readRDS(lat_model_path)

# 1. Infection Model Digitization Sensitivity
n_sim <- SENSITIVITY_B
cat(sprintf("[1/3] Running Infection Submodel Digitization Sensitivity (B = %d iterations)...\n", n_sim))

grid_temp <- expand.grid(
  T = seq(15, 32, length.out = 50),
  dew_h = c(8, 12, 24)
)

inf_sim_preds <- matrix(NA, nrow = nrow(grid_temp), ncol = n_sim)
inf_sim_params <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(inf_sim_params) <- c("Topt", "sigma_Wbar", "gamma", "lambda", "k")

primary_formula <- I_obs_pert ~ exp(-0.5 * ((T - Topt) / pmax(1.0, sigma_Wbar + gamma * (dew_h - W_BAR)))^2) * (1 - exp(-(dew_h / lambda)^k))

successful_inf_sims <- 0

for (b in 1:n_sim) {
  noise <- runif(nrow(inf_data), min = -0.05, max = 0.05)
  inf_pert <- inf_data %>%
    mutate(
      pustules_pert = pmax(0.01, pustules_cm2 + noise),
      I_obs_pert = pustules_pert / max(pustules_pert)
    )
  
  fit_pert <- try(nlsLM(
    primary_formula,
    data = inf_pert,
    start = as.list(coef_prim),
    lower = c(Topt = 15, sigma_Wbar = 0.1, gamma = 0, lambda = 1.0, k = 1.0),
    upper = c(Topt = 35, sigma_Wbar = 10, gamma = 5, lambda = 50, k = 20),
    control = nls.lm.control(maxiter = 300)
  ), silent = TRUE)
  
  if (!inherits(fit_pert, "try-error")) {
    successful_inf_sims <- successful_inf_sims + 1
    inf_sim_params[b, ] <- coef(fit_pert)
    inf_sim_preds[, b] <- predict(fit_pert, newdata = grid_temp)
  }
}

cat("  - Successful infection fits:", successful_inf_sims, "/", n_sim, "\n")

grid_temp$q05 <- apply(inf_sim_preds, 1, quantile, probs = 0.05, na.rm = TRUE)
grid_temp$q50 <- apply(inf_sim_preds, 1, quantile, probs = 0.50, na.rm = TRUE)
grid_temp$q95 <- apply(inf_sim_preds, 1, quantile, probs = 0.95, na.rm = TRUE)

inf_param_quantiles <- apply(inf_sim_params, 2, quantile, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)

# 2. Latent Period Submodel Digitization Sensitivity
cat(sprintf("[2/3] Running Latent Period Digitization Sensitivity (B = %d iterations)...\n", n_sim))

grid_lat <- tibble(T = seq(18, 27, length.out = 50))
lat_sim_preds <- matrix(NA, nrow = nrow(grid_lat), ncol = n_sim)
lat_sim_params <- matrix(NA, nrow = n_sim, ncol = 2)
colnames(lat_sim_params) <- c("a", "b")

for (b in 1:n_sim) {
  noise_lp <- runif(nrow(latent_data), min = -0.5, max = 0.5)
  lat_pert <- latent_data %>%
    mutate(
      LP_pert = pmax(5, LP_days + noise_lp),
      r_pert = 1 / LP_pert
    )
  
  fit_lat_pert <- lm(r_pert ~ T, data = lat_pert)
  lat_sim_params[b, ] <- coef(fit_lat_pert)
  r_pred <- predict(fit_lat_pert, newdata = grid_lat)
  lat_sim_preds[, b] <- 1 / r_pred
}

grid_lat$q05 <- apply(lat_sim_preds, 1, quantile, probs = 0.05, na.rm = TRUE)
grid_lat$q50 <- apply(lat_sim_preds, 1, quantile, probs = 0.50, na.rm = TRUE)
grid_lat$q95 <- apply(lat_sim_preds, 1, quantile, probs = 0.95, na.rm = TRUE)

lat_param_quantiles <- apply(lat_sim_params, 2, quantile, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)

# 3. Export Supplementary Figure S4 (Publication Standard)
cat("[3/3] Exporting Manuscript Figure S4 (FigS4_digitization_sensitivity)...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )

fig4a <- ggplot() +
  geom_ribbon(data = grid_temp, aes(x = T, ymin = q05, ymax = q95, fill = factor(dew_h)), alpha = 0.25) +
  geom_line(data = grid_temp, aes(x = T, y = q50, color = factor(dew_h)), linewidth = 0.9) +
  geom_point(data = inf_data %>% filter(dew_h %in% c(8, 12, 24)), aes(x = T, y = I_obs, color = factor(dew_h)), size = 2.2) +
  scale_fill_viridis_d(name = "Wetness (h)") +
  scale_color_viridis_d(name = "Wetness (h)") +
  labs(
    x = "Temperature (°C)",
    y = "Relative infection I(T,W)"
  ) +
  theme_clean

fig4b <- ggplot() +
  geom_ribbon(data = grid_lat, aes(x = T, ymin = q05, ymax = q95), fill = "blue", alpha = 0.2) +
  geom_line(data = grid_lat, aes(x = T, y = q50), color = "blue", linewidth = 1) +
  geom_point(data = latent_data, aes(x = T, y = LP_days), size = 2.5, color = "black") +
  labs(
    x = "Temperature (°C)",
    y = "Latent period (days)"
  ) +
  theme_clean

fig4_combined <- (fig4a + fig4b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig4_combined <- fig4_combined + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS4_digitization_sensitivity.png"), fig4_combined, width = 10, height = 4.2, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "FigS4_digitization_sensitivity.pdf"), fig4_combined, width = 10, height = 4.2)

sens_table <- tibble(
  Submodel = c("Infection - Topt (°C)", "Infection - sigma_Wbar (°C)", "Infection - gamma (°C h^-1)", "Infection - lambda (h)", "Infection - k",
               "Latent - a (intercept)", "Latent - b (slope)"),
  Baseline_Estimate = c(coef_prim["Topt"], coef_prim["sigma_Wbar"], coef_prim["gamma"], coef_prim["lambda"], coef_prim["k"],
                        coef(fit_lin_rate)[1], coef(fit_lin_rate)[2]),
  P05_Sensitivity = c(inf_param_quantiles[1, "Topt"], inf_param_quantiles[1, "sigma_Wbar"], inf_param_quantiles[1, "gamma"], inf_param_quantiles[1, "lambda"], inf_param_quantiles[1, "k"],
                      lat_param_quantiles[1, "a"], lat_param_quantiles[1, "b"]),
  P50_Sensitivity = c(inf_param_quantiles[2, "Topt"], inf_param_quantiles[2, "sigma_Wbar"], inf_param_quantiles[2, "gamma"], inf_param_quantiles[2, "lambda"], inf_param_quantiles[2, "k"],
                      lat_param_quantiles[2, "a"], lat_param_quantiles[2, "b"]),
  P95_Sensitivity = c(inf_param_quantiles[3, "Topt"], inf_param_quantiles[3, "sigma_Wbar"], inf_param_quantiles[3, "gamma"], inf_param_quantiles[3, "lambda"], inf_param_quantiles[3, "k"],
                      lat_param_quantiles[3, "a"], lat_param_quantiles[3, "b"])
)

write_csv(sens_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "digitization_sensitivity_params.csv"))
cat("  -> Saved FigS4_digitization_sensitivity.png, .pdf and tables/digitization_sensitivity_params.csv\n\n")

cat("Step 6 completed successfully.\n")
cat("==============================================================================\n")
