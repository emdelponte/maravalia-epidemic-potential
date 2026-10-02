# ==============================================================================
# Script: R/07_wetness_episodes.R
# Purpose: Step 7 Extraction of Microclimatic Wetness Episodes and Quantification
#          of Infection Potential using Revised Centered Model 3 across 65 Sites (2001-2025)
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS),
            file.path(PROJ_ROOT, DIR_DATA_DERIVED))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 7: EXTRACTION OF HOURLY WETNESS EPISODES (2001-2025)\n")
cat("==============================================================================\n\n")

# 1. Load Primary Infection Submodel
inf_model_path <- file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_primary.rds")
if (!file.exists(inf_model_path)) inf_model_path <- file.path(PROJ_ROOT, "models_v2/model_infection_primary.rds")
fit_primary <- readRDS(inf_model_path)
coefs <- coef(fit_primary)

cat("Loaded Model 3 Parameters:\n")
print(round(coefs, 4))
cat("\n")

weather_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_maravalia/hourly")
weather_files <- list.files(weather_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(weather_files) == 0) {
  # Fallback to example subset if full raw weather is not downloaded
  subset_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_subset_example")
  weather_files <- list.files(subset_dir, pattern = "\\.rds$", full.names = TRUE)
  weather_files <- weather_files[!grepl("aus_site", basename(weather_files))]
  if (length(weather_files) > 0) {
    warning("Full weather record not found; running on subset example site.")
  } else {
    stop("No weather files found in data/raw/nasapower_maravalia/hourly. Run scripts/get_data.R first.")
  }
}

cat("Found", length(weather_files), "hourly meteorological data files.\n")

calc_infection_potential <- function(T_wet, W) calc_I(T_wet, W, coefs)

extract_episodes_fast <- function(df, rh_thresh = RH_THRESH, include_precip = TRUE) {
  is_wet <- wet_hours(df, rh = rh_thresh, p = if (include_precip) PRECIP_THRESH else NA)
  if (!any(is_wet)) return(NULL)
  ep <- wet_runs(is_wet)
  cs <- c(0, cumsum(df$T2M))
  mean_T <- (cs[ep$end_idx + 1] - cs[ep$start_idx]) / ep$wet_h
  tibble(
    cell_id = df$cell_id[1],
    lat = df$LAT[1],
    lon = df$LON[1],
    start_year = df$YEAR[ep$start_idx],
    start_month = df$MO[ep$start_idx],
    start_day = df$DY[ep$start_idx],
    start_hour = df$HR[ep$start_idx],
    end_hour = df$HR[ep$end_idx],
    start_idx = ep$start_idx,
    end_idx = ep$end_idx,
    wet_h = ep$wet_h,
    T_wet = mean_T,
    I_pot = calc_infection_potential(mean_T, ep$wet_h)
  )
}

# 1. Sensitivity Analysis across RH Thresholds (85%, 90%, 95%)
cat("[1/2] Extracting episodes across RH thresholds (85%, 90%, 95%)...\n")
rh_thresholds <- c(85, 90, 95)
rh_results <- list()
n_cells <- length(weather_files)

for (rh in rh_thresholds) {
  cat("  - Processing RH threshold:", rh, "%\n")
  site_episodes <- purrr::map(weather_files, function(f) {
    df <- readRDS(f)
    extract_episodes_fast(df, rh_thresh = rh, include_precip = TRUE)
  })
  all_episodes <- bind_rows(site_episodes)
  
  summary_rh <- all_episodes %>%
    summarise(
      RH_Threshold = paste0(rh, "%"),
      Total_Episodes = n(),
      Episodes_ge_6h = sum(wet_h >= 6),
      Episodes_ge_8h = sum(wet_h >= 8),
      Episodes_ge_12h = sum(wet_h >= 12),
      Mean_Wetness_Duration_h = mean(wet_h),
      Mean_Wetness_Temp_C = mean(T_wet),
      Total_Infection_Potential = sum(I_pot),
      Mean_Annual_Infection_Potential_per_Site = sum(I_pot) / (n_cells * N_YEARS)
    )
  
  rh_results[[as.character(rh)]] <- list(summary = summary_rh, data = all_episodes)
}

rh_comp_table <- bind_rows(purrr::map(rh_results, "summary"))
write_csv(rh_comp_table, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "rh_sensitivity_comparison.csv"))

episodes_rh90 <- rh_results[["90"]]$data
saveRDS(episodes_rh90, file.path(PROJ_ROOT, DIR_DATA_DERIVED, "phase2_episodes_rh90.rds"))
cat("  -> Saved data/derived/phase2_episodes_rh90.rds with N =", nrow(episodes_rh90), "total episodes.\n")
cat("  -> Saved outputs/tables/rh_sensitivity_comparison.csv\n\n")

# 2. Export Diagnostic Figure
cat("[2/2] Exporting Diagnostic Wetness Figure...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 12)
  )

fig5a <- ggplot(rh_comp_table, aes(x = RH_Threshold, y = Mean_Annual_Infection_Potential_per_Site, fill = RH_Threshold)) +
  geom_col(width = 0.45, color = "black") +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  labs(
    x = "Relative humidity threshold",
    y = expression("Mean annual cumulative " * EPI[inf])
  ) +
  theme_clean

fig5b <- ggplot(episodes_rh90 %>% filter(wet_h >= 6), aes(x = wet_h, y = I_pot, color = T_wet)) +
  geom_point(alpha = 0.3, size = 1.1) +
  scale_color_viridis_c(name = expression("Mean " * T[wet] * " (°C)")) +
  labs(
    x = "Wetness duration (h)",
    y = expression("Relative infection " * I(T[wet], W))
  ) +
  theme_clean

fig5_combined <- (fig5a + fig5b) + plot_annotation(tag_levels = 'A')

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "diagnostic_wetness_episodes.png"), fig5_combined, width = 10, height = 4.2, dpi = 300)

cat("Step 7 completed successfully.\n")
cat("==============================================================================\n")
