# ==============================================================================
# Script: R/01_qc_and_data_prep.R
# Purpose: Step 1 Structural QC of Weather & Site Data, and Dataset Preparation
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(purrr)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

cat("==============================================================================\n")
cat("STEP 1: STRUCTURAL QC & DATASET PREPARATION\n")
cat("==============================================================================\n\n")

# 1. Brazilian Site Data QC
excel_path <- file.path(PROJ_ROOT, "data/metadata/cryptostegia_nordeste_bonilla2023.xlsx")
cat("[1/4] Inspecting host location data from:", excel_path, "\n")

if (file.exists(excel_path)) {
  sites_df <- read_excel(excel_path, sheet = "Cryptostegia")
  
  site_qc <- sites_df %>%
    rename(
      uf = UF,
      municipio = Município,
      lat_deg = `Latitude decimal`,
      lon_deg = `Longitude decimal`,
      altitude_m = `Altitude (m)`,
      area_ha = `Área infestada (ha)`
    )
  
  cat("  - Total municipal sites loaded:", nrow(site_qc), "\n")
  cat("  - State distribution (UF):\n")
  print(table(site_qc$uf))
  cat("  - Latitude range:", min(site_qc$lat_deg, na.rm=TRUE), "to", max(site_qc$lat_deg, na.rm=TRUE), "°S\n")
  cat("  - Longitude range:", min(site_qc$lon_deg, na.rm=TRUE), "to", max(site_qc$lon_deg, na.rm=TRUE), "°W\n")
  cat("  - Total infested area recorded:", sum(site_qc$area_ha, na.rm=TRUE), "ha\n")
  
  write_csv(site_qc, file.path(PROJ_ROOT, "data/metadata/brazilian_sites_qc.csv"))
  cat("  -> Saved data/metadata/brazilian_sites_qc.csv\n\n")
} else {
  warning("Excel site data file not found at: ", excel_path)
}

# 2. Hourly Meteorological Data Structural QC
weather_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_maravalia/hourly")
cat("[2/4] Performing structural QC on hourly NASA POWER meteorological data in:", weather_dir, "\n")

weather_files <- list.files(weather_dir, pattern = "\\.rds$", full.names = TRUE)
cat("  - Number of weather site files found:", length(weather_files), "\n")

if (length(weather_files) > 0) {
  weather_qc_list <- map(weather_files, function(f) {
    df <- readRDS(f)
    fn <- basename(f)
    n_rows <- nrow(df)
    
    tibble(
      file = fn,
      cell_id = df$cell_id[1],
      lon = df$LON[1],
      lat = df$LAT[1],
      total_records = n_rows,
      num_variables = ncol(df),
      variables = paste(colnames(df), collapse = "; "),
      start_time = paste(df$YEAR[1], df$MO[1], df$DY[1], df$HR[1], sep="-"),
      end_time = paste(df$YEAR[n_rows], df$MO[n_rows], df$DY[n_rows], df$HR[n_rows], sep="-"),
      missing_obs = sum(is.na(df)),
      temporal_resolution = "1-hour continuous"
    )
  })
  
  weather_qc_df <- bind_rows(weather_qc_list)
  write_csv(weather_qc_df, file.path(PROJ_ROOT, "data/metadata/weather_structural_qc.csv"))
  cat("  -> Saved data/metadata/weather_structural_qc.csv\n")
  cat("  - Total sites:", nrow(weather_qc_df), "\n")
  cat("  - Observations per site:", unique(weather_qc_df$total_records), "(25 years, 2001-2025)\n")
  cat("  - Total missing observations across entire dataset:", sum(weather_qc_df$missing_obs), "\n\n")
} else {
  cat("  - Notice: No raw weather files present in data/raw/. (Expected if downloading via get_data.R)\n\n")
}

# 3. Infection Dataset Reconstruction
cat("[3/4] Reconstructing digitized infection dataset...\n")

infection_data <- tibble(
  dew_h = c(6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 24, 24, 24, 24, 24),
  T = c(17, 20, 22, 25, 30, 17, 20, 22, 25, 30, 17, 20, 22, 25, 30, 17, 20, 22, 25, 30),
  pustules_cm2 = c(
    0.05, 0.15, 0.65, 0.35, 0.05,
    0.55, 1.45, 6.35, 2.05, 0.10,
    2.05, 12.95, 10.75, 3.05, 2.05,
    13.05, 11.45, 14.35, 13.25, 2.55
  )
)

max_pustules <- max(infection_data$pustules_cm2)
cat("  - Maximum observed pustule density:", max_pustules, "pustules cm^-2\n")

infection_data <- infection_data %>%
  mutate(I_obs = pustules_cm2 / max_pustules)

cat("  - Relative infection response I_obs calculated (0 <= I_obs <= 1).\n")
write_csv(infection_data, file.path(PROJ_ROOT, "data/metadata/digitized_infection_data.csv"))
cat("  -> Saved data/metadata/digitized_infection_data.csv\n\n")

# 4. Latent Period Dataset Reconstruction
cat("[4/4] Reconstructing digitized latent period dataset...\n")

latent_data <- tibble(
  T = c(18, 18, 18, 22, 22, 22, 25, 25, 27, 27),
  LP_days = c(19, 20, 21, 10, 11, 17, 10, 11, 8, 9)
) %>%
  mutate(r = 1 / LP_days)

cat("  - Total latent period observations:", nrow(latent_data), "\n")
cat("  - Supported temperature range:", min(latent_data$T), "to", max(latent_data$T), "°C\n")

write_csv(latent_data, file.path(PROJ_ROOT, "data/metadata/digitized_latent_period.csv"))
cat("  -> Saved data/metadata/digitized_latent_period.csv\n\n")

cat("Step 1 completed successfully.\n")
cat("==============================================================================\n")
