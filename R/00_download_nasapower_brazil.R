# ==============================================================================
# Script: R/00_download_nasapower_brazil.R
# Purpose: Download hourly NASA POWER meteorological data (2001-2025) for 65
#          grid cells in Northeastern Brazil covering Cryptostegia madagascariensis
#          infestation sites (Bonilla et al. 2023).
# Status: NEWLY WRITTEN for reproducible repository build (Step 2).
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(nasapower)
})

cat("==============================================================================\n")
cat("BRAZIL WEATHER DATA DOWNLOAD (NASA POWER HOURLY 2001-2025)\n")
cat("[NOTE: Newly written automated retrieval script]\n")
cat("==============================================================================\n\n")

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

out_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_maravalia/hourly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load target grid cell coordinates
qc_meta_file <- file.path(PROJ_ROOT, "data/metadata/weather_structural_qc.csv")
site_meta_file <- file.path(PROJ_ROOT, "data/metadata/brazilian_sites_qc.csv")

if (file.exists(qc_meta_file)) {
  cells <- read_csv(qc_meta_file, show_col_types = FALSE) %>%
    distinct(cell_id, lon, lat, file) %>%
    arrange(lat, lon)
} else if (file.exists(site_meta_file)) {
  pts <- read_csv(site_meta_file, show_col_types = FALSE)
  cells <- pts %>%
    filter(!is.na(lat_deg), !is.na(lon_deg)) %>%
    mutate(
      lat = round(lat_deg, 1),
      lon = round(lon_deg, 1),
      cell_id = sprintf("lat_%.1f_lon_%.1f", lat, lon),
      file = paste0(cell_id, ".rds")
    ) %>%
    distinct(cell_id, lon, lat, file) %>%
    arrange(lat, lon)
} else {
  stop("Grid cell coordinate metadata file not found in data/metadata/.")
}

cat("Target weather grid cells to download:", nrow(cells), "\n")
years <- 2001:2025
power_pars <- c("T2M", "RH2M", "T2MDEW", "PRECTOTCORR")

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  outfile <- file.path(out_dir, cell$file)
  
  if (file.exists(outfile)) {
    cat(sprintf("[%02d/%02d] Skipping %s — already exists.\n", i, nrow(cells), cell$cell_id))
    next
  }
  
  cat(sprintf("[%02d/%02d] Downloading %s (Lat: %.2f, Lon: %.2f)...\n", i, nrow(cells), cell$cell_id, cell$lat, cell$lon))
  
  out <- map_dfr(years, function(y) {
    cat("   Year:", y, "\n")
    tryCatch({
      x <- get_power(
        community = "AG",
        lonlat = c(cell$lon, cell$lat),
        pars = power_pars,
        dates = c(paste0(y, "-01-01"), paste0(y, "-12-31")),
        temporal_api = "hourly",
        time_standard = "UTC"
      )
      Sys.sleep(1.2)
      x
    }, error = function(e) {
      cat("   Warning: Download failed for year", y, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  })
  
  if (!is.null(out) && nrow(out) > 0) {
    out <- out %>%
      mutate(cell_id = cell$cell_id)
    saveRDS(out, outfile)
    cat(sprintf("   Saved: %s | %d hourly records.\n\n", outfile, nrow(out)))
  }
  Sys.sleep(1.5)
}

cat("Brazilian NASA POWER download pipeline ready.\n")
cat("==============================================================================\n")
