# ==============================================================================
# Script: R/11_download_australian_nasapower.R
# Purpose: Download Hourly NASA POWER Weather Data (2001-2025) for 13 Australian
#          Validation Sites from Tomley & Evans (2004) Plant Pathology
# Project: Maravalia cryptostegiae Epidemiological Submodels (Australian Validation)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(nasapower)
  library(tibble)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

cat("==============================================================================\n")
cat("AUSTRALIAN VALIDATION — NASA POWER HOURLY WEATHER DATA DOWNLOAD (2001-2025)\n")
cat("==============================================================================\n\n")

# 1. Define / Load Australian Release & Monitoring Sites (Tomley & Evans 2004)
metadata_file <- file.path(PROJ_ROOT, "data/metadata/australian_sites_metadata.csv")

if (file.exists(metadata_file)) {
  aus_sites <- read_csv(metadata_file, show_col_types = FALSE)
} else {
  aus_sites <- tribble(
    ~site_id, ~site_name, ~region, ~lat, ~lon, ~inoculation_date, ~environment_type, ~observed_field_impact,
    "aus_site_01", "Chillagoe NP", "Tropical North QLD", -17.1517, 144.5244, "18/01/1995", "Inland / Karst", "Spread to 80 km by Aug 1996",
    "aus_site_02", "Lakefield NP", "Cape York Peninsula", -14.9356, 144.2047, "19/01/1995", "Moist Coastal Airflow", "High defoliation & 20-75% weed mortality",
    "aus_site_03", "Wrotham Park", "Gulf Region", -16.6342, 144.0042, "06/02/1995", "Gulf Savanna", "Rapid spread to 50 km in 8 months",
    "aus_site_04", "Rutland Plains", "Gulf Coast", -15.5458, 141.7258, "06/02/1995", "Gulf Lowlands", "Satellite spread >20 km detected",
    "aus_site_05", "Inkerman", "Burdekin Coast", -19.7431, 147.4647, "06/02/1995", "Sand Ridge / Clay Floodplain", "Severe rust damage; up to 75% mortality",
    "aus_site_06", "Delta Downs", "Gulf of Carpentaria", -17.1008, 140.9167, "07/02/1995", "Coastal Gulf", "Spread to 70 km in 18 months",
    "aus_site_07", "Georgetown", "Etheridge Inland", -18.2917, 143.5475, "07/02/1995", "Elevated Dry Inland", "Moderate rust; shorter wet season duration",
    "aus_site_08", "McLeod River", "Cooktown Hinterland", -16.4833, 145.2333, "08/02/1995", "Riparian / Watercourse", "Heavy rust (47 pustules cm^-2); total defoliation",
    "aus_site_09", "Strathmore", "Gulf / Croydon", -17.8500, 142.3667, "08/03/1995", "Gulf Savanna", "Spread to 45 km by Aug 1996",
    "aus_site_10", "Rockhampton", "Central QLD Coast", -23.3750, 150.5117, "18/01/1995", "Subtropical Coast", "Rapid spread over 100 km in 14 months",
    "aus_site_11", "Charters Towers", "Inland Central QLD", -20.0764, 146.2611, "1995-1996", "Inland Research Hub (TWRC)", "Key monitoring hub; heavy rust in wet years",
    "aus_site_12", "Hughenden", "Dry Inland QLD", -20.8433, 144.2008, "1995-1996", "Semi-Arid (<600mm rain)", "Inadequate control in most years due to dry air",
    "aus_site_13", "Laura Lea", "Cape York Peninsula", -15.5583, 144.4447, "1995-1996", "Moist South-East Airflow", "Enhanced rust intensity; 20% mortality by 1999"
  ) %>%
    mutate(point_id = sprintf("%s_lat_%.4f_lon_%.4f", site_id, lat, lon))
  
  write_csv(aus_sites, metadata_file)
}

hourly_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_australia/hourly")
dir.create(hourly_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loaded metadata for", nrow(aus_sites), "Australian field sites from Tomley & Evans (2004).\n")

# 2. Download Configuration (2001-2025, Hourly)
years <- 2001:2025
power_pars <- c("T2M", "RH2M", "T2MDEW", "PRECTOTCORR")

cat("[Downloading Hourly Weather Data for 13 Australian Locations (2001-2025)]\n")

for (i in seq_len(nrow(aus_sites))) {
  cell <- aus_sites[i, ]
  outfile <- file.path(hourly_dir, paste0(cell$point_id, ".rds"))
  
  if (file.exists(outfile)) {
    cat(sprintf("[%02d/%02d] Skipping %s (%s) — file already exists.\n", i, nrow(aus_sites), cell$site_name, cell$point_id))
    next
  }
  
  cat(sprintf("[%02d/%02d] Downloading %s (Lat: %.4f, Lon: %.4f)...\n", i, nrow(aus_sites), cell$site_name, cell$lat, cell$lon))
  
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
      mutate(
        site_id = cell$site_id,
        site_name = cell$site_name,
        region = cell$region,
        environment_type = cell$environment_type,
        original_lat = cell$lat,
        original_lon = cell$lon
      )
    saveRDS(out, outfile)
    cat(sprintf("   Saved: %s | %d hourly records.\n\n", outfile, nrow(out)))
  }
  Sys.sleep(1.5)
}

cat("Australian NASA POWER weather data download pipeline ready.\n")
cat("==============================================================================\n")
