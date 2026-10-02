# ==============================================================================
# Script: R/08_latent_tracking.R
# Purpose: Step 8 Hourly Latent Development Cohort Tracking across 65 Sites
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_DATA_DERIVED),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 8: HOURLY LATENT DEVELOPMENT COHORT TRACKING\n")
cat("==============================================================================\n\n")

episodes_file <- file.path(PROJ_ROOT, DIR_DATA_DERIVED, "phase2_episodes_rh90.rds")
if (!file.exists(episodes_file)) {
  stop("phase2_episodes_rh90.rds not found in data/derived/. Run R/07_wetness_episodes.R first.")
}

episodes <- readRDS(episodes_file) %>%
  filter(I_pot > 0)

cat("Loaded", nrow(episodes), "infection cohorts with I_pot > 0 for latent tracking.\n")

weather_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_maravalia/hourly")
weather_files <- list.files(weather_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(weather_files) == 0) {
  subset_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_subset_example")
  weather_files <- list.files(subset_dir, pattern = "\\.rds$", full.names = TRUE)
  weather_files <- weather_files[!grepl("aus_site", basename(weather_files))]
}

track_site_cohorts <- function(site_file, site_episodes) {
  if (nrow(site_episodes) == 0) return(NULL)
  
  df <- readRDS(site_file)
  n_hours <- nrow(df)
  
  delta_D <- calc_rate_hourly(df$T2M)
  cum_D   <- c(0, cumsum(delta_D))
  
  init_indices <- site_episodes$end_idx + 1
  valid_mask <- init_indices <= (n_hours - 24)
  if (!any(valid_mask)) return(NULL)
  
  sub_episodes <- site_episodes[valid_mask, ]
  inits <- init_indices[valid_mask]
  
  n_cohorts <- length(inits)
  spore_idx <- integer(n_cohorts)
  completed <- logical(n_cohorts)
  
  for (c_idx in seq_len(n_cohorts)) {
    i_start <- inits[c_idx]
    target_cum <- cum_D[i_start] + 1.0
    match_pos <- findInterval(target_cum, cum_D)
    
    if (match_pos <= n_hours) {
      spore_idx[c_idx] <- match_pos
      completed[c_idx] <- TRUE
    } else {
      spore_idx[c_idx] <- NA_integer_
      completed[c_idx] <- FALSE
    }
  }
  
  sub_episodes %>%
    mutate(
      init_idx = inits,
      spore_idx = spore_idx,
      completed = completed,
      lp_hours = ifelse(completed, spore_idx - init_idx + 1, NA_real_),
      lp_days = lp_hours / 24,
      spore_year = ifelse(completed, df$YEAR[spore_idx], NA_integer_),
      spore_month = ifelse(completed, df$MO[spore_idx], NA_integer_),
      spore_day = ifelse(completed, df$DY[spore_idx], NA_integer_)
    )
}

cat("[1/2] Tracking latent development for all cohorts...\n")

site_split <- split(episodes, episodes$cell_id)

tracked_cohorts_list <- purrr::map(weather_files, function(f) {
  cell_id_val <- gsub("\\.rds$", "", basename(f))
  s_episodes <- site_split[[cell_id_val]]
  if (is.null(s_episodes) || nrow(s_episodes) == 0) return(NULL)
  track_site_cohorts(f, s_episodes)
})

tracked_cohorts <- bind_rows(tracked_cohorts_list)
cohorts_out_file <- file.path(PROJ_ROOT, DIR_DATA_DERIVED, "phase2_latent_cohorts.rds")
saveRDS(tracked_cohorts, cohorts_out_file)

cat("  -> Saved data/derived/phase2_latent_cohorts.rds\n")
cat("  - Total cohorts tracked:", nrow(tracked_cohorts), "\n")
cat("  - Completed cohorts (D >= 1.0):", sum(tracked_cohorts$completed), sprintf("(%.1f%%)\n", 100 * mean(tracked_cohorts$completed)))

comp_cohorts <- tracked_cohorts %>% filter(completed)
cat("  - Mean Latent Period (completed cohorts):", round(mean(comp_cohorts$lp_days), 2), "days\n")
cat("  - Latent Period Range:", round(min(comp_cohorts$lp_days), 2), "to", round(max(comp_cohorts$lp_days), 2), "days\n\n")

cat("[2/2] Exporting Latent Period Climatology Figure...\n")

fig6a <- ggplot(comp_cohorts, aes(x = factor(start_month), y = lp_days)) +
  geom_boxplot(fill = "lightblue", outlier.alpha = 0.2, outlier.size = 0.8) +
  labs(
    x = "Infection initiation month",
    y = "Latent period duration (days)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "diagnostic_latent_period_seasonal.png"), fig6a, width = 7, height = 4.5, dpi = 300)

cat("Step 8 completed successfully.\n")
cat("==============================================================================\n")
