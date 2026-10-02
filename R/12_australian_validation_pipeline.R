# ==============================================================================
# Script: R/12_australian_validation_pipeline.R
# Purpose: Step 12 Model Validation under Australian Field Conditions (Tomley & Evans 2004)
#          Processes 25-Year Hourly Weather (2001-2025) across 13 Queensland Sites
#          using Revised Centered Model 3 (Dynamic Breadth)
#          Exports Validation Summary Table and Manuscript Figure 5 (Fig5_australia)
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(knitr)
  library(maps)
  library(ggrepel)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 12: AUSTRALIAN FIELD VALIDATION PIPELINE (TOMLEY & EVANS 2004)\n")
cat("==============================================================================\n\n")

# 1. Load Metadata & Parameters
metadata_path <- file.path(PROJ_ROOT, "data/metadata/australian_sites_metadata.csv")
if (!file.exists(metadata_path)) metadata_path <- file.path(PROJ_ROOT, "nasapower_australia/australian_sites_metadata.csv")

aus_metadata <- read_csv(metadata_path, show_col_types = FALSE)
cat("Loaded metadata for", nrow(aus_metadata), "Australian field locations.\n")

inf_model_path <- file.path(PROJ_ROOT, DIR_OUTPUT_MODELS, "model_infection_primary.rds")
if (!file.exists(inf_model_path)) inf_model_path <- file.path(PROJ_ROOT, "models_v2/model_infection_primary.rds")
fit_primary <- readRDS(inf_model_path)
coefs <- coef(fit_primary)

calc_I_pot <- function(T_val, W_val) calc_I(T_val, W_val, coefs)
calc_rate  <- function(T_val) calc_rate_hourly(T_val)

# 2. Extract Wetness Episodes & Calculate Infection Potential for Australian Sites
cat("[1/4] Processing hourly weather & extracting wetness episodes...\n")

aus_weather_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_australia/hourly")
aus_weather_files <- list.files(aus_weather_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(aus_weather_files) == 0) {
  subset_dir <- file.path(PROJ_ROOT, "data/raw/nasapower_subset_example")
  aus_weather_files <- list.files(subset_dir, pattern = "\\.rds$", full.names = TRUE)
  aus_weather_files <- aus_weather_files[grepl("aus_site", basename(aus_weather_files))]
}

cat("Found", length(aus_weather_files), "hourly weather files in", aus_weather_dir, "\n")

site_episodes_list <- purrr::map(aus_weather_files, function(f) {
  df <- readRDS(f)
  n_rows <- nrow(df)
  if (n_rows == 0) return(NULL)
  
  cell_id_val <- df$site_id[1]
  site_name_val <- df$site_name[1]
  region_val <- df$region[1]
  env_type_val <- df$environment_type[1]
  
  is_wet <- wet_hours(df)
  ep <- wet_runs(is_wet)
  episodes_df <- tibble(start_idx = ep$start_idx, end_idx = ep$end_idx, duration_h = ep$wet_h) %>%
    filter(duration_h >= W_MIN)
  
  if (nrow(episodes_df) == 0) return(NULL)
  
  episodes_df %>%
    mutate(
      site_id = cell_id_val,
      site_name = site_name_val,
      region = region_val,
      environment_type = env_type_val,
      start_year = df$YEAR[start_idx],
      start_month= df$MO[start_idx],
      start_day  = df$DY[start_idx],
      start_hour = df$HR[start_idx],
      T_mean = purrr::map2_dbl(start_idx, end_idx, ~mean(df$T2M[.x:.y], na.rm = TRUE)),
      RH_mean = purrr::map2_dbl(start_idx, end_idx, ~mean(df$RH2M[.x:.y], na.rm = TRUE)),
      I_pot = purrr::map2_dbl(T_mean, duration_h, calc_I_pot)
    )
})

aus_episodes <- bind_rows(site_episodes_list)
cat("  - Total wetness episodes extracted (W >= 6h):", nrow(aus_episodes), "\n")
cat("  - Mean episodes per site per year:", round(nrow(aus_episodes) / (length(aus_weather_files) * N_YEARS), 1), "\n\n")

# 3. Track Latent Cohort Development across Australian Sites
cat("[2/4] Tracking latent period development for all cohorts...\n")

site_split <- split(aus_episodes, aus_episodes$site_id)

tracked_list <- purrr::map(aus_weather_files, function(f) {
  df <- readRDS(f)
  if (nrow(df) == 0) return(NULL)
  
  site_id_val <- df$site_id[1]
  s_episodes <- site_split[[site_id_val]]
  if (is.null(s_episodes) || nrow(s_episodes) == 0) return(NULL)
  
  n_hours <- nrow(df)
  delta_D <- calc_rate(df$T2M)
  cum_D   <- c(0, cumsum(delta_D))
  
  inits <- s_episodes$end_idx + 1
  valid_mask <- inits <= (n_hours - 24)
  if (!any(valid_mask)) return(NULL)
  
  sub_ep <- s_episodes[valid_mask, ]
  inits <- inits[valid_mask]
  
  n_c <- length(inits)
  spore_idx <- integer(n_c)
  completed <- logical(n_c)
  
  for (c_idx in seq_len(n_c)) {
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
  
  sub_ep %>%
    mutate(
      init_idx = inits,
      spore_idx = spore_idx,
      completed = completed,
      lp_hours = ifelse(completed, spore_idx - init_idx + 1, NA_real_),
      lp_days = lp_hours / 24,
      spore_year = ifelse(completed, df$YEAR[spore_idx], NA_integer_),
      spore_month = ifelse(completed, df$MO[spore_idx], NA_integer_)
    )
})

aus_cohorts <- bind_rows(tracked_list)
comp_aus <- aus_cohorts %>% filter(completed)
cat("  - Total cohorts tracked:", nrow(aus_cohorts), "\n")
cat("  - Completed cohorts:", sum(aus_cohorts$completed), sprintf("(%.2f%%)\n", 100 * mean(aus_cohorts$completed)))
cat("  - Mean Latent Period:", round(mean(comp_aus$lp_days), 2), "days\n\n")

# 4. Compute Site-Level Epidemic Pressure Indices & Validation Summaries
cat("[3/4] Computing site-level annual & monthly epidemic pressure metrics...\n")

aus_annual_inf <- aus_cohorts %>%
  group_by(site_id, site_name, region, environment_type, year = start_year) %>%
  summarise(
    episodes_count = n(),
    high_inf_episodes = sum(I_pot >= 0.5),
    EPI_inf = sum(I_pot),
    mean_lp_days = mean(lp_days, na.rm = TRUE),
    .groups = "drop"
  )

aus_annual_spore <- aus_cohorts %>%
  filter(completed) %>%
  group_by(site_id, site_name, year = spore_year) %>%
  summarise(
    EPI_spore = sum(I_pot),
    completed_count = n(),
    .groups = "drop"
  )

aus_annual_metrics <- full_join(aus_annual_inf, aus_annual_spore, by = c("site_id", "site_name", "year")) %>%
  replace_na(list(EPI_spore = 0, completed_count = 0)) %>%
  rename(start_year = year)

# Divide by N_YEARS (years without infection = 0)
site_info <- aus_episodes %>% distinct(site_id, site_name, region, environment_type)
aus_site_summary <- aus_annual_metrics %>%
  group_by(site_id) %>%
  summarise(
    Mean_Annual_EPI_inf = sum(EPI_inf, na.rm = TRUE) / N_YEARS,
    Mean_Annual_EPI_spore = sum(EPI_spore, na.rm = TRUE) / N_YEARS,
    Mean_Annual_Episodes = sum(episodes_count, na.rm = TRUE) / N_YEARS,
    Mean_Annual_High_Inf_Episodes = sum(high_inf_episodes, na.rm = TRUE) / N_YEARS,
    Mean_LP_days = mean(mean_lp_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(site_info, by = "site_id") %>%
  relocate(site_name, region, environment_type, .after = site_id) %>%
  left_join(aus_metadata %>% select(site_id, lat, lon, observed_field_impact), by = "site_id") %>%
  arrange(desc(Mean_Annual_EPI_inf))

write_csv(aus_site_summary, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "australian_model_validation_summary.csv"))
cat("  -> Saved outputs/tables/australian_model_validation_summary.csv\n")

# Monthly climatology
monthly_inf <- aus_cohorts %>%
  group_by(site_name, month = start_month) %>%
  summarise(monthly_EPI_inf = sum(I_pot) / N_YEARS, .groups = "drop")

monthly_spore <- aus_cohorts %>%
  filter(completed) %>%
  group_by(site_name, month = spore_month) %>%
  summarise(monthly_EPI_spore = sum(I_pot) / N_YEARS, .groups = "drop")

aus_monthly_clim <- full_join(monthly_inf, monthly_spore, by = c("site_name", "month")) %>%
  replace_na(list(monthly_EPI_spore = 0, monthly_EPI_inf = 0)) %>%
  complete(site_name, month = 1:12, fill = list(monthly_EPI_inf = 0, monthly_EPI_spore = 0))

# 5. Export Diagnostic Validation Figures
cat("[4/4] Generating & exporting Australian validation figures (Fig5_australia)...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 12)
  )

# Panel A (Top): Monthly Climatology across Australian sites (EPI_inf bars only)
overall_monthly <- aus_monthly_clim %>%
  group_by(month) %>%
  summarise(
    mean_inf = mean(monthly_EPI_inf),
    mean_spore = mean(monthly_EPI_spore),
    .groups = "drop"
  )

fig5a_monthly <- ggplot(overall_monthly, aes(x = factor(month), y = mean_inf)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.6) +
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  labs(
    x = "Month",
    y = expression("Mean monthly " * EPI[inf])
  ) +
  theme_clean

# Panel B (Middle): Site Risk Ranking & Field Impact ("High impact" vs "Other sites")
aus_site_summary <- aus_site_summary %>%
  mutate(
    impact_cat = if_else(
      grepl("mortality|defoliation|Heavy rust", observed_field_impact, ignore.case = TRUE),
      "High impact",
      "Other sites"
    )
  )

fig5b_ranking <- ggplot(aus_site_summary, aes(x = reorder(site_name, Mean_Annual_EPI_inf), y = Mean_Annual_EPI_inf, fill = impact_cat)) +
  geom_col(width = 0.65) +
  coord_flip() +
  scale_fill_manual(name = "Field impact", values = c("High impact" = "darkgreen", "Other sites" = "goldenrod3")) +
  labs(
    x = "Field release site",
    y = expression("Mean annual " * EPI[inf])
  ) +
  theme_clean +
  theme(legend.position = "bottom")

# Panel C (Bottom): Spatial map of Queensland release locations
world_map <- map_data("world", region = "Australia")

fig5c_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey93", color = "grey65", linewidth = 0.4) +
  geom_point(data = aus_site_summary, aes(x = lon, y = lat, fill = Mean_Annual_EPI_inf, size = Mean_Annual_EPI_inf), shape = 21, color = "black", stroke = 0.7) +
  geom_text_repel(data = aus_site_summary, aes(x = lon, y = lat, label = site_name), size = 3.2, fontface = "bold", box.padding = 0.35, point.padding = 0.3, max.overlaps = 20) +
  scale_fill_viridis_c(name = expression("Mean annual " * EPI[inf]), option = "plasma") +
  scale_size_continuous(guide = "none", range = c(3, 8)) +
  coord_quickmap(xlim = c(139, 153), ylim = c(-24.5, -13.5)) +
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°S)"
  ) +
  theme_clean +
  theme(
    legend.position = "right",
    legend.box.margin = margin(l = 5)
  )

# Multi-panel composite (3x1 vertical stack: Monthly -> Ranking -> Map at bottom)
fig5_combined <- (fig5a_monthly / fig5b_ranking / fig5c_map) +
  plot_layout(heights = c(1, 1.2, 1.4)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig5_australia.png"), fig5_combined, width = 8.5, height = 14, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig5_australia.pdf"), fig5_combined, width = 8.5, height = 14)

cat("  -> Saved Fig5_australia.png and .pdf\n\n")

high_impact_EPI <- aus_site_summary %>% filter(impact_cat == "High impact") %>% pull(Mean_Annual_EPI_inf)
mod_impact_EPI  <- aus_site_summary %>% filter(impact_cat == "Other sites") %>% pull(Mean_Annual_EPI_inf)
t_res <- t.test(high_impact_EPI, mod_impact_EPI)
w_res <- wilcox.test(high_impact_EPI, mod_impact_EPI, exact = TRUE)

cat("--- MODEL ACCURACY & VALIDATION STATISTICAL SUMMARY ---\n")
cat("Mean Annual EPI_inf for High Impact Sites:", round(mean(high_impact_EPI), 2), "\n")
cat("Mean Annual EPI_inf for Other (non-high-impact) sites:", round(mean(mod_impact_EPI), 2), "\n")
cat("Difference in Predicted Infection Potential:", round(mean(high_impact_EPI) - mean(mod_impact_EPI), 2), "EPI units\n")
cat(sprintf("Welch t = %.2f, df = %.1f, P = %.3f\n", t_res$statistic, t_res$parameter, t_res$p.value))
cat(sprintf("Mann-Whitney W = %.0f, P = %.3f (n = %d vs %d)\n\n", w_res$statistic, w_res$p.value, length(high_impact_EPI), length(mod_impact_EPI)))

cat("Step 12 completed successfully.\n")
cat("==============================================================================\n")
