# ==============================================================================
# Script: R/09_epidemic_pressure_indices.R
# Purpose: Step 9 Calculate Weather-Driven Epidemic Pressure Index (EPI_inf and EPI_spore),
#          Export Site Rankings and Climatology Summary Tables, and Generate
#          Publication Manuscript Figure 3 (Fig3_seasonal_climatology)
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

for (d in c(file.path(PROJ_ROOT, DIR_OUTPUT_TABLES),
            file.path(PROJ_ROOT, DIR_OUTPUT_FIGS),
            file.path(PROJ_ROOT, DIR_DATA_DERIVED))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("==============================================================================\n")
cat("STEP 9: WEATHER-DRIVEN EPIDEMIC PRESSURE INDEX (EPI)\n")
cat("==============================================================================\n\n")

cohorts_file <- file.path(PROJ_ROOT, DIR_DATA_DERIVED, "phase2_latent_cohorts.rds")
if (!file.exists(cohorts_file)) {
  stop("phase2_latent_cohorts.rds not found in data/derived/. Run R/08_latent_tracking.R first.")
}

cohorts <- readRDS(cohorts_file)
cat("Loaded", nrow(cohorts), "tracked infection cohorts.\n")

cat("[1/2] Calculating site-level annual & monthly epidemic metrics...\n")

annual_site_inf <- cohorts %>%
  group_by(cell_id, lat, lon, year = start_year) %>%
  summarise(
    episodes_count = n(),
    high_inf_episodes = sum(I_pot >= 0.5),
    EPI_inf = sum(I_pot),
    mean_lp_days = mean(lp_days, na.rm = TRUE),
    .groups = "drop"
  )

annual_site_spore <- cohorts %>%
  filter(completed) %>%
  group_by(cell_id, lat, lon, year = spore_year) %>%
  summarise(
    completed_count = n(),
    EPI_spore = sum(I_pot),
    .groups = "drop"
  )

annual_site_metrics <- full_join(
  annual_site_inf, annual_site_spore,
  by = c("cell_id", "lat", "lon", "year")
) %>%
  replace_na(list(EPI_spore = 0, completed_count = 0)) %>%
  rename(start_year = year)

write_csv(annual_site_metrics, file.path(PROJ_ROOT, DIR_DATA_DERIVED, "phase2_annual_site_metrics.csv"))

# Divide by N_YEARS (years without infection = 0)
cell_coords <- cohorts %>% distinct(cell_id, lat, lon)
site_summary <- annual_site_metrics %>%
  group_by(cell_id) %>%
  summarise(
    Mean_Annual_EPI_inf = sum(EPI_inf) / N_YEARS,
    Mean_Annual_EPI_spore = sum(EPI_spore) / N_YEARS,
    Mean_Annual_Episodes = sum(episodes_count) / N_YEARS,
    Mean_Annual_High_Inf_Episodes = sum(high_inf_episodes) / N_YEARS,
    Mean_Annual_Completed_Cohorts = sum(completed_count) / N_YEARS,
    Mean_LP_days = mean(mean_lp_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(cell_coords, by = "cell_id") %>%
  relocate(lat, lon, .after = cell_id) %>%
  arrange(desc(Mean_Annual_EPI_inf))

write_csv(site_summary, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "site_epidemic_rankings.csv"))
cat("  -> Saved outputs/tables/site_epidemic_rankings.csv (Mean Brazil EPI_inf =", 
    round(mean(site_summary$Mean_Annual_EPI_inf), 2), ")\n")

monthly_inf <- cohorts %>%
  group_by(cell_id, start_month) %>%
  summarise(monthly_EPI_inf = sum(I_pot) / N_YEARS, .groups = "drop")

monthly_spore <- cohorts %>%
  filter(completed) %>%
  group_by(cell_id, start_month = spore_month) %>%
  summarise(monthly_EPI_spore = sum(I_pot) / N_YEARS, .groups = "drop")

monthly_metrics <- full_join(monthly_inf, monthly_spore, by = c("cell_id", "start_month")) %>%
  replace_na(list(monthly_EPI_spore = 0, monthly_EPI_inf = 0)) %>%
  complete(cell_id, start_month = 1:12, fill = list(monthly_EPI_inf = 0, monthly_EPI_spore = 0))

monthly_clim_summary <- monthly_metrics %>%
  group_by(start_month) %>%
  summarise(
    Mean_Monthly_EPI_inf = mean(monthly_EPI_inf),
    SD_Monthly_EPI_inf = sd(monthly_EPI_inf),
    Mean_Monthly_EPI_spore = mean(monthly_EPI_spore),
    SD_Monthly_EPI_spore = sd(monthly_EPI_spore),
    .groups = "drop"
  )

lp_monthly_summary <- cohorts %>%
  filter(completed) %>%
  group_by(start_month) %>%
  summarise(Mean_Latent_Period_days = mean(lp_days), .groups = "drop")

monthly_clim_summary <- left_join(monthly_clim_summary, lp_monthly_summary, by = "start_month")
write_csv(monthly_clim_summary, file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "monthly_climatology_summary.csv"))
cat("  -> Saved outputs/tables/monthly_climatology_summary.csv\n\n")

# 2. Export Manuscript Figure 3 (Publication Standard: EPI_inf only, no sporulation)
cat("[2/2] Exporting Manuscript Figure 3 (Fig3_seasonal_climatology)...\n")

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )

# Panel A: Mean monthly EPI_inf across 65 cells with ±SD error bars (no legend, no sporulation line)
fig3a <- ggplot(monthly_clim_summary, aes(x = factor(start_month), y = Mean_Monthly_EPI_inf)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.6) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean_Monthly_EPI_inf - SD_Monthly_EPI_inf),
        ymax = Mean_Monthly_EPI_inf + SD_Monthly_EPI_inf),
    width = 0.25, color = "black", linewidth = 0.5
  ) +
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  labs(
    x = "Month",
    y = expression("Mean monthly " * EPI[inf])
  ) +
  theme_clean

# Panel B: Latent period duration (days) by month of infection initiation
comp_cohorts <- cohorts %>% filter(completed)

fig3b <- ggplot(comp_cohorts, aes(x = factor(start_month), y = lp_days)) +
  geom_boxplot(fill = "lightblue", outlier.alpha = 0.2, outlier.size = 0.8) +
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  labs(
    x = "Infection initiation month",
    y = "Latent period duration (days)"
  ) +
  theme_clean

# Composite 2x1 vertical stack (Panel A over Panel B)
fig3_combined <- (fig3a / fig3b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig3_seasonal_climatology.png"), fig3_combined, width = 8, height = 7.5, dpi = 300)
ggsave(file.path(PROJ_ROOT, DIR_OUTPUT_FIGS, "Fig3_seasonal_climatology.pdf"), fig3_combined, width = 8, height = 7.5)

cat("  -> Saved Fig3_seasonal_climatology.png and .pdf\n\n")

cat("Step 9 completed successfully.\n")
cat("==============================================================================\n")
