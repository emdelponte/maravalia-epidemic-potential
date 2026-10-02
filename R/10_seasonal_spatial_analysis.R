# ==============================================================================
# Script: R/10_seasonal_spatial_analysis.R
# Purpose: Step 10 Seasonal Climatology & Spatial Risk Mapping across Northeastern Brazil
#          Exports Publication Manuscript Figure 4 (Fig4_spatial_maps) to outputs/figures/
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(rnaturalearth)
  library(sf)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."
source(file.path(PROJ_ROOT, "R/maravalia_functions.R"))

fig_dir <- file.path(PROJ_ROOT, DIR_OUTPUT_FIGS)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("STEP 10: SEASONAL CLIMATOLOGY & SPATIAL RISK MAPPING\n")
cat("==============================================================================\n\n")

rankings_file <- file.path(PROJ_ROOT, DIR_OUTPUT_TABLES, "site_epidemic_rankings.csv")
if (!file.exists(rankings_file)) {
  stop("site_epidemic_rankings.csv not found in outputs/tables/. Run R/09_epidemic_pressure_indices.R first.")
}

site_rankings <- read_csv(rankings_file, show_col_types = FALSE)

host_sites_file <- file.path(PROJ_ROOT, "data/metadata/brazilian_sites_qc.csv")
if (!file.exists(host_sites_file)) host_sites_file <- file.path(PROJ_ROOT, "data/brazilian_sites_qc.csv")
host_sites <- read_csv(host_sites_file, show_col_types = FALSE)

host_localities_scen_path <- file.path(PROJ_ROOT, "data/metadata/host_localities_epi_scenarios.csv")
if (file.exists(host_localities_scen_path)) {
  host_scen <- read_csv(host_localities_scen_path, show_col_types = FALSE)
} else {
  host_scen <- NULL
}

cat("Loaded", nrow(site_rankings), "weather grid cell sites and", nrow(host_sites), "host infestation records.\n\n")

# 1. Spatial Map of Northeastern Brazil (EPI_inf only)
cat("[1/2] Loading Brazilian state boundaries and generating spatial risk maps...\n")

br_states <- ne_states(country = "Brazil", returnclass = "sf")
visible_states_list <- c("CE", "RN", "PB", "PI")
br_ne <- br_states %>% filter(postal %in% visible_states_list)

theme_map <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 12)
  )

host_sites_clean <- host_sites %>% filter(!is.na(lat_deg) & !is.na(lon_deg))

fig4a <- ggplot() +
  geom_sf(data = br_states, fill = "gray96", color = "gray55", linewidth = 0.5) +
  geom_point(
    data = site_rankings,
    aes(x = lon, y = lat, color = Mean_Annual_EPI_inf, size = Mean_Annual_EPI_inf),
    alpha = 0.85
  ) +
  geom_point(
    data = host_sites_clean,
    aes(x = lon_deg, y = lat_deg),
    shape = 17, color = "red", size = 2, alpha = 0.75
  ) +
  geom_sf_text(data = br_ne, aes(label = postal), fontface = "bold", color = "gray25", size = 3.8) +
  scale_color_viridis_c(
    name = expression("Mean annual " * EPI[inf]),
    option = "viridis",
    limits = c(0, 120)
  ) +
  scale_size_continuous(range = c(2.5, 6), guide = "none") +
  coord_sf(xlim = c(-43.2, -35.8), ylim = c(-8.2, -2.5), expand = FALSE) +
  labs(
    x = "Longitude (°W)",
    y = "Latitude (°S)"
  ) +
  theme_map +
  theme(
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 10)
  )

# 2. Fig 4B: log10(infested area) vs EPI_inf with linear trend & statistical annotations
cat("[2/2] Generating Infested Area vs EPI_inf Correlation Plot...\n")

if (!is.null(host_scen)) {
  df_corr <- host_scen %>%
    rename(area_ha = `Infested Area (ha)`, EPI_val = S2_RH90_cap24, uf = State) %>%
    filter(!is.na(area_ha) & !is.na(EPI_val)) %>%
    mutate(log10_area = log10(area_ha))

  # Locality level (n = 73)
  r_loc <- cor.test(df_corr$log10_area, df_corr$EPI_val, method = "pearson")
  rho_loc <- cor.test(df_corr$log10_area, df_corr$EPI_val, method = "spearman")

  # Grid cell level (n = 28 unique cells)
  cell_agg <- df_corr %>%
    group_by(cell) %>%
    summarise(
      mean_log_area = mean(log10_area),
      EPI_val = first(EPI_val),
      .groups = "drop"
    )
  r_cell <- cor.test(cell_agg$mean_log_area, cell_agg$EPI_val, method = "pearson")
  rho_cell <- cor.test(cell_agg$mean_log_area, cell_agg$EPI_val, method = "spearman")

  n_loc <- nrow(df_corr)
  n_cell <- nrow(cell_agg)
} else {
  find_nearest_epi <- function(h_lat, h_lon, site_df) {
    if (is.na(h_lat) || is.na(h_lon)) return(NA_real_)
    dist <- sqrt((site_df$lat - h_lat)^2 + (site_df$lon - h_lon)^2)
    site_df$Mean_Annual_EPI_inf[which.min(dist)]
  }
  find_nearest_cell <- function(h_lat, h_lon, site_df) {
    if (is.na(h_lat) || is.na(h_lon)) return(NA_character_)
    dist <- sqrt((site_df$lat - h_lat)^2 + (site_df$lon - h_lon)^2)
    site_df$cell_id[which.min(dist)]
  }

  df_corr <- host_sites_clean %>%
    rowwise() %>%
    mutate(
      EPI_val = find_nearest_epi(lat_deg, lon_deg, site_rankings),
      cell = find_nearest_cell(lat_deg, lon_deg, site_rankings),
      log10_area = log10(area_ha)
    ) %>%
    ungroup() %>%
    filter(!is.na(log10_area) & !is.na(EPI_val))

  r_loc <- cor.test(df_corr$log10_area, df_corr$EPI_val, method = "pearson")
  rho_loc <- cor.test(df_corr$log10_area, df_corr$EPI_val, method = "spearman")

  cell_agg <- df_corr %>%
    group_by(cell) %>%
    summarise(mean_log_area = mean(log10_area), EPI_val = first(EPI_val), .groups = "drop")
  r_cell <- cor.test(cell_agg$mean_log_area, cell_agg$EPI_val, method = "pearson")
  rho_cell <- cor.test(cell_agg$mean_log_area, cell_agg$EPI_val, method = "spearman")

  n_loc <- nrow(df_corr)
  n_cell <- nrow(cell_agg)
}

cat(sprintf("Locality correlation (n = %d): Pearson r = %.3f (P = %.4f), Spearman rho = %.3f (P = %.4f)\n",
            n_loc, r_loc$estimate, r_loc$p.value, rho_loc$estimate, rho_loc$p.value))
cat(sprintf("Grid cell correlation (n = %d): Pearson r = %.3f (P = %.4f), Spearman rho = %.3f (P = %.4f)\n\n",
            n_cell, r_cell$estimate, r_cell$p.value, rho_cell$estimate, rho_cell$p.value))

p_format <- function(p) if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)

ann_text <- paste0(
  "Locality (n = ", n_loc, "):\n",
  "  Pearson r = ", sprintf("%.3f", r_loc$estimate), " (P ", p_format(r_loc$p.value), ")\n",
  "  Spearman ρ = ", sprintf("%.3f", rho_loc$estimate), " (P ", p_format(rho_loc$p.value), ")\n\n",
  "Grid cell (n = ", n_cell, "):\n",
  "  Pearson r = ", sprintf("%.3f", r_cell$estimate), " (P ", p_format(r_cell$p.value), ")\n",
  "  Spearman ρ = ", sprintf("%.3f", rho_cell$estimate), " (P ", p_format(rho_cell$p.value), ")"
)

fig4b <- ggplot(df_corr, aes(x = EPI_val, y = log10_area)) +
  geom_point(aes(color = uf), size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate(
    "text", x = 118, y = 2.85, label = ann_text,
    hjust = 1, vjust = 1, size = 3.2, lineheight = 1.05,
    fontface = "plain", color = "gray20"
  ) +
  scale_color_brewer(palette = "Set1", name = "State") +
  scale_x_continuous(limits = c(0, 125)) +
  labs(
    x = expression("Nearest grid cell annual " * EPI[inf]),
    y = expression("Log"[10] * " infested area (ha)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Composite Figure 4: Side-by-Side (Panel A | Panel B)
fig4_combined <- (fig4a | fig4b) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Fig4_spatial_maps.png"), fig4_combined, width = 11.5, height = 5.5, dpi = 300)

# Cairo PDF for Greek rho rendering
tryCatch({
  ggsave(file.path(fig_dir, "Fig4_spatial_maps.pdf"), fig4_combined, width = 11.5, height = 5.5, device = cairo_pdf)
}, error = function(e) {
  ggsave(file.path(fig_dir, "Fig4_spatial_maps.pdf"), fig4_combined, width = 11.5, height = 5.5)
})

cat("  -> Saved Fig4_spatial_maps.png and .pdf\n\n")

cat("Step 10 completed successfully.\n")
cat("==============================================================================\n")
