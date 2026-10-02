# =============================================================================
# R/maravalia_functions.R
# Shared Epidemiological Functions & Submodels for Brazil and Australia Pipelines
# Sourced by all analysis scripts. Thresholds are defined in config/settings.R.
# =============================================================================

# Locate and source global configuration
find_proj_root <- function() {
  if (exists("PROJ_ROOT", envir = .GlobalEnv) && !is.null(get("PROJ_ROOT", envir = .GlobalEnv))) {
    return(get("PROJ_ROOT", envir = .GlobalEnv))
  }
  if (file.exists("config/settings.R")) return(".")
  if (file.exists("../config/settings.R")) return("..")
  "."
}

PROJ_ROOT <- find_proj_root()
settings_path <- file.path(PROJ_ROOT, "config/settings.R")
if (!file.exists(settings_path)) {
  # Fallback if working directory is inside R/
  if (file.exists("../config/settings.R")) {
    settings_path <- "../config/settings.R"
    PROJ_ROOT <- ".."
  }
}

if (file.exists(settings_path)) {
  source(settings_path)
} else {
  stop("config/settings.R not found. Scripts must be run from the repository root directory.", call. = FALSE)
}

# ---- Wetness Classification -------------------------------------------------
wet_hours <- function(df, rh = RH_THRESH, p = PRECIP_THRESH) {
  w <- df$RH2M >= rh
  if (!is.na(p)) w <- w | (df$PRECTOTCORR >= p)
  w[is.na(w)] <- FALSE
  w
}

# Returns start/end indices and duration of continuous wet runs (optionally split)
wet_runs <- function(is_wet, max_h = MAX_EPISODE_H) {
  r   <- rle(is_wet)
  end <- cumsum(r$lengths); st <- end - r$lengths + 1
  st  <- st[r$values]; end <- end[r$values]
  if (!is.na(max_h) && length(st) > 0) {
    s2 <- integer(0); e2 <- integer(0)
    for (i in seq_along(st)) {
      b <- seq(st[i], end[i], by = max_h)
      s2 <- c(s2, b); e2 <- c(e2, pmin(b + max_h - 1, end[i]))
    }
    st <- s2; end <- e2
  }
  data.frame(start_idx = st, end_idx = end, wet_h = end - st + 1)
}

# ---- Infection Submodel (Model 3, sigma capped at W_CAP) ---------------------
calc_I <- function(T_wet, W, coefs, w_cap = W_CAP, w_min = W_MIN) {
  Ws  <- pmin(W, w_cap)
  sig <- pmax(1.0, coefs[["sigma_Wbar"]] + coefs[["gamma"]] * (Ws - W_BAR))
  I   <- exp(-0.5 * ((T_wet - coefs[["Topt"]]) / sig)^2) *
         (1 - exp(-(W / coefs[["lambda"]])^coefs[["k"]]))
  ifelse(W < w_min, 0, pmin(1, pmax(0, I)))
}

# ---- Latent Development Rate (per hour), plateau above LAT_T_PLATEAU --------
A_LAT <- -0.07814868; B_LAT <- 0.00718076
calc_rate_hourly <- function(T_vals, t_plateau = LAT_T_PLATEAU) {
  pmax(0, (A_LAT + B_LAT * pmin(T_vals, t_plateau)) / 24)
}
