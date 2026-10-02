# =============================================================================
# config/settings.R
# Global Model Parameters and Epidemiological Thresholds
# Manuscript: Epidemic Potential of Maravalia cryptostegiae in Northeastern Brazil
# =============================================================================

# ---- Repository Root Guard --------------------------------------------------
if (!file.exists("config/settings.R") || !file.exists("run_all.R")) {
  stop("Execution halted: Scripts must be executed from the repository root directory containing 'config/settings.R' and 'run_all.R'. Current working directory: ", getwd(), call. = FALSE)
}

# ---- Wetness Parameters -----------------------------------------------------
if (!exists("RH_THRESH"))      RH_THRESH      <- 90    # Hour is wet if RH2M >= RH_THRESH (%)
if (!exists("PRECIP_THRESH"))  PRECIP_THRESH  <- NA    # Precipitation threshold (mm/h). NA = disabled (pure RH threshold).
                                                       # In sensitivity runs, set to 1.0 mm/h.
if (!exists("W_CAP"))          W_CAP          <- 24    # Thermal breadth sigma(W) capped at W = 24 h (experimental range)
if (!exists("W_MIN"))          W_MIN          <- 6     # Minimum wetness episode duration (h); episodes < W_MIN yield I = 0
if (!exists("MAX_EPISODE_H"))  MAX_EPISODE_H  <- NA    # Maximum episode duration before splitting (h). NA = disabled.
                                                       # In sensitivity runs, set to 24 h to split into 24-h blocks.
if (!exists("LAT_T_PLATEAU"))  LAT_T_PLATEAU  <- 27    # Latent development rate held constant at r(27 C) above 27 C
if (!exists("W_BAR"))          W_BAR          <- 12.5  # Centering constant for wetness duration in Model 3 (h)
if (!exists("N_YEARS"))        N_YEARS        <- 25    # Total years in weather record (2001-2025)

# ---- Stochastic Seeds -------------------------------------------------------
if (!exists("SEED_OPTIM_CV"))   SEED_OPTIM_CV   <- 123  # Seed for LOWO/LOTO multi-start cross-validation
if (!exists("SEED_BOOTSTRAP"))  SEED_BOOTSTRAP  <- 2026 # Seed for parametric bootstrap of dynamic breadth parameter gamma
if (!exists("SEED_H1_FIT"))     SEED_H1_FIT     <- 1    # Seed for H1 multistart optimization
if (!exists("SEED_SENSITIVITY")) SEED_SENSITIVITY <- 123 # Seed for digitization sensitivity perturbations

# ---- Bootstrap Configurations -----------------------------------------------
if (!exists("BOOTSTRAP_B"))     BOOTSTRAP_B     <- 2000 # Number of parametric bootstrap replicates for gamma test
if (!exists("SENSITIVITY_B"))   SENSITIVITY_B   <- 500  # Number of stochastic perturbation replicates for digitization test

# ---- Directory Conventions --------------------------------------------------
if (!exists("DIR_DATA_RAW"))      DIR_DATA_RAW      <- "data/raw"
if (!exists("DIR_DATA_DERIVED"))  DIR_DATA_DERIVED  <- "data/derived"
if (!exists("DIR_DATA_METADATA")) DIR_DATA_METADATA <- "data/metadata"
if (!exists("DIR_OUTPUT_TABLES")) DIR_OUTPUT_TABLES <- "outputs/tables"
if (!exists("DIR_OUTPUT_FIGS"))   DIR_OUTPUT_FIGS   <- "outputs/figures"
if (!exists("DIR_OUTPUT_MODELS")) DIR_OUTPUT_MODELS <- "outputs/models"
if (!exists("DIR_OUTPUT_LOGS"))   DIR_OUTPUT_LOGS   <- "outputs/logs"
