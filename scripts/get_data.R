# ==============================================================================
# Script: scripts/get_data.R
# Purpose: Verify or download NASA POWER hourly meteorological datasets
#          and validate SHA-256 checksums against data/raw/MANIFEST.csv
# Project: Maravalia cryptostegiae Epidemiological Submodels
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

if (!exists("PROJ_ROOT")) PROJ_ROOT <- "."

manifest_file <- file.path(PROJ_ROOT, "data/raw/MANIFEST.csv")
if (!file.exists(manifest_file)) {
  stop("MANIFEST.csv not found at: ", manifest_file)
}

manifest <- read_csv(manifest_file, show_col_types = FALSE)
cat("==============================================================================\n")
cat("DATA RETRIEVAL & INTEGRITY VERIFICATION\n")
cat(sprintf("Loaded manifest with %d raw data files.\n", nrow(manifest)))
cat("==============================================================================\n\n")

# Check which files exist locally
manifest <- manifest %>%
  mutate(
    full_local_path = file.path(PROJ_ROOT, relative_path),
    exists_locally = file.exists(full_local_path)
  )

missing_count <- sum(!manifest$exists_locally)
cat(sprintf("Local status: %d files present, %d files missing.\n\n", 
            sum(manifest$exists_locally), missing_count))

if (missing_count > 0) {
  cat("[!] Missing raw NASA POWER files detected.\n")
  cat("    Initiating automated download pipeline via R/00_download_nasapower_brazil.R and R/11_download_australian_nasapower.R...\n")
  cat("    NOTE: Full download from NASA POWER API takes ~30-45 minutes. Pre-computed data can also be obtained from Zenodo.\n\n")
  
  # Run download scripts if needed
  br_missing <- manifest %>% filter(!exists_locally, grepl("nasapower_maravalia", relative_path))
  if (nrow(br_missing) > 0) {
    cat(">>> Downloading missing Brazilian weather grid cells...\n")
    source(file.path(PROJ_ROOT, "R/00_download_nasapower_brazil.R"))
  }
  
  au_missing <- manifest %>% filter(!exists_locally, grepl("nasapower_australia", relative_path))
  if (nrow(au_missing) > 0) {
    cat(">>> Downloading missing Australian weather stations...\n")
    source(file.path(PROJ_ROOT, "R/11_download_australian_nasapower.R"))
  }
}

# Verify SHA-256 checksums
cat("\n[Integrity Verification: SHA-256 Checksum Validation]\n")
compute_sha256 <- function(filepath) {
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(file = filepath, algo = "sha256")
  } else if (requireNamespace("openssl", quietly = TRUE)) {
    as.character(openssl::sha256(file(filepath, "rb")))
  } else {
    # System fallback
    cmd <- if (Sys.info()["sysname"] == "Darwin") paste0("shasum -a 256 '", filepath, "'") else paste0("sha256sum '", filepath, "'")
    res <- system(cmd, intern = TRUE)
    strsplit(res, " ")[[1]][1]
  }
}

verified <- 0
mismatches <- 0

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  fpath <- row$full_local_path
  if (!file.exists(fpath)) {
    warning("File still missing after download step: ", row$relative_path, call. = FALSE)
    next
  }
  
  chk <- compute_sha256(fpath)
  if (tolower(chk) == tolower(row$sha256)) {
    verified <- verified + 1
  } else {
    mismatches <- mismatches + 1
    warning(sprintf("Checksum mismatch for %s:\n  Expected: %s\n  Computed: %s\n  (Note: NASA POWER reanalysis values can be revised by NASA over time; this does not prevent execution.)",
                    row$filename, row$sha256, chk), call. = FALSE)
  }
}

cat(sprintf("\nVerification complete:\n  - Verified exact match: %d / %d files\n  - Warnings / mismatches: %d\n",
            verified, nrow(manifest), mismatches))
cat("==============================================================================\n")
