# ==============================================================================
# Makefile: Automated Pipeline for Maravalia cryptostegiae Epidemic Potential
# ==============================================================================

.PHONY: all data analysis figures verify quick sensitivity clean help

help:
	@echo "Available targets:"
	@echo "  make all         - Run full workflow: data -> analysis -> figures -> verify"
	@echo "  make data        - Check and retrieve raw meteorological data"
	@echo "  make analysis    - Execute full epidemiological pipeline (B=2000)"
	@echo "  make figures     - Regenerate publication figures (Fig 1-5, S1-S4)"
	@echo "  make verify      - Numerically verify regenerated tables against reference outputs"
	@echo "  make quick       - Fast test run with B=100 bootstrap replicates"
	@echo "  make sensitivity - Run sensitivity scenarios (PRECIP_THRESH=1.0, MAX_EPISODE_H=24)"
	@echo "  make clean       - Remove generated tables, figures, models, and logs"

all: data analysis verify

data:
	@echo ">>> [Target: data] Verifying raw NASA POWER datasets..."
	Rscript scripts/get_data.R

analysis:
	@echo ">>> [Target: analysis] Running full epidemiological modeling pipeline..."
	Rscript run_all.R

figures:
	@echo ">>> [Target: figures] Generating publication figures..."
	Rscript R/02_infection_submodel.R
	Rscript R/04_bootstrap_gamma.R
	Rscript R/05_latent_submodel.R
	Rscript R/06_digitization_sensitivity.R
	Rscript R/09_epidemic_pressure_indices.R
	Rscript R/10_seasonal_spatial_analysis.R
	Rscript R/12_australian_validation_pipeline.R

verify:
	@echo ">>> [Target: verify] Verifying outputs against reference benchmarks..."
	Rscript scripts/verify_outputs.R

quick:
	@echo ">>> [Target: quick] Running quick verification pipeline (B=100)..."
	Rscript run_all.R --quick
	Rscript scripts/verify_outputs.R

sensitivity:
	@echo ">>> [Target: sensitivity] Running sensitivity pipeline..."
	Rscript run_all.R --sensitivity

clean:
	@echo ">>> [Target: clean] Cleaning output directories..."
	rm -rf outputs/tables/* outputs/figures/* outputs/models/* outputs/logs/*
