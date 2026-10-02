# Reproducibility Summary

**Study**: Epidemic Potential of *Maravalia cryptostegiae*, a Classical Biological Control Agent of Rubber Vine, in Northeastern Brazil  
**Status**: All analytical results, tables, and figures fully reproduced and verified.

---

## 1. Pipeline Verification Overview

The modeling pipeline executes end-to-end via `run_all.R`. Execution from scratch runs all 11 stages and sensitivity analyses:

- **Master Driver**: `Rscript run_all.R --sensitivity`
- **Verification Suite**: `Rscript scripts/verify_outputs.R`
- **Verification Outcome**: **27 / 27 tables PASSED**; **7 / 7 scientific benchmarks PASSED**.

### Key Scientific Benchmarks Verified
| Metric | Benchmark Target | Obtained Value | Status |
| :--- | :---: | :---: | :---: |
| Brazil Mean Annual EPI_inf | 62.83 | 62.83 | PASS |
| Parametric Bootstrap Delta RSS | 0.2601 | 0.2601 | PASS |
| Bootstrap Empirical P-value | 0.0020 | 0.0020 | PASS |
| Australia McLeod River Mean EPI | 161.00 | 161.00 | PASS |
| Australia Inkerman Mean EPI | 5.74 | 5.74 | PASS |
| Australia Delta Downs Mean EPI | 0.72 | 0.72 | PASS |
| Australia Welch t-test P-value | 0.140 | 0.140 | PASS |

---

## 2. Table Verification Summary

| Suite / Scenario | Target Directory | Evaluated Tables | Status |
| :--- | :--- | :---: | :---: |
| **Default Analysis** | `outputs/tables/` vs. `reference_outputs/tables_v2/` | 11 | **11 / 11 PASS** |
| **Sensitivity Run (a)** ($P \ge 1.0\text{ mm/h}$) | `outputs/tables_v2_precip1/` | 8 | **8 / 8 PASS** |
| **Sensitivity Run (b)** ($W \le 24\text{ h}$) | `outputs/tables_v2_maxh24/` | 8 | **8 / 8 PASS** |
| **Total** | | **27** | **27 / 27 PASS** |

---

## 3. Quick Start

```bash
# 1. Verify or retrieve meteorological data
Rscript scripts/get_data.R

# 2. Run full modeling pipeline
Rscript run_all.R --sensitivity

# 3. Verify numerical identity against reference benchmarks
Rscript scripts/verify_outputs.R
```
