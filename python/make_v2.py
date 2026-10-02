import re,io,os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
def rd(f): return open('R/'+f,encoding='utf-8').read()
def wr(f,s): open('R/'+f,'w',encoding='utf-8').write(s); print('wrote R/'+f)
def rep(s,a,b,n=1,regex=False):
    c = len(re.findall(a,s,flags=re.S)) if regex else s.count(a)
    assert c>=n, f'pattern not found: {a[:60]}'
    return re.sub(a,b,s,flags=re.S) if regex else s.replace(a,b)
HDR='# >>> v2 revision (2026-09-28): see R/maravalia_functions_v2.R and review_2026-09-28/CHANGES_v2.md\nsource("R/maravalia_functions_v2.R")\nfor (d in c("data_v2","tables_v2","figures_v2","models_v2")) if (!dir.exists(d)) dir.create(d)\n'
def paths(s):
    for a,b in [('"data/phase2_','"data_v2/phase2_'),('"tables/rh_','"tables_v2/rh_'),('"tables/site_epidemic','"tables_v2/site_epidemic'),
                ('"tables/monthly_clim','"tables_v2/monthly_clim'),('"tables/australian_model_validation','"tables_v2/australian_model_validation'),
                ('"tables/bootstrap','"tables_v2/bootstrap'),('"data/bootstrap_gamma_samples','"data_v2/bootstrap_gamma_samples'),
                ('"report/figures/','"figures_v2/'),('"figures/','"figures_v2/'),('Saved figures/','Saved figures_v2/'),
                ('"tables/infection_model','"tables_v2/infection_model'),('"tables/all_infection','"tables_v2/all_infection'),('"models/model_infection_','"models_v2/model_infection_')]:
        s=s.replace(a,b)
    return s
first_code=lambda s: s.index('cat("====')
ins=lambda s,t: s[:first_code(s)]+t+'\n'+s[first_code(s):]

# ---- 07: wetness episodes ----
s=paths(rd('07_wetness_episodes.R'))
s=s.replace('"models_v2/model_infection_primary.rds"','"models/model_infection_primary.rds"')
over='''
# >>> v2 overrides: shared wetness rule (no P>0 clause), sigma capped at W_CAP
calc_infection_potential <- function(T_wet, W) calc_I(T_wet, W, coefs)
extract_episodes_fast <- function(df, rh_thresh = RH_THRESH, include_precip = TRUE) {
  is_wet <- wet_hours(df, rh = rh_thresh, p = if (include_precip) PRECIP_THRESH else NA)
  if (!any(is_wet)) return(NULL)
  ep <- wet_runs(is_wet)
  cs <- c(0, cumsum(df$T2M))
  mean_T <- (cs[ep$end_idx + 1] - cs[ep$start_idx]) / ep$wet_h
  tibble(cell_id = df$cell_id[1], lat = df$LAT[1], lon = df$LON[1],
         start_year = df$YEAR[ep$start_idx], start_month = df$MO[ep$start_idx],
         start_day = df$DY[ep$start_idx], start_hour = df$HR[ep$start_idx],
         end_hour = df$HR[ep$end_idx], start_idx = ep$start_idx, end_idx = ep$end_idx,
         wet_h = ep$wet_h, T_wet = mean_T, I_pot = calc_infection_potential(mean_T, ep$wet_h))
}
'''
s=rep(s,'cat("[1/2] Extracting',over+'\ncat("[1/2] Extracting')
s=ins(s,HDR); wr('07_wetness_episodes_v2.R',s)

# ---- 08: latent tracking ----
s=paths(rd('08_latent_tracking.R'))
s=rep(s,'track_site_cohorts <- function','# >>> v2 override: latent rate held constant above LAT_T_PLATEAU\ncalc_rate_primary <- function(T_vals) calc_rate_hourly(T_vals)\n\ntrack_site_cohorts <- function')
s=ins(s,HDR); wr('08_latent_tracking_v2.R',s)

# ---- 09, 10: paths only ----
for f in ['09_epidemic_pressure_indices.R','10_seasonal_spatial_analysis.R']:
    s=paths(rd(f)); s=s.replace('read_csv("tables/site_epidemic','read_csv("tables_v2/site_epidemic'); s=ins(s,HDR); wr(f.replace('.R','_v2.R'),s)

# ---- 12: Australia ----
s=paths(rd('12_australian_validation_pipeline.R'))
s=s.replace('"models_v2/model_infection_primary.rds"','"models/model_infection_primary.rds"')
s=rep(s,'is_wet <- df$RH2M >= 90','is_wet <- wet_hours(df)   # v2: same rule as Brazil')
s=rep(s,'# 2. Extract Wetness Episodes','# >>> v2 overrides: identical infection/latent functions to Brazil\ncalc_I_pot <- function(T_val, W_val) calc_I(T_val, W_val, coefs)\ncalc_rate  <- function(T_val) calc_rate_hourly(T_val) * 24   # original used /24 inside\n\n# 2. Extract Wetness Episodes')
s=rep(s,'Dry Inland Sites:','Other (non-high-impact) sites:')
s=rep(s,'t_res <- t.test(high_impact_EPI, mod_impact_EPI)','t_res <- t.test(high_impact_EPI, mod_impact_EPI)\nw_res <- wilcox.test(high_impact_EPI, mod_impact_EPI, exact = TRUE)')
s=rep(s,'cat("Welch Two-Sample t-test p-value:", format.pval(t_res$p.value, digits = 4), "\\n\\n")',
 'cat(sprintf("Welch t = %.2f, df = %.1f, P = %.3f\\n", t_res$statistic, t_res$parameter, t_res$p.value))\ncat(sprintf("Mann-Whitney W = %.0f, P = %.3f (n = %d vs %d)\\n\\n", w_res$statistic, w_res$p.value, length(high_impact_EPI), length(mod_impact_EPI)))')
s=ins(s,HDR); wr('12_australian_validation_pipeline_v2.R',s)

# ---- 04: bootstrap ----
s=paths(rd('04_bootstrap_gamma.R'))
s=rep(s,'lower = c(Topt = 15, sigma_Wbar = 0.1, lambda = 0.1, k = 0.1)','lower = c(Topt = 15, sigma_Wbar = 0.1, lambda = 1.0, k = 1.0)',n=2)
helper='''
# >>> v2: H1 refit started from the H0 solution (gamma = 0) plus fixed/random starts,
#     keeping the lowest RSS. Guarantees RSS_H1 <= RSS_H0 (nested), so Delta RSS >= 0.
fit_h1_best <- function(df, h0) {
  starts <- list(st_h1)
  if (!inherits(h0, "try-error")) {
    c0 <- coef(h0)
    for (g in c(0, 0.1, 0.3)) starts[[length(starts) + 1]] <- list(Topt = c0[["Topt"]],
      sigma_Wbar = c0[["sigma_Wbar"]], gamma = g, lambda = max(1, c0[["lambda"]]), k = max(1, c0[["k"]]))
  }
  for (j in 1:4) starts[[length(starts) + 1]] <- lapply(st_h1, function(v) v * exp(runif(1, -0.4, 0.4)))
  best <- NULL
  for (st in starts) {
    f <- try(nlsLM(m_h1_formula, data = df, start = st,
                   lower = c(Topt = 15, sigma_Wbar = 0.1, gamma = 0, lambda = 1.0, k = 1.0),
                   upper = c(Topt = 35, sigma_Wbar = 10, gamma = 5, lambda = 50, k = 20),
                   control = nls.lm.control(maxiter = 300)), silent = TRUE)
    if (!inherits(f, "try-error") && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  if (is.null(best)) structure("fail", class = "try-error") else best
}
fit_h1 <- fit_h1_best(inf_data, fit_h0)
'''
s=rep(s,'rss_h0_obs <- deviance(fit_h0)',helper+'\nrss_h0_obs <- deviance(fit_h0)')
s=rep(s,r'fit_h1_b <- try\(nlsLM\(m_h1_formula, data = df_sim.*?silent = TRUE\)','fit_h1_b <- fit_h1_best(df_sim, fit_h0_b)',regex=True)
s=rep(s,'p_value <- mean(valid_delta >= delta_rss_obs)','p_value <- (sum(valid_delta >= delta_rss_obs) + 1) / (length(valid_delta) + 1)   # v2')
s=rep(s,'(P < 0.001)','(P = %.4f)'); s=rep(s,'label = sprintf("Observed ΔRSS\\n= %.4f\\n(P = %.4f)", delta_rss_obs)','label = sprintf("Observed ΔRSS\\n= %.4f\\n(P = %.4f)", delta_rss_obs, p_value)')
s=ins(s,HDR); wr('04_bootstrap_gamma_v2.R',s)

# ---- 02: widen Model 2 c0 bound; one AIC convention ----
s=paths(rd('02_infection_submodel.R'))
s=rep(s,'aic_val  <- N * log(rss_val / N) + 2 * p_val','aic_val  <- -2 * loglik + 2 * (p_val + 1)   # v2: same likelihood as LogLik; +1 for residual variance')
s=rep(s,'aicc_val <- aic_val + (2 * p_val * (p_val + 1)) / (N - p_val - 1)','aicc_val <- aic_val + (2 * (p_val + 1) * (p_val + 2)) / (N - p_val - 2)')
s=rep(s,'upper = c(a = 100, b = 20, c0 = 20, c1 = 1.0, lambda = 50, k = 20)','upper = c(a = 100, b = 20, c0 = 200, c1 = 5.0, lambda = 50, k = 20)')
s=rep(s,'c0 = c(1.0, 5.0, 10.0),','c0 = c(1.0, 5.0, 10.0, 25, 50),')
s=ins(s,HDR); wr('02_infection_submodel_v2.R',s)

open('R/run_phase2_v2.R','w').write('''# Runs the v2 (review) pipeline. Outputs go to data_v2/, tables_v2/, figures_v2/, models_v2/.
# Originals in data/, tables/, figures/, models/ are not touched.
# Run from the project root (setwd to the folder containing R/).
t0 <- Sys.time()
source("R/02_infection_submodel_v2.R")   # Model 2 bound + AIC convention (Model 3 unchanged)
source("R/04_bootstrap_gamma_v2.R")      # corrected H1 refits; ~2000 x 8 fits, may take a while
source("R/07_wetness_episodes_v2.R")
source("R/08_latent_tracking_v2.R")
source("R/09_epidemic_pressure_indices_v2.R")
source("R/10_seasonal_spatial_analysis_v2.R")
source("R/12_australian_validation_pipeline_v2.R")
cat("v2 pipeline finished in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\\n")
''')
print('wrote R/run_phase2_v2.R')
