# Epidemiological Methods & Mathematical Formulations

This document provides a concise summary of the mathematical models, parameterizations, microclimatic wetness definitions, information criteria conventions, and pipeline corrections implemented in this study.

---

## 1. Infection Potential Submodels

The infection potential submodel predicts relative infection response $I(T, W) \in [0, 1]$ as a function of mean incubation temperature $T$ (°C) and surface wetness duration $W$ (hours), normalized to the maximum observed pustule density ($14.35\text{ pustules cm}^{-2}$).

### Candidate Submodels

#### Model 1: Separable Generalized Beta $\times$ Weibull (5 parameters)
$$I_1(T, W) = a (T - T_{\min})^d (T_{\max} - T)^e \left[1 - \exp\left(-\left(\frac{W}{\lambda}\right)^k\right)\right]$$
where $T_{\min} = 10\text{ °C}$ and $T_{\max} = 35\text{ °C}$ are biological cardinal temperature bounds, and $\lambda, k$ are Weibull wetness threshold and shape parameters.

#### Model 2: Wetness-Expanded Analytis $\times$ Weibull (6 parameters)
$$I_2(T, W) = \left[a \left(\frac{T - 10}{25}\right)^b \left(1 - \frac{T - 10}{25}\right)\right]^{c_0 + c_1 W} \left[1 - \exp\left(-\left(\frac{W}{\lambda}\right)^k\right)\right]$$
where $c_0 + c_1 W$ expands thermal breadth as a function of wetness duration.

#### Model 3: Centered Gaussian Dynamic Breadth $\times$ Weibull (**Primary Model**, 5 parameters)
$$I_3(T, W) = \exp\left(-0.5 \left(\frac{T - T_{\text{opt}}}{\sigma(W)}\right)^2\right) \left[1 - \exp\left(-\left(\frac{W}{\lambda}\right)^k\right)\right]$$
with dynamic thermal breadth:
$$\sigma(W) = \max\left(1.0, \, \sigma_{\bar{W}} + \gamma (W_{\text{capped}} - \bar{W})\right)$$
where:
- $\bar{W} = 12.5\text{ h}$ is the empirical centering constant (mean of experimental durations: 6, 8, 12, 24 h).
- $W_{\text{capped}} = \min(W, 24\text{ h})$ prevents non-biological extrapolation beyond experimental limits.
- $\sigma_{\bar{W}}$ is the thermal breadth at $\bar{W} = 12.5\text{ h}$.
- $\gamma$ is the dynamic breadth expansion rate (°C h⁻¹).
- If $W < 6\text{ h}$, $I_3(T, W) = 0$ (biological minimum wetness threshold).

#### Supplementary Models S1–S3
- **Model S1 (Simple Analytis $\times$ Weibull, 5 parameters)**: Static power exponent $c$ replacing $c_0 + c_1 W$.
- **Model S2 (Duthie Hyperbolic Secant $\times$ Weibull, 4 parameters)**:
  $$I_{\text{S2}}(T, W) = \frac{1 - \exp(-(B \cdot W)^D)}{\cosh((T - F) G / 2)}$$
- **Model S3 (Fixed-Breadth Gaussian $\times$ Weibull, Null $H_0$, 4 parameters)**:
  Static $\sigma(W) = \sigma_{\bar{W}}$ ($\gamma = 0$).

---

## 2. Information Criteria & Goodness-of-Fit Conventions

For nonlinear least-squares regression with $N$ observations and $p$ estimated parameters, the residual variance $\sigma^2$ is an additional estimated parameter ($K = p + 1$). The log-likelihood under Gaussian errors is:
$$\log L = -\frac{N}{2} \ln(2\pi) - \frac{N}{2} \ln\left(\frac{\text{RSS}}{N}\right) - \frac{N}{2}$$

Information criteria are defined consistently with standard likelihood principles:
$$\text{AIC} = -2 \log L + 2(p + 1)$$
$$\text{AICc} = \text{AIC} + \frac{2(p + 1)(p + 2)}{N - (p + 1) - 1} = \text{AIC} + \frac{2(p + 1)(p + 2)}{N - p - 2}$$

---

## 3. Microclimatic Wetness Definition

In NASA POWER hourly reanalysis:
- **Wet hour definition**: An hour is classified as wet if and only if relative humidity at 2 m meets or exceeds the critical threshold:
  $$\text{Wet Hour} \iff \text{RH2M} \ge 90\%$$
- **Precipitation clause**: The precipitation clause (`PRECTOTCORR > 0`) is **disabled** in the primary model (`PRECIP_THRESH = NA`) because trace rain occurs in 50–74% of hours in Northeastern Brazil, which unrealistically merged distinct dew events into continuous multi-day wetness episodes.
- **Episode construction**: Continuous runs of wet hours form a wetness episode of duration $W$ hours and mean temperature $T_{\text{wet}} = \frac{1}{W} \sum T_h$.

---

## 4. Latent Period Submodel & Cohort Tracking

The latent development rate $r(T)$ ($\text{day}^{-1}$) was parameterized from Evans & Fleureau (1993) data across 18–27 °C:
$$r(T) = a + b T = -0.07814868 + 0.00718076 \cdot T \quad (R^2 = 0.942, \, P < 0.0001)$$
- **Hourly rate**: $r_{\text{hourly}}(T) = \max\left(0, \frac{r(T)}{24}\right)$
- **Thermal plateau**: Above $T = 27\text{ °C}$ (the experimental limit), development rate is held constant:
  $$r_{\text{hourly}}(T) = \max\left(0, \frac{a + b \min(T, 27)}{24}\right)$$
- **Cohort development**: For each infection episode initiating at hour $t_0$, accumulated physiological development accumulates hourly:
  $$D(t) = \sum_{h=t_0}^t r_{\text{hourly}}(T_h)$$
  Sporulation completion occurs at the earliest hour where $D(t) \ge 1.0$.

---

## 5. Epidemic Pressure Indices (EPI)

- **Infection Pressure ($EPI_{\text{inf}}$)**:
  Cumulative potential infection opportunities initiated within a given time interval:
  $$EPI_{\text{inf}} = \sum_{i \in \text{episodes}} I_3(T_i, W_i)$$
- **Annual Normalization**:
  Annual site means are computed strictly over the 25-year record ($N_{\text{years}} = 25$):
  $$\overline{EPI}_{\text{inf}} = \frac{1}{25} \sum_{y=2001}^{2025} EPI_{\text{inf}, y}$$
  ensuring that dry years without infection events correctly count as zero.

---

## 6. Summary of Corrections in Final Version (v2 / Final)

1. **Wetness Definition**: Removed `PRECTOTCORR > 0` clause, restoring mean episode duration to 9.4 h (from 31.2 h).
2. **Model 2 Boundary**: Expanded $c_0$ upper bound from 20 to 200, finding an interior unconstrained optimum ($c_0 = 26.20$).
3. **Parametric Bootstrap**: Initialized $H_1$ optimizations directly from the nested $H_0$ solution, guaranteeing $\Delta\text{RSS} \ge 0$ ($P = 0.0020$).
4. **Zero-Event Averaging**: Normalised annual means by $N_{\text{years}} = 25$ rather than non-zero years, correcting dry-site inflation in Australia (Delta Downs: 0.72 vs 0.82).
5. **Figure De-cluttering (v3 figures)**: Removed all theoretical sporulation pressure lines ($EPI_{\text{spore}}$) from manuscript figures, focusing exclusively on realized weather-driven infection potential ($EPI_{\text{inf}}$).
