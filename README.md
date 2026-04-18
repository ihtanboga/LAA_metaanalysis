# LAAC vs OAC: Domain-Based Meta-Analysis of Clinical Trade-Offs

A domain-based meta-analysis with harmonised endpoint reconstruction evaluating **percutaneous left atrial appendage closure (LAAC) versus oral anticoagulation (OAC)** in atrial fibrillation. All pooled estimates are Mantel–Haenszel random-effects risk ratios with REML τ² and Hartung–Knapp confidence intervals.

## Key Findings

| Endpoint | Contemporary Pool (k=3, n=5,002) | Expanded Pool (k=4, n=5,890) | All Trials (k=5, n=7,004) |
|---|---|---|---|
| Ischaemic stroke + systemic embolism | RR 1.29 (0.83–1.99); I²=0% | RR 1.22 (0.93–1.60); I²=0% | **RR 1.30 (1.01–1.67); P=0.040**; I²=0% |
| Non-procedural clinically relevant bleeding | **RR 0.51 (0.37–0.70)**; I²=0% | RR 0.60 (0.41–0.88); I²=55% | **RR 0.56 (0.40–0.76)**; I²=53% |
| All major bleeding (incl. procedural) | RR 0.87 (0.62–1.21); I²=0% | RR 0.95 (0.72–1.26); I²=0% | RR 0.93 (0.77–1.13); I²=0% |
| All-cause death | RR 0.91 (0.60–1.38); I²=0% | RR 1.00 (0.79–1.27); I²=1% | RR 0.92 (0.74–1.15); I²=38% |
| Hemorrhagic stroke | RR 0.90 (0.43–1.89) | — | RR 0.58 (0.24–1.41); I²=30% |
| Any stroke + SE | RR 1.19 (0.62–2.27) | RR 1.15 (0.81–1.63) | RR 1.09 (0.84–1.41); I²=0% |

**Historical Warfarin Benchmark (PROTECT-AF/PREVAIL, 5-year pooled, k=1, n=1,114):** IS+SE RR 1.68 (0.93–3.02); non-procedural bleeding RR 0.48 (0.32–0.71); hemorrhagic stroke **RR 0.20 (0.07–0.56)**.

**Headline.** In the contemporary DOAC-era pool, LAAC reduces non-procedural bleeding by ≈49% versus OAC without a statistically significant difference in thromboembolic protection; a modest IS+SE excess reaches significance only when warfarin-era evidence is pooled. The hemorrhagic stroke benefit is a warfarin-era artefact that does not persist in contemporary practice.

## Included Trials

| Trial | N | Comparator | Device | Follow-up |
|---|---|---|---|---|
| CHAMPION-AF | 3,000 | DOAC | Watchman FLX | 3.0 yr |
| OPTION | 1,600 | OAC (95% DOAC) | Watchman FLX | 3.0 yr |
| PRAGUE-17 | 402 | DOAC (96% apixaban) | Amulet (61%)/Watchman (36%) | 3.5 yr (median) |
| CLOSURE-AF | 888 | Best medical care (85% DOAC) | Watchman/FLX/Amulet | 3.0 yr (median) |
| PROTECT-AF / PREVAIL | 707 / 407 (1,114 pooled) | Warfarin | Watchman Gen 1 | 5 yr (pooled) |

**Total: 7,004 patients across 6 RCTs (5 analysis units — PROTECT-AF + PREVAIL analysed as a single pooled unit per Reddy et al. 2017, JACC).**

## Analysis Pools

- **Primary Contemporary (k=3; n=5,002)** — CHAMPION-AF + OPTION + PRAGUE-17 (primary inferential frame, all DOAC-era comparators + contemporary device platforms)
- **Expanded Contemporary (k=4; n=5,890)** — + CLOSURE-AF (adds a high-risk, phenotype-expansion population on best medical care)
- **Historical Warfarin Benchmark (k=1; n=1,114)** — PROTECT-AF/PREVAIL pooled 5-year dataset
- **All Trials (k=5; n=7,004)** — all of the above combined, for endpoints where pooling is defensible

## Methodology

- **Domain-based endpoint harmonisation** — trial-specific composites deconstructed into clinically interpretable domains (primary efficacy: IS+SE; key safety: non-procedural clinically relevant bleeding; supportive safety: all major bleeding; key harm: procedural complications; supportive efficacy: all-cause / CV death, hemorrhagic stroke, any stroke + SE)
- **KM curve reconstruction** — Kaplan–Meier digitisation and IPD reconstruction via the Guyot algorithm (`reconstructKM`) for endpoints not directly reported
- **Compositional endpoint derivation** — Monte Carlo event thinning with temporal weighting (e.g., CLOSURE-AF Table S11 for peri-procedural vs chronic bleeding timing) where the published Kaplan–Meier did not separate a harmonized domain
- **Comparator-era stratification** — warfarin-era vs DOAC-era data analysed separately; DOAC-era Contemporary pool is the primary inferential frame
- **Statistical framework** — two-stage random-effects meta-analysis (Mantel–Haenszel RR, REML τ², Hartung–Knapp–Sidik–Jonkman CI); HR-based pooling as secondary; prediction intervals emphasised over I² for heterogeneity interpretation
- **Time-varying incidence modelling (Figure 4)** — Poisson GAM with arm-specific thin-plate splines on time and a trial-level random intercept on 1-stage pooled IPD (`mgcv`, REML)
- **Benefit–risk assessment** — common-denominator Kaul framework: absolute event differences per 1,000 patient-years and per 1,000 patients over trial-specific horizons, without severity weighting
- **Data provenance tiers (T1–T4)** — tracking the directness and reliability of each data point; primary inference restricted to endpoints with ≥3 Tier 1–2 trials

## Repository Structure

```
R/                          # Analysis scripts
  00_setup.R                # Package loading & global settings
  01_data.R                 # All trial-level data input
  02_primary_analysis.R     # Primary efficacy + key safety meta-analysis
  03_secondary_analysis.R   # Secondary endpoints, landmark, procedural
  04_sensitivity.R          # Sensitivity analyses & benefit-risk
  06_main_figures.R         # Figures 1–3, 5–6
  07_suppl_figures.R        # 6 supplementary figures
  08_tables.R               # Main + supplementary tables
  run_all.R                 # Master pipeline script

  # Figure 4: time-varying Poisson GAM incidence (added 2026-04)
  10_figure4_stroke_incidence.R          # Figure 4A — stroke
  11_figure4_bleeding_incidence.R        # Figure 4B — all bleeding + PRAGUE-sensitivity
  12_figure4_nonproc_bleeding_incidence.R # Figure 4C — non-procedural bleeding

  # 1-stage meta-analysis helpers (IPD forest by time band)
  forest_stroke.R              # Period-specific Cox HR forest plot
  forest_stroke_check.R        # Overall (time-band-free) HR sanity check

  # CLOSURE-AF-specific IPD reconstruction & endpoint cascade
  closure_af_km.R              # CLOSURE-AF primary composite KM reconstruction
  closure_af_bleeding_free.R   # Bleeding-free composite (Table S13 approximation)
  closure_af_death_free.R      # Model 1 (bleeding-free) + Model 2 (thromboembolic only)
  closure_af_sensitivity.R     # Tipping point, crossover, exclusion sensitivity

  # Per-trial endpoint IPD reconstruction + calibration (added 2026-04)
  reconstruction/
    reconstruct_all_endpoints.R   # Guyot reconstruction for 13 trial×endpoint curves
    champion_stroke_calibrate.R   # NAR calibration to paper event counts (CHAMPION)
    prague_stroke_calibrate.R     # Arm-swap correction + NAR calibration (PRAGUE-17)
    option_stroke_se.R            # Derive OPTION stroke+SE from death+stroke+SE composite
    closure_af_bleeding_only.R    # Derive pure major bleeding from CLOSURE-AF primary
    closure_af_stroke_se_only.R   # Derive stroke+SE from CLOSURE-AF primary composite

data/
  km_reconstructed/         # 46 digitised KM curve CSVs (all trials)
  closure_af_ipd/           # Reconstructed pseudo-IPD for CLOSURE-AF primary composite
  ipd_reconstructed/        # Per-trial × per-endpoint reconstructed IPD (added 2026-04)
                            # 17 (arm, time, status) files spanning CHAMPION-AF, OPTION,
                            # PRAGUE-17, CLOSURE-AF × stroke/bleeding/composite endpoints
```

## How to Run

```r
# Full pipeline (data -> analysis -> tables -> figures)
source("R/run_all.R")

# Or step by step:
source("R/01_data.R")             # Load all data
source("R/02_primary_analysis.R") # Primary efficacy & safety
source("R/03_secondary_analysis.R") # Secondary endpoints
source("R/04_sensitivity.R")      # Sensitivity & benefit-risk
source("R/08_tables.R")           # Generate tables
source("R/06_main_figures.R")     # Generate Figures 1–3, 5–6
source("R/07_suppl_figures.R")    # Generate supplementary figures

# Figure 4 — time-varying Poisson GAM (independent pipeline)
source("R/10_figure4_stroke_incidence.R")
source("R/11_figure4_bleeding_incidence.R")
source("R/12_figure4_nonproc_bleeding_incidence.R")
```

### Requirements

R ≥ 4.3 with packages: `meta`, `metafor`, `survival`, `survRM2`, `mgcv`, `tidyverse`, `patchwork`, `gt`, `flextable`, `officer`, `forestploter`, `ggtext`, `scales`, `RColorBrewer`, `gridExtra`, `checkmate`, `reconstructKM`, `jsonlite`

## Figure 4: Time-Varying Incidence (GAM)

Figure 4 (Panels A–C) was generated by fitting Poisson generalized additive models (GAM) to time-split reconstructed IPD, with separate arm-specific thin-plate splines on time, a trial-level random intercept (`s(trial, bs = "re")`), and a log person-time offset. The prediction excludes the random-effect contribution to give a "typical trial" population-average curve.

| Panel | Endpoint | Cohort(s) | Script |
|---|---|---|---|
| 4A | Stroke | Contemporary (k=3), Expanded (+ CLOSURE-AF) | `R/10_figure4_stroke_incidence.R` |
| 4B | All bleeding | Contemporary, Expanded + PRAGUE-sensitivity | `R/11_figure4_bleeding_incidence.R` |
| 4C | Non-procedural bleeding | Contemporary only (CLOSURE-AF lacks this endpoint) | `R/12_figure4_nonproc_bleeding_incidence.R` |

## Per-Trial Endpoint Reconstruction

`R/reconstruction/reconstruct_all_endpoints.R` reconstructs 13 trial × endpoint Kaplan–Meier curves in one pass from the JSON/CSV digitization files, producing harmonized `(arm, time, status)` IPD under `data/ipd_reconstructed/`. Two per-trial calibration scripts address known digitization limitations:

- `champion_stroke_calibrate.R` — CHAMPION-AF stroke curves are sparsely digitized (12–14 clicks over 3 years). Iterative linear interpolation + NAR-scaling binary search reproduces the paper event counts (device 50, AC 33).
- `prague_stroke_calibrate.R` — PRAGUE-17 stroke/TIA digitization had its arm labels swapped; this script applies the swap and calibrates to paper counts (AC 13, device 14).

Compositional endpoint derivations (where the published Kaplan–Meier did not separate the harmonized domain) use Monte Carlo event thinning with arm-specific removal targets and, where available, temporal weighting from supplementary tables (e.g., CLOSURE-AF Supplementary Table S11 for peri-procedural vs chronic bleeding timing). See `data/ipd_reconstructed/README_methods.md` for full methodological detail.

## KM Curve Digitisation

The KM-reconstructed datasets in `data/km_reconstructed/` were generated using our open-source KM digitisation tool:

**[KM Reconstruction Tool (Demo)](https://ihtanboga.github.io/KMR/)**

This web-based tool enables extraction of individual patient-level time-to-event data from published Kaplan-Meier survival curves using the Guyot algorithm.

## Benefit–Risk (Kaul Common-Denominator Framework)

Integrating non-procedural bleeding (benefit), hemorrhagic stroke (benefit), all-cause death (neutral), IS+SE (harm), and procedural complications (harm) on a per-1,000-patient-years scale:

| Pool | Benefit–Risk Ratio | Net benefit (events per 1,000 patient-years) |
|---|---|---|
| Contemporary DOAC-era (k=3; n=5,002) | **3.17** | **≈ +18.5** (favouring LAAC) |
| All Trials (k=5; n=7,004) | 1.53 | attenuated — higher early procedural harm and a greater relative contribution of IS+SE |

Procedural harm is concentrated in the first months; the non-procedural bleeding benefit accrues progressively over follow-up. Full weighting assumptions and sensitivity analyses are in the Supplementary Appendix.

## Citation

Preprint / in press: Martins Filho E, Filby SJ, Suruagy Motta RFO, Palma Dallan LA, Tebet MA, Tanboga IH. *Left Atrial Appendage Closure Versus Oral Anticoagulation in Atrial Fibrillation: A Domain-Based Meta-analysis of Clinical Trade-Offs.*

- PROSPERO: **CRD420261361187** (registered 06 April 2026)
- Statistical Analysis Plan (OSF): [10.17605/OSF.IO/ZVRDN](https://doi.org/10.17605/OSF.IO/ZVRDN)

## License

This repository contains analysis code and reconstructed datasets for academic research purposes. Individual trial data remain the property of the original investigators.
