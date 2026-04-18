# LAA Closure vs Oral Anticoagulation: Domain-Based Meta-Analysis

A domain-based meta-analysis with harmonised endpoint reconstruction evaluating **percutaneous left atrial appendage closure (LAAC) versus oral anticoagulation (OAC)** in atrial fibrillation.

## Key Findings

| Endpoint | Contemporary Pool (k=3) | All Trials (k=5) | Direction |
|---|---|---|---|
| Ischaemic stroke + SE | RR 1.35 (0.77–2.35) | RR 1.41 (1.11–1.78)* | **Significantly higher with LAAC (All Trials)** |
| Non-procedural bleeding | HR 0.51 (0.37–0.70)* | HR 0.56 (0.40–0.76)* | **Significantly lower with LAAC** |
| All major bleeding | HR 0.85 (0.64–1.13) | HR 0.93 (0.76–1.14) | Neutral |
| All-cause death | RR 0.91 (0.60–1.38) | RR 0.92 (0.74–1.15) | Neutral |
| Haemorrhagic stroke | HR 0.97 (0.85–1.10) | HR 0.58 (0.17–1.97) | Neutral (DOAC-era: no difference) |

\* P < 0.05

## Included Trials

| Trial | N | Comparator | Device | Follow-up |
|---|---|---|---|---|
| CHAMPION-AF | 3,000 | NOAC | Watchman FLX | 3.0 yr |
| OPTION | 1,600 | OAC (95% DOAC) | Watchman FLX | 3.0 yr |
| PRAGUE-17 | 402 | DOAC (96% apixaban) | Amulet (61%)/Watchman (36%) | 3.5 yr (median) |
| CLOSURE-AF | 888 | Best medical care (85% DOAC) | Watchman/FLX/Amulet | 3.0 yr (median) |
| PROTECT-AF/PREVAIL | 1,114 | Warfarin | Watchman Gen 1 | 5.0 yr (max) |

**Total: 7,004 patients across 6 RCTs (5 analysis units)**

## Analysis Pools

- **Primary Contemporary**: CHAMPION-AF + OPTION + PRAGUE-17 (n=5,002)
- **Expanded Contemporary**: + CLOSURE-AF (n=5,890)
- **Historical Warfarin Benchmark**: PROTECT-AF/PREVAIL pooled (n=1,114)

## Methodology

- **Domain-based endpoint harmonisation**: trial-specific composites deconstructed into clinically interpretable domains (efficacy, safety, harm) rather than pooling incomparable composites
- **KM curve reconstruction**: Kaplan-Meier digitisation and IPD reconstruction via the Guyot method for endpoints not directly reported
- **Comparator-era stratification**: warfarin-era vs DOAC-era data analysed separately
- **Statistical framework**: Two-stage random-effects (REML + HKSJ), Mantel-Haenszel RR as primary estimand, HR-based pooling as secondary
- **Data provenance tiers (T1-T4)**: tracking the directness and reliability of each data point

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

## License

This repository contains analysis code and reconstructed datasets for academic research purposes. Individual trial data remain the property of the original investigators.
