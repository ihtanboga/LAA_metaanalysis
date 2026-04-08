# LAA Closure vs Oral Anticoagulation: Domain-Based Meta-Analysis

A domain-based meta-analysis with harmonised endpoint reconstruction evaluating **percutaneous left atrial appendage closure (LAAC) versus oral anticoagulation (OAC)** in atrial fibrillation.

## Key Findings

| Endpoint | Contemporary Pool (k=3) | Direction |
|---|---|---|
| Ischaemic stroke + SE | RR 1.33 (0.88-2.02) | Numerically higher with LAAC |
| Non-procedural bleeding | HR 0.52 (0.44-0.62) | **Significantly lower with LAAC** |
| All major bleeding | HR 0.85 (0.64-1.13) | Numerically lower with LAAC |
| All-cause death | RR 0.95 (0.78-1.16) | Neutral |
| Haemorrhagic stroke | HR 0.57 (0.27-1.19) | Numerically lower with LAAC |

## Included Trials

| Trial | N | Comparator | Device | Follow-up |
|---|---|---|---|---|
| CHAMPION-AF | 3,000 | NOAC | Watchman FLX | 3.0 yr |
| OPTION | 1,600 | OAC (95% DOAC) | Watchman FLX | 3.0 yr |
| PRAGUE-17 | 402 | DOAC (96% apixaban) | Amulet/Watchman | 3.5 yr |
| CLOSURE-AF | 888 | Best medical care | Mixed | 3.0 yr |
| PROTECT-AF/PREVAIL | 1,114 | Warfarin | Watchman Gen 1 | 5.0 yr |

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
  06_main_figures.R         # 5 main figures
  07_suppl_figures.R        # 6 supplementary figures
  08_tables.R               # Main + supplementary tables
  run_all.R                 # Master pipeline script
  closure_af_km.R           # CLOSURE-AF KM reconstruction
  closure_af_bleeding_free.R    # CLOSURE-AF endpoint cascade
  closure_af_death_free.R       # CLOSURE-AF death-free analysis
  closure_af_sensitivity.R      # CLOSURE-AF sensitivity analyses

data/
  km_reconstructed/         # 46 digitised KM curve CSVs (all trials)
  closure_af_ipd/           # Reconstructed pseudo-IPD for CLOSURE-AF
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
source("R/06_main_figures.R")     # Generate main figures
source("R/07_suppl_figures.R")    # Generate supplementary figures
```

### Requirements

R >= 4.3 with packages: `meta`, `metafor`, `survival`, `survRM2`, `tidyverse`, `patchwork`, `gt`, `flextable`, `officer`, `forestploter`, `ggtext`, `scales`, `RColorBrewer`, `gridExtra`, `checkmate`

## KM Curve Digitisation

The KM-reconstructed datasets in `data/km_reconstructed/` were generated using our open-source KM digitisation tool:

**[KM Reconstruction Tool (Demo)](https://ihtanboga.github.io/KMR/)**

This web-based tool enables extraction of individual patient-level time-to-event data from published Kaplan-Meier survival curves using the Guyot algorithm.

## Citation

> Tanboga IH, et al. Percutaneous Left Atrial Appendage Closure Versus Oral Anticoagulation in Atrial Fibrillation: A Domain-Based Meta-Analysis With Harmonised Endpoint Reconstruction. *European Heart Journal* (submitted).

## License

This repository contains analysis code and reconstructed datasets for academic research purposes. Individual trial data remain the property of the original investigators.
