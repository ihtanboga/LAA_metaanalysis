# Supplementary Methods: Individual Patient Data Reconstruction and Endpoint Derivation

## Overview

Because individual patient data (IPD) were not available from the source trials, we reconstructed time-to-event data from published Kaplan–Meier (KM) curves for every reported bleeding and thromboembolic endpoint across the included randomized controlled trials (CHAMPION-AF, OPTION, PRAGUE-17, and CLOSURE-AF). Where the published endpoint definitions did not align across trials (e.g., composite endpoints mixing death, stroke, and systemic embolism), we applied principled Monte Carlo event-thinning procedures to derive harmonized endpoints suitable for meta-analysis. All analyses were performed in R (v4.5.2) using the `reconstructKM`, `survival`, and `survRM2` packages.

## 1. KM Curve Digitization

Published KM curves (cumulative incidence format) were digitized from the primary manuscripts and supplementary appendices of each trial. Each curve was manually traced using WebPlotDigitizer, producing (time, cumulative incidence %) coordinate pairs. Number-at-risk (NAR) tables were extracted directly from the figures at their reported time points (typically yearly intervals). For each endpoint, both arms (anticoagulation [AC] and device) were digitized separately and stored as arm-specific CSV files with accompanying JSON metadata containing the NAR table and calibration parameters.

## 2. IPD Reconstruction (Guyot Algorithm)

Time-to-event IPD were reconstructed from each digitized KM curve using the algorithm of **Guyot et al. (BMC Med Res Methodol 2012;12:9)**, as implemented in the `reconstructKM` R package. The pipeline for each arm was:

1. Convert the digitized cumulative incidence values to survival probabilities: S(t) = 1 − CIF(t)/100.
2. Enforce monotonically non-increasing survival and a fixed origin at (t = 0, S = 1); noise points at t ≤ 0 from digitization were discarded before re-anchoring the origin.
3. Trim any trailing horizontal segment beyond the last observed drop (required by `format_raw_tabs`).
4. Pass the cleaned (time, survival) sequence together with the NAR table to `format_raw_tabs()` followed by `KM_reconstruct()`, which iteratively allocates events and censoring times within each NAR interval consistent with the published KM step pattern.

The reconstructed IPD for each trial × endpoint pair were then combined across arms into a single file with three columns (`arm`, `time`, `status`) and saved as a per-endpoint CSV. A summary table listing N and event counts for every reconstructed cohort is provided (`_reconstruction_summary.csv`). Typical digitization loss was 3–5% of events relative to the published totals, consistent with prior benchmarking of the Guyot algorithm.

## 3. Derivation of Harmonized Endpoints by Monte Carlo Event Thinning

Several included trials reported only composite endpoints that mixed the event types of primary interest (bleeding, stroke/systemic embolism, and death). Because the reconstructed composite IPD does not retain the identity of each patient's first event, we used Monte Carlo (MC) event thinning to derive harmonized single-component endpoints from the composite IPD. The general approach was:

- **Targets:** number of events to remove in each arm, taken from the published event totals (or, where not reported, estimated from arm-specific KM cumulative incidence values at a reference time point multiplied by the enrolled N).
- **Scaling:** targets were multiplied by an arm-specific scale ratio (reconstructed events ÷ published events) to correct for the ~5% digitization undercount.
- **Stochastic rounding:** non-integer removal targets were rounded stochastically, i.e., `⌊x⌋ + Bernoulli(x − ⌊x⌋)`, to preserve the expected count across simulations.
- **Extended follow-up:** each "removed" event patient had their status set to 0 and their observation time re-sampled from the arm-specific distribution of administratively censored times occurring after the original event time. This preserves the study's observed censoring pattern.
- **Simulation:** for each derivation, 1,000 MC replicates were run; the representative IPD saved to disk was the replicate whose restricted mean survival time difference (Δ-RMST, tau = 6 years for CLOSURE-AF; tau = 3 years for OPTION) was closest to the across-simulation median.

Two selection rules were used depending on whether published temporal information was available for the event type being removed (or retained):

### 3.1 Temporal-weighted thinning

When the source publication reported a stratified event count by pre-specified time intervals (e.g., < 3 months, 3–6 months, ≥ 6 months), removals were drawn independently from the composite events falling into each interval, in proportion to that interval's share of the total. For CLOSURE-AF, this approach was used for major-bleeding events, stratified according to Supplementary Table S11 (LAA Closure: 57.1% < 3 months, 10.0% 3–6 months, 32.9% ≥ 6 months; Medical Therapy: 27.9% / 9.8% / 62.3%). This preserves the strong early-vs-late asymmetry of peri-procedural versus chronic-anticoagulation bleeding.

### 3.2 Proportional (uniform) thinning

When no temporal breakdown was published for the event type in question, removals were drawn uniformly at random from the full set of composite events in each arm. This was used for all death removals (no trial provided a temporal breakdown of all-cause or CV death), as well as for any event type whose temporal distribution could not be reasonably approximated from the source reporting.

## 4. Endpoint-Specific Derivations

The following harmonized IPD files were derived by applying the thinning procedure above:

| File | Source composite | Removal rule | Removal target (device / AC) |
|---|---|---|---|
| `CLOSURE_majorbleeding.csv` | CLOSURE-AF primary composite (stroke + SE + major bleeding + CV/unexplained death) | Temporal-weighted (S11) to *keep* bleeding; all non-bleeding events censored with extended follow-up | Keep ~66.4 / 58.1 bleeding events (Table 2: 70 / 61, scaled) |
| `CLOSURE_stroke_se.csv` | CLOSURE-AF primary composite | Two-phase: remove bleeding (temporal, S11) then remove CV/unexplained death (proportional) | Remove bleeding 66.4 / 58.1 + death 53.1 / 36.2 (scaled from Table 2) |
| `OPTION_stroke_se.csv` | OPTION composite (death from any cause + stroke + SE) | Proportional removal of deaths | Remove 28.3 / 32.6 death events (published 36-month KM: 29 / 34, scaled) |

The CLOSURE-AF bleeding extraction is the conceptual inverse of the "bleeding-free composite" sensitivity analysis: rather than removing bleeding events to isolate the non-bleeding components, we retained only the events probabilistically attributed to bleeding and censored all others.

## 5. Limitations and Assumptions

1. **Event identity is unobserved.** Reconstructed IPD contain only composite event indicators; the mapping from a specific IPD event to its underlying component (bleeding, stroke, death) is probabilistic, not deterministic. All derived single-component endpoints therefore inherit simulation-level uncertainty, which we quantified as the 2.5th–97.5th percentile range of the effect estimate across the 1,000 MC replicates.

2. **Temporal distribution assumed for uncontrolled components.** For CLOSURE-AF, Supplementary Table S11 provides a temporal breakdown only for major bleeding. When constructing the stroke + SE endpoint, death removals used proportional (uniform) thinning because no analogous temporal stratification was published for CV or unexplained death. If deaths were strongly concentrated in a specific time window, the derived endpoint could shift modestly; however, the near-null result for stroke + SE was robust to this assumption in prior sensitivity analyses.

3. **KM vs. cumulative incidence function.** Our reconstruction is based on 1 − KM curves; when the source trial reported competing-risk-adjusted cumulative incidence (Aalen-Johansen), the reconstructed event rate may slightly exceed the true component-specific incidence in the presence of competing events. The magnitude of this bias is small (~0.02 yr on Δ-RMST in CLOSURE-AF) and is acknowledged as a limitation.

4. **Endpoint definition heterogeneity is preserved, not harmonized at the definition level.** The thinning procedures harmonize composite structure but cannot reconcile differences in adjudication standard (e.g., ISTH major vs. BARC 3–5 vs. "clinically relevant bleeding"), inclusion/exclusion of procedure-related events, or inclusion of TIA in stroke endpoints. These definitional differences are summarized in Supplementary Table [X] and are addressed in the meta-analysis through pre-specified subgroup and sensitivity analyses.

5. **Digitization error.** Guyot reconstruction introduces a typical 3–5% event undercount relative to published totals, scaled proportionally across arms. All removal targets were rescaled to reconstructed event counts to preserve effect-estimate consistency.

## 6. Reproducibility

All reconstruction and derivation scripts, together with the input JSON/CSV digitization files, are archived in the project repository under `json_csv/` (digitized curves and metadata), `ipd_data/` (reconstructed IPD per trial × endpoint), and the corresponding R scripts (`reconstruct_all_endpoints.R`, `closure_af_bleeding_only.R`, `closure_af_stroke_se_only.R`, `option_stroke_se.R`). Random seeds were fixed (set.seed(2026)) for exact reproducibility of the representative IPD.

## Software

R 4.5.2; R packages: `reconstructKM`, `survival`, `survRM2`, `survminer`, `dplyr`, `readr`, `jsonlite`, `ggplot2`, `ggpubr`.

## Key reference

Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan–Meier survival curves. *BMC Med Res Methodol.* 2012;12:9.
