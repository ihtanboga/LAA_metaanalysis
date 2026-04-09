###############################################################################
# 01_data.R
# LAA Closure vs OAC Meta-Analysis — Comprehensive Data Input
#
# Sources: CHAMPION-AF (NEJM 2026), OPTION (NEJM 2024), PRAGUE-17 (JACC 2022),
#          CLOSURE-AF (NEJM 2026), PROTECT-AF/PREVAIL 5yr pooled (JACC 2017),
#          + sub-analyses & supplement data
#
# Data provenance tiers:
#   T1 = trial-reported exact harmonized
#   T2 = trial-reported near-harmonized
#   T3 = KM-reconstructed fixed-horizon
#   T4 = derived / model-based
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "00_setup.R"))

# =============================================================================
# 1. STUDY DESIGN TABLE
# =============================================================================

study_design <- tibble::tribble(
  ~study,              ~year, ~journal,   ~n_total, ~n_device, ~n_control,
  ~device,             ~comparator,        ~fu_years, ~ni_margin,
  ~ni_type,            ~ni_met,  ~pool,

  "CHAMPION-AF",       2026, "NEJM",      3000, 1499, 1501,
  "Watchman FLX",      "NOAC",             3.0,  "4.8 ppt",
  "absolute",          TRUE,     "contemporary",

 "OPTION",            2024, "NEJM",      1600,  803,  797,
 "Watchman FLX",      "OAC (95% DOAC)",   3.0,  "5.0 ppt",
 "absolute",          TRUE,     "contemporary",

  "PRAGUE-17",         2022, "JACC",       402,  201,  201,
  "Amulet/Watchman",   "DOAC (96% apixaban)", 3.5, "sHR 1.47",
  "relative",          TRUE,     "contemporary",

  "CLOSURE-AF",        2026, "NEJM",       888,  446,  442,
  "Mixed (FLX/Amulet)","Best medical care", 3.0, "HR 1.30",
  "relative",          FALSE,    "expanded",

  "PROTECT-AF/PREVAIL",2017, "JACC",      1114,  732,  382,
  "Watchman",          "Warfarin",          5.0, "RR 2.0/1.75",
  "relative",          TRUE,     "historical"
)

# =============================================================================
# 2. BASELINE CHARACTERISTICS
# =============================================================================

baseline <- tibble::tribble(
  ~study,              ~age_mean, ~age_sd, ~female_pct, ~chadsvasc_mean,
  ~chadsvasc_sd, ~hasbled_mean, ~hasbled_sd, ~prior_stroke_pct,
  ~prior_bleed_pct, ~ckd_severe_pct, ~paroxysmal_af_pct,

  "CHAMPION-AF",       71.7, 7.5,  31.9, 3.5, 1.3, 1.3, 0.8,  7.8,  1.0,  NA,  68.9,
  "OPTION",            69.6, 7.7,  34.1, 3.5, 1.3, 1.2, 0.8,  NA,   1.0,  NA,  NA,
  "PRAGUE-17",         73.3, 7.0,  34.3, 4.7, 1.5, 3.1, 0.9, 35.3, 47.8,  32.0, NA,
  "CLOSURE-AF",        77.9, 7.1,  38.6, 5.2, 1.5, 3.0, 0.9, 33.0, 30.4,  24.1, NA,
  "PROTECT-AF/PREVAIL",72.5, 8.5,  30.0, 3.7, 1.4, NA,  NA,  18.0,  NA,   NA,  NA
)

# =============================================================================
# 3. ENDPOINT DATA — Events & Effect Sizes per Study
# =============================================================================

# Helper: derive log(HR) and SE from HR and 95% CI
derive_loghr <- function(hr, lo, hi) {
  loghr <- log(hr)
  se    <- (log(hi) - log(lo)) / (2 * qnorm(0.975))
  list(loghr = loghr, se = se)
}

# -------------------------------------------------------------------------
# A. PRIMARY EFFICACY: Ischemic Stroke + Systemic Embolism
# -------------------------------------------------------------------------
# Tier: T1 for OPTION; T2 for CHAMPION (IS+SE derived from components);
#       T2 for PRAGUE-17 (IS+SE from component counts, sHR for stroke+TIA);
#       T2 for CLOSURE-AF (components from supplement);
#       T1 for PP-5yr (pooled patient-level HR)

eff_is_se <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~effect_type, ~tier,

  # CHAMPION-AF: IS=45 vs 27; SE=0 vs 2 → IS+SE = 45 vs 29
  # HR from IS alone: 1.61 (1.00-2.59); SE ctrl=2 (Table S17)
  "CHAMPION-AF",         45, 1499,   29, 1501,
   1.61,  1.00,   2.59,   "HR",         "T2",

  # OPTION: IS+SE = 11 vs 11; HR 0.97 (0.42-2.25) (Table S17)
  "OPTION",              11,  803,   11,  797,
   0.97,  0.42,   2.25,   "HR",         "T1",

  # PRAGUE-17 4yr: IS=13, SE=0 → 13 vs IS=10, SE=1 → 11
  # sHR for all stroke excl TIA = 1.38 (0.63-3.03); proxy for IS+SE
  "PRAGUE-17",           13,  201,   11,  201,
   1.38,  0.63,   3.03,   "sHR",        "T2",

  # CLOSURE-AF: IS=18+SE=3 =21 vs IS=15+SE=1 =16
  # HR derived from KM-reconstructed IPD (trial reports dRMST only)
  "CLOSURE-AF",          21,  446,   16,  442,
   1.30,  0.69,   2.45,   "HR(derived)","T3",

  # PROTECT-AF/PREVAIL 5yr pooled: IS+SE HR 1.71 (0.94-3.11)
  # Events: PROTECT IS 24+SE 3 + PREVAIL IS 17+SE 1 = 45 vs 14 (Reddy 2017 Table 3/4)
  "PROTECT-AF/PREVAIL",  45,  732,   14,  382,
   1.71,  0.94,   3.11,   "HR",         "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# B. KEY SAFETY: Non-procedural Clinically Relevant Bleeding
#    (ISTH Major + CRNMB, excluding procedural)
# -------------------------------------------------------------------------

saf_nonproc_bleed <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~effect_type, ~tier,

  # CHAMPION-AF: non-proc ISTH major+CRNMB = 154 vs 260; HR 0.55 (0.45-0.67)
  "CHAMPION-AF",        154, 1499,  260, 1501,
   0.55,  0.45,   0.67,   "HR",         "T1",

  # OPTION: non-proc major+CRNMB = 65 vs 137; HR 0.44 (0.33-0.59)
  "OPTION",              65,  803,  137,  797,
   0.44,  0.33,   0.59,   "HR",         "T1",

  # PRAGUE-17 4yr: non-proc major+CRNMB sHR 0.55 (0.31-0.97)
  # 18 patients (23 events) vs 32 patients (40 events) — Table 2, 4yr
  "PRAGUE-17",           18,  201,   32,  201,
   0.55,  0.31,   0.97,   "sHR",        "T2",

  # CLOSURE-AF: non-proc major bleeding = 70-18proc = 52 vs 61
  # CRNMB not separately reported; non-proc major only
  # Direction: device slightly fewer non-proc major ≥6mo
  "CLOSURE-AF",          52,  446,   61,  442,
   0.86,  0.60,   1.24,   "HR(derived)","T3",

  # PP-5yr: post-procedure (>7d) major bleeding HR 0.48 (0.32-0.71)
  # Events: 48 vs 51 (Reddy 2017 Table 4)
  "PROTECT-AF/PREVAIL",  48,  732,   51,  382,
   0.48,  0.32,   0.71,   "HR",         "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# C. KEY HARM: Procedural Major Complications (device arm only)
# -------------------------------------------------------------------------

proc_complications <- tibble::tribble(
  ~study,              ~events, ~n_attempted, ~rate_pct, ~deaths,
  ~pericardial, ~embolization, ~stroke,

  "CHAMPION-AF",        14,  1408,  1.0,  0,    8,  0,  1,
  "OPTION",              8,   762,  1.0,  0,    2,  0,  0,
  "PRAGUE-17",           9,   187,  4.8,  2,    3,  1,  0,
  "CLOSURE-AF",         24,   421,  5.7,  2,    5,  1,  0,
  "PROTECT-AF",         36,   449,  7.7,  0,   22,  3,  5,
  "PREVAIL",            11,   265,  4.2,  0,    1,  2,  1
)

# -------------------------------------------------------------------------
# D. SUPPORTIVE SAFETY: All Major Bleeding (incl procedural)
# -------------------------------------------------------------------------

saf_all_major_bleed <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # CHAMPION-AF: ISTH major (proc+non-proc) 83 vs 87; HR 0.92 (0.68-1.24)
  "CHAMPION-AF",         83, 1499,   87, 1501,
   0.92,  0.68,   1.24,   "T1",

  # OPTION: ISTH major (proc incl) 30 vs 38; HR 0.77 (0.48-1.24)
  "OPTION",              30,  803,   38,  797,
   0.77,  0.48,   1.24,   "T1",

  # PRAGUE-17: clinically relevant bleeding (all) sHR 0.75 (0.44-1.27)
  # 24 patients vs 32 patients (Table 2, 4yr)
  "PRAGUE-17",           24,  201,   32,  201,
   0.75,  0.44,   1.27,   "T2",

  # CLOSURE-AF: major bleeding 70 vs 61
  "CLOSURE-AF",          70,  446,   61,  442,
   1.16,  0.83,   1.62,   "T2",

  # PP-5yr: all major bleeding HR 0.91 (0.64-1.29) — Reddy 2017 Table 4
  "PROTECT-AF/PREVAIL",  85,  732,   50,  382,
   0.91,  0.64,   1.29,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# E. ALL-CAUSE DEATH
# -------------------------------------------------------------------------

sec_acm <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # CHAMPION-AF: 71 vs 67; HR 1.01 (0.72-1.41)
  "CHAMPION-AF",         71, 1499,   67, 1501,
   1.01,  0.72,   1.41,   "T1",

  # OPTION: 29 vs 34; HR 0.83 (0.51-1.36)
  "OPTION",              29,  803,   34,  797,
   0.83,  0.51,   1.36,   "T1",

  # PRAGUE-17 4yr: HR 0.81 (0.54-1.22) — 42 vs 53 deaths (4yr text)
  "PRAGUE-17",           42,  201,   53,  201,
   0.81,  0.54,   1.22,   "T1",

  # CLOSURE-AF: 155 vs 141
  "CLOSURE-AF",         155,  446,  141,  442,
   1.11,  0.89,   1.39,   "T2",

  # PP-5yr: HR 0.73 (0.54-0.98) — 106 vs 73 deaths (Reddy 2017 Table 4)
  "PROTECT-AF/PREVAIL", 106,  732,   73,  382,
   0.73,  0.54,   0.98,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# F. CV / UNEXPLAINED DEATH
# -------------------------------------------------------------------------

sec_cvdeath <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # CHAMPION-AF: CV/unexpl death 38 vs 36; HR 1.01 (0.64-1.59) — Table S17
  "CHAMPION-AF",         38, 1499,   36, 1501,
   1.01,  0.64,   1.59,   "T1",

  # OPTION: 15 vs 15; HR 0.97 (0.48-1.99) — Table S17
  "OPTION",              15,  803,   15,  797,
   0.97,  0.48,   1.99,   "T1",

  # PRAGUE-17: sHR 0.68 (0.39-1.20) — 20 vs 30 (4yr Table 2)
  "PRAGUE-17",           20,  201,   30,  201,
   0.68,  0.39,   1.20,   "T1",

  # CLOSURE-AF: CV/unexplained death 99 vs 81
  "CLOSURE-AF",          99,  446,   81,  442,
   1.24,  0.93,   1.66,   "T2",

  # PP-5yr: CV/unexplained death HR 0.59 (0.37-0.94) — 39 vs 33 (Reddy 2017 Table 4)
  "PROTECT-AF/PREVAIL",  39,  732,   33,  382,
   0.59,  0.37,   0.94,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# G. HEMORRHAGIC STROKE
# -------------------------------------------------------------------------

sec_hemstroke <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # CHAMPION-AF: 5 vs 5; HR 0.96 (0.28-3.32)
  "CHAMPION-AF",          5, 1499,    5, 1501,
   0.96,  0.28,   3.32,   "T1",

  # OPTION: 3 vs 3; HR 0.98 (0.20-4.88)
  "OPTION",               3,  803,    3,  797,
   0.98,  0.20,   4.88,   "T1",

  # PRAGUE-17: 1 vs 2 (4yr text: 1 hem LAAC, 2 hem DOAC)
  "PRAGUE-17",            1,  201,    2,  201,
   NA,    NA,     NA,      "T1",

  # CLOSURE-AF: 10 vs 13
  "CLOSURE-AF",          10,  446,   13,  442,
   0.77,  0.34,   1.75,   "T2",

  # PP-5yr: HR 0.20 (0.07-0.56) — 5 vs 13 (Reddy 2017 Table 3)
  "PROTECT-AF/PREVAIL",   5,  732,   13,  382,
   0.20,  0.07,   0.56,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# H. ANY STROKE + SE (includes hemorrhagic)
# -------------------------------------------------------------------------

sec_anystroke_se <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # CHAMPION-AF: all stroke 50 + SE 0 = 50 vs stroke 33 + SE 2 = 35
  # HR for all strokes = 1.46 (0.94-2.27) — Table S17
  "CHAMPION-AF",         50, 1499,   35, 1501,
   1.46,  0.94,   2.27,   "T1",

  # OPTION: stroke+SE 14 vs 16; HR 0.85 (0.42-1.75)
  "OPTION",              14,  803,   16,  797,
   0.85,  0.42,   1.75,   "T1",

  # PRAGUE-17 4yr: all stroke excl TIA 14 + SE 0 = 14 vs 12 + SE 1 = 13
  # sHR 1.14 (0.56-2.30) for stroke/TIA — proxy
  "PRAGUE-17",           14,  201,   13,  201,
   1.14,  0.56,   2.30,   "T2",

  # CLOSURE-AF: stroke 27 + SE 3 = first-event composite 29 vs stroke 27 + SE 1 = 28
  "CLOSURE-AF",          29,  446,   28,  442,
   1.08,  0.64,   1.81,   "T2",

  # PP-5yr: all stroke/SE HR 0.96 (0.60-1.54) — 49 vs 27 (Reddy 2017 Table 4)
  "PROTECT-AF/PREVAIL",  49,  732,   27,  382,
   0.96,  0.60,   1.54,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# I. ISCHEMIC STROKE ALONE
# -------------------------------------------------------------------------

sec_is_alone <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  "CHAMPION-AF",         45, 1499,   27, 1501,
   1.61,  1.00,   2.59,   "T1",
  "OPTION",               9,  803,   10,  797,
   0.88,  0.36,   2.17,   "T1",
  # PRAGUE-17 4yr: IS 13 vs 10 (excl TIA)
  "PRAGUE-17",           13,  201,   10,  201,
   1.13,  0.44,   2.93,   "T2",
  "CLOSURE-AF",          18,  446,   15,  442,
   1.21,  0.61,   2.41,   "T2",
  # PP-5yr: IS 41 vs 14 (Reddy 2017 Table 3)
  "PROTECT-AF/PREVAIL",  41,  732,   14,  382,
   1.40,  0.76,   2.59,   "T1"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# J. DISABLING / FATAL STROKE
# -------------------------------------------------------------------------

sec_disabling_stroke <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~tier,

  # PP-5yr: HR 0.45 (0.21-0.94)
  "PROTECT-AF/PREVAIL",  11,  732,   13,  382,
   0.45,  0.21,   0.94,   "T1"
  # Other trials: data not available for this specific endpoint
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# K. TRIAL-REPORTED PRIMARY COMPOSITES (for reference/NACE)
# -------------------------------------------------------------------------

trial_primary <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,  ~primary_def,

  # CHAMPION-AF: CV death+stroke+SE: 81 vs 65; HR 1.20 (0.87-1.66)
  "CHAMPION-AF",         81, 1499,   65, 1501,
   1.20,  0.87,   1.66,  "CV death + stroke + SE",

  # OPTION: death+stroke+SE: 41 vs 44; HR 0.91 (0.59-1.39)
  "OPTION",              41,  803,   44,  797,
   0.91,  0.59,   1.39,  "Death + stroke + SE",

  # PRAGUE-17: composite sHR 0.81 (0.56-1.18)
  "PRAGUE-17",           57,  201,   68,  201,
   0.81,  0.56,   1.18,  "Stroke/TIA + SE + CV death + bleeding + complications",

  # CLOSURE-AF: NACE 155 vs 127; dRMST -0.36; HR(reconstruct) 1.26
  "CLOSURE-AF",         155,  446,  127,  442,
   1.26,  0.99,   1.60,  "Stroke + SE + major bleeding + CV/unexplained death",

  # PP composite: stroke+SE+CV death HR 0.82 (0.58-1.17)
  "PROTECT-AF/PREVAIL",  76,  732,   49,  382,
   0.82,  0.58,   1.17,  "Stroke + SE + CV death"
) %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

# -------------------------------------------------------------------------
# L. NET CLINICAL BENEFIT (trial-reported)
# -------------------------------------------------------------------------

ncb <- tibble::tribble(
  ~study,               ~ev_device, ~n_device, ~ev_control, ~n_control,
  ~hr,    ~hr_lo, ~hr_hi,

  # CHAMPION-AF: CV death+stroke+SE+non-proc bleeding: 215 vs 300; HR 0.66
  "CHAMPION-AF",        215, 1499,  300, 1501,
   0.66,  0.56,   0.79,

  # PP-5yr NCB: 1.42%/yr, RR 0.74
  "PROTECT-AF/PREVAIL",  NA,  732,   NA,  382,
   0.74,  0.56,   1.00
)

# =============================================================================
# 4. LANDMARK DATA (0-6 months vs >6 months)
# =============================================================================

# CLOSURE-AF landmark from plan_old
landmark_closure <- tibble::tribble(
  ~period,      ~ev_device, ~ev_control, ~hr,  ~hr_lo, ~hr_hi,
  "0-6 months",  68,         43,          1.63, 1.11,   2.39,
  ">6 months",   79,         78,          1.05, 0.77,   1.44
)

# =============================================================================
# 5. SUBGROUP INTERACTION DATA
# =============================================================================

# CHAMPION-AF subgroup interactions (safety endpoint)
subgroup_champion_safety <- tibble::tribble(
  ~subgroup,              ~category,    ~hr,  ~hr_lo, ~hr_hi,
  "Age",                  ">=75",       0.57, 0.42,   0.78,
  "Age",                  "<75",        0.54, 0.42,   0.69,
  "Sex",                  "Female",     0.71, 0.49,   1.02,
  "Sex",                  "Male",       0.48, 0.38,   0.61,
  "HAS-BLED",            ">=3",        0.39, 0.21,   0.71,
  "HAS-BLED",            "<3",         0.58, 0.47,   0.72,
  "CHA2DS2-VASc",        ">=4",        0.53, 0.39,   0.73,
  "CHA2DS2-VASc",        "<4",         0.55, 0.43,   0.71,
  "Prior embolic",        "Yes",        0.72, 0.37,   1.41,
  "Prior embolic",        "No",         0.53, 0.43,   0.66
)

# OPTION subgroup interactions (safety endpoint)
subgroup_option_safety <- tibble::tribble(
  ~subgroup,              ~category,    ~hr,  ~hr_lo, ~hr_hi,
  "Age",                  ">=75",       0.53, 0.28,   1.00,
  "Age",                  "<75",        0.41, 0.29,   0.57,
  "Sex",                  "Female",     0.53, 0.31,   0.90,
  "Sex",                  "Male",       0.40, 0.28,   0.57,
  "HAS-BLED",            ">=3",        0.40, 0.13,   1.26,
  "HAS-BLED",            "<3",         0.45, 0.33,   0.61,
  "CHA2DS2-VASc",        ">=4",        0.47, 0.30,   0.74,
  "CHA2DS2-VASc",        "<4",         0.42, 0.29,   0.60
)

# =============================================================================
# 6. ANTITHROMBOTIC REGIMEN DATA
# =============================================================================

antithrombotic <- tibble::tribble(
  ~study,             ~device_oac_start, ~device_dapt_start, ~device_ap_start,
  ~device_none_start, ~ctrl_oac_start,   ~ctrl_none_start,

  "CHAMPION-AF",       100,  0,   0,   0,  100,   0,
  "OPTION",              0, 88,  12,   0,  100,   0,
  "PRAGUE-17",           0, 65,  35,   0,  100,   0,
  "CLOSURE-AF",        9.1,79.6, 4.9, 4.0, 88.3,  7.5,
  "PROTECT-AF",        100,  0,   0,   0,  100,   0,
  "PREVAIL",           100,  0,   0,   0,  100,   0
)

# =============================================================================
# 7. DEVICE-RELATED THROMBUS (DRT) DATA
# =============================================================================

drt_data <- tibble::tribble(
  ~study,              ~drt_events, ~drt_n_screened, ~drt_pct, ~timepoint,

  "CHAMPION-AF",        63, 1320, 4.8, "4 months",
  "OPTION",             14,  738, 1.9, "12 months",
  "PRAGUE-17",          NA,   NA,  NA, NA,
  "CLOSURE-AF",         14,  320, 4.4, "3 months",
  "PP (Dukkipati 2018)", 65, 1739, 3.7, "pooled"
)

# =============================================================================
# 8. PROCEDURAL PERFORMANCE (single-arm meta-analysis data)
# =============================================================================

procedural_perf <- tibble::tribble(
  ~study,              ~implant_success, ~n_attempted,  ~pericardial_eff,
  ~device_embol, ~proc_stroke, ~proc_death,

  "CHAMPION-AF",       1386, 1403,  10, 0, 0, 0,
  "OPTION",             753,  762,   2, 0, 0, 0,
  "PRAGUE-17",          181,  187,   4, 1, 0, 2,
  "CLOSURE-AF",         414,  421,   5, 1, 0, 2,
  "PROTECT-AF",         408,  449,  22, 3, 5, 0,
  "PREVAIL",            252,  265,   1, 0, 1, 0
)

# =============================================================================
# 9. CLOSURE-AF ENDPOINT CASCADE (reconstructed IPD)
# =============================================================================

closure_cascade <- tibble::tribble(
  ~model,                ~endpoint_def,                          ~ev_device, ~ev_ctrl,
  ~hr_mc,   ~drmst,

  "Published NACE",      "Stroke+SE+bleeding+CV/unexp death",    155, 127,
   1.26,    -0.36,
  "Reconstructed NACE",  "Same (KM-reconstructed)",              147, 121,
   1.26,    -0.343,
  "Model 1: MACE",       "Stroke+SE+CV/unexp death",              81,  63,
   1.29,    -0.245,
  "Model 2: Thromboembolic", "Stroke+SE only",                    27,  27,
   1.02,    -0.016
)

# =============================================================================
# 10. POOL ASSIGNMENTS (for grouped analysis)
# =============================================================================

pool_assign <- tibble::tribble(
  ~study,               ~pool_primary, ~pool_expanded, ~pool_historical,
  ~pool_flx,

  "CHAMPION-AF",         TRUE,  TRUE,   FALSE, TRUE,
  "OPTION",              TRUE,  TRUE,   FALSE, TRUE,
  "PRAGUE-17",           TRUE,  TRUE,   FALSE, FALSE,
  "CLOSURE-AF",          FALSE, TRUE,   FALSE, FALSE,
  "PROTECT-AF/PREVAIL",  FALSE, FALSE,  TRUE,  FALSE
)

# =============================================================================
# 11. NI MARGIN RE-EXPRESSION DATA
# =============================================================================

ni_margin_data <- tibble::tribble(
  ~study,              ~rr_115, ~rr_125, ~rr_140, ~abs_05, ~abs_10, ~abs_20,

  # rr_115 = "is upper 95% CI of pooled RR < 1.15?"
  # TRUE = compatible with NI at that threshold
  "CHAMPION-AF",        FALSE,  FALSE,   TRUE,    FALSE,   TRUE,    TRUE,
  "OPTION",             TRUE,   TRUE,    TRUE,    TRUE,    TRUE,    TRUE,
  "PRAGUE-17",          TRUE,   TRUE,    TRUE,    TRUE,    TRUE,    TRUE,
  "CLOSURE-AF",         FALSE,  FALSE,   FALSE,   FALSE,   FALSE,   TRUE,
  "PROTECT-AF/PREVAIL", FALSE,  FALSE,   FALSE,   FALSE,   FALSE,   FALSE
)

# =============================================================================
# 12. LAAOS-III BENCHMARK DATA (for Discussion only)
# =============================================================================

laaos <- list(
  n_total     = 4770,
  n_occl      = 2379,
  n_no_occl   = 2391,
  fu_years    = 3.8,
  primary_hr  = 0.67,
  primary_lo  = 0.53,
  primary_hi  = 0.85,
  acm_hr      = 1.00,
  acm_lo      = 0.89,
  acm_hi      = 1.13,
  bleed_hr    = 0.93,
  note        = "Surgical LAA occlusion during cardiac surgery; patients continued OAC"
)

# =============================================================================
# SAVE ALL DATA
# =============================================================================

save(
  study_design, baseline, eff_is_se, saf_nonproc_bleed,
  proc_complications, saf_all_major_bleed,
  sec_acm, sec_cvdeath, sec_hemstroke,
  sec_anystroke_se, sec_is_alone, sec_disabling_stroke,
  trial_primary, ncb, landmark_closure,
  subgroup_champion_safety, subgroup_option_safety,
  antithrombotic, drt_data, procedural_perf,
  closure_cascade, pool_assign, ni_margin_data, laaos,
  file = file.path(DIR_D, "meta_data.RData")
)

cat(">> 01_data.R loaded and saved. Datasets:\n")
cat("   - study_design          (5 trials)\n")
cat("   - baseline              (5 trials)\n")
cat("   - eff_is_se             (primary efficacy)\n")
cat("   - saf_nonproc_bleed     (key safety)\n")
cat("   - saf_all_major_bleed   (supportive safety)\n")
cat("   - sec_acm, sec_cvdeath, sec_hemstroke, etc.\n")
cat("   - proc_complications    (procedural harm)\n")
cat("   - pool_assign           (pool membership)\n")
