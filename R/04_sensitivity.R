###############################################################################
# 04_sensitivity.R
# LAA Closure vs OAC Meta-Analysis — Sensitivity & Robustness Analyses
#
# Contents:
#   1. Leave-one-out (expanded pool)
#   2. Fixed-effect vs random-effects comparison
#   3. REML vs Paule-Mandel tau-squared estimator
#   4. HKSJ vs DerSimonian-Laird CI comparison
#   5. Tier 1-2 only (exclude T3/T4 derived data)
#   6. FLX-only subset (CHAMPION + OPTION)
#   7. Contemporary-only (exclude CLOSURE and PP)
#   8. Historical benchmark (PROTECT-AF/PREVAIL alone)
#   9. Benefit-risk per 1000 patients over 3 years
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "01_data.R"))

cat("\n========================================\n")
cat("  04_sensitivity.R\n")
cat("========================================\n\n")

# =============================================================================
# 0. HELPER FUNCTIONS
# =============================================================================

get_pool_studies <- function(pool_name) {
  switch(pool_name,
    "contemporary" = pool_assign %>% filter(pool_primary)   %>% pull(study),
    "expanded"     = pool_assign %>% filter(pool_expanded)  %>% pull(study),
    "all"          = pool_assign$study,
    "flx"          = pool_assign %>% filter(pool_flx)       %>% pull(study),
    stop("Unknown pool: ", pool_name)
  )
}

# Generic RR meta-analysis with configurable options
run_rr <- function(data, studies, label,
                   method.tau = "REML", hakn = TRUE, comb.random = TRUE) {
  d <- data %>% filter(study %in% studies)
  if (nrow(d) < 2) return(NULL)
  metabin(
    event.e = ev_device, n.e = n_device,
    event.c = ev_control, n.c = n_control,
    studlab = study, data = d, sm = "RR",
    method = "MH", method.tau = method.tau,
    hakn = hakn, prediction = TRUE,
    incr = 0.5, allincr = FALSE,
    comb.random = comb.random,
    title = label
  )
}

# Generic HR meta-analysis with configurable options
run_hr <- function(data, studies, label,
                   method.tau = "REML", hakn = TRUE, comb.random = TRUE) {
  d <- data %>%
    filter(study %in% studies) %>%
    filter(!is.na(loghr) & !is.na(se_loghr))
  if (nrow(d) < 2) return(NULL)
  metagen(
    TE = loghr, seTE = se_loghr,
    studlab = study, data = d, sm = "HR",
    method.tau = method.tau, hakn = hakn,
    prediction = TRUE, comb.random = comb.random,
    title = label
  )
}

# Extract key metrics from meta object
extract_metrics <- function(m, type = "RE") {
  if (is.null(m)) return(NULL)
  if (type == "RE") {
    tibble(
      estimate = exp(m$TE.random),
      ci_lo    = exp(m$lower.random),
      ci_hi    = exp(m$upper.random),
      p_value  = m$pval.random,
      pi_lo    = exp(m$lower.predict),
      pi_hi    = exp(m$upper.predict),
      tau2     = m$tau2,
      I2       = m$I2,
      Q        = m$Q,
      Q_pval   = m$pval.Q,
      k        = m$k
    )
  } else {
    tibble(
      estimate = exp(m$TE.fixed),
      ci_lo    = exp(m$lower.fixed),
      ci_hi    = exp(m$upper.fixed),
      p_value  = m$pval.fixed,
      pi_lo    = NA_real_,
      pi_hi    = NA_real_,
      tau2     = m$tau2,
      I2       = m$I2,
      Q        = m$Q,
      Q_pval   = m$pval.Q,
      k        = m$k
    )
  }
}

# =============================================================================
# 1. LEAVE-ONE-OUT SENSITIVITY (Expanded Pool)
# =============================================================================

cat("--- 1. LEAVE-ONE-OUT SENSITIVITY (Expanded Pool) ---\n\n")

expanded_studies <- get_pool_studies("expanded")

# Focus on the two primary endpoints: IS+SE (efficacy) and NonProc Bleed (safety)
loo_targets <- list(
  list(data = eff_is_se,        name = "IS+SE",          sm = "RR"),
  list(data = eff_is_se,        name = "IS+SE",          sm = "HR"),
  list(data = saf_nonproc_bleed, name = "NonProc Bleed",  sm = "RR"),
  list(data = saf_nonproc_bleed, name = "NonProc Bleed",  sm = "HR")
)

loo_results <- list()

for (target in loo_targets) {
  cat(sprintf("  Leave-one-out: %s (%s)\n", target$name, target$sm))

  loo_rows <- list()

  for (drop_study in expanded_studies) {
    remaining <- setdiff(expanded_studies, drop_study)

    if (target$sm == "RR") {
      m <- run_rr(target$data, remaining,
                  sprintf("%s RR excl %s", target$name, drop_study))
    } else {
      m <- run_hr(target$data, remaining,
                  sprintf("%s HR excl %s", target$name, drop_study))
    }

    if (!is.null(m)) {
      metrics <- extract_metrics(m, "RE")
      metrics$dropped <- drop_study
      loo_rows[[drop_study]] <- metrics
    }
  }

  key <- sprintf("%s_%s", target$name, target$sm)
  loo_df <- bind_rows(loo_rows)
  loo_results[[key]] <- loo_df

  # Print
  for (i in 1:nrow(loo_df)) {
    r <- loo_df[i, ]
    cat(sprintf("    Excl %-22s: %s = %.2f (%.2f-%.2f), I2 = %.0f%%\n",
                r$dropped, target$sm, r$estimate, r$ci_lo, r$ci_hi, r$I2 * 100))
  }
  cat("\n")
}

# =============================================================================
# 2. FIXED-EFFECT vs RANDOM-EFFECTS COMPARISON
# =============================================================================

cat("--- 2. FIXED-EFFECT vs RANDOM-EFFECTS ---\n\n")

fe_vs_re_list <- list()

for (ep_info in list(
  list(data = eff_is_se,          name = "IS+SE"),
  list(data = saf_nonproc_bleed,  name = "NonProc Bleed"),
  list(data = saf_all_major_bleed, name = "All Major Bleed"),
  list(data = sec_acm,            name = "All-Cause Mortality")
)) {
  m_rr <- run_rr(ep_info$data, expanded_studies,
                  sprintf("%s RR expanded", ep_info$name))
  m_hr <- run_hr(ep_info$data, expanded_studies,
                  sprintf("%s HR expanded", ep_info$name))

  if (!is.null(m_rr)) {
    re <- extract_metrics(m_rr, "RE") %>% mutate(model = "RE", sm = "RR")
    fe <- extract_metrics(m_rr, "FE") %>% mutate(model = "FE", sm = "RR")
    fe_vs_re_list[[sprintf("%s_RR", ep_info$name)]] <- bind_rows(re, fe) %>%
      mutate(endpoint = ep_info$name)
  }

  if (!is.null(m_hr)) {
    re <- extract_metrics(m_hr, "RE") %>% mutate(model = "RE", sm = "HR")
    fe <- extract_metrics(m_hr, "FE") %>% mutate(model = "FE", sm = "HR")
    fe_vs_re_list[[sprintf("%s_HR", ep_info$name)]] <- bind_rows(re, fe) %>%
      mutate(endpoint = ep_info$name)
  }
}

fe_vs_re <- bind_rows(fe_vs_re_list)

cat("  Endpoint             | SM | Model | Estimate (95% CI)       | tau2    | I2\n")
cat("  ---------------------------------------------------------------------------\n")
for (i in 1:nrow(fe_vs_re)) {
  r <- fe_vs_re[i, ]
  cat(sprintf("  %-22s | %s | %-5s | %.2f (%.2f-%.2f)    | %.4f  | %.0f%%\n",
              r$endpoint, r$sm, r$model,
              r$estimate, r$ci_lo, r$ci_hi,
              r$tau2, r$I2 * 100))
}

# =============================================================================
# 3. REML vs PAULE-MANDEL TAU-SQUARED ESTIMATOR
# =============================================================================

cat("\n--- 3. REML vs PAULE-MANDEL TAU-SQUARED ---\n\n")

tau_comparison <- list()

for (ep_info in list(
  list(data = eff_is_se,         name = "IS+SE"),
  list(data = saf_nonproc_bleed, name = "NonProc Bleed")
)) {
  m_reml <- run_rr(ep_info$data, expanded_studies,
                    sprintf("%s REML", ep_info$name), method.tau = "REML")
  m_pm   <- run_rr(ep_info$data, expanded_studies,
                    sprintf("%s PM", ep_info$name), method.tau = "PM")

  if (!is.null(m_reml) && !is.null(m_pm)) {
    tau_comparison[[ep_info$name]] <- tibble(
      endpoint    = ep_info$name,
      method      = c("REML", "Paule-Mandel"),
      estimate    = c(exp(m_reml$TE.random), exp(m_pm$TE.random)),
      ci_lo       = c(exp(m_reml$lower.random), exp(m_pm$lower.random)),
      ci_hi       = c(exp(m_reml$upper.random), exp(m_pm$upper.random)),
      tau2        = c(m_reml$tau2, m_pm$tau2),
      I2          = c(m_reml$I2, m_pm$I2)
    )
  }
}

tau_comp_df <- bind_rows(tau_comparison)
cat("  Endpoint          | Method       | RR (95% CI)       | tau2   | I2\n")
cat("  -------------------------------------------------------------------\n")
for (i in 1:nrow(tau_comp_df)) {
  r <- tau_comp_df[i, ]
  cat(sprintf("  %-19s | %-12s | %.2f (%.2f-%.2f) | %.4f | %.0f%%\n",
              r$endpoint, r$method, r$estimate, r$ci_lo, r$ci_hi,
              r$tau2, r$I2 * 100))
}

# =============================================================================
# 4. HKSJ vs DERSIMONIAN-LAIRD CI
# =============================================================================

cat("\n--- 4. HKSJ vs DerSimonian-LAIRD ---\n\n")

hksj_vs_dl <- list()

for (ep_info in list(
  list(data = eff_is_se,         name = "IS+SE"),
  list(data = saf_nonproc_bleed, name = "NonProc Bleed")
)) {
  m_hksj <- run_rr(ep_info$data, expanded_studies,
                    sprintf("%s HKSJ", ep_info$name), hakn = TRUE)
  m_dl   <- run_rr(ep_info$data, expanded_studies,
                    sprintf("%s DL", ep_info$name), hakn = FALSE)

  if (!is.null(m_hksj) && !is.null(m_dl)) {
    hksj_vs_dl[[ep_info$name]] <- tibble(
      endpoint    = ep_info$name,
      ci_method   = c("HKSJ", "DerSimonian-Laird"),
      estimate    = c(exp(m_hksj$TE.random), exp(m_dl$TE.random)),
      ci_lo       = c(exp(m_hksj$lower.random), exp(m_dl$lower.random)),
      ci_hi       = c(exp(m_hksj$upper.random), exp(m_dl$upper.random)),
      ci_width    = c(exp(m_hksj$upper.random) - exp(m_hksj$lower.random),
                      exp(m_dl$upper.random) - exp(m_dl$lower.random))
    )
  }
}

hksj_dl_df <- bind_rows(hksj_vs_dl)
cat("  Endpoint          | CI Method         | RR (95% CI)       | CI width\n")
cat("  ---------------------------------------------------------------------\n")
for (i in 1:nrow(hksj_dl_df)) {
  r <- hksj_dl_df[i, ]
  cat(sprintf("  %-19s | %-17s | %.2f (%.2f-%.2f) | %.3f\n",
              r$endpoint, r$ci_method, r$estimate, r$ci_lo, r$ci_hi, r$ci_width))
}

# =============================================================================
# 5. TIER 1-2 ONLY (exclude T3/T4 derived data)
# =============================================================================

cat("\n--- 5. TIER 1-2 ONLY (exclude T3/T4) ---\n\n")

# Filter to tier T1 or T2 only
tier12_results <- list()

for (ep_info in list(
  list(data = eff_is_se,          name = "IS+SE"),
  list(data = saf_nonproc_bleed,  name = "NonProc Bleed"),
  list(data = saf_all_major_bleed, name = "All Major Bleed")
)) {
  d_tier <- ep_info$data %>%
    filter(study %in% expanded_studies) %>%
    filter(tier %in% c("T1", "T2"))

  studies_tier <- d_tier$study

  if (length(studies_tier) >= 2) {
    m_rr <- run_rr(ep_info$data %>% filter(tier %in% c("T1", "T2")),
                    studies_tier,
                    sprintf("%s RR T1-T2", ep_info$name))
    if (!is.null(m_rr)) {
      metrics <- extract_metrics(m_rr, "RE")
      metrics$endpoint <- ep_info$name
      metrics$analysis <- "Tier 1-2 only"
      tier12_results[[ep_info$name]] <- metrics
      cat(sprintf("  %-20s: RR = %.2f (%.2f-%.2f), PI (%.2f-%.2f), I2 = %.0f%%, k = %d\n",
                  ep_info$name, metrics$estimate, metrics$ci_lo, metrics$ci_hi,
                  metrics$pi_lo, metrics$pi_hi, metrics$I2 * 100, metrics$k))
    }
  } else {
    cat(sprintf("  %-20s: < 2 T1/T2 studies; skipping.\n", ep_info$name))
  }
}

tier12_df <- bind_rows(tier12_results)

# =============================================================================
# 6. FLX-ONLY SUBSET (CHAMPION + OPTION)
# =============================================================================

cat("\n--- 6. FLX-ONLY (CHAMPION + OPTION) ---\n\n")

flx_studies <- get_pool_studies("flx")
flx_results <- list()

for (ep_info in list(
  list(data = eff_is_se,          name = "IS+SE"),
  list(data = saf_nonproc_bleed,  name = "NonProc Bleed"),
  list(data = saf_all_major_bleed, name = "All Major Bleed"),
  list(data = sec_acm,            name = "All-Cause Mortality")
)) {
  m_rr <- run_rr(ep_info$data, flx_studies,
                  sprintf("%s RR FLX", ep_info$name))
  m_hr <- run_hr(ep_info$data, flx_studies,
                  sprintf("%s HR FLX", ep_info$name))

  if (!is.null(m_rr)) {
    metrics <- extract_metrics(m_rr, "RE")
    metrics$endpoint <- ep_info$name
    metrics$estimand <- "RR"
    flx_results[[paste0(ep_info$name, "_RR")]] <- metrics
    cat(sprintf("  %-20s RR = %.2f (%.2f-%.2f)\n",
                ep_info$name, metrics$estimate, metrics$ci_lo, metrics$ci_hi))
  }
  if (!is.null(m_hr)) {
    metrics <- extract_metrics(m_hr, "RE")
    metrics$endpoint <- ep_info$name
    metrics$estimand <- "HR"
    flx_results[[paste0(ep_info$name, "_HR")]] <- metrics
    cat(sprintf("  %-20s HR = %.2f (%.2f-%.2f)\n",
                ep_info$name, metrics$estimate, metrics$ci_lo, metrics$ci_hi))
  }
}

flx_df <- bind_rows(flx_results)

# =============================================================================
# 7. CONTEMPORARY-ONLY (exclude CLOSURE and PP)
# =============================================================================

cat("\n--- 7. CONTEMPORARY-ONLY (CHAMPION + OPTION + PRAGUE) ---\n\n")

contemp_studies <- get_pool_studies("contemporary")
contemp_results <- list()

for (ep_info in list(
  list(data = eff_is_se,          name = "IS+SE"),
  list(data = saf_nonproc_bleed,  name = "NonProc Bleed"),
  list(data = saf_all_major_bleed, name = "All Major Bleed"),
  list(data = sec_acm,            name = "All-Cause Mortality")
)) {
  m_rr <- run_rr(ep_info$data, contemp_studies,
                  sprintf("%s RR Contemporary", ep_info$name))
  m_hr <- run_hr(ep_info$data, contemp_studies,
                  sprintf("%s HR Contemporary", ep_info$name))

  if (!is.null(m_rr)) {
    metrics <- extract_metrics(m_rr, "RE")
    metrics$endpoint <- ep_info$name
    metrics$estimand <- "RR"
    contemp_results[[paste0(ep_info$name, "_RR")]] <- metrics
    cat(sprintf("  %-20s RR = %.2f (%.2f-%.2f), PI (%.2f-%.2f), I2 = %.0f%%\n",
                ep_info$name, metrics$estimate, metrics$ci_lo, metrics$ci_hi,
                metrics$pi_lo, metrics$pi_hi, metrics$I2 * 100))
  }
  if (!is.null(m_hr)) {
    metrics <- extract_metrics(m_hr, "RE")
    metrics$endpoint <- ep_info$name
    metrics$estimand <- "HR"
    contemp_results[[paste0(ep_info$name, "_HR")]] <- metrics
    cat(sprintf("  %-20s HR = %.2f (%.2f-%.2f), PI (%.2f-%.2f), I2 = %.0f%%\n",
                ep_info$name, metrics$estimate, metrics$ci_lo, metrics$ci_hi,
                metrics$pi_lo, metrics$pi_hi, metrics$I2 * 100))
  }
}

contemp_df <- bind_rows(contemp_results)

# =============================================================================
# 8. HISTORICAL BENCHMARK (PROTECT-AF/PREVAIL alone — descriptive)
# =============================================================================

cat("\n--- 8. HISTORICAL BENCHMARK (PROTECT-AF/PREVAIL) ---\n\n")

hist_study <- "PROTECT-AF/PREVAIL"
hist_results <- list()

for (ep_info in list(
  list(data = eff_is_se,          name = "IS+SE"),
  list(data = saf_nonproc_bleed,  name = "NonProc Bleed"),
  list(data = sec_acm,            name = "All-Cause Mortality"),
  list(data = sec_hemstroke,      name = "Hemorrhagic Stroke")
)) {
  d <- ep_info$data %>% filter(study == hist_study)
  if (nrow(d) == 1) {
    rr_val <- (d$ev_device / d$n_device) / (d$ev_control / d$n_control)
    rr_se  <- sqrt(1/d$ev_device - 1/d$n_device + 1/d$ev_control - 1/d$n_control)
    rr_lo  <- exp(log(rr_val) - 1.96 * rr_se)
    rr_hi  <- exp(log(rr_val) + 1.96 * rr_se)

    hist_results[[ep_info$name]] <- tibble(
      endpoint = ep_info$name,
      RR       = rr_val,
      RR_lo    = rr_lo,
      RR_hi    = rr_hi,
      HR       = d$hr,
      HR_lo    = d$hr_lo,
      HR_hi    = d$hr_hi
    )

    cat(sprintf("  %-20s: RR = %.2f (%.2f-%.2f) | HR = %.2f (%.2f-%.2f)\n",
                ep_info$name, rr_val, rr_lo, rr_hi,
                d$hr, d$hr_lo, d$hr_hi))
  }
}

hist_df <- bind_rows(hist_results)
cat("\n  Note: PP results reflect warfarin comparator and first-gen Watchman.\n")
cat("  These serve as historical context, not direct comparison.\n")

# =============================================================================
# 9. BENEFIT-RISK ASSESSMENT (per 1000 patients over 3 years)
# =============================================================================

cat("\n--- 9. BENEFIT-RISK: Per 1000 Patients / 3 Years ---\n\n")

# Use expanded pool RR estimates for the calculation
# We re-run expanded pool analyses to extract point estimates

m_eff     <- run_rr(eff_is_se, expanded_studies, "IS+SE — benefit-risk")
m_bleed   <- run_rr(saf_nonproc_bleed, expanded_studies, "NonProc Bleed — benefit-risk")
m_hem     <- run_rr(sec_hemstroke, expanded_studies, "Hem Stroke — benefit-risk")
m_acm     <- run_rr(sec_acm, expanded_studies, "ACM — benefit-risk")

# Control arm event rates (3-year, from expanded pool, weighted average)
calc_ctrl_rate <- function(data, studies) {
  d <- data %>% filter(study %in% studies)
  sum(d$ev_control) / sum(d$n_control)
}

ctrl_is_se   <- calc_ctrl_rate(eff_is_se, expanded_studies)
ctrl_bleed   <- calc_ctrl_rate(saf_nonproc_bleed, expanded_studies)
ctrl_hem     <- calc_ctrl_rate(sec_hemstroke, expanded_studies)
ctrl_acm     <- calc_ctrl_rate(sec_acm, expanded_studies)

# Procedural complication rate (from contemporary FLX-era)
proc_flx_d   <- proc_complications %>% filter(study %in% c("CHAMPION-AF", "OPTION"))
proc_rate    <- sum(proc_flx_d$events) / sum(proc_flx_d$n_attempted)

# Per 1000 patients over 3 years
N <- 1000

# Device arm expected events = control_rate * RR * N
# Difference = device_events - control_events = control_rate * (RR - 1) * N

# For procedural: one-time at implant
proc_events <- proc_rate * N

get_rr <- function(m) if (!is.null(m)) exp(m$TE.random) else NA_real_

benefit_risk <- tibble::tribble(
  ~component,                    ~direction,  ~ctrl_rate_3yr, ~pooled_RR,
  ~ctrl_per_1000, ~device_per_1000, ~difference,

  "Ischemic stroke + SE",       "harm",
    ctrl_is_se, get_rr(m_eff),
    ctrl_is_se * N, ctrl_is_se * get_rr(m_eff) * N,
    ctrl_is_se * (get_rr(m_eff) - 1) * N,

  "Non-proc bleeding",          "benefit",
    ctrl_bleed, get_rr(m_bleed),
    ctrl_bleed * N, ctrl_bleed * get_rr(m_bleed) * N,
    ctrl_bleed * (1 - get_rr(m_bleed)) * N,

  "Hemorrhagic stroke",         "benefit",
    ctrl_hem, get_rr(m_hem),
    ctrl_hem * N, ctrl_hem * get_rr(m_hem) * N,
    ctrl_hem * (1 - get_rr(m_hem)) * N,

  "Procedural complications",   "harm_onetime",
    proc_rate, NA_real_,
    0, proc_events,
    proc_events,

  "All-cause mortality",        "neutral",
    ctrl_acm, get_rr(m_acm),
    ctrl_acm * N, ctrl_acm * get_rr(m_acm) * N,
    ctrl_acm * (get_rr(m_acm) - 1) * N
)

cat("  Component                   | Ctrl/1000 | Device/1000 | Diff/1000 | Direction\n")
cat("  -------------------------------------------------------------------------------\n")
for (i in 1:nrow(benefit_risk)) {
  r <- benefit_risk[i, ]
  cat(sprintf("  %-29s | %7.1f   | %7.1f     | %+7.1f   | %s\n",
              r$component, r$ctrl_per_1000, r$device_per_1000,
              r$difference, r$direction))
}

cat("\n  Interpretation:\n")
cat("    (+) difference = more events with device (harm of device)\n")
cat("    (-) difference = fewer events with device (benefit of device)\n")
cat("    Procedural complications are one-time at implant.\n")

# Net calculation
net_is_extra     <- benefit_risk %>% filter(component == "Ischemic stroke + SE") %>% pull(difference)
net_bleed_saved  <- benefit_risk %>% filter(component == "Non-proc bleeding") %>% pull(difference)
net_hem_saved    <- benefit_risk %>% filter(component == "Hemorrhagic stroke") %>% pull(difference)
net_proc         <- benefit_risk %>% filter(component == "Procedural complications") %>% pull(difference)
net_death_diff   <- benefit_risk %>% filter(component == "All-cause mortality") %>% pull(difference)

cat(sprintf("\n  NET SUMMARY (per 1000 patients, 3 years):\n"))
cat(sprintf("    Extra IS+SE with device:        %+.1f\n", net_is_extra))
cat(sprintf("    Bleeds prevented by device:     %+.1f (negative = prevented)\n", -net_bleed_saved))
cat(sprintf("    Hem strokes prevented:          %+.1f (negative = prevented)\n", -net_hem_saved))
cat(sprintf("    Procedural complications:       %+.1f (one-time cost)\n", net_proc))
cat(sprintf("    Death difference:               %+.1f\n", net_death_diff))

# =============================================================================
# 10. SAVE ALL SENSITIVITY RESULTS
# =============================================================================

save(
  # Leave-one-out
  loo_results,
  # FE vs RE
  fe_vs_re,
  # Tau estimator comparison
  tau_comp_df,
  # HKSJ vs DL
  hksj_dl_df,
  # Tier 1-2
  tier12_df,
  # FLX only
  flx_df,
  # Contemporary only
  contemp_df,
  # Historical benchmark
  hist_df,
  # Benefit-risk
  benefit_risk,
  file = file.path(DIR_O, "sensitivity_results.RData")
)

cat(sprintf("\n>> Results saved to %s\n", file.path(DIR_O, "sensitivity_results.RData")))
cat(">> 04_sensitivity.R complete.\n")
