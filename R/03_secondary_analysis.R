###############################################################################
# 03_secondary_analysis.R
# LAA Closure vs OAC Meta-Analysis — Secondary Endpoints, Landmark,
# Procedural Complications, CLOSURE-AF Cascade, NI Margin Heatmap
#
# Method: same two-stage RE (REML + HKSJ) as 02_primary_analysis.R
# Pool: Expanded (primary), with all-trial sensitivity where relevant
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "01_data.R"))

cat("\n========================================\n")
cat("  03_secondary_analysis.R\n")
cat("========================================\n\n")

# =============================================================================
# 0. HELPER FUNCTIONS (re-defined here for standalone use)
# =============================================================================

get_pool_studies <- function(pool_name) {
  switch(pool_name,
    "contemporary" = pool_assign %>% filter(pool_primary)   %>% pull(study),
    "expanded"     = pool_assign %>% filter(pool_expanded)  %>% pull(study),
    "all"          = pool_assign$study,
    stop("Unknown pool: ", pool_name)
  )
}

run_meta_rr <- function(data, pool_filter, label) {
  cat(sprintf("  [RR] %s\n", label))
  studies <- get_pool_studies(pool_filter)
  d <- data %>% filter(study %in% studies)
  if (nrow(d) < 2) {
    cat("    >> Fewer than 2 studies; skipping.\n")
    return(NULL)
  }
  m <- metabin(
    event.e = ev_device, n.e = n_device,
    event.c = ev_control, n.c = n_control,
    studlab = study, data = d, sm = "RR",
    method = "MH", method.tau = "REML",
    hakn = TRUE, prediction = TRUE,
    incr = 0.5, allincr = FALSE, title = label
  )
  list(
    label = label, pool = pool_filter, estimand = "RR", k = m$k,
    RE_est = exp(m$TE.random), RE_lo = exp(m$lower.random),
    RE_hi = exp(m$upper.random), RE_pval = m$pval.random,
    PI_lo = exp(m$lower.predict), PI_hi = exp(m$upper.predict),
    tau2 = m$tau2, Q = m$Q, Q_pval = m$pval.Q, I2 = m$I2,
    FE_est = exp(m$TE.fixed), FE_lo = exp(m$lower.fixed),
    FE_hi = exp(m$upper.fixed), FE_pval = m$pval.fixed,
    meta_obj = m
  )
}

run_meta_hr <- function(data, pool_filter, label) {
  cat(sprintf("  [HR] %s\n", label))
  studies <- get_pool_studies(pool_filter)
  d <- data %>%
    filter(study %in% studies) %>%
    filter(!is.na(loghr) & !is.na(se_loghr))
  if (nrow(d) < 2) {
    cat("    >> Fewer than 2 studies with HR; skipping.\n")
    return(NULL)
  }
  m <- metagen(
    TE = loghr, seTE = se_loghr, studlab = study,
    data = d, sm = "HR", method.tau = "REML",
    hakn = TRUE, prediction = TRUE, title = label
  )
  list(
    label = label, pool = pool_filter, estimand = "HR", k = m$k,
    RE_est = exp(m$TE.random), RE_lo = exp(m$lower.random),
    RE_hi = exp(m$upper.random), RE_pval = m$pval.random,
    PI_lo = exp(m$lower.predict), PI_hi = exp(m$upper.predict),
    tau2 = m$tau2, Q = m$Q, Q_pval = m$pval.Q, I2 = m$I2,
    FE_est = exp(m$TE.fixed), FE_lo = exp(m$lower.fixed),
    FE_hi = exp(m$upper.fixed), FE_pval = m$pval.fixed,
    meta_obj = m
  )
}

# Compact printer
print_result <- function(r) {
  if (is.null(r)) return(invisible(NULL))
  cat(sprintf("    %s = %.2f (%.2f-%.2f), PI (%.2f-%.2f), tau2=%.4f, I2=%.0f%%, k=%d\n",
              r$estimand, r$RE_est, r$RE_lo, r$RE_hi,
              r$PI_lo, r$PI_hi, r$tau2, r$I2 * 100, r$k))
}

# =============================================================================
# 1. SECONDARY ENDPOINTS — Expanded Pool (RR + HR for each)
# =============================================================================

# Define endpoints to iterate over
sec_endpoints <- list(
  list(data = sec_acm,          name = "All-Cause Mortality"),
  list(data = sec_cvdeath,      name = "CV/Unexplained Death"),
  list(data = sec_hemstroke,    name = "Hemorrhagic Stroke"),
  list(data = sec_anystroke_se, name = "Any Stroke + SE"),
  list(data = sec_is_alone,     name = "Ischemic Stroke Alone"),
  list(data = trial_primary,    name = "Trial-Reported Primary Composite")
)

sec_results_rr <- list()
sec_results_hr <- list()

cat("--- SECONDARY ENDPOINTS (Expanded Pool) ---\n\n")

for (ep in sec_endpoints) {
  cat(sprintf(">> %s\n", ep$name))

  rr_label <- sprintf("%s RR — Expanded", ep$name)
  hr_label <- sprintf("%s HR — Expanded", ep$name)

  rr_res <- run_meta_rr(ep$data, "expanded", rr_label)
  hr_res <- run_meta_hr(ep$data, "expanded", hr_label)

  print_result(rr_res)
  print_result(hr_res)
  cat("\n")

  sec_results_rr[[ep$name]] <- rr_res
  sec_results_hr[[ep$name]] <- hr_res
}

# Also run all-trial pool for key secondary (ACM)
cat("--- ALL-CAUSE MORTALITY: All-Trial Pool ---\n")
sec_acm_rr_all <- run_meta_rr(sec_acm, "all", "All-Cause Mortality RR — All trials")
sec_acm_hr_all <- run_meta_hr(sec_acm, "all", "All-Cause Mortality HR — All trials")
print_result(sec_acm_rr_all)
print_result(sec_acm_hr_all)

# =============================================================================
# 2. LANDMARK ANALYSIS (CLOSURE-AF: 0-6 months vs >6 months)
# =============================================================================

cat("\n--- LANDMARK ANALYSIS: CLOSURE-AF ---\n")

# The landmark_closure data has period, ev_device, ev_control, hr, hr_lo, hr_hi
# This is a within-trial comparison, so we report descriptively
# and test interaction (ratio of HRs)

landmark_results <- landmark_closure %>%
  mutate(
    loghr    = log(hr),
    se_loghr = (log(hr_hi) - log(hr_lo)) / (2 * qnorm(0.975))
  )

cat("  Period         | Device ev | Control ev | HR (95% CI)\n")
cat("  --------------------------------------------------------\n")
for (i in 1:nrow(landmark_results)) {
  r <- landmark_results[i, ]
  cat(sprintf("  %-14s | %5d     | %5d      | %.2f (%.2f-%.2f)\n",
              r$period, r$ev_device, r$ev_control,
              r$hr, r$hr_lo, r$hr_hi))
}

# Interaction test: ratio of HRs across periods
# log(HR_early / HR_late) with pooled SE
lr <- landmark_results
log_ratio  <- lr$loghr[1] - lr$loghr[2]
se_ratio   <- sqrt(lr$se_loghr[1]^2 + lr$se_loghr[2]^2)
z_interact <- log_ratio / se_ratio
p_interact <- 2 * pnorm(-abs(z_interact))

landmark_interaction <- list(
  ratio_HR    = exp(log_ratio),
  ratio_lo    = exp(log_ratio - 1.96 * se_ratio),
  ratio_hi    = exp(log_ratio + 1.96 * se_ratio),
  z           = z_interact,
  p_value     = p_interact
)

cat(sprintf("\n  Interaction test (early vs late):\n"))
cat(sprintf("    HR ratio = %.2f (%.2f-%.2f), p = %.3f\n",
            landmark_interaction$ratio_HR,
            landmark_interaction$ratio_lo,
            landmark_interaction$ratio_hi,
            landmark_interaction$p_value))

# =============================================================================
# 3. PROCEDURAL COMPLICATIONS — Single-Arm Meta-Analysis (GLMM)
# =============================================================================

cat("\n--- PROCEDURAL COMPLICATIONS: Single-Arm Pooling (metaprop) ---\n\n")

# Overall major procedural complication rate
proc_overall <- metaprop(
  event  = proc_complications$events,
  n      = proc_complications$n_attempted,
  studlab = proc_complications$study,
  sm     = "PLOGIT",          # logit transform for proportions
  method.tau = "ML",
  method.random.ci = "HK",
  prediction = TRUE,
  title  = "Overall Procedural Major Complications"
)

cat(sprintf("  Overall complication rate (GLMM):\n"))
cat(sprintf("    Pooled = %.1f%% (%.1f%%–%.1f%%), PI (%.1f%%–%.1f%%)\n",
            proc_overall$TE.random * 100,  # Note: metaprop back-transforms
            proc_overall$lower.random * 100,
            proc_overall$upper.random * 100,
            proc_overall$lower.predict * 100,
            proc_overall$upper.predict * 100))
cat(sprintf("    tau2 = %.4f, I2 = %.0f%%, Q = %.2f (p = %.3f), k = %d\n",
            proc_overall$tau2, proc_overall$I2 * 100,
            proc_overall$Q, proc_overall$pval.Q, proc_overall$k))

# Pericardial effusion/tamponade rate
proc_pericard <- metaprop(
  event   = proc_complications$pericardial,
  n       = proc_complications$n_attempted,
  studlab = proc_complications$study,
  sm      = "PLOGIT",
  method.tau = "ML",
  method.random.ci = "HK",
  prediction = TRUE,
  title   = "Pericardial Effusion/Tamponade"
)

cat(sprintf("\n  Pericardial effusion/tamponade:\n"))
cat(sprintf("    Pooled = %.1f%% (%.1f%%–%.1f%%)\n",
            proc_pericard$TE.random * 100,
            proc_pericard$lower.random * 100,
            proc_pericard$upper.random * 100))

# Device embolization rate
proc_embol <- metaprop(
  event   = proc_complications$embolization,
  n       = proc_complications$n_attempted,
  studlab = proc_complications$study,
  sm      = "PLOGIT",
  method.tau = "ML",
  method.random.ci = "HK",
  prediction = TRUE,
  title   = "Device Embolization"
)

cat(sprintf("\n  Device embolization:\n"))
cat(sprintf("    Pooled = %.2f%% (%.2f%%–%.2f%%)\n",
            proc_embol$TE.random * 100,
            proc_embol$lower.random * 100,
            proc_embol$upper.random * 100))

# Procedural death rate
proc_death <- metaprop(
  event   = proc_complications$deaths,
  n       = proc_complications$n_attempted,
  studlab = proc_complications$study,
  sm      = "PLOGIT",
  method.tau = "ML",
  method.random.ci = "HK",
  prediction = TRUE,
  title   = "Procedural Death"
)

cat(sprintf("\n  Procedural death:\n"))
cat(sprintf("    Pooled = %.2f%% (%.2f%%–%.2f%%)\n",
            proc_death$TE.random * 100,
            proc_death$lower.random * 100,
            proc_death$upper.random * 100))

proc_results <- list(
  overall     = proc_overall,
  pericardial = proc_pericard,
  embolization = proc_embol,
  death       = proc_death
)

# Contemporary-only procedural rate (FLX-era: CHAMPION + OPTION)
cat("\n  Contemporary FLX-era (CHAMPION + OPTION) complication rate:\n")
proc_flx <- proc_complications %>%
  filter(study %in% c("CHAMPION-AF", "OPTION"))

proc_flx_meta <- metaprop(
  event   = proc_flx$events,
  n       = proc_flx$n_attempted,
  studlab = proc_flx$study,
  sm      = "PLOGIT",
  method.tau = "ML",
  method.random.ci = "HK",
  title   = "FLX-era Procedural Complications"
)

cat(sprintf("    Pooled = %.1f%% (%.1f%%–%.1f%%), k = %d\n",
            proc_flx_meta$TE.random * 100,
            proc_flx_meta$lower.random * 100,
            proc_flx_meta$upper.random * 100,
            proc_flx_meta$k))

proc_results$flx_era <- proc_flx_meta

# =============================================================================
# 4. CLOSURE-AF ENDPOINT CASCADE ANALYSIS
# =============================================================================

cat("\n--- CLOSURE-AF ENDPOINT CASCADE ---\n\n")

# For each cascade model, compute RR from events and assess direction
# Total N: device=446, control=442
cascade_analysis <- closure_cascade %>%
  mutate(
    n_device  = 446,
    n_control = 442,
    rr        = (ev_device / n_device) / (ev_ctrl / n_control),
    rr_se     = sqrt(1/ev_device - 1/n_device + 1/ev_ctrl - 1/n_control),
    rr_lo     = exp(log(rr) - 1.96 * rr_se),
    rr_hi     = exp(log(rr) + 1.96 * rr_se),
    direction = ifelse(rr > 1, "Favors OAC", "Favors Device")
  )

cat("  Model                          | ev_dev | ev_ctrl | RR (95% CI)       | HR(mc) | dRMST  | Direction\n")
cat("  ---------------------------------------------------------------------------------------------------------\n")
for (i in 1:nrow(cascade_analysis)) {
  r <- cascade_analysis[i, ]
  cat(sprintf("  %-32s | %5d  | %5d   | %.2f (%.2f-%.2f) | %.2f   | %+.3f | %s\n",
              r$model, r$ev_device, r$ev_ctrl,
              r$rr, r$rr_lo, r$rr_hi,
              r$hr_mc, r$drmst,
              r$direction))
}

cat("\n  Key insight: As the composite is disaggregated from NACE to\n")
cat("  thromboembolic events alone, the device-control difference attenuates.\n")

# =============================================================================
# 5. NI MARGIN HEATMAP DATA PREPARATION
# =============================================================================

cat("\n--- NI MARGIN RE-EXPRESSION HEATMAP ---\n\n")

# Reshape ni_margin_data for heatmap plotting
ni_heatmap_long <- ni_margin_data %>%
  pivot_longer(
    cols      = c(rr_115, rr_125, rr_140, abs_05, abs_10, abs_20),
    names_to  = "threshold",
    values_to = "compatible"
  ) %>%
  mutate(
    # Human-readable labels
    margin_type = case_when(
      str_starts(threshold, "rr")  ~ "Relative (RR)",
      str_starts(threshold, "abs") ~ "Absolute (ppt)"
    ),
    margin_label = case_when(
      threshold == "rr_115" ~ "RR < 1.15",
      threshold == "rr_125" ~ "RR < 1.25",
      threshold == "rr_140" ~ "RR < 1.40",
      threshold == "abs_05" ~ "< 0.5 ppt",
      threshold == "abs_10" ~ "< 1.0 ppt",
      threshold == "abs_20" ~ "< 2.0 ppt"
    ),
    # Order factors for plotting
    margin_label = factor(margin_label,
                          levels = c("RR < 1.15", "RR < 1.25", "RR < 1.40",
                                     "< 0.5 ppt", "< 1.0 ppt", "< 2.0 ppt")),
    study = factor(study,
                   levels = c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                              "CLOSURE-AF", "PROTECT-AF/PREVAIL"))
  )

# Summary: how many trials compatible at each threshold?
ni_summary <- ni_heatmap_long %>%
  group_by(threshold, margin_label) %>%
  summarise(
    n_compatible    = sum(compatible),
    n_total         = n(),
    pct_compatible  = n_compatible / n_total * 100,
    .groups = "drop"
  )

cat("  NI Threshold    | Compatible / Total | Pct\n")
cat("  -------------------------------------------------\n")
for (i in 1:nrow(ni_summary)) {
  r <- ni_summary[i, ]
  cat(sprintf("  %-17s | %d / %d              | %.0f%%\n",
              as.character(r$margin_label), r$n_compatible, r$n_total, r$pct_compatible))
}

# =============================================================================
# 6. COMPILE SECONDARY SUMMARY TABLE
# =============================================================================

cat("\n--- SECONDARY ANALYSIS SUMMARY ---\n\n")

all_sec_results <- c(sec_results_rr, sec_results_hr)
all_sec_results <- all_sec_results[!sapply(all_sec_results, is.null)]

secondary_summary <- bind_rows(lapply(all_sec_results, function(r) {
  tibble(
    Endpoint  = r$label,
    Pool      = r$pool,
    Estimand  = r$estimand,
    k         = r$k,
    Estimate  = r$RE_est,
    CI_lo     = r$RE_lo,
    CI_hi     = r$RE_hi,
    P_value   = r$RE_pval,
    PI_lo     = r$PI_lo,
    PI_hi     = r$PI_hi,
    tau2      = r$tau2,
    I2_pct    = r$I2 * 100,
    Q         = r$Q,
    Q_pval    = r$Q_pval,
    FE_est    = r$FE_est,
    FE_lo     = r$FE_lo,
    FE_hi     = r$FE_hi
  )
}))

cat("  Endpoint                             | Est    | 95% CI          | PI              | I2\n")
cat("  ---------------------------------------------------------------------------------------\n")
for (i in 1:nrow(secondary_summary)) {
  r <- secondary_summary[i, ]
  cat(sprintf("  %-38s | %5.2f  | (%.2f-%.2f)  | (%.2f-%.2f)  | %.0f%%\n",
              r$Endpoint, r$Estimate, r$CI_lo, r$CI_hi,
              r$PI_lo, r$PI_hi, r$I2_pct))
}

# =============================================================================
# 7. SAVE ALL RESULTS
# =============================================================================

save(
  # Secondary endpoint results
  sec_results_rr, sec_results_hr,
  secondary_summary,
  # All-trial ACM
  sec_acm_rr_all, sec_acm_hr_all,
  # Landmark
  landmark_results, landmark_interaction,
  # Procedural
  proc_results,
  # CLOSURE cascade
  cascade_analysis,
  # NI margin
  ni_heatmap_long, ni_summary,
  file = file.path(DIR_O, "secondary_results.RData")
)

cat(sprintf("\n>> Results saved to %s\n", file.path(DIR_O, "secondary_results.RData")))
cat(">> 03_secondary_analysis.R complete.\n")
