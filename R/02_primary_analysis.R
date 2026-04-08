###############################################################################
# 02_primary_analysis.R
# LAA Closure vs OAC Meta-Analysis — Primary & Key Safety Analyses
#
# Method: Two-stage RE meta-analysis
#   - Primary estimand: RR (Mantel-Haenszel for rare events, REML + HKSJ for RE)
#   - Secondary estimand: HR (generic inverse-variance via metagen, REML + HKSJ)
#   - Fixed-effect as built-in sensitivity
#   - Prediction interval mandatory
#
# Pool hierarchy:
#   Contemporary = CHAMPION + OPTION + PRAGUE-17
#   Expanded     = Contemporary + CLOSURE-AF
#   All          = Expanded + PROTECT-AF/PREVAIL (historical)
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "01_data.R"))

cat("\n========================================\n")
cat("  02_primary_analysis.R\n")
cat("========================================\n\n")

# =============================================================================
# 1. HELPER: get studies by pool
# =============================================================================

get_pool_studies <- function(pool_name) {
  switch(pool_name,
    "contemporary" = pool_assign %>% filter(pool_primary)   %>% pull(study),
    "expanded"     = pool_assign %>% filter(pool_expanded)  %>% pull(study),
    "all"          = pool_assign$study,
    stop("Unknown pool: ", pool_name)
  )
}

# =============================================================================
# 2. CORE META-ANALYSIS FUNCTIONS
# =============================================================================

# ---------------------------------------------------------------------------
# run_meta_rr: RR-based meta-analysis from 2x2 event counts
#   Uses metabin() with Mantel-Haenszel (appropriate for rare events)
#   RE model with REML tau-squared and Hartung-Knapp-Sidik-Jonkman correction
# ---------------------------------------------------------------------------

run_meta_rr <- function(data, pool_filter, label) {

  cat(sprintf("  [RR] %s\n", label))

  studies <- get_pool_studies(pool_filter)
  d <- data %>% filter(study %in% studies)

  if (nrow(d) < 2) {
    cat("    >> Fewer than 2 studies; skipping RE.\n")
    return(NULL)
  }

  # Mantel-Haenszel fixed-effect + REML random-effects with HKSJ

  m <- metabin(
    event.e  = ev_device,
    n.e      = n_device,
    event.c  = ev_control,
    n.c      = n_control,
    studlab  = study,
    data     = d,
    sm       = "RR",
    method   = "MH",                     # Mantel-Haenszel for fixed
    method.tau = "REML",                 # REML tau-squared estimator
    hakn     = TRUE,                     # HKSJ correction for RE CI
    prediction = TRUE,                   # prediction interval
    incr     = 0.5,                      # continuity correction for zero cells
    allincr  = FALSE,                    # only add 0.5 where zero cells exist
    title    = label
  )

  # Build tidy result list
  res <- list(
    label       = label,
    pool        = pool_filter,
    estimand    = "RR",
    k           = m$k,
    # Random-effects (primary)
    RE_est      = exp(m$TE.random),
    RE_lo       = exp(m$lower.random),
    RE_hi       = exp(m$upper.random),
    RE_pval     = m$pval.random,
    # Prediction interval
    PI_lo       = exp(m$lower.predict),
    PI_hi       = exp(m$upper.predict),
    # Heterogeneity
    tau2        = m$tau2,
    Q           = m$Q,
    Q_pval      = m$pval.Q,
    I2          = m$I2,
    # Fixed-effect (sensitivity)
    FE_est      = exp(m$TE.fixed),
    FE_lo       = exp(m$lower.fixed),
    FE_hi       = exp(m$upper.fixed),
    FE_pval     = m$pval.fixed,
    # Full meta object
    meta_obj    = m
  )

  cat(sprintf("    RE: RR = %.2f (%.2f-%.2f), PI = (%.2f-%.2f), tau2 = %.4f, I2 = %.0f%%, k = %d\n",
              res$RE_est, res$RE_lo, res$RE_hi,
              res$PI_lo, res$PI_hi,
              res$tau2, res$I2 * 100, res$k))

  return(res)
}

# ---------------------------------------------------------------------------
# run_meta_hr: HR-based meta-analysis from log(HR) and SE
#   Uses metagen() with inverse-variance weighting
#   CAUTION: never naively pool HR and sHR (noted in data tier column)
# ---------------------------------------------------------------------------

run_meta_hr <- function(data, pool_filter, label) {

  cat(sprintf("  [HR] %s\n", label))

  studies <- get_pool_studies(pool_filter)
  d <- data %>%
    filter(study %in% studies) %>%
    filter(!is.na(loghr) & !is.na(se_loghr))

  if (nrow(d) < 2) {
    cat("    >> Fewer than 2 studies with HR data; skipping RE.\n")
    return(NULL)
  }

  m <- metagen(
    TE       = loghr,
    seTE     = se_loghr,
    studlab  = study,
    data     = d,
    sm       = "HR",
    method.tau = "REML",
    hakn     = TRUE,
    prediction = TRUE,
    title    = label
  )

  res <- list(
    label       = label,
    pool        = pool_filter,
    estimand    = "HR",
    k           = m$k,
    RE_est      = exp(m$TE.random),
    RE_lo       = exp(m$lower.random),
    RE_hi       = exp(m$upper.random),
    RE_pval     = m$pval.random,
    PI_lo       = exp(m$lower.predict),
    PI_hi       = exp(m$upper.predict),
    tau2        = m$tau2,
    Q           = m$Q,
    Q_pval      = m$pval.Q,
    I2          = m$I2,
    FE_est      = exp(m$TE.fixed),
    FE_lo       = exp(m$lower.fixed),
    FE_hi       = exp(m$upper.fixed),
    FE_pval     = m$pval.fixed,
    meta_obj    = m
  )

  cat(sprintf("    RE: HR = %.2f (%.2f-%.2f), PI = (%.2f-%.2f), tau2 = %.4f, I2 = %.0f%%, k = %d\n",
              res$RE_est, res$RE_lo, res$RE_hi,
              res$PI_lo, res$PI_hi,
              res$tau2, res$I2 * 100, res$k))

  return(res)
}

# =============================================================================
# 3. PRIMARY EFFICACY: Ischemic Stroke + Systemic Embolism
# =============================================================================

cat("\n--- PRIMARY EFFICACY: Ischemic Stroke + SE ---\n")

eff_rr_contemp   <- run_meta_rr(eff_is_se, "contemporary", "IS+SE RR — Contemporary")
eff_rr_expanded  <- run_meta_rr(eff_is_se, "expanded",     "IS+SE RR — Expanded")
eff_rr_all       <- run_meta_rr(eff_is_se, "all",          "IS+SE RR — All trials")

eff_hr_contemp   <- run_meta_hr(eff_is_se, "contemporary", "IS+SE HR — Contemporary")
eff_hr_expanded  <- run_meta_hr(eff_is_se, "expanded",     "IS+SE HR — Expanded")
eff_hr_all       <- run_meta_hr(eff_is_se, "all",          "IS+SE HR — All trials")

# =============================================================================
# 4. KEY SAFETY: Non-Procedural Bleeding
# =============================================================================

cat("\n--- KEY SAFETY: Non-Procedural Clinically Relevant Bleeding ---\n")

saf_np_rr_contemp   <- run_meta_rr(saf_nonproc_bleed, "contemporary", "NonProc Bleed RR — Contemporary")
saf_np_rr_expanded  <- run_meta_rr(saf_nonproc_bleed, "expanded",     "NonProc Bleed RR — Expanded")
saf_np_rr_all       <- run_meta_rr(saf_nonproc_bleed, "all",          "NonProc Bleed RR — All trials")

saf_np_hr_contemp   <- run_meta_hr(saf_nonproc_bleed, "contemporary", "NonProc Bleed HR — Contemporary")
saf_np_hr_expanded  <- run_meta_hr(saf_nonproc_bleed, "expanded",     "NonProc Bleed HR — Expanded")
saf_np_hr_all       <- run_meta_hr(saf_nonproc_bleed, "all",          "NonProc Bleed HR — All trials")

# =============================================================================
# 5. SUPPORTIVE SAFETY: All Major Bleeding (including procedural)
# =============================================================================

cat("\n--- SUPPORTIVE SAFETY: All Major Bleeding ---\n")

saf_all_rr_contemp   <- run_meta_rr(saf_all_major_bleed, "contemporary", "All Major Bleed RR — Contemporary")
saf_all_rr_expanded  <- run_meta_rr(saf_all_major_bleed, "expanded",     "All Major Bleed RR — Expanded")
saf_all_rr_all       <- run_meta_rr(saf_all_major_bleed, "all",          "All Major Bleed RR — All trials")

saf_all_hr_contemp   <- run_meta_hr(saf_all_major_bleed, "contemporary", "All Major Bleed HR — Contemporary")
saf_all_hr_expanded  <- run_meta_hr(saf_all_major_bleed, "expanded",     "All Major Bleed HR — Expanded")
saf_all_hr_all       <- run_meta_hr(saf_all_major_bleed, "all",          "All Major Bleed HR — All trials")

# =============================================================================
# 6. COMPILE SUMMARY TABLE
# =============================================================================

cat("\n--- Building summary table ---\n")

# Collect all non-NULL results into a list
all_primary_results <- list(
  # Efficacy RR
  eff_rr_contemp, eff_rr_expanded, eff_rr_all,
  # Efficacy HR
  eff_hr_contemp, eff_hr_expanded, eff_hr_all,
  # NonProc Bleed RR
  saf_np_rr_contemp, saf_np_rr_expanded, saf_np_rr_all,
  # NonProc Bleed HR
  saf_np_hr_contemp, saf_np_hr_expanded, saf_np_hr_all,
  # All Major Bleed RR
  saf_all_rr_contemp, saf_all_rr_expanded, saf_all_rr_all,
  # All Major Bleed HR
  saf_all_hr_contemp, saf_all_hr_expanded, saf_all_hr_all
)

# Remove NULLs
all_primary_results <- all_primary_results[!sapply(all_primary_results, is.null)]

# Build tidy summary data.frame
primary_summary <- bind_rows(lapply(all_primary_results, function(r) {
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

# =============================================================================
# 7. PRINT CLEAN SUMMARY
# =============================================================================

cat("\n")
cat("==========================================================================\n")
cat("  PRIMARY ANALYSIS SUMMARY TABLE (Random-Effects, REML + HKSJ)\n")
cat("==========================================================================\n\n")

print_row <- function(r) {
  cat(sprintf("  %-42s | %s | k=%d | %s = %5.2f (%5.2f-%5.2f) | PI (%5.2f-%5.2f) | tau2=%.4f | I2=%4.1f%%\n",
              r$label, r$pool, r$k, r$estimand,
              r$RE_est, r$RE_lo, r$RE_hi,
              r$PI_lo, r$PI_hi,
              r$tau2, r$I2 * 100))
}

cat("--- ISCHEMIC STROKE + SE ---\n")
for (r in list(eff_rr_contemp, eff_rr_expanded, eff_rr_all,
               eff_hr_contemp, eff_hr_expanded, eff_hr_all)) {
  if (!is.null(r)) print_row(r)
}

cat("\n--- NON-PROCEDURAL BLEEDING ---\n")
for (r in list(saf_np_rr_contemp, saf_np_rr_expanded, saf_np_rr_all,
               saf_np_hr_contemp, saf_np_hr_expanded, saf_np_hr_all)) {
  if (!is.null(r)) print_row(r)
}

cat("\n--- ALL MAJOR BLEEDING ---\n")
for (r in list(saf_all_rr_contemp, saf_all_rr_expanded, saf_all_rr_all,
               saf_all_hr_contemp, saf_all_hr_expanded, saf_all_hr_all)) {
  if (!is.null(r)) print_row(r)
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

save(
  # Result lists
  all_primary_results, primary_summary,
  # Individual RR results - Efficacy
  eff_rr_contemp, eff_rr_expanded, eff_rr_all,
  # Individual HR results - Efficacy
  eff_hr_contemp, eff_hr_expanded, eff_hr_all,
  # Individual RR results - NonProc Bleed
  saf_np_rr_contemp, saf_np_rr_expanded, saf_np_rr_all,
  # Individual HR results - NonProc Bleed
  saf_np_hr_contemp, saf_np_hr_expanded, saf_np_hr_all,
  # Individual RR results - All Major Bleed
  saf_all_rr_contemp, saf_all_rr_expanded, saf_all_rr_all,
  # Individual HR results - All Major Bleed
  saf_all_hr_contemp, saf_all_hr_expanded, saf_all_hr_all,
  # Helper functions (for use in downstream scripts)
  file = file.path(DIR_O, "primary_results.RData")
)

cat(sprintf("\n>> Results saved to %s\n", file.path(DIR_O, "primary_results.RData")))
cat(">> 02_primary_analysis.R complete.\n")
