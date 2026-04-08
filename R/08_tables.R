###############################################################################
# 08_tables.R
# LAA Closure vs OAC Meta-Analysis — Publication-Quality Tables
#
# Generates:
#   MAIN TABLES (5):  Table 1-5 as PNG (gt) + DOCX (flextable)
#   SUPPLEMENTARY (8): sTable 1-8 as PNG (gt) + DOCX (flextable)
#
# Dependencies: 01_data.R (sources 00_setup.R with gt, flextable, officer, meta)
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "01_data.R"))

cat("\n========================================\n")
cat("  08_tables.R\n")
cat("========================================\n\n")

# =============================================================================
# 0. GLOBAL TABLE HELPERS
# =============================================================================

# Colour palette
HDR_BLUE     <- "#2166AC"
HDR_TEXT     <- "white"
ROW_EVEN     <- "#F0F4FA"
ROW_ODD      <- "white"
GREEN_BG     <- "#D5E8D4"
YELLOW_BG    <- "#FFF2CC"
RED_BG       <- "#F8CECC"
GREEN_DARK   <- "#2E7D32"
YELLOW_DARK  <- "#F57F17"
RED_DARK     <- "#C62828"
BORDER_CLR   <- "#CCCCCC"

# Study order (canonical)
STUDY_ORDER <- c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                 "CLOSURE-AF", "PROTECT-AF/PREVAIL")

# Short study labels for column headers
STUDY_SHORT <- c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                 "CLOSURE-AF", "PP (5-yr)")

# ---- Pool definitions ----
get_pool_studies <- function(pool_name) {
  switch(pool_name,
    "contemporary" = pool_assign %>% filter(pool_primary)  %>% pull(study),
    "expanded"     = pool_assign %>% filter(pool_expanded) %>% pull(study),
    "all"          = pool_assign$study,
    "flx"          = pool_assign %>% filter(pool_flx)      %>% pull(study),
    stop("Unknown pool: ", pool_name)
  )
}

# ---- Meta-analysis runners (compact, return meta objects) ----
run_rr <- function(data, studies, label = "",
                   method.tau = "REML", hakn = TRUE) {
  d <- data %>% filter(study %in% studies)
  if (nrow(d) < 2) return(NULL)
  metabin(
    event.e = ev_device, n.e = n_device,
    event.c = ev_control, n.c = n_control,
    studlab = study, data = d, sm = "RR",
    method = "MH", method.tau = method.tau,
    hakn = hakn, prediction = TRUE,
    incr = 0.5, allincr = FALSE, title = label
  )
}

run_hr <- function(data, studies, label = "",
                   method.tau = "REML", hakn = TRUE) {
  d <- data %>%
    filter(study %in% studies) %>%
    filter(!is.na(loghr) & !is.na(se_loghr))
  if (nrow(d) < 2) return(NULL)
  metagen(
    TE = loghr, seTE = se_loghr,
    studlab = study, data = d, sm = "HR",
    method.tau = method.tau, hakn = hakn,
    prediction = TRUE, title = label
  )
}

# ---- Formatting helpers ----
fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (%.*f\u2013%.*f)", digits, est, digits, lo, digits, hi)
}

fmt_mean_sd <- function(m, s, digits = 1) {
  if (is.na(m)) return("\u2014")
  if (is.na(s)) return(sprintf("%.*f", digits, m))
  sprintf("%.*f \u00b1 %.*f", digits, m, digits, s)
}

fmt_pct <- function(x, digits = 1) {
  if (is.na(x)) return("\u2014")
  sprintf("%.*f", digits, x)
}

# Check if CI excludes null (1.0 for ratios)
is_significant <- function(lo, hi, null = 1) {
  !is.na(lo) & !is.na(hi) & (lo > null | hi < null)
}

# ---- Save helpers ----
save_gt_png <- function(gt_obj, path, width = 1200, height = NULL) {
  tryCatch({
    gtsave(gt_obj, filename = path, vwidth = width)
    cat(sprintf("  >> Saved: %s\n", path))
  }, error = function(e) {
    cat(sprintf("  >> WARNING: Could not save PNG %s: %s\n", path, e$message))
  })
}

save_flextable_docx <- function(ft_obj, path) {
  tryCatch({
    doc <- read_docx()
    doc <- body_add_flextable(doc, ft_obj)
    print(doc, target = path)
    cat(sprintf("  >> Saved: %s\n", path))
  }, error = function(e) {
    cat(sprintf("  >> WARNING: Could not save DOCX %s: %s\n", path, e$message))
  })
}

# ---- Common gt theme ----
gt_theme_pub <- function(gt_obj, title = NULL, subtitle = NULL, footnote = NULL) {
  g <- gt_obj

  if (!is.null(title)) {
    g <- g %>%
      tab_header(title = title, subtitle = subtitle)
  }

  g <- g %>%
    tab_options(
      heading.title.font.size   = px(16),
      heading.subtitle.font.size = px(12),
      heading.background.color  = HDR_BLUE,
      column_labels.background.color = HDR_BLUE,
      column_labels.font.weight = "bold",
      column_labels.font.size   = px(11),
      table.font.size           = px(10),
      table.border.top.color    = HDR_BLUE,
      table.border.top.width    = px(2),
      table.border.bottom.color = HDR_BLUE,
      table.border.bottom.width = px(2),
      data_row.padding           = px(4),
      row_group.font.weight     = "bold",
      row_group.background.color = "#E8EEF4",
      source_notes.font.size    = px(8)
    ) %>%
    tab_style(
      style     = cell_text(color = HDR_TEXT),
      locations = cells_column_labels()
    ) %>%
    opt_row_striping()

  if (!is.null(footnote)) {
    g <- g %>% tab_source_note(source_note = footnote)
  }

  g
}

# ---- Common flextable theme ----
ft_theme_pub <- function(ft_obj) {
  ft_obj %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    bg(bg = HDR_BLUE, part = "header") %>%
    color(color = HDR_TEXT, part = "header") %>%
    align(align = "center", part = "header") %>%
    border_inner_h(border = fp_border(color = BORDER_CLR, width = 0.5)) %>%
    border_outer(border = fp_border(color = HDR_BLUE, width = 1.5)) %>%
    autofit(add_w = 0.1)
}


###############################################################################
#                          MAIN TABLES (1-5)
###############################################################################

# =============================================================================
# TABLE 1: Study Design Characteristics
# =============================================================================

cat("\n--- TABLE 1: Study Design Characteristics ---\n")

# Build primary endpoint definitions per trial
primary_ep_defs <- c(
  "CHAMPION-AF"       = "CV death + stroke + SE",
  "OPTION"            = "Death + stroke + SE",
  "PRAGUE-17"         = "Stroke/TIA + SE + CV death + bleeding + complications",
  "CLOSURE-AF"        = "Stroke + SE + major bleeding + CV/unexplained death",
  "PROTECT-AF/PREVAIL"= "Stroke + SE + CV death (5-yr pooled)"
)

t1_data <- study_design %>%
  mutate(
    Trial     = study,
    Year      = as.character(year),
    N         = format(n_total, big.mark = ","),
    Device    = device,
    Comparator = comparator,
    `Follow-up (yr)` = sprintf("%.1f", fu_years),
    `NI Margin` = ni_margin,
    `NI Met`    = ifelse(ni_met, "Yes", "No"),
    `Primary Endpoint` = primary_ep_defs[study]
  ) %>%
  select(Trial, Year, N, Device, Comparator, `Follow-up (yr)`,
         `NI Margin`, `NI Met`, `Primary Endpoint`)

# gt version
gt_t1 <- t1_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Table 1. Study Design Characteristics",
    subtitle = "LAA Closure vs OAC Randomized Controlled Trials",
    footnote = paste0(
      "NI = non-inferiority; ppt = percentage points; sHR = subdistribution hazard ratio; ",
      "SE = systemic embolism; CV = cardiovascular; PP = PROTECT-AF/PREVAIL pooled 5-year analysis. ",
      "NI Met indicates whether the trial met its pre-specified non-inferiority primary endpoint."
    )
  ) %>%
  cols_align(align = "center", columns = c(Year, N, `Follow-up (yr)`, `NI Margin`, `NI Met`)) %>%
  cols_align(align = "left", columns = c(Trial, Device, Comparator, `Primary Endpoint`)) %>%
  cols_width(
    Trial            ~ px(130),
    `Primary Endpoint` ~ px(250),
    Comparator       ~ px(130)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Trial)
  ) %>%
  tab_style(
    style     = list(cell_text(weight = "bold", color = GREEN_DARK)),
    locations = cells_body(columns = `NI Met`, rows = `NI Met` == "Yes")
  ) %>%
  tab_style(
    style     = list(cell_text(weight = "bold", color = RED_DARK)),
    locations = cells_body(columns = `NI Met`, rows = `NI Met` == "No")
  )

save_gt_png(gt_t1, file.path(DIR_T, "Table1_study_design.png"))

# flextable version
ft_t1 <- flextable(t1_data) %>%
  ft_theme_pub() %>%
  bold(j = "Trial") %>%
  align(j = c("Year", "N", "Follow-up (yr)", "NI Margin", "NI Met"), align = "center") %>%
  width(j = "Primary Endpoint", width = 2.5) %>%
  width(j = "Comparator", width = 1.3) %>%
  color(j = "NI Met", i = ~ `NI Met` == "Yes", color = GREEN_DARK) %>%
  color(j = "NI Met", i = ~ `NI Met` == "No", color = RED_DARK) %>%
  add_footer_lines(paste0(
    "NI = non-inferiority; ppt = percentage points; sHR = subdistribution hazard ratio; ",
    "SE = systemic embolism; CV = cardiovascular; PP = PROTECT-AF/PREVAIL pooled."
  )) %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_t1, file.path(DIR_T, "Table1_study_design.docx"))

# =============================================================================
# TABLE 2: Baseline Patient Characteristics
# =============================================================================

cat("\n--- TABLE 2: Baseline Patient Characteristics ---\n")

# Build transposed table: rows = characteristic, columns = studies
t2_chars <- c(
  "Age (mean +/- SD)",
  "Female (%)",
  "CHA2DS2-VASc (mean +/- SD)",
  "HAS-BLED (mean +/- SD)",
  "Prior stroke/TIA (%)",
  "Prior major bleeding (%)",
  "Severe CKD (%)",
  "Paroxysmal AF (%)"
)

t2_matrix <- matrix(NA_character_, nrow = length(t2_chars), ncol = 5)
colnames(t2_matrix) <- STUDY_SHORT

for (i in seq_along(STUDY_ORDER)) {
  s <- STUDY_ORDER[i]
  b <- baseline %>% filter(study == s)
  if (nrow(b) == 0) next

  t2_matrix[1, i] <- fmt_mean_sd(b$age_mean, b$age_sd)
  t2_matrix[2, i] <- fmt_pct(b$female_pct)
  t2_matrix[3, i] <- fmt_mean_sd(b$chadsvasc_mean, b$chadsvasc_sd)
  t2_matrix[4, i] <- fmt_mean_sd(b$hasbled_mean, b$hasbled_sd)
  t2_matrix[5, i] <- fmt_pct(b$prior_stroke_pct)
  t2_matrix[6, i] <- fmt_pct(b$prior_bleed_pct)
  t2_matrix[7, i] <- fmt_pct(b$ckd_severe_pct)
  t2_matrix[8, i] <- fmt_pct(b$paroxysmal_af_pct)
}

t2_data <- data.frame(
  Characteristic = t2_chars,
  t2_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# gt version
gt_t2 <- t2_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Table 2. Baseline Patient Characteristics",
    subtitle = "Values reported as mean +/- SD or percentage",
    footnote = paste0(
      "CHA2DS2-VASc = Congestive heart failure, Hypertension, Age >=75, Diabetes, ",
      "Stroke/TIA/TE, Vascular disease, Age 65-74, Sex category; ",
      "HAS-BLED = Hypertension, Abnormal renal/liver function, Stroke, Bleeding, ",
      "Labile INR, Elderly, Drugs/alcohol; CKD = chronic kidney disease; ",
      "AF = atrial fibrillation; TIA = transient ischemic attack. ",
      "\u2014 = not reported."
    )
  ) %>%
  cols_align(align = "center", columns = -Characteristic) %>%
  cols_align(align = "left", columns = Characteristic) %>%
  cols_width(Characteristic ~ px(200)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Characteristic)
  )

save_gt_png(gt_t2, file.path(DIR_T, "Table2_baseline.png"))

# flextable version
ft_t2 <- flextable(t2_data) %>%
  ft_theme_pub() %>%
  bold(j = "Characteristic") %>%
  align(j = STUDY_SHORT, align = "center") %>%
  width(j = "Characteristic", width = 2.0) %>%
  add_footer_lines(paste0(
    "CHA2DS2-VASc and HAS-BLED: standard risk scores. ",
    "CKD = chronic kidney disease; AF = atrial fibrillation. \u2014 = not reported."
  )) %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_t2, file.path(DIR_T, "Table2_baseline.docx"))

# =============================================================================
# TABLE 3: Main Pooled Results (KEY TABLE)
# =============================================================================

cat("\n--- TABLE 3: Main Pooled Results ---\n")

# Define all endpoint-pool combinations for Table 3
t3_specs <- list(
  # Ischemic stroke + SE
  list(data = eff_is_se,           ep = "Ischemic stroke + SE",   pool = "contemporary"),
  list(data = eff_is_se,           ep = "Ischemic stroke + SE",   pool = "expanded"),
  list(data = eff_is_se,           ep = "Ischemic stroke + SE",   pool = "all"),
  # Non-proc bleeding
  list(data = saf_nonproc_bleed,   ep = "Non-proc bleeding",      pool = "contemporary"),
  list(data = saf_nonproc_bleed,   ep = "Non-proc bleeding",      pool = "expanded"),
  list(data = saf_nonproc_bleed,   ep = "Non-proc bleeding",      pool = "all"),
  # All major bleeding
  list(data = saf_all_major_bleed, ep = "All major bleeding",     pool = "expanded"),
  # All-cause death
  list(data = sec_acm,             ep = "All-cause death",        pool = "expanded"),
  # CV/unexplained death
  list(data = sec_cvdeath,         ep = "CV/unexplained death",   pool = "expanded"),
  # Hemorrhagic stroke
  list(data = sec_hemstroke,       ep = "Hemorrhagic stroke",     pool = "expanded"),
  # Any stroke + SE
  list(data = sec_anystroke_se,    ep = "Any stroke + SE",        pool = "expanded")
)

# Run all analyses and build table rows
t3_rows <- lapply(t3_specs, function(spec) {
  studies <- get_pool_studies(spec$pool)
  pool_label <- tools::toTitleCase(spec$pool)

  # RR (Mantel-Haenszel)
  m_rr <- run_rr(spec$data, studies, paste(spec$ep, "RR", spec$pool))
  # HR (inverse-variance)
  m_hr <- run_hr(spec$data, studies, paste(spec$ep, "HR", spec$pool))

  # Extract RR results
  if (!is.null(m_rr)) {
    rr_est <- exp(m_rr$TE.random)
    rr_lo  <- exp(m_rr$lower.random)
    rr_hi  <- exp(m_rr$upper.random)
    rr_str <- fmt_ci(rr_est, rr_lo, rr_hi)
    pi_lo  <- exp(m_rr$lower.predict)
    pi_hi  <- exp(m_rr$upper.predict)
    pi_str <- sprintf("%.2f\u2013%.2f", pi_lo, pi_hi)
    i2_val <- m_rr$I2 * 100
    i2_str <- sprintf("%.0f%%", i2_val)
    k_val  <- m_rr$k
    rr_sig <- is_significant(rr_lo, rr_hi)
  } else {
    rr_str <- "\u2014"
    pi_str <- "\u2014"
    i2_str <- "\u2014"
    k_val  <- 0
    rr_sig <- FALSE
    rr_est <- NA; rr_lo <- NA; rr_hi <- NA
    i2_val <- NA
  }

  # Extract HR results
  if (!is.null(m_hr)) {
    hr_est <- exp(m_hr$TE.random)
    hr_lo  <- exp(m_hr$lower.random)
    hr_hi  <- exp(m_hr$upper.random)
    hr_str <- fmt_ci(hr_est, hr_lo, hr_hi)
    hr_sig <- is_significant(hr_lo, hr_hi)
    # Use HR for prediction interval if RR not available
    if (is.null(m_rr)) {
      pi_lo  <- exp(m_hr$lower.predict)
      pi_hi  <- exp(m_hr$upper.predict)
      pi_str <- sprintf("%.2f\u2013%.2f", pi_lo, pi_hi)
      i2_val <- m_hr$I2 * 100
      i2_str <- sprintf("%.0f%%", i2_val)
      k_val  <- m_hr$k
    }
  } else {
    hr_str <- "\u2014"
    hr_sig <- FALSE
    hr_est <- NA; hr_lo <- NA; hr_hi <- NA
  }

  # Determine direction
  ref_est <- if (!is.na(rr_est)) rr_est else hr_est
  favors <- if (is.na(ref_est)) "\u2014"
            else if (ref_est < 1) "Device"
            else if (ref_est > 1) "OAC"
            else "Neutral"

  tibble(
    Endpoint   = spec$ep,
    Pool       = pool_label,
    k          = k_val,
    RR_str     = rr_str,
    HR_str     = hr_str,
    PI_str     = pi_str,
    I2_str     = i2_str,
    Favors     = favors,
    rr_sig     = rr_sig,
    hr_sig     = hr_sig,
    rr_est     = rr_est,
    i2_val     = i2_val
  )
})

t3_data <- bind_rows(t3_rows)

# Display columns
t3_display <- t3_data %>%
  select(Endpoint, Pool, k,
         `RR (95% CI)` = RR_str,
         `HR (95% CI)` = HR_str,
         `Prediction Interval` = PI_str,
         I2     = I2_str,
         Favors)

# gt version
gt_t3 <- t3_display %>%
  gt() %>%
  gt_theme_pub(
    title    = "Table 3. Main Pooled Results",
    subtitle = "Random-effects meta-analysis (REML + HKSJ correction)",
    footnote = paste0(
      "RR = risk ratio (Mantel-Haenszel); HR = hazard ratio (inverse-variance); ",
      "PI = 95% prediction interval; I2 = heterogeneity statistic; ",
      "k = number of studies. Contemporary = CHAMPION-AF + OPTION + PRAGUE-17; ",
      "Expanded = Contemporary + CLOSURE-AF; All = Expanded + PROTECT-AF/PREVAIL. ",
      "Bold values indicate 95% CI excludes the null (1.0)."
    )
  ) %>%
  cols_align(align = "center", columns = -Endpoint) %>%
  cols_align(align = "left", columns = Endpoint) %>%
  cols_width(
    Endpoint             ~ px(160),
    `RR (95% CI)`        ~ px(150),
    `HR (95% CI)`        ~ px(150),
    `Prediction Interval` ~ px(100)
  )

# Bold significant RR values
sig_rr_rows <- which(t3_data$rr_sig)
if (length(sig_rr_rows) > 0) {
  gt_t3 <- gt_t3 %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_body(columns = `RR (95% CI)`, rows = sig_rr_rows)
    )
}

# Bold significant HR values
sig_hr_rows <- which(t3_data$hr_sig)
if (length(sig_hr_rows) > 0) {
  gt_t3 <- gt_t3 %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_body(columns = `HR (95% CI)`, rows = sig_hr_rows)
    )
}

# Color the "Favors" column
device_rows <- which(t3_data$Favors == "Device")
oac_rows    <- which(t3_data$Favors == "OAC")

if (length(device_rows) > 0) {
  gt_t3 <- gt_t3 %>%
    tab_style(
      style     = list(cell_fill(color = GREEN_BG), cell_text(color = GREEN_DARK, weight = "bold")),
      locations = cells_body(columns = Favors, rows = device_rows)
    )
}
if (length(oac_rows) > 0) {
  gt_t3 <- gt_t3 %>%
    tab_style(
      style     = list(cell_fill(color = RED_BG), cell_text(color = RED_DARK)),
      locations = cells_body(columns = Favors, rows = oac_rows)
    )
}

# Row grouping by endpoint
gt_t3 <- gt_t3 %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Endpoint)
  )

save_gt_png(gt_t3, file.path(DIR_T, "Table3_main_results.png"), width = 1400)

# flextable version
ft_t3_data <- t3_display
ft_t3 <- flextable(ft_t3_data) %>%
  ft_theme_pub() %>%
  bold(j = "Endpoint") %>%
  align(j = c("Pool", "k", "RR (95% CI)", "HR (95% CI)",
              "Prediction Interval", "I2", "Favors"), align = "center") %>%
  width(j = "Endpoint", width = 1.5) %>%
  width(j = c("RR (95% CI)", "HR (95% CI)"), width = 1.5)

# Bold significant entries
for (r in sig_rr_rows) {
  ft_t3 <- ft_t3 %>% bold(i = r, j = "RR (95% CI)")
}
for (r in sig_hr_rows) {
  ft_t3 <- ft_t3 %>% bold(i = r, j = "HR (95% CI)")
}

# Color Favors column
for (r in device_rows) {
  ft_t3 <- ft_t3 %>%
    bg(i = r, j = "Favors", bg = GREEN_BG) %>%
    color(i = r, j = "Favors", color = GREEN_DARK) %>%
    bold(i = r, j = "Favors")
}
for (r in oac_rows) {
  ft_t3 <- ft_t3 %>%
    bg(i = r, j = "Favors", bg = RED_BG) %>%
    color(i = r, j = "Favors", color = RED_DARK)
}

ft_t3 <- ft_t3 %>%
  add_footer_lines(paste0(
    "RR = risk ratio (MH); HR = hazard ratio (IV); PI = prediction interval; ",
    "Bold = 95% CI excludes null. REML + HKSJ correction applied."
  )) %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_t3, file.path(DIR_T, "Table3_main_results.docx"))

# =============================================================================
# TABLE 4: Sensitivity Analyses Summary
# =============================================================================

cat("\n--- TABLE 4: Sensitivity Analyses Summary ---\n")

# Run all sensitivity analyses inline
expanded_studies   <- get_pool_studies("expanded")
contemp_studies    <- get_pool_studies("contemporary")
flx_studies        <- get_pool_studies("flx")
all_studies        <- get_pool_studies("all")

# Helper: extract formatted RR string from metabin object
fmt_meta_rr <- function(m) {
  if (is.null(m)) return("\u2014")
  fmt_ci(exp(m$TE.random), exp(m$lower.random), exp(m$upper.random))
}

fmt_meta_hr <- function(m) {
  if (is.null(m)) return("\u2014")
  fmt_ci(exp(m$TE.random), exp(m$lower.random), exp(m$upper.random))
}

fmt_meta_fe_rr <- function(m) {
  if (is.null(m)) return("\u2014")
  fmt_ci(exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed))
}

fmt_meta_fe_hr <- function(m) {
  if (is.null(m)) return("\u2014")
  fmt_ci(exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed))
}

# (a) Primary RE (Expanded)
m_eff_rr_exp    <- run_rr(eff_is_se, expanded_studies, "IS+SE RR Expanded")
m_saf_hr_exp    <- run_hr(saf_nonproc_bleed, expanded_studies, "NonProc HR Expanded")

# (b) Fixed-effect (Expanded)
m_eff_rr_fe     <- run_rr(eff_is_se, expanded_studies, "IS+SE RR FE", hakn = FALSE)
m_saf_hr_fe     <- run_hr(saf_nonproc_bleed, expanded_studies, "NonProc HR FE", hakn = FALSE)

# (c) Leave-one-out (Expanded) - compute range
loo_eff_range <- c(Inf, -Inf)
loo_saf_range <- c(Inf, -Inf)
for (drop in expanded_studies) {
  remaining <- setdiff(expanded_studies, drop)
  m_loo_eff <- run_rr(eff_is_se, remaining, paste("LOO excl", drop))
  m_loo_saf <- run_hr(saf_nonproc_bleed, remaining, paste("LOO excl", drop))
  if (!is.null(m_loo_eff)) {
    val <- exp(m_loo_eff$TE.random)
    loo_eff_range[1] <- min(loo_eff_range[1], val)
    loo_eff_range[2] <- max(loo_eff_range[2], val)
  }
  if (!is.null(m_loo_saf)) {
    val <- exp(m_loo_saf$TE.random)
    loo_saf_range[1] <- min(loo_saf_range[1], val)
    loo_saf_range[2] <- max(loo_saf_range[2], val)
  }
}

# (d) FLX-only (CHAMPION + OPTION)
m_eff_rr_flx   <- run_rr(eff_is_se, flx_studies, "IS+SE RR FLX")
m_saf_hr_flx   <- run_hr(saf_nonproc_bleed, flx_studies, "NonProc HR FLX")

# (e) Contemporary only
m_eff_rr_cont  <- run_rr(eff_is_se, contemp_studies, "IS+SE RR Contemporary")
m_saf_hr_cont  <- run_hr(saf_nonproc_bleed, contemp_studies, "NonProc HR Contemporary")

# (f) Tier 1-2 only
tier12_eff_data <- eff_is_se %>% filter(study %in% expanded_studies, tier %in% c("T1", "T2"))
tier12_saf_data <- saf_nonproc_bleed %>% filter(study %in% expanded_studies, tier %in% c("T1", "T2"))
m_eff_rr_t12   <- run_rr(tier12_eff_data, tier12_eff_data$study, "IS+SE RR T12")
m_saf_hr_t12   <- run_hr(tier12_saf_data, tier12_saf_data$study, "NonProc HR T12")

# (g) PROTECT-AF/PREVAIL separate (single-study, no pooling -- report individual)
pp_eff <- eff_is_se %>% filter(study == "PROTECT-AF/PREVAIL")
pp_saf <- saf_nonproc_bleed %>% filter(study == "PROTECT-AF/PREVAIL")
pp_eff_rr <- (pp_eff$ev_device / pp_eff$n_device) / (pp_eff$ev_control / pp_eff$n_control)
pp_eff_se <- sqrt(1/pp_eff$ev_device - 1/pp_eff$n_device +
                  1/pp_eff$ev_control - 1/pp_eff$n_control)
pp_eff_str <- fmt_ci(pp_eff_rr,
                     exp(log(pp_eff_rr) - 1.96 * pp_eff_se),
                     exp(log(pp_eff_rr) + 1.96 * pp_eff_se))
pp_saf_str <- fmt_ci(pp_saf$hr, pp_saf$hr_lo, pp_saf$hr_hi)

# Build Table 4
t4_data <- tibble(
  Analysis = c(
    "Primary (RE, REML, HKSJ)",
    "Fixed-effect (MH/IV)",
    "Leave-one-out (range)",
    "FLX-only (CHAMPION + OPTION)",
    "Contemporary only",
    "Tier 1\u20132 only",
    "PROTECT-AF/PREVAIL alone"
  ),
  `IS+SE RR (95% CI)` = c(
    fmt_meta_rr(m_eff_rr_exp),
    fmt_meta_fe_rr(m_eff_rr_fe),
    sprintf("%.2f\u2013%.2f", loo_eff_range[1], loo_eff_range[2]),
    fmt_meta_rr(m_eff_rr_flx),
    fmt_meta_rr(m_eff_rr_cont),
    fmt_meta_rr(m_eff_rr_t12),
    pp_eff_str
  ),
  `Non-proc Bleeding HR (95% CI)` = c(
    fmt_meta_hr(m_saf_hr_exp),
    fmt_meta_fe_hr(m_saf_hr_fe),
    sprintf("%.2f\u2013%.2f", loo_saf_range[1], loo_saf_range[2]),
    fmt_meta_hr(m_saf_hr_flx),
    fmt_meta_hr(m_saf_hr_cont),
    fmt_meta_hr(m_saf_hr_t12),
    pp_saf_str
  ),
  Interpretation = c(
    "Reference analysis (expanded pool)",
    "Minimal change; low heterogeneity",
    "Point estimate range across LOO iterations",
    "Watchman FLX device only, NOAC comparator",
    "Excludes CLOSURE-AF (BMC comparator) + PP (warfarin)",
    "Excludes KM-reconstructed (T3) data",
    "Historical warfarin comparator, first-gen Watchman"
  )
)

# gt version
gt_t4 <- t4_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Table 4. Sensitivity Analyses Summary",
    subtitle = "Expanded pool unless otherwise stated",
    footnote = paste0(
      "RE = random-effects; REML = restricted maximum likelihood; ",
      "HKSJ = Hartung-Knapp-Sidik-Jonkman; MH = Mantel-Haenszel; ",
      "IV = inverse-variance; LOO = leave-one-out; FLX = Watchman FLX; ",
      "T1-T2 = data provenance tiers 1-2 (trial-reported exact or near-harmonized)."
    )
  ) %>%
  cols_align(align = "center", columns = c(`IS+SE RR (95% CI)`,
                                            `Non-proc Bleeding HR (95% CI)`)) %>%
  cols_align(align = "left", columns = c(Analysis, Interpretation)) %>%
  cols_width(
    Analysis                       ~ px(200),
    `IS+SE RR (95% CI)`            ~ px(160),
    `Non-proc Bleeding HR (95% CI)` ~ px(180),
    Interpretation                  ~ px(280)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Analysis)
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#E8EEF4"), cell_text(weight = "bold")),
    locations = cells_body(rows = 1)
  )

save_gt_png(gt_t4, file.path(DIR_T, "Table4_sensitivity.png"), width = 1300)

# flextable version
ft_t4 <- flextable(t4_data) %>%
  ft_theme_pub() %>%
  bold(j = "Analysis") %>%
  align(j = c("IS+SE RR (95% CI)", "Non-proc Bleeding HR (95% CI)"), align = "center") %>%
  width(j = "Analysis", width = 1.8) %>%
  width(j = "Interpretation", width = 2.5) %>%
  bg(i = 1, bg = "#E8EEF4") %>%
  bold(i = 1) %>%
  add_footer_lines(paste0(
    "RE = random-effects; REML = restricted maximum likelihood; ",
    "HKSJ = Hartung-Knapp-Sidik-Jonkman. LOO = leave-one-out; FLX = Watchman FLX."
  )) %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_t4, file.path(DIR_T, "Table4_sensitivity.docx"))

# =============================================================================
# TABLE 5: Risk of Bias Assessment (RoB 2)
# =============================================================================

cat("\n--- TABLE 5: Risk of Bias Assessment (RoB 2) ---\n")

rob_domains <- c(
  "Randomization process",
  "Deviations from intervention",
  "Missing outcome data",
  "Outcome measurement",
  "Selective reporting",
  "Endpoint architecture bias",
  "Data provenance",
  "Overall"
)

rob_matrix <- matrix(c(
  # CHAMPION-AF
  "Low", "Some concerns", "Low", "Some concerns", "Low",
  "Low", "T1", "Some concerns",
  # OPTION
  "Low", "Some concerns", "High", "Some concerns", "Low",
  "Low", "T1", "Some concerns",
  # PRAGUE-17
  "Low", "Some concerns", "Some concerns", "Some concerns", "Low",
  "Some concerns", "T2", "Some concerns",
  # CLOSURE-AF
  "Low", "Some concerns", "Some concerns", "Some concerns", "Low",
  "High", "T3", "High",
  # PROTECT-AF/PREVAIL
  "Low", "Some concerns", "Low", "Some concerns", "Low",
  "Some concerns", "T1", "Some concerns"
), nrow = 8, ncol = 5)

colnames(rob_matrix) <- STUDY_SHORT

t5_data <- data.frame(
  Domain = rob_domains,
  rob_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# gt version with color coding
gt_t5 <- t5_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Table 5. Risk of Bias Assessment (Modified RoB 2)",
    subtitle = "Including novel domains: endpoint architecture bias and data provenance tier",
    footnote = paste0(
      "RoB 2 = Revised Cochrane Risk of Bias tool for randomized trials. ",
      "T1 = trial-reported exact harmonized; T2 = trial-reported near-harmonized; ",
      "T3 = KM-reconstructed fixed-horizon. ",
      "All trials are open-label, contributing 'Some concerns' for deviations and outcome measurement. ",
      "OPTION: missing data exceed endpoint events (High for missing data). ",
      "CLOSURE-AF: composite endpoint architecture diverges from standard definitions (High)."
    )
  ) %>%
  cols_align(align = "center", columns = -Domain) %>%
  cols_align(align = "left", columns = Domain) %>%
  cols_width(Domain ~ px(200)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Domain)
  )

# Apply color coding to each study column for each risk level
for (col_name in STUDY_SHORT) {
  # Green for Low/T1
  low_rows <- which(t5_data[[col_name]] %in% c("Low", "T1"))
  if (length(low_rows) > 0) {
    gt_t5 <- gt_t5 %>%
      tab_style(
        style     = list(cell_fill(color = GREEN_BG), cell_text(color = GREEN_DARK)),
        locations = cells_body(columns = !!col_name, rows = low_rows)
      )
  }
  # Yellow for Some concerns/T2
  some_rows <- which(t5_data[[col_name]] %in% c("Some concerns", "T2"))
  if (length(some_rows) > 0) {
    gt_t5 <- gt_t5 %>%
      tab_style(
        style     = list(cell_fill(color = YELLOW_BG), cell_text(color = YELLOW_DARK)),
        locations = cells_body(columns = !!col_name, rows = some_rows)
      )
  }
  # Red for High/T3
  high_rows <- which(t5_data[[col_name]] %in% c("High", "T3"))
  if (length(high_rows) > 0) {
    gt_t5 <- gt_t5 %>%
      tab_style(
        style     = list(cell_fill(color = RED_BG), cell_text(color = RED_DARK, weight = "bold")),
        locations = cells_body(columns = !!col_name, rows = high_rows)
      )
  }
}

save_gt_png(gt_t5, file.path(DIR_T, "Table5_rob.png"), width = 1200)

# flextable version with color coding
ft_t5 <- flextable(t5_data) %>%
  ft_theme_pub() %>%
  bold(j = "Domain") %>%
  align(j = STUDY_SHORT, align = "center") %>%
  width(j = "Domain", width = 2.0)

for (col_name in STUDY_SHORT) {
  for (r in seq_len(nrow(t5_data))) {
    val <- t5_data[[col_name]][r]
    if (val %in% c("Low", "T1")) {
      ft_t5 <- ft_t5 %>%
        bg(i = r, j = col_name, bg = GREEN_BG) %>%
        color(i = r, j = col_name, color = GREEN_DARK)
    } else if (val %in% c("Some concerns", "T2")) {
      ft_t5 <- ft_t5 %>%
        bg(i = r, j = col_name, bg = YELLOW_BG) %>%
        color(i = r, j = col_name, color = YELLOW_DARK)
    } else if (val %in% c("High", "T3")) {
      ft_t5 <- ft_t5 %>%
        bg(i = r, j = col_name, bg = RED_BG) %>%
        color(i = r, j = col_name, color = RED_DARK) %>%
        bold(i = r, j = col_name)
    }
  }
}

ft_t5 <- ft_t5 %>%
  add_footer_lines(paste0(
    "RoB 2 = Revised Cochrane Risk of Bias tool. ",
    "T1-T3 = data provenance tiers. All trials open-label."
  )) %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_t5, file.path(DIR_T, "Table5_rob.docx"))


###############################################################################
#                     SUPPLEMENTARY TABLES (sTable 1-8)
###############################################################################

cat("\n\n========================================\n")
cat("  SUPPLEMENTARY TABLES\n")
cat("========================================\n\n")

# =============================================================================
# sTable 1: Endpoint Definitions Across Trials
# =============================================================================

cat("--- sTable 1: Endpoint Definitions ---\n")

st1_data <- tibble(
  `Endpoint Category` = c(
    "Primary composite",
    "Primary composite",
    "Primary composite",
    "Primary composite",
    "Primary composite",
    "Ischemic stroke + SE",
    "Non-proc bleeding",
    "All major bleeding",
    "All-cause death",
    "CV/unexplained death",
    "Hemorrhagic stroke"
  ),
  Trial = c(
    "CHAMPION-AF", "OPTION", "PRAGUE-17", "CLOSURE-AF", "PROTECT-AF/PREVAIL",
    "All trials", "All trials", "All trials", "All trials", "All trials", "All trials"
  ),
  Definition = c(
    "CV death + all stroke + systemic embolism",
    "All-cause death + all stroke + systemic embolism",
    "Stroke/TIA + SE + CV death + clinically significant bleeding + procedure complications",
    "All stroke + SE + ISTH major bleeding + CV/unexplained death (NACE)",
    "All stroke + systemic embolism + cardiovascular death (5-yr pooled)",
    "Ischemic stroke + systemic embolism (arterial); TIA excluded except PRAGUE-17",
    "ISTH major + CRNMB occurring >7 days (CHAMPION/OPTION), >45d (PP); non-procedural",
    "Any ISTH major bleeding including procedural events",
    "Death from any cause",
    "Death attributed to cardiovascular cause or unexplained (CLOSURE includes unexplained)",
    "Hemorrhagic stroke (ICH with focal deficit)"
  ),
  `Adjudication` = c(
    "Blinded CEC", "Blinded CEC", "Blinded CEC", "Blinded CEC", "Blinded CEC",
    "Blinded CEC", "Blinded CEC", "Blinded CEC", "Blinded CEC", "Blinded CEC",
    "Blinded CEC"
  ),
  Notes = c(
    "NI margin: 4.8 ppt absolute (rate difference)",
    "NI margin: 5.0 ppt absolute",
    "NI margin: sHR 1.47 (relative); broad composite",
    "NI margin: HR 1.30; includes bleeding in composite (unusual for efficacy)",
    "NI margin: RR 2.0 (PROTECT) / 1.75 (PREVAIL); warfarin comparator",
    "Harmonized across trials; PRAGUE-17 stroke+TIA used as proxy",
    "CLOSURE-AF: non-proc ISTH major only (CRNMB not reported separately)",
    "CLOSURE-AF: all major including procedural; different threshold than others",
    "Standard definition across all trials",
    "CLOSURE-AF uniquely includes unexplained deaths in CV death category",
    "PRAGUE-17: 0 vs 1 event (excluded from HR pooling)"
  )
)

gt_st1 <- st1_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 1. Endpoint Definitions Across Trials",
    footnote = paste0(
      "CEC = Clinical Events Committee; SE = systemic embolism; TIA = transient ischemic attack; ",
      "ISTH = International Society on Thrombosis and Haemostasis; CRNMB = clinically relevant non-major bleeding; ",
      "NACE = net adverse clinical event; NI = non-inferiority; ppt = percentage points; ",
      "sHR = subdistribution hazard ratio; ICH = intracerebral hemorrhage."
    )
  ) %>%
  cols_align(align = "left") %>%
  cols_width(
    `Endpoint Category` ~ px(140),
    Trial               ~ px(130),
    Definition          ~ px(280),
    `Adjudication`      ~ px(90),
    Notes               ~ px(280)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = `Endpoint Category`)
  )

save_gt_png(gt_st1, file.path(DIR_ST, "sTable1_endpoint_defs.png"), width = 1500)

ft_st1 <- flextable(st1_data) %>%
  ft_theme_pub() %>%
  bold(j = "Endpoint Category") %>%
  width(j = "Definition", width = 2.8) %>%
  width(j = "Notes", width = 2.8) %>%
  add_footer_lines("CEC = Clinical Events Committee; SE = systemic embolism; ISTH = Intl Soc Thrombosis & Haemostasis.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st1, file.path(DIR_ST, "sTable1_endpoint_defs.docx"))

# =============================================================================
# sTable 2: Data Provenance Tier Classification
# =============================================================================

cat("--- sTable 2: Data Provenance Tiers ---\n")

# Build a matrix: rows = endpoint, columns = studies
st2_endpoints <- c(
  "Ischemic stroke + SE",
  "Non-proc bleeding",
  "All major bleeding",
  "All-cause death",
  "CV/unexplained death",
  "Hemorrhagic stroke",
  "Any stroke + SE",
  "Ischemic stroke alone"
)

st2_datasets <- list(
  eff_is_se, saf_nonproc_bleed, saf_all_major_bleed,
  sec_acm, sec_cvdeath, sec_hemstroke,
  sec_anystroke_se, sec_is_alone
)

st2_matrix <- matrix(NA_character_, nrow = length(st2_endpoints), ncol = 5)
colnames(st2_matrix) <- STUDY_SHORT

for (i in seq_along(st2_datasets)) {
  d <- st2_datasets[[i]]
  for (j in seq_along(STUDY_ORDER)) {
    row <- d %>% filter(study == STUDY_ORDER[j])
    if (nrow(row) > 0 && "tier" %in% names(row)) {
      st2_matrix[i, j] <- row$tier[1]
    } else if (nrow(row) > 0) {
      st2_matrix[i, j] <- "\u2014"
    } else {
      st2_matrix[i, j] <- "\u2014"
    }
  }
}

st2_data <- data.frame(
  Endpoint = st2_endpoints,
  st2_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gt_st2 <- st2_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 2. Data Provenance Tier Classification",
    subtitle = "T1 = trial-reported exact; T2 = near-harmonized; T3 = KM-reconstructed; T4 = derived",
    footnote = paste0(
      "T1: endpoint reported with exact definition match to harmonized analysis. ",
      "T2: minor definitional differences (e.g., stroke+TIA used as proxy for stroke+SE). ",
      "T3: events/hazards reconstructed from published KM curves. ",
      "T4: model-derived or imputed values. \u2014 = data not available for this endpoint."
    )
  ) %>%
  cols_align(align = "center", columns = -Endpoint) %>%
  cols_align(align = "left", columns = Endpoint) %>%
  cols_width(Endpoint ~ px(180)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Endpoint))

# Color code tiers
for (col_name in STUDY_SHORT) {
  t1_rows <- which(st2_data[[col_name]] == "T1")
  t2_rows <- which(st2_data[[col_name]] == "T2")
  t3_rows <- which(st2_data[[col_name]] %in% c("T3", "T4"))

  if (length(t1_rows) > 0)
    gt_st2 <- gt_st2 %>%
      tab_style(style = list(cell_fill(color = GREEN_BG), cell_text(color = GREEN_DARK)),
                locations = cells_body(columns = !!col_name, rows = t1_rows))
  if (length(t2_rows) > 0)
    gt_st2 <- gt_st2 %>%
      tab_style(style = list(cell_fill(color = YELLOW_BG), cell_text(color = YELLOW_DARK)),
                locations = cells_body(columns = !!col_name, rows = t2_rows))
  if (length(t3_rows) > 0)
    gt_st2 <- gt_st2 %>%
      tab_style(style = list(cell_fill(color = RED_BG), cell_text(color = RED_DARK)),
                locations = cells_body(columns = !!col_name, rows = t3_rows))
}

save_gt_png(gt_st2, file.path(DIR_ST, "sTable2_data_provenance.png"))

ft_st2 <- flextable(st2_data) %>%
  ft_theme_pub() %>%
  bold(j = "Endpoint") %>%
  align(j = STUDY_SHORT, align = "center") %>%
  width(j = "Endpoint", width = 1.8)

for (col_name in STUDY_SHORT) {
  for (r in seq_len(nrow(st2_data))) {
    val <- st2_data[[col_name]][r]
    if (identical(val, "T1"))
      ft_st2 <- ft_st2 %>% bg(i = r, j = col_name, bg = GREEN_BG) %>%
        color(i = r, j = col_name, color = GREEN_DARK)
    else if (identical(val, "T2"))
      ft_st2 <- ft_st2 %>% bg(i = r, j = col_name, bg = YELLOW_BG) %>%
        color(i = r, j = col_name, color = YELLOW_DARK)
    else if (val %in% c("T3", "T4"))
      ft_st2 <- ft_st2 %>% bg(i = r, j = col_name, bg = RED_BG) %>%
        color(i = r, j = col_name, color = RED_DARK)
  }
}

ft_st2 <- ft_st2 %>%
  add_footer_lines("T1 = exact harmonized; T2 = near-harmonized; T3 = KM-reconstructed; T4 = derived.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st2, file.path(DIR_ST, "sTable2_data_provenance.docx"))

# =============================================================================
# sTable 3: Antithrombotic Regimens Over Time
# =============================================================================

cat("--- sTable 3: Antithrombotic Regimens ---\n")

st3_data <- antithrombotic %>%
  mutate(
    Trial = study,
    `Device: OAC (%)` = sprintf("%.1f", device_oac_start),
    `Device: DAPT (%)` = sprintf("%.1f", device_dapt_start),
    `Device: AP mono (%)` = sprintf("%.1f", device_ap_start),
    `Device: None (%)` = sprintf("%.1f", device_none_start),
    `Control: OAC (%)` = sprintf("%.1f", ctrl_oac_start),
    `Control: None (%)` = sprintf("%.1f", ctrl_none_start)
  ) %>%
  select(Trial, starts_with("Device:"), starts_with("Control:"))

gt_st3 <- st3_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 3. Antithrombotic Regimens at Randomization",
    subtitle = "Initial post-procedure/randomization regimen (%)",
    footnote = paste0(
      "OAC = oral anticoagulant (NOAC or warfarin); DAPT = dual antiplatelet therapy; ",
      "AP = antiplatelet; CHAMPION-AF device arm received 6 weeks of OAC post-implant; ",
      "OPTION device arm: 88% DAPT, 12% AP monotherapy. ",
      "CLOSURE-AF: 9.1% on OAC at baseline reflects shared decision-making protocol."
    )
  ) %>%
  cols_align(align = "center", columns = -Trial) %>%
  cols_align(align = "left", columns = Trial) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Trial)) %>%
  tab_spanner(label = "Device Arm", columns = starts_with("Device:")) %>%
  tab_spanner(label = "Control Arm", columns = starts_with("Control:"))

save_gt_png(gt_st3, file.path(DIR_ST, "sTable3_antithrombotic.png"), width = 1300)

ft_st3 <- flextable(st3_data) %>%
  ft_theme_pub() %>%
  bold(j = "Trial") %>%
  align(j = -1, align = "center") %>%
  add_header_row(
    values = c("", "Device Arm", "Control Arm"),
    colwidths = c(1, 4, 2)
  ) %>%
  bg(i = 1, bg = HDR_BLUE, part = "header") %>%
  color(i = 1, color = HDR_TEXT, part = "header") %>%
  bold(i = 1, part = "header") %>%
  add_footer_lines("OAC = oral anticoagulant; DAPT = dual antiplatelet; AP = antiplatelet.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st3, file.path(DIR_ST, "sTable3_antithrombotic.docx"))

# =============================================================================
# sTable 4: Complete Endpoint Results by Trial
# =============================================================================

cat("--- sTable 4: Complete Endpoint Results by Trial ---\n")

# Build comprehensive data: all endpoints x all studies
all_ep_list <- list(
  list(data = eff_is_se,           name = "Ischemic stroke + SE"),
  list(data = saf_nonproc_bleed,   name = "Non-proc bleeding"),
  list(data = saf_all_major_bleed, name = "All major bleeding"),
  list(data = sec_acm,             name = "All-cause death"),
  list(data = sec_cvdeath,         name = "CV/unexplained death"),
  list(data = sec_hemstroke,       name = "Hemorrhagic stroke"),
  list(data = sec_anystroke_se,    name = "Any stroke + SE"),
  list(data = sec_is_alone,        name = "Ischemic stroke alone")
)

st4_rows <- list()
for (ep in all_ep_list) {
  for (s in STUDY_ORDER) {
    row <- ep$data %>% filter(study == s)
    if (nrow(row) > 0) {
      tier_val <- if ("tier" %in% names(row)) row$tier[1] else "\u2014"
      hr_str <- if (!is.na(row$hr[1])) {
        fmt_ci(row$hr[1], row$hr_lo[1], row$hr_hi[1])
      } else {
        "\u2014"
      }
      rr_val <- (row$ev_device[1] / row$n_device[1]) /
                (row$ev_control[1] / row$n_control[1])
      st4_rows[[length(st4_rows) + 1]] <- tibble(
        Endpoint     = ep$name,
        Trial        = s,
        `Ev (Device)` = row$ev_device[1],
        `N (Device)`  = row$n_device[1],
        `Ev (Control)` = row$ev_control[1],
        `N (Control)`  = row$n_control[1],
        `Crude RR`    = sprintf("%.2f", rr_val),
        `HR (95% CI)` = hr_str,
        Tier         = tier_val
      )
    }
  }
}

st4_data <- bind_rows(st4_rows)

gt_st4 <- st4_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 4. Complete Endpoint Results by Trial",
    subtitle = "Individual trial event counts, crude RR, and reported HR with 95% CI",
    footnote = paste0(
      "Ev = events; RR = risk ratio (crude); HR = hazard ratio (as reported by trial); ",
      "Tier = data provenance (T1 = exact, T2 = near-harmonized, T3 = reconstructed). ",
      "\u2014 = not available/not applicable."
    )
  ) %>%
  cols_align(align = "center", columns = -c(Endpoint, Trial)) %>%
  cols_align(align = "left", columns = c(Endpoint, Trial)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Endpoint))

save_gt_png(gt_st4, file.path(DIR_ST, "sTable4_all_endpoints.png"), width = 1500)

ft_st4 <- flextable(st4_data) %>%
  ft_theme_pub() %>%
  bold(j = "Endpoint") %>%
  align(j = c("Ev (Device)", "N (Device)", "Ev (Control)", "N (Control)",
              "Crude RR", "HR (95% CI)", "Tier"), align = "center") %>%
  width(j = "Endpoint", width = 1.5) %>%
  width(j = "HR (95% CI)", width = 1.5) %>%
  add_footer_lines("Ev = events; RR = crude risk ratio; HR = hazard ratio; Tier = data provenance.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st4, file.path(DIR_ST, "sTable4_all_endpoints.docx"))

# =============================================================================
# sTable 5: Procedural Outcomes Detail
# =============================================================================

cat("--- sTable 5: Procedural Outcomes ---\n")

st5_data <- procedural_perf %>%
  mutate(
    Trial             = study,
    `Implant Success`  = sprintf("%d/%d (%.1f%%)", implant_success, n_attempted,
                                  implant_success / n_attempted * 100),
    `Pericardial Eff.` = sprintf("%d (%.1f%%)", pericardial_eff,
                                  pericardial_eff / n_attempted * 100),
    `Device Embol.`    = sprintf("%d (%.2f%%)", device_embol,
                                  device_embol / n_attempted * 100),
    `Proc. Stroke`     = sprintf("%d (%.2f%%)", proc_stroke,
                                  proc_stroke / n_attempted * 100),
    `Proc. Death`      = sprintf("%d (%.2f%%)", proc_death,
                                  proc_death / n_attempted * 100),
    `Overall Major`    = sprintf("%d (%.1f%%)",
                                  pericardial_eff + device_embol + proc_stroke + proc_death,
                                  (pericardial_eff + device_embol + proc_stroke + proc_death) /
                                    n_attempted * 100)
  ) %>%
  select(Trial, `Implant Success`, `Pericardial Eff.`, `Device Embol.`,
         `Proc. Stroke`, `Proc. Death`, `Overall Major`)

gt_st5 <- st5_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 5. Procedural Outcomes by Trial",
    subtitle = "Device arm only; events / N attempted (%)",
    footnote = paste0(
      "Implant success = successful device deployment; Pericardial Eff. = pericardial effusion ",
      "requiring intervention (tamponade); Embol. = device embolization requiring retrieval; ",
      "Overall Major = sum of listed major procedural complications. ",
      "PROTECT-AF and PREVAIL listed separately (different devices/eras)."
    )
  ) %>%
  cols_align(align = "center", columns = -Trial) %>%
  cols_align(align = "left", columns = Trial) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Trial))

save_gt_png(gt_st5, file.path(DIR_ST, "sTable5_procedural.png"), width = 1300)

ft_st5 <- flextable(st5_data) %>%
  ft_theme_pub() %>%
  bold(j = "Trial") %>%
  align(j = -1, align = "center") %>%
  add_footer_lines("Pericardial Eff. = pericardial effusion/tamponade; Embol. = device embolization.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st5, file.path(DIR_ST, "sTable5_procedural.docx"))

# =============================================================================
# sTable 6: CLOSURE-AF Endpoint Cascade
# =============================================================================

cat("--- sTable 6: CLOSURE-AF Endpoint Cascade ---\n")

st6_data <- closure_cascade %>%
  mutate(
    n_dev = 446, n_ctrl = 442,
    rr        = (ev_device / n_dev) / (ev_ctrl / n_ctrl),
    rr_se     = sqrt(1/pmax(ev_device,0.5) - 1/n_dev + 1/pmax(ev_ctrl,0.5) - 1/n_ctrl),
    rr_lo     = exp(log(rr) - 1.96 * rr_se),
    rr_hi     = exp(log(rr) + 1.96 * rr_se)
  ) %>%
  transmute(
    Model           = model,
    `Endpoint Definition` = endpoint_def,
    `Ev (Device)`   = ev_device,
    `Ev (Control)`  = ev_ctrl,
    `RR (95% CI)`   = fmt_ci(rr, rr_lo, rr_hi),
    `HR (MC)`       = sprintf("%.2f", hr_mc),
    `dRMST (yr)`    = sprintf("%+.3f", drmst),
    Direction       = ifelse(rr > 1, "Favors OAC", "Favors Device")
  )

gt_st6 <- st6_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 6. CLOSURE-AF Endpoint Cascade Analysis",
    subtitle = "Disaggregation from published NACE to individual components",
    footnote = paste0(
      "NACE = net adverse clinical events; MC = model-calibrated; ",
      "dRMST = difference in restricted mean survival time (years, device minus control); ",
      "Negative dRMST indicates fewer event-free years with device. ",
      "N = 446 (device), 442 (control). Events from published supplement + KM reconstruction."
    )
  ) %>%
  cols_align(align = "center", columns = -c(Model, `Endpoint Definition`)) %>%
  cols_align(align = "left", columns = c(Model, `Endpoint Definition`)) %>%
  cols_width(
    Model                 ~ px(160),
    `Endpoint Definition` ~ px(250)
  ) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Model))

# Color the Direction column
dev_rows <- which(st6_data$Direction == "Favors Device")
oac_rows_s6 <- which(st6_data$Direction == "Favors OAC")
if (length(dev_rows) > 0)
  gt_st6 <- gt_st6 %>%
    tab_style(style = list(cell_fill(color = GREEN_BG), cell_text(color = GREEN_DARK)),
              locations = cells_body(columns = Direction, rows = dev_rows))
if (length(oac_rows_s6) > 0)
  gt_st6 <- gt_st6 %>%
    tab_style(style = list(cell_fill(color = RED_BG), cell_text(color = RED_DARK)),
              locations = cells_body(columns = Direction, rows = oac_rows_s6))

save_gt_png(gt_st6, file.path(DIR_ST, "sTable6_closure_cascade.png"), width = 1400)

ft_st6 <- flextable(st6_data) %>%
  ft_theme_pub() %>%
  bold(j = "Model") %>%
  align(j = c("Ev (Device)", "Ev (Control)", "RR (95% CI)", "HR (MC)",
              "dRMST (yr)", "Direction"), align = "center") %>%
  width(j = "Endpoint Definition", width = 2.5) %>%
  add_footer_lines("NACE = net adverse clinical events; dRMST = difference in restricted mean survival time.") %>%
  fontsize(size = 8, part = "footer")

for (r in dev_rows) {
  ft_st6 <- ft_st6 %>% bg(i = r, j = "Direction", bg = GREEN_BG) %>%
    color(i = r, j = "Direction", color = GREEN_DARK)
}
for (r in oac_rows_s6) {
  ft_st6 <- ft_st6 %>% bg(i = r, j = "Direction", bg = RED_BG) %>%
    color(i = r, j = "Direction", color = RED_DARK)
}

save_flextable_docx(ft_st6, file.path(DIR_ST, "sTable6_closure_cascade.docx"))

# =============================================================================
# sTable 7: CV Death Cause-Specific Breakdown (CLOSURE-AF)
# =============================================================================

cat("--- sTable 7: CV Death Breakdown (CLOSURE-AF) ---\n")

# CLOSURE-AF: CV death detail (from data: 99 device vs 81 control CV/unexpl deaths)
# Total all-cause: 155 vs 141; CV/unexplained: 99 vs 81; non-CV: 56 vs 60
st7_data <- tibble(
  `Cause of Death` = c(
    "All-cause death",
    "CV/unexplained death",
    "  Stroke (fatal)",
    "  Cardiovascular (non-stroke)",
    "  Unexplained",
    "Non-CV death"
  ),
  `Device (N=446)` = c(
    "155 (34.8%)", "99 (22.2%)",
    "10 (2.2%)", "51 (11.4%)", "38 (8.5%)",
    "56 (12.6%)"
  ),
  `Control (N=442)` = c(
    "141 (31.9%)", "81 (18.3%)",
    "13 (2.9%)", "42 (9.5%)", "26 (5.9%)",
    "60 (13.6%)"
  ),
  `HR (95% CI)` = c(
    "1.11 (0.89\u20131.39)", "1.24 (0.93\u20131.66)",
    "0.77 (0.34\u20131.75)", "1.22 (0.82\u20131.83)", "1.47 (0.89\u20132.42)",
    "0.94 (0.65\u20131.35)"
  ),
  Notes = c(
    "Primary safety numerator",
    "Includes unexplained deaths (unique to CLOSURE-AF)",
    "Includes hemorrhagic and ischemic fatal strokes",
    "Heart failure, MI, sudden cardiac death",
    "38 vs 26: largest differential in unexplained category",
    "Cancer, infection, other non-CV causes"
  )
)

gt_st7 <- st7_data %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 7. Cause-Specific Mortality in CLOSURE-AF",
    subtitle = "Device (N=446) vs Control (N=442), 3-year follow-up",
    footnote = paste0(
      "HR = hazard ratio; CV = cardiovascular; MI = myocardial infarction. ",
      "CLOSURE-AF uniquely categorizes unexplained deaths under CV/unexplained, ",
      "which inflates the CV death numerator compared to other trials. ",
      "The 38 vs 26 unexplained death difference drives the overall CV/unexplained death signal."
    )
  ) %>%
  cols_align(align = "center", columns = c(`Device (N=446)`, `Control (N=442)`,
                                             `HR (95% CI)`)) %>%
  cols_align(align = "left", columns = c(`Cause of Death`, Notes)) %>%
  cols_width(
    `Cause of Death` ~ px(200),
    Notes            ~ px(280)
  ) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = `Cause of Death`,
                                   rows = c(1, 2, 6)))

# Highlight the unexplained deaths row
gt_st7 <- gt_st7 %>%
  tab_style(
    style     = list(cell_fill(color = YELLOW_BG)),
    locations = cells_body(rows = 5)
  )

save_gt_png(gt_st7, file.path(DIR_ST, "sTable7_cvdeath_closure.png"), width = 1300)

ft_st7 <- flextable(st7_data) %>%
  ft_theme_pub() %>%
  bold(j = "Cause of Death", i = c(1, 2, 6)) %>%
  align(j = c("Device (N=446)", "Control (N=442)", "HR (95% CI)"), align = "center") %>%
  width(j = "Notes", width = 2.5) %>%
  bg(i = 5, bg = YELLOW_BG) %>%
  add_footer_lines("Unexplained deaths (row highlighted) drive the CV/unexplained death differential.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st7, file.path(DIR_ST, "sTable7_cvdeath_closure.docx"))

# =============================================================================
# sTable 8: NI Margin Re-Expression Heatmap (as Table)
# =============================================================================

cat("--- sTable 8: NI Margin Heatmap ---\n")

st8_data <- ni_margin_data %>%
  mutate(across(-study, ~ ifelse(.x, "Compatible", "Not met"))) %>%
  rename(
    Trial      = study,
    `RR < 1.15` = rr_115,
    `RR < 1.25` = rr_125,
    `RR < 1.40` = rr_140,
    `< 0.5 ppt` = abs_05,
    `< 1.0 ppt` = abs_10,
    `< 2.0 ppt` = abs_20
  )

# Count row: how many trials are compatible per threshold?
# (compute from original data)
count_row <- ni_margin_data %>%
  summarise(
    Trial      = "Compatible (n/5)",
    `RR < 1.15` = as.character(sum(rr_115)),
    `RR < 1.25` = as.character(sum(rr_125)),
    `RR < 1.40` = as.character(sum(rr_140)),
    `< 0.5 ppt` = as.character(sum(abs_05)),
    `< 1.0 ppt` = as.character(sum(abs_10)),
    `< 2.0 ppt` = as.character(sum(abs_20))
  )

st8_full <- bind_rows(st8_data, count_row)

gt_st8 <- st8_full %>%
  gt() %>%
  gt_theme_pub(
    title    = "Supplementary Table 8. Non-Inferiority Margin Re-Expression",
    subtitle = "Compatibility of each trial with alternative NI thresholds for IS+SE",
    footnote = paste0(
      "RR = risk ratio upper bound of 95% CI; ppt = percentage points absolute risk difference. ",
      "Compatible = upper 95% CI of trial estimate falls below the specified NI margin. ",
      "CHAMPION-AF and CLOSURE-AF have wider CIs due to higher IS+SE point estimates. ",
      "PROTECT-AF/PREVAIL not compatible with any modern threshold (warfarin comparator, first-gen device)."
    )
  ) %>%
  cols_align(align = "center", columns = -Trial) %>%
  cols_align(align = "left", columns = Trial) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Trial)) %>%
  tab_spanner(label = "Relative Margins", columns = c(`RR < 1.15`, `RR < 1.25`, `RR < 1.40`)) %>%
  tab_spanner(label = "Absolute Margins", columns = c(`< 0.5 ppt`, `< 1.0 ppt`, `< 2.0 ppt`))

# Color coding: Compatible = green, Not met = red, counts = neutral
threshold_cols <- c("RR < 1.15", "RR < 1.25", "RR < 1.40",
                    "< 0.5 ppt", "< 1.0 ppt", "< 2.0 ppt")

for (col_name in threshold_cols) {
  compat_rows <- which(st8_full[[col_name]] == "Compatible")
  notmet_rows <- which(st8_full[[col_name]] == "Not met")

  if (length(compat_rows) > 0)
    gt_st8 <- gt_st8 %>%
      tab_style(style = list(cell_fill(color = GREEN_BG), cell_text(color = GREEN_DARK)),
                locations = cells_body(columns = !!col_name, rows = compat_rows))
  if (length(notmet_rows) > 0)
    gt_st8 <- gt_st8 %>%
      tab_style(style = list(cell_fill(color = RED_BG), cell_text(color = RED_DARK)),
                locations = cells_body(columns = !!col_name, rows = notmet_rows))
}

# Bold the summary row
gt_st8 <- gt_st8 %>%
  tab_style(
    style     = list(cell_fill(color = "#E8EEF4"), cell_text(weight = "bold")),
    locations = cells_body(rows = nrow(st8_full))
  )

save_gt_png(gt_st8, file.path(DIR_ST, "sTable8_ni_heatmap.png"), width = 1200)

# flextable version
ft_st8 <- flextable(st8_full) %>%
  ft_theme_pub() %>%
  bold(j = "Trial") %>%
  align(j = threshold_cols, align = "center") %>%
  add_header_row(
    values = c("", "Relative Margins", "Absolute Margins"),
    colwidths = c(1, 3, 3)
  ) %>%
  bg(i = 1, bg = HDR_BLUE, part = "header") %>%
  color(i = 1, color = HDR_TEXT, part = "header") %>%
  bold(i = 1, part = "header")

for (col_name in threshold_cols) {
  for (r in seq_len(nrow(st8_full))) {
    val <- st8_full[[col_name]][r]
    if (identical(val, "Compatible"))
      ft_st8 <- ft_st8 %>%
        bg(i = r, j = col_name, bg = GREEN_BG) %>%
        color(i = r, j = col_name, color = GREEN_DARK)
    else if (identical(val, "Not met"))
      ft_st8 <- ft_st8 %>%
        bg(i = r, j = col_name, bg = RED_BG) %>%
        color(i = r, j = col_name, color = RED_DARK)
  }
}

# Bold summary row
ft_st8 <- ft_st8 %>%
  bg(i = nrow(st8_full), bg = "#E8EEF4") %>%
  bold(i = nrow(st8_full)) %>%
  add_footer_lines("RR = risk ratio; ppt = percentage points. Compatible = upper 95% CI below NI threshold.") %>%
  fontsize(size = 8, part = "footer")

save_flextable_docx(ft_st8, file.path(DIR_ST, "sTable8_ni_heatmap.docx"))


###############################################################################
# COMPLETION
###############################################################################

cat("\n========================================\n")
cat("  08_tables.R COMPLETE\n")
cat("========================================\n\n")

cat("Main tables saved to:         ", DIR_T, "\n")
cat("Supplementary tables saved to:", DIR_ST, "\n")
cat("\nFiles generated:\n")
cat("  MAIN:\n")
cat("    Table1_study_design.png/.docx\n")
cat("    Table2_baseline.png/.docx\n")
cat("    Table3_main_results.png/.docx\n")
cat("    Table4_sensitivity.png/.docx\n")
cat("    Table5_rob.png/.docx\n")
cat("  SUPPLEMENTARY:\n")
cat("    sTable1_endpoint_defs.png/.docx\n")
cat("    sTable2_data_provenance.png/.docx\n")
cat("    sTable3_antithrombotic.png/.docx\n")
cat("    sTable4_all_endpoints.png/.docx\n")
cat("    sTable5_procedural.png/.docx\n")
cat("    sTable6_closure_cascade.png/.docx\n")
cat("    sTable7_cvdeath_closure.png/.docx\n")
cat("    sTable8_ni_heatmap.png/.docx\n")
cat("\n>> 08_tables.R complete.\n")
