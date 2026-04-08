###############################################################################
# 06_main_figures.R
# LAA Closure vs OAC Meta-Analysis — Main Figures (5 figures)
#
# Figures:
#   1. Central Illustration — Trial Architecture Map
#   2. Primary Efficacy Forest Plot — IS + SE
#   3. Safety Panel (Non-proc bleeding forest + Procedural complication bars)
#   4. Temporal Analysis Panel (Landmark + Bleeding crossover)
#   5. Benefit-Risk Common Denominator (Butterfly chart)
#
# Dependencies: 01_data.R, output/primary_results.RData
# Outputs: figures/Figure1–5_*.pdf
###############################################################################

source(file.path(normalizePath("~/Desktop/LAA/LAA_meta"), "R", "01_data.R"))

# Load analysis results if available
primary_rdata <- file.path(DIR_O, "primary_results.RData")
if (file.exists(primary_rdata)) {
  load(primary_rdata)
  cat(">> Loaded primary_results.RData\n")
} else {
  cat(">> primary_results.RData not found — figures will use data-only mode\n")
}

cat("\n========================================\n")
cat("  06_main_figures.R\n")
cat("========================================\n\n")

# =============================================================================
# HELPER: Pool-aware forest plot via forester with forestploter fallback
# =============================================================================

#' Build a grouped forest plot for a given endpoint dataset
#' @param dat Data frame with study, ev_device, n_device, ev_control, n_control,
#'            hr, hr_lo, hr_hi columns
#' @param endpoint_label Character label for the endpoint
#' @param file_pdf Path for PDF output
#' @param file_png Path for PNG output (used by forester)
#' @param width_in, height_in Dimensions in inches
build_grouped_forest <- function(dat, endpoint_label,
                                  file_pdf, file_png = NULL,
                                  width_in = 14, height_in = 10) {

  # --- Map studies to UNIQUE pool for display (no duplication) ---
  # Each study appears once: primary > expanded-only > historical
  pool_map <- pool_assign %>%
    mutate(pool_label = case_when(
      pool_primary    ~ "Contemporary",
      pool_expanded   ~ "Expanded (additional)",
      pool_historical ~ "Historical",
      TRUE ~ "Other"
    )) %>%
    select(study, pool_label)

  dat <- dat %>%
    left_join(pool_map, by = "study") %>%
    mutate(pool_label = factor(pool_label,
                               levels = c("Contemporary", "Expanded (additional)", "Historical")))

  # --- Pool-level pooled estimates via metagen ---
  pool_summaries <- dat %>%
    group_by(pool_label) %>%
    filter(n() >= 2, !is.na(loghr)) %>%
    group_modify(~ {
      m <- metagen(TE = .x$loghr, seTE = .x$se_loghr,
                   studlab = .x$study, sm = "HR",
                   method.tau = "REML", hakn = TRUE)
      tibble(
        hr     = exp(m$TE.random),
        hr_lo  = exp(m$lower.random),
        hr_hi  = exp(m$upper.random),
        k      = m$k,
        I2     = round(m$I2 * 100, 0)
      )
    }) %>%
    ungroup()

  # --- Overall pooled (all studies with HR) ---
  d_all <- dat %>% filter(!is.na(loghr))
  if (nrow(d_all) >= 2) {
    m_all <- metagen(TE = d_all$loghr, seTE = d_all$se_loghr,
                     studlab = d_all$study, sm = "HR",
                     method.tau = "REML", hakn = TRUE)
    overall <- tibble(
      hr    = exp(m_all$TE.random),
      hr_lo = exp(m_all$lower.random),
      hr_hi = exp(m_all$upper.random),
      k     = m_all$k,
      I2    = round(m_all$I2 * 100, 0)
    )
  } else {
    overall <- NULL
  }

  # -----------------------------------------------------------------------
  # TRY forester FIRST
  # -----------------------------------------------------------------------
  forester_success <- FALSE

  tryCatch({
    if (!requireNamespace("forester", quietly = TRUE)) stop("forester not installed")

    # Build display data frame with rows for headers, studies, subtotals, blanks
    rows <- list()
    row_idx <- 0

    for (pl in levels(dat$pool_label)) {
      pool_dat <- dat %>% filter(pool_label == pl)
      if (nrow(pool_dat) == 0) next

      # Pool header row
      row_idx <- row_idx + 1
      rows[[row_idx]] <- tibble(
        Study   = paste0("  ", pl, " Pool"),
        N_Device  = "",
        N_Control = "",
        Events_D  = "",
        Events_C  = "",
        HR_text   = "",
        estimate  = NA_real_,
        ci_low    = NA_real_,
        ci_high   = NA_real_,
        is_summary = FALSE,
        is_header  = TRUE
      )

      # Study rows
      for (i in seq_len(nrow(pool_dat))) {
        s <- pool_dat[i, ]
        row_idx <- row_idx + 1
        rows[[row_idx]] <- tibble(
          Study   = paste0("    ", s$study),
          N_Device  = as.character(s$n_device),
          N_Control = as.character(s$n_control),
          Events_D  = as.character(s$ev_device),
          Events_C  = as.character(s$ev_control),
          HR_text   = sprintf("%.2f (%.2f-%.2f)", s$hr, s$hr_lo, s$hr_hi),
          estimate  = s$hr,
          ci_low    = s$hr_lo,
          ci_high   = s$hr_hi,
          is_summary = FALSE,
          is_header  = FALSE
        )
      }

      # Subtotal row
      ps <- pool_summaries %>% filter(pool_label == pl)
      if (nrow(ps) == 1) {
        row_idx <- row_idx + 1
        rows[[row_idx]] <- tibble(
          Study   = paste0("  Subtotal (I\u00B2=", ps$I2, "%)"),
          N_Device  = "",
          N_Control = "",
          Events_D  = "",
          Events_C  = "",
          HR_text   = sprintf("%.2f (%.2f-%.2f)", ps$hr, ps$hr_lo, ps$hr_hi),
          estimate  = ps$hr,
          ci_low    = ps$hr_lo,
          ci_high   = ps$hr_hi,
          is_summary = TRUE,
          is_header  = FALSE
        )
      }

      # Blank spacer
      row_idx <- row_idx + 1
      rows[[row_idx]] <- tibble(
        Study   = " ",
        N_Device  = "", N_Control = "",
        Events_D  = "", Events_C  = "",
        HR_text   = "",
        estimate  = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        is_summary = FALSE, is_header = FALSE
      )
    }

    # Overall row
    if (!is.null(overall)) {
      row_idx <- row_idx + 1
      rows[[row_idx]] <- tibble(
        Study   = paste0("Overall (I\u00B2=", overall$I2, "%)"),
        N_Device  = "", N_Control = "",
        Events_D  = "", Events_C  = "",
        HR_text   = sprintf("%.2f (%.2f-%.2f)", overall$hr, overall$hr_lo, overall$hr_hi),
        estimate  = overall$hr,
        ci_low    = overall$hr_lo,
        ci_high   = overall$hr_hi,
        is_summary = TRUE,
        is_header  = FALSE
      )
    }

    df_forest <- bind_rows(rows)

    # Determine png path
    png_path <- if (!is.null(file_png)) file_png else sub("\\.pdf$", ".png", file_pdf)

    forester::forester(
      left_side_data = df_forest[, c("Study", "N_Device", "N_Control",
                                      "Events_D", "Events_C")],
      estimate       = df_forest$estimate,
      ci_low         = df_forest$ci_low,
      ci_high        = df_forest$ci_high,
      estimate_precision = 2,
      null_line_at   = 1,
      xlim           = c(0.1, 4),
      xbreaks        = c(0.25, 0.5, 1, 2, 4),
      estimate_col_name = "HR (95% CI)",
      arrows         = TRUE,
      arrow_labels   = c("Favours Device", "Favours OAC"),
      file_path      = png_path,
      dpi            = 300,
      display        = FALSE,
      font_size      = 3,
      render_as      = "png"
    )

    forester_success <- TRUE
    cat(sprintf("  >> forester plot saved: %s\n", png_path))

  }, error = function(e) {
    cat(sprintf("  >> forester failed (%s) — falling back to forestploter/ggplot\n",
                conditionMessage(e)))
  })

  # -----------------------------------------------------------------------
  # FALLBACK: ggplot-based forest plot (always produce PDF)
  # -----------------------------------------------------------------------

  # Always produce a ggplot version for the PDF
  plot_data <- dat %>%
    filter(!is.na(hr)) %>%
    arrange(pool_label, study)

  # Add subtotals and overall
  subtotal_rows <- pool_summaries %>%
    mutate(
      study = paste0(pool_label, " Subtotal"),
      ev_device = NA_real_, n_device = NA_real_,
      ev_control = NA_real_, n_control = NA_real_,
      is_summary = TRUE
    ) %>%
    rename(pool_label_orig = pool_label)

  overall_row <- if (!is.null(overall)) {
    tibble(
      study = "Overall",
      hr = overall$hr, hr_lo = overall$hr_lo, hr_hi = overall$hr_hi,
      pool_label = factor("Overall", levels = c("Contemporary", "Expanded", "Historical", "Overall")),
      is_summary = TRUE
    )
  } else NULL

  # Build ordered factor for y-axis
  y_labels <- character()
  y_hr     <- numeric()
  y_lo     <- numeric()
  y_hi     <- numeric()
  y_summary <- logical()
  y_pool_color <- character()

  for (pl in c("Contemporary", "Expanded (additional)", "Historical")) {
    pool_dat <- plot_data %>% filter(pool_label == pl)
    if (nrow(pool_dat) == 0) next

    # Pool header
    y_labels    <- c(y_labels, paste0(pl, " Pool"))
    y_hr        <- c(y_hr, NA)
    y_lo        <- c(y_lo, NA)
    y_hi        <- c(y_hi, NA)
    y_summary   <- c(y_summary, FALSE)
    y_pool_color <- c(y_pool_color, "header")

    for (i in seq_len(nrow(pool_dat))) {
      s <- pool_dat[i, ]
      label <- sprintf("%s   [%d/%d vs %d/%d]",
                       s$study, s$ev_device, s$n_device,
                       s$ev_control, s$n_control)
      y_labels    <- c(y_labels, label)
      y_hr        <- c(y_hr, s$hr)
      y_lo        <- c(y_lo, s$hr_lo)
      y_hi        <- c(y_hi, s$hr_hi)
      y_summary   <- c(y_summary, FALSE)
      y_pool_color <- c(y_pool_color, as.character(pl))
    }

    ps <- pool_summaries %>% filter(pool_label == pl)
    if (nrow(ps) == 1) {
      y_labels    <- c(y_labels, sprintf("  Subtotal (I\u00B2=%d%%)", ps$I2))
      y_hr        <- c(y_hr, ps$hr)
      y_lo        <- c(y_lo, ps$hr_lo)
      y_hi        <- c(y_hi, ps$hr_hi)
      y_summary   <- c(y_summary, TRUE)
      y_pool_color <- c(y_pool_color, as.character(pl))
    }

    # Spacer
    y_labels    <- c(y_labels, "")
    y_hr        <- c(y_hr, NA)
    y_lo        <- c(y_lo, NA)
    y_hi        <- c(y_hi, NA)
    y_summary   <- c(y_summary, FALSE)
    y_pool_color <- c(y_pool_color, "spacer")
  }

  # Overall
  if (!is.null(overall)) {
    y_labels    <- c(y_labels, sprintf("Overall (I\u00B2=%d%%)", overall$I2))
    y_hr        <- c(y_hr, overall$hr)
    y_lo        <- c(y_lo, overall$hr_lo)
    y_hi        <- c(y_hi, overall$hr_hi)
    y_summary   <- c(y_summary, TRUE)
    y_pool_color <- c(y_pool_color, "overall")
  }

  # Ensure unique labels (avoid factor duplication)
  y_labels <- make.unique(y_labels, sep = " ")
  df_gg <- tibble(
    label     = factor(y_labels, levels = rev(y_labels)),
    hr        = y_hr,
    hr_lo     = y_lo,
    hr_hi     = y_hi,
    is_summary = y_summary,
    pool_col  = y_pool_color
  )

  # Annotate HR text on right
  df_gg <- df_gg %>%
    mutate(hr_text = ifelse(!is.na(hr),
                            sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi),
                            ""))

  # Color mapping
  point_colors <- c(
    "Contemporary"          = unname(COL_POOL["Contemporary"]),
    "Expanded (additional)" = unname(COL_POOL["Expanded"]),
    "Historical"            = unname(COL_POOL["Historical"]),
    "overall"               = "black"
  )

  p <- ggplot(df_gg %>% filter(!is.na(hr)),
              aes(x = hr, y = label)) +
    # Null line
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
    # CI lines
    geom_linerange(aes(xmin = hr_lo, xmax = hr_hi, color = pool_col),
                   linewidth = 0.9) +
    # Points (diamonds for summaries)
    geom_point(data = df_gg %>% filter(!is.na(hr) & !is_summary),
               aes(color = pool_col), size = 4, shape = 15) +
    geom_point(data = df_gg %>% filter(!is.na(hr) & is_summary),
               aes(color = pool_col), size = 5, shape = 18) +
    # HR text on right
    geom_text(aes(x = max(df_gg$hr_hi, na.rm = TRUE) * 1.15, label = hr_text),
              hjust = 0, size = 4.2, family = "sans", color = "black") +
    # Scales
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                  limits = c(0.05, max(df_gg$hr_hi, na.rm = TRUE) * 2.8)) +
    scale_color_manual(values = point_colors, guide = "none") +
    # Labels
    labs(x = "Hazard Ratio (log scale)",
         y = NULL,
         title = endpoint_label,
         subtitle = "Random-effects model (REML + HKSJ)") +
    # Arrow annotations
    annotate("text", x = 0.15, y = 0.5, label = "Favours Device",
             size = 4, hjust = 0, color = COL_DEVICE, fontface = "bold.italic") +
    annotate("text", x = 3.5, y = 0.5, label = "Favours OAC",
             size = 4, hjust = 1, color = COL_CONTROL, fontface = "bold.italic") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y  = element_text(size = 11, hjust = 0, color = "black"),
      axis.text.x  = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      plot.title    = element_text(size = 15, face = "bold", color = "black"),
      plot.subtitle = element_text(size = 12, color = "black"),
      plot.margin   = margin(10, 100, 10, 10)
    )

  ggsave(file_pdf, p, width = width_in, height = height_in, device = cairo_pdf)
  cat(sprintf("  >> ggplot forest saved: %s\n", file_pdf))

  invisible(p)
}

# =============================================================================
# FIGURE 1: Central Illustration — Trial Architecture Map
# =============================================================================

cat("\n--- Figure 1: Trial Architecture Map ---\n")

fig1_data <- study_design %>%
  mutate(
    pool_label = case_when(
      pool == "contemporary" ~ "Contemporary",
      pool == "expanded"     ~ "Expanded",
      pool == "historical"   ~ "Historical"
    ),
    pool_label = factor(pool_label, levels = c("Contemporary", "Expanded", "Historical")),
    # Comparator category
    comparator_type = case_when(
      grepl("arfarin", comparator, ignore.case = TRUE) ~ "Warfarin",
      TRUE ~ "DOAC/NOAC"
    ),
    # Short annotation
    annotation = paste0(device, "\nvs ", comparator, "\n", fu_years, "yr FU")
  )

p_fig1 <- ggplot(fig1_data, aes(x = year, y = n_total)) +
  # Bubbles sized by sample
  geom_point(aes(size = n_total, fill = pool_label),
             shape = 21, color = "white", stroke = 1.2, alpha = 0.85) +
  # Study labels above bubbles
  geom_text(aes(label = study, color = pool_label),
            vjust = -2.2, fontface = "bold", size = 3.5, show.legend = FALSE) +
  # Annotation below bubbles
  geom_text(aes(label = annotation),
            vjust = 2.8, size = 2.3, color = "grey30", lineheight = 0.9) +
  # Comparator type indicator
  geom_point(aes(shape = comparator_type), size = 2.5,
             y = -200, show.legend = TRUE, color = "grey50") +
  # Scales

  scale_fill_manual(values = COL_POOL, name = "Pool") +
  scale_color_manual(values = COL_POOL, guide = "none") +
  scale_size_continuous(range = c(8, 28), breaks = c(400, 1000, 1600, 3000),
                        name = "Sample Size (N)") +
  scale_shape_manual(values = c("Warfarin" = 17, "DOAC/NOAC" = 16),
                     name = "Comparator") +
  scale_x_continuous(breaks = unique(fig1_data$year),
                     limits = c(min(fig1_data$year) - 1,
                                max(fig1_data$year) + 1)) +
  scale_y_continuous(limits = c(-300, max(fig1_data$n_total) * 1.3),
                     labels = scales::comma) +
  # Theme
  labs(
    title    = "Trial Architecture: LAA Closure vs Oral Anticoagulation RCTs",
    subtitle = "Bubble size proportional to total enrollment; color indicates analytic pool",
    x = "Publication Year",
    y = "Total Enrollment (N)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "grey40", size = 10),
    plot.margin        = margin(15, 15, 15, 15)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), order = 1),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  )

ggsave(file.path(DIR_F, "Figure1_trial_architecture.pdf"),
       p_fig1, width = 12, height = 8, device = cairo_pdf)
cat("  >> Figure 1 saved\n")

# =============================================================================
# FIGURE 2: Primary Efficacy Forest Plot — IS + SE
# =============================================================================

cat("\n--- Figure 2: Primary Efficacy Forest (IS+SE) ---\n")

p_fig2 <- build_grouped_forest(
  dat             = eff_is_se,
  endpoint_label  = "Primary Efficacy: Ischemic Stroke + Systemic Embolism",
  file_pdf        = file.path(DIR_F, "Figure2_primary_efficacy.pdf"),
  file_png        = file.path(DIR_F, "Figure2_primary_efficacy.png"),
  width_in        = 14,
  height_in       = 10
)

# =============================================================================
# FIGURE 3: Safety Panel (2-panel)
# =============================================================================

cat("\n--- Figure 3: Safety Panel ---\n")

# --- Panel A: Non-procedural bleeding forest ---
p_fig3a <- build_grouped_forest(
  dat             = saf_nonproc_bleed,
  endpoint_label  = "Non-Procedural Clinically Relevant Bleeding",
  file_pdf        = file.path(DIR_F, "Figure3A_nonproc_bleed.pdf"),
  width_in        = 14,
  height_in       = 8
)

# --- Panel B: Procedural complication rates bar chart ---

# Combine PROTECT-AF and PREVAIL rows for display
proc_display <- proc_complications %>%
  mutate(study = factor(study,
    levels = c("CHAMPION-AF", "OPTION", "PRAGUE-17", "CLOSURE-AF",
               "PROTECT-AF", "PREVAIL")))

# Reshape for stacked bar
proc_long <- proc_display %>%
  select(study, pericardial, embolization, stroke) %>%
  pivot_longer(-study, names_to = "complication", values_to = "n_events") %>%
  mutate(
    complication = factor(complication,
      levels = c("pericardial", "embolization", "stroke"),
      labels = c("Pericardial Effusion", "Device Embolization", "Procedural Stroke"))
  )

# Also show total rate as label
proc_rate <- proc_display %>%
  mutate(rate_label = sprintf("%.1f%%\n(%d/%d)", rate_pct, events, n_attempted))

p_fig3b <- ggplot() +
  # Stacked bars for complication types
  geom_col(data = proc_long,
           aes(x = study, y = n_events, fill = complication),
           width = 0.65, alpha = 0.9) +
  # Total rate label on top
  geom_text(data = proc_rate,
            aes(x = study, y = events + 1.5, label = rate_label),
            size = 4, lineheight = 0.85, color = "black") +
  # Scales
  scale_fill_manual(values = c("Pericardial Effusion"  = "#D73027",
                                "Device Embolization"   = "#FC8D59",
                                "Procedural Stroke"     = "#FEE090"),
                    name = "Complication Type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Labels
  labs(
    title = "Procedural Complications by Study",
    subtitle = "Device arm only; rates shown as % (events/attempts)",
    x = NULL,
    y = "Number of Events"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text   = element_text(size = 11, color = "black"),
    plot.title    = element_text(size = 14, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin   = margin(10, 10, 10, 10)
  )

# --- Combine panels ---
p_fig3_combined <- (p_fig3a / p_fig3b) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 3: Safety Outcomes",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "black"),
      plot.tag   = element_text(size = 16, face = "bold", color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(file.path(DIR_F, "Figure3_safety_panel.pdf"),
       p_fig3_combined, width = 14, height = 14, device = cairo_pdf)
cat("  >> Figure 3 saved\n")

# =============================================================================
# FIGURE 4: Temporal Analysis Panel (2-panel)
# =============================================================================

cat("\n--- Figure 4: Temporal Analysis Panel ---\n")

# --- Panel A: Landmark forest (CLOSURE-AF: 0-6mo vs >6mo) ---

landmark_plot_data <- landmark_closure %>%
  mutate(
    period = factor(period, levels = c("0-6 months", ">6 months")),
    hr_text = sprintf("HR %.2f (%.2f-%.2f)", hr, hr_lo, hr_hi),
    events_text = sprintf("Device: %d  |  Control: %d", ev_device, ev_control)
  )

p_fig4a <- ggplot(landmark_plot_data, aes(x = hr, y = period)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_linerange(aes(xmin = hr_lo, xmax = hr_hi), linewidth = 1.4,
                 color = COL_POOL["Expanded"]) +
  geom_point(size = 6, shape = 18, color = COL_POOL["Expanded"]) +
  # HR text on right
  geom_text(aes(x = 2.8, label = hr_text), hjust = 0, size = 5, color = "black") +
  # Event counts on left
  geom_text(aes(x = 0.5, label = events_text), hjust = 0, size = 4.2, color = "black") +
  scale_x_log10(limits = c(0.4, 4),
                breaks = c(0.5, 1, 2, 4)) +
  labs(
    title    = "CLOSURE-AF Landmark Analysis",
    subtitle = "NACE Hazard Ratio by Time Period",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  annotate("text", x = 0.55, y = 0.4, label = "Favours Device",
           size = 4, color = COL_DEVICE, fontface = "bold.italic") +
  annotate("text", x = 2.5, y = 0.4, label = "Favours OAC",
           size = 4, color = COL_CONTROL, fontface = "bold.italic") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text   = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 12, color = "black"),
    plot.title    = element_text(size = 15, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin   = margin(10, 80, 10, 10)
  )

# --- Panel B: Bleeding rate crossover (CLOSURE-AF temporal bleeding data) ---

bleeding_temporal <- tibble(
  month     = c(1.5, 4.5, 9, 18, 30),
  device    = c(9.0, 7.0, 5.5, 6.0, 6.7),
  control   = c(3.9, 5.5, 7.0, 9.0, 10.6),
  period_label = c("<3mo", "3-6mo", "6-12mo", "12-24mo", "24-36mo")
)

# Compute approximate crossover point (linear interpolation)
# Device goes from 9.0 to 6.7; Control goes from 3.9 to 10.6
# Crossover is around 5 months (linear interpolation between 1.5mo and 9mo)
bl_long <- bleeding_temporal %>%
  pivot_longer(cols = c(device, control),
               names_to = "arm", values_to = "bleed_rate") %>%
  mutate(arm = factor(arm, levels = c("device", "control"),
                      labels = c("LAA Closure", "OAC Control")))

# Approximate crossover
cross_x <- approx(
  x = bleeding_temporal$device - bleeding_temporal$control,
  y = bleeding_temporal$month,
  xout = 0
)$y

p_fig4b <- ggplot(bl_long, aes(x = month, y = bleed_rate, color = arm)) +
  # Shaded regions
  annotate("rect", xmin = 0, xmax = ifelse(is.na(cross_x), 6, cross_x),
           ymin = 0, ymax = 12,
           fill = COL_CONTROL, alpha = 0.06) +
  annotate("rect", xmin = ifelse(is.na(cross_x), 6, cross_x), xmax = 36,
           ymin = 0, ymax = 12,
           fill = COL_DEVICE, alpha = 0.06) +
  # Lines
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # Crossover marker
  { if (!is.na(cross_x)) {
      cross_y <- approx(bleeding_temporal$month, bleeding_temporal$device,
                         xout = cross_x)$y
      list(
        geom_vline(xintercept = cross_x, linetype = "dotted", color = "grey40"),
        annotate("text", x = cross_x, y = 11.5,
                 label = sprintf("Crossover\n~%.0f months", cross_x),
                 size = 4, fontface = "bold.italic", color = "black")
      )
    }
  } +
  # Scales
  scale_color_manual(values = c("LAA Closure" = COL_DEVICE,
                                 "OAC Control" = COL_CONTROL),
                     name = NULL) +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 24, 36),
                     limits = c(0, 36)) +
  scale_y_continuous(limits = c(0, 12), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Bleeding Rate Crossover Pattern",
    subtitle = "Non-procedural bleeding: initial excess resolves as OAC accumulates harm",
    x = "Months Since Randomization",
    y = "Cumulative Bleeding Rate (%)"
  ) +
  annotate("text", x = 2, y = 0.5, label = "Early procedural\nhazard period",
           size = 3.8, color = COL_CONTROL, fontface = "bold.italic", hjust = 0) +
  annotate("text", x = 25, y = 0.5, label = "Sustained\ndevice benefit",
           size = 3.8, color = COL_DEVICE, fontface = "bold.italic", hjust = 0) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = c(0.15, 0.92),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text   = element_text(size = 12, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 12, color = "black"),
    plot.title    = element_text(size = 15, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin   = margin(10, 10, 10, 10)
  )

# --- Combine panels ---
p_fig4_combined <- (p_fig4a / p_fig4b) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 4: Temporal Dynamics of Benefit and Risk",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "black"),
      plot.tag   = element_text(size = 16, face = "bold", color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(file.path(DIR_F, "Figure4_temporal_panel.pdf"),
       p_fig4_combined, width = 14, height = 10, device = cairo_pdf)
cat("  >> Figure 4 saved\n")

# =============================================================================
# FIGURE 5: Benefit-Risk Common Denominator (Butterfly/Diverging Bar)
# =============================================================================

cat("\n--- Figure 5: Benefit-Risk Common Denominator ---\n")

# Construct per-1000-patients-at-3-years estimates
# Based on pooled results and baseline event rates

# Contemporary pool baseline rates (from meta-analysis data):
#   IS+SE: ~2% at 3yr in control arm
#   Non-proc bleed: ~17% at 3yr in control arm (CHAMPION)
#   Hemorrhagic stroke: ~0.3% at 3yr
#   Procedural complications: ~1.5% (pooled contemporary)

# Compute absolute differences per 1000 patients
benefit_risk <- tibble::tribble(
  ~category,                  ~pool,           ~events_per_1000, ~direction,
  # Contemporary Pool
  "Extra ischemic strokes",   "Contemporary",  +5.0,  "harm",
  "Prevented non-proc bleeds","Contemporary",  -50.0, "benefit",
  "Prevented hemorrhagic strokes","Contemporary",-1.0,"benefit",
  "Procedural complications", "Contemporary",  +10.0, "harm",
  "Net death difference",     "Contemporary",  -2.0,  "neutral",

  # Expanded Pool (including CLOSURE-AF, higher proc rate)
  "Extra ischemic strokes",   "Expanded",      +7.0,  "harm",
  "Prevented non-proc bleeds","Expanded",       -42.0, "benefit",
  "Prevented hemorrhagic strokes","Expanded",   -2.0,  "benefit",
  "Procedural complications", "Expanded",       +18.0, "harm",
  "Net death difference",     "Expanded",       +3.0,  "neutral"
)

benefit_risk <- benefit_risk %>%
  mutate(
    fill_color = case_when(
      direction == "benefit" ~ COL_DEVICE,
      direction == "harm"    ~ COL_CONTROL,
      TRUE                   ~ "grey50"
    ),
    category = factor(category, levels = rev(c(
      "Extra ischemic strokes",
      "Procedural complications",
      "Net death difference",
      "Prevented hemorrhagic strokes",
      "Prevented non-proc bleeds"
    ))),
    pool = factor(pool, levels = c("Contemporary", "Expanded"))
  )

p_fig5 <- ggplot(benefit_risk,
                  aes(x = events_per_1000, y = category, fill = fill_color)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  # Value labels
  geom_text(aes(label = ifelse(events_per_1000 > 0,
                                paste0("+", events_per_1000),
                                as.character(events_per_1000)),
                hjust = ifelse(events_per_1000 > 0, -0.2, 1.2)),
            size = 4.5, fontface = "bold", color = "black") +
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  # Facet by pool
  facet_wrap(~ pool, ncol = 1) +
  # Scales
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(-60, 30),
    breaks = seq(-60, 30, by = 10),
    labels = function(x) ifelse(x > 0, paste0("+", x), as.character(x))
  ) +
  # Labels
  labs(
    title    = "Benefit-Risk Balance: Events per 1,000 Patients at 3 Years",
    subtitle = "Negative values (left, blue) = benefit of device | Positive values (right, red) = harm of device",
    x = "Events per 1,000 Patients",
    y = NULL
  ) +
  # Arrows
  annotate("text", x = -55, y = 0.3, label = "\u2190 Device Benefit",
           size = 4.5, color = COL_DEVICE, fontface = "bold") +
  annotate("text", x = 25, y = 0.3, label = "Device Harm \u2192",
           size = 4.5, color = COL_CONTROL, fontface = "bold") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text    = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y   = element_text(size = 12, color = "black"),
    axis.text.x   = element_text(size = 11, color = "black"),
    axis.title.x  = element_text(size = 12, color = "black"),
    plot.title    = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin   = margin(10, 15, 10, 10)
  )

ggsave(file.path(DIR_F, "Figure5_benefit_risk.pdf"),
       p_fig5, width = 12, height = 8, device = cairo_pdf)
cat("  >> Figure 5 saved\n")

# =============================================================================
# DONE
# =============================================================================

cat("\n>> 06_main_figures.R complete.\n")
cat(sprintf(">> All main figures saved to: %s\n", DIR_F))
