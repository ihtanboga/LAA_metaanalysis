###############################################################################
# 07_suppl_figures.R
# LAA Closure vs OAC Meta-Analysis — Supplementary Figures (6 figures)
#
# sFigure 1: PRISMA Flow Diagram
# sFigure 2: Secondary Endpoints Forest Plots (6-panel)
# sFigure 3: Sensitivity Analysis (Leave-One-Out)
# sFigure 4: NI Margin Heatmap
# sFigure 5: CLOSURE-AF Endpoint Cascade
# sFigure 6: Subgroup Interaction Forest
#
# Dependencies: 01_data.R, output/primary_results.RData
# Outputs: suppl_figures/sFigure1–6_*.pdf
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
cat("  07_suppl_figures.R\n")
cat("========================================\n\n")

# =============================================================================
# HELPER: Minimal ggplot-based forest plot for a single endpoint
# =============================================================================

#' Generate a compact forest plot panel for one secondary endpoint
#' @param dat Data frame with study, hr, hr_lo, hr_hi, ev_device, n_device,
#'            ev_control, n_control, loghr, se_loghr
#' @param title_text Title for the panel
#' @return A ggplot object
mini_forest <- function(dat, title_text) {

  # Filter to studies with non-missing HR
  d <- dat %>% filter(!is.na(hr) & !is.na(hr_lo) & !is.na(hr_hi))

  if (nrow(d) == 0) {
    # Return a placeholder plot
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(title_text, "\n(No data)"),
                 size = 4, color = "grey50") +
        theme_void() +
        labs(title = title_text)
    )
  }

  # Pool with metagen if >= 2 studies
  pooled <- NULL
  if (nrow(d %>% filter(!is.na(loghr) & !is.na(se_loghr))) >= 2) {
    d_m <- d %>% filter(!is.na(loghr) & !is.na(se_loghr))
    m <- metagen(TE = d_m$loghr, seTE = d_m$se_loghr,
                 studlab = d_m$study, sm = "HR",
                 method.tau = "REML", hakn = TRUE)
    pooled <- tibble(
      study = sprintf("Pooled (I\u00B2=%.0f%%)", m$I2 * 100),
      hr    = exp(m$TE.random),
      hr_lo = exp(m$lower.random),
      hr_hi = exp(m$upper.random),
      is_pooled = TRUE
    )
  }

  # Build plot data
  plot_d <- d %>%
    mutate(
      label = sprintf("%s  [%d/%d vs %d/%d]",
                      study, ev_device, n_device, ev_control, n_control),
      is_pooled = FALSE
    ) %>%
    select(study, label, hr, hr_lo, hr_hi, is_pooled)

  if (!is.null(pooled)) {
    pooled$label <- pooled$study
    plot_d <- bind_rows(plot_d, pooled)
  }

  # Spacer between studies and pooled
  if (!is.null(pooled)) {
    spacer <- tibble(study = " ", label = " ", hr = NA, hr_lo = NA, hr_hi = NA,
                     is_pooled = FALSE)
    n_studies <- sum(!plot_d$is_pooled & plot_d$study != " ")
    plot_d <- bind_rows(
      plot_d %>% filter(!is_pooled),
      spacer,
      plot_d %>% filter(is_pooled)
    )
  }

  plot_d$label <- factor(plot_d$label, levels = rev(plot_d$label))

  # HR text
  plot_d <- plot_d %>%
    mutate(hr_text = ifelse(!is.na(hr),
                            sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi), ""))

  # Compute x range — use log scale so we need generous right space for text
  x_max <- max(c(plot_d$hr_hi, 4), na.rm = TRUE)

  p <- ggplot(plot_d %>% filter(!is.na(hr)),
              aes(x = hr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_linerange(aes(xmin = hr_lo, xmax = hr_hi),
                   linewidth = 0.8, color = "black") +
    # Study points
    geom_point(data = . %>% filter(!is_pooled),
               size = 3.5, shape = 15, color = COL_DEVICE) +
    # Pooled diamond
    geom_point(data = . %>% filter(is_pooled),
               size = 4.5, shape = 18, color = "black") +
    # HR text — placed outside plot area via clip="off"
    geom_text(aes(x = x_max * 1.15, label = hr_text),
              hjust = 0, size = 3.5, family = "sans", color = "black") +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 4),
      limits = c(0.05, x_max * 1.3)
    ) +
    coord_cartesian(clip = "off") +
    labs(title = title_text, x = "HR (log scale)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y  = element_text(size = 9, hjust = 0, color = "black"),
      axis.text.x  = element_text(size = 9, color = "black"),
      axis.title   = element_text(size = 10, color = "black"),
      plot.title   = element_text(size = 11, face = "bold", color = "black"),
      plot.margin  = margin(5, 120, 5, 5)
    )

  return(p)
}

# =============================================================================
# sFIGURE 1: PRISMA Flow Diagram
# =============================================================================

cat("\n--- sFigure 1: PRISMA Flow Diagram ---\n")

# Create a simplified PRISMA-style flowchart using ggplot

prisma_data <- tibble(
  box = c(
    "Records identified\nthrough database searching\n(n = 2,847)",
    "Additional records from\ntrial registries & references\n(n = 42)",
    "Records after\nduplicates removed\n(n = 2,104)",
    "Records screened\n(n = 2,104)",
    "Records excluded\n(n = 2,041)",
    "Full-text articles\nassessed for eligibility\n(n = 63)",
    "Full-text articles excluded (n = 57)\n- Non-randomized (n = 23)\n- Not LAA closure vs OAC (n = 18)\n- Duplicate cohort (n = 9)\n- Ongoing, no results (n = 7)",
    "Studies included in\nquantitative synthesis\n(n = 6 [5 entries])"
  ),
  x = c(2, 6, 4, 4, 7, 4, 7.5, 4),
  y = c(9, 9, 7.5, 6, 6, 4, 4, 2),
  w = c(3.0, 3.0, 3.0, 3.0, 2.5, 3.0, 3.5, 3.0),
  h = c(1.2, 1.2, 0.9, 0.9, 0.7, 0.9, 1.5, 0.9),
  fill_col = c("white", "white", "#E8F4FD", "#E8F4FD", "#FEE0D2",
                "#E8F4FD", "#FEE0D2", "#D4EDDA")
)

# Arrows
arrows <- tibble(
  x_start = c(2, 6, 4, 4, 4, 4),
  y_start = c(8.4, 8.4, 7.05, 5.55, 5.55, 3.55),
  x_end   = c(4, 4, 4, 4, 6.5, 4),
  y_end   = c(7.95, 7.95, 6.45, 4.45, 4.0, 2.45)
)

# Exclusion arrows (horizontal)
arrows_excl <- tibble(
  x_start = c(5.5, 5.5),
  y_start = c(6, 4),
  x_end   = c(5.8, 5.8),
  y_end   = c(6, 4)
)

p_sfig1 <- ggplot() +
  # Boxes
  geom_rect(data = prisma_data,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2,
                fill = fill_col),
            color = "grey40", linewidth = 0.5) +
  scale_fill_identity() +
  # Text in boxes
  geom_text(data = prisma_data,
            aes(x = x, y = y, label = box),
            size = 2.6, lineheight = 0.85) +
  # Arrows
  geom_segment(data = arrows,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.4) +
  # Horizontal exclusion connectors
  geom_segment(aes(x = 5.5, y = 6, xend = 5.75, yend = 6),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.4) +
  geom_segment(aes(x = 5.5, y = 4, xend = 5.75, yend = 4),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.4) +
  # Title
  labs(title = "sFigure 1: PRISMA Flow Diagram") +
  coord_fixed(ratio = 0.6, xlim = c(0, 10), ylim = c(1, 10.5)) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(DIR_SF, "sFigure1_prisma_flow.pdf"),
       p_sfig1, width = 10, height = 12, device = cairo_pdf)
cat("  >> sFigure 1 saved\n")

# =============================================================================
# sFIGURE 2: Secondary Endpoints Forest Plots (6-panel, 3x2 grid)
# =============================================================================

cat("\n--- sFigure 2: Secondary Endpoints Forest (6-panel) ---\n")

p_acm       <- mini_forest(sec_acm,          "All-Cause Mortality")
p_cvdeath   <- mini_forest(sec_cvdeath,      "Cardiovascular / Unexplained Death")
p_hemstroke <- mini_forest(sec_hemstroke,     "Hemorrhagic Stroke")
p_anystroke <- mini_forest(sec_anystroke_se,  "Any Stroke + Systemic Embolism")
p_is_alone  <- mini_forest(sec_is_alone,     "Ischemic Stroke Alone")
p_all_bleed <- mini_forest(saf_all_major_bleed, "All Major Bleeding (incl. procedural)")

p_sfig2 <- (p_acm | p_cvdeath) /
            (p_hemstroke | p_anystroke) /
            (p_is_alone | p_all_bleed) +
  plot_annotation(
    title = "sFigure 2: Secondary Endpoint Forest Plots",
    subtitle = "Random-effects (REML + HKSJ) — all available trials",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "black"),
      plot.subtitle = element_text(size = 12, color = "black"),
      plot.tag = element_text(size = 14, face = "bold", color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(file.path(DIR_SF, "sFigure2_secondary_forest.pdf"),
       p_sfig2, width = 16, height = 18, device = cairo_pdf)
cat("  >> sFigure 2 saved\n")

# =============================================================================
# sFIGURE 3: Sensitivity Analysis — Leave-One-Out
# =============================================================================

cat("\n--- sFigure 3: Sensitivity Analysis (Leave-One-Out) ---\n")

# Leave-one-out for primary efficacy (IS+SE) using HR
loo_data <- eff_is_se %>% filter(!is.na(loghr) & !is.na(se_loghr))
n_studies <- nrow(loo_data)

loo_results <- map_dfr(seq_len(n_studies), function(i) {
  d_loo <- loo_data[-i, ]
  if (nrow(d_loo) < 2) return(NULL)

  m <- metagen(TE = d_loo$loghr, seTE = d_loo$se_loghr,
               studlab = d_loo$study, sm = "HR",
               method.tau = "REML", hakn = TRUE)
  tibble(
    dropped  = loo_data$study[i],
    hr       = exp(m$TE.random),
    hr_lo    = exp(m$lower.random),
    hr_hi    = exp(m$upper.random),
    I2       = m$I2 * 100,
    tau2     = m$tau2
  )
})

# Also add the overall (no study dropped)
if (n_studies >= 2) {
  m_full <- metagen(TE = loo_data$loghr, seTE = loo_data$se_loghr,
                    studlab = loo_data$study, sm = "HR",
                    method.tau = "REML", hakn = TRUE)
  loo_overall <- tibble(
    dropped = "None (Overall)",
    hr      = exp(m_full$TE.random),
    hr_lo   = exp(m_full$lower.random),
    hr_hi   = exp(m_full$upper.random),
    I2      = m_full$I2 * 100,
    tau2    = m_full$tau2
  )
  loo_results <- bind_rows(loo_results, loo_overall)
}

loo_results <- loo_results %>%
  mutate(
    is_overall = (dropped == "None (Overall)"),
    label = factor(dropped, levels = rev(dropped)),
    hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi),
    I2_text = sprintf("I\u00B2=%.0f%%", I2)
  )

# Leave-one-out for non-proc bleeding
loo_bleed_data <- saf_nonproc_bleed %>% filter(!is.na(loghr) & !is.na(se_loghr))

loo_bleed <- map_dfr(seq_len(nrow(loo_bleed_data)), function(i) {
  d_loo <- loo_bleed_data[-i, ]
  if (nrow(d_loo) < 2) return(NULL)
  m <- metagen(TE = d_loo$loghr, seTE = d_loo$se_loghr,
               studlab = d_loo$study, sm = "HR",
               method.tau = "REML", hakn = TRUE)
  tibble(
    dropped = loo_bleed_data$study[i],
    hr      = exp(m$TE.random),
    hr_lo   = exp(m$lower.random),
    hr_hi   = exp(m$upper.random),
    I2      = m$I2 * 100
  )
})

if (nrow(loo_bleed_data) >= 2) {
  m_bl <- metagen(TE = loo_bleed_data$loghr, seTE = loo_bleed_data$se_loghr,
                  studlab = loo_bleed_data$study, sm = "HR",
                  method.tau = "REML", hakn = TRUE)
  loo_bleed <- bind_rows(loo_bleed, tibble(
    dropped = "None (Overall)", hr = exp(m_bl$TE.random),
    hr_lo = exp(m_bl$lower.random), hr_hi = exp(m_bl$upper.random),
    I2 = m_bl$I2 * 100
  ))
}

loo_bleed <- loo_bleed %>%
  mutate(
    is_overall = (dropped == "None (Overall)"),
    label = factor(dropped, levels = rev(dropped)),
    hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi)
  )

# --- Panel A: IS+SE LOO ---
p_loo_eff <- ggplot(loo_results %>% filter(!is.na(hr)),
                     aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_linerange(aes(xmin = hr_lo, xmax = hr_hi,
                      color = is_overall), linewidth = 0.9) +
  geom_point(aes(shape = is_overall, color = is_overall), size = 4) +
  geom_text(aes(x = max(loo_results$hr_hi, na.rm = TRUE) * 1.15,
                label = paste(hr_text, I2_text, sep = "  ")),
            hjust = 0, size = 4, color = "black") +
  scale_x_log10(breaks = c(0.5, 1, 2), limits = c(0.3, 5)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = COL_POOL["Expanded"]),
                     guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 15), guide = "none") +
  labs(
    title = "Ischemic Stroke + SE: Leave-One-Out Analysis",
    subtitle = "Study dropped shown on y-axis; pooled HR of remaining studies",
    x = "Pooled HR (log scale)", y = "Study Dropped"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    plot.margin = margin(10, 90, 10, 10)
  )

# --- Panel B: Non-proc bleeding LOO ---
p_loo_bleed <- ggplot(loo_bleed %>% filter(!is.na(hr)),
                       aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_linerange(aes(xmin = hr_lo, xmax = hr_hi,
                      color = is_overall), linewidth = 0.9) +
  geom_point(aes(shape = is_overall, color = is_overall), size = 4) +
  geom_text(aes(x = max(loo_bleed$hr_hi, na.rm = TRUE) * 1.15,
                label = hr_text),
            hjust = 0, size = 4, color = "black") +
  scale_x_log10(breaks = c(0.25, 0.5, 1), limits = c(0.2, 1.8)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = COL_DEVICE),
                     guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 15), guide = "none") +
  labs(
    title = "Non-Procedural Bleeding: Leave-One-Out Analysis",
    x = "Pooled HR (log scale)", y = "Study Dropped"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    plot.margin = margin(10, 90, 10, 10)
  )

p_sfig3 <- (p_loo_eff / p_loo_bleed) +
  plot_annotation(
    tag_levels = "A",
    title = "sFigure 3: Sensitivity Analysis — Leave-One-Out",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "black"),
      plot.tag = element_text(size = 16, face = "bold", color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(file.path(DIR_SF, "sFigure3_sensitivity.pdf"),
       p_sfig3, width = 12, height = 12, device = cairo_pdf)
cat("  >> sFigure 3 saved\n")

# =============================================================================
# sFIGURE 4: NI Margin Heatmap
# =============================================================================

cat("\n--- sFigure 4: NI Margin Heatmap ---\n")

# Reshape ni_margin_data for heatmap
ni_long <- ni_margin_data %>%
  pivot_longer(-study, names_to = "threshold", values_to = "compatible") %>%
  mutate(
    threshold_label = case_when(
      threshold == "rr_115" ~ "RR < 1.15",
      threshold == "rr_125" ~ "RR < 1.25",
      threshold == "rr_140" ~ "RR < 1.40",
      threshold == "abs_05" ~ "Abs < 0.5%",
      threshold == "abs_10" ~ "Abs < 1.0%",
      threshold == "abs_20" ~ "Abs < 2.0%"
    ),
    threshold_label = factor(threshold_label,
      levels = c("RR < 1.15", "RR < 1.25", "RR < 1.40",
                 "Abs < 0.5%", "Abs < 1.0%", "Abs < 2.0%")),
    study = factor(study, levels = c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                                     "CLOSURE-AF", "PROTECT-AF/PREVAIL")),
    fill_val = case_when(
      compatible  ~ "Compatible",
      !compatible ~ "Not Compatible"
    )
  )

# Add pooled row estimates (approximate)
# Contemporary pool: HR ~1.28 (0.82-2.00) for IS+SE → check each threshold
pooled_ni <- tibble::tribble(
  ~study,    ~threshold,  ~compatible,
  "Pooled (Contemporary)", "rr_115", FALSE,
  "Pooled (Contemporary)", "rr_125", FALSE,
  "Pooled (Contemporary)", "rr_140", TRUE,
  "Pooled (Contemporary)", "abs_05", FALSE,
  "Pooled (Contemporary)", "abs_10", TRUE,
  "Pooled (Contemporary)", "abs_20", TRUE,
  "Pooled (Expanded)",     "rr_115", FALSE,
  "Pooled (Expanded)",     "rr_125", FALSE,
  "Pooled (Expanded)",     "rr_140", TRUE,
  "Pooled (Expanded)",     "abs_05", FALSE,
  "Pooled (Expanded)",     "abs_10", FALSE,
  "Pooled (Expanded)",     "abs_20", TRUE
) %>%
  mutate(
    threshold_label = case_when(
      threshold == "rr_115" ~ "RR < 1.15",
      threshold == "rr_125" ~ "RR < 1.25",
      threshold == "rr_140" ~ "RR < 1.40",
      threshold == "abs_05" ~ "Abs < 0.5%",
      threshold == "abs_10" ~ "Abs < 1.0%",
      threshold == "abs_20" ~ "Abs < 2.0%"
    ),
    threshold_label = factor(threshold_label,
      levels = c("RR < 1.15", "RR < 1.25", "RR < 1.40",
                 "Abs < 0.5%", "Abs < 1.0%", "Abs < 2.0%")),
    study = factor(study, levels = c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                                     "CLOSURE-AF", "PROTECT-AF/PREVAIL",
                                     "Pooled (Contemporary)", "Pooled (Expanded)")),
    fill_val = ifelse(compatible, "Compatible", "Not Compatible")
  )

ni_all <- bind_rows(
  ni_long %>% select(study, threshold_label, fill_val, compatible),
  pooled_ni %>% select(study, threshold_label, fill_val, compatible)
) %>%
  mutate(study = factor(study,
    levels = rev(c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                   "CLOSURE-AF", "PROTECT-AF/PREVAIL",
                   "Pooled (Contemporary)", "Pooled (Expanded)"))))

p_sfig4 <- ggplot(ni_all, aes(x = threshold_label, y = study, fill = fill_val)) +
  geom_tile(color = "white", linewidth = 1.5) +
  # Tick/cross symbols
  geom_text(aes(label = ifelse(compatible, "\u2713", "\u2717")),
            size = 5, fontface = "bold",
            color = ifelse(ni_all$compatible, "darkgreen", "darkred")) +
  scale_fill_manual(
    values = c("Compatible" = "#B8E186", "Not Compatible" = "#F4A582"),
    name = "NI Status"
  ) +
  # Separator line between individual and pooled
  geom_hline(yintercept = 2.5, linewidth = 1, color = "grey30") +
  labs(
    title = "sFigure 4: Non-Inferiority Margin Compatibility",
    subtitle = "Ischemic stroke + SE: can the 95% CI upper bound exclude each NI threshold?",
    x = "Non-Inferiority Threshold",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(DIR_SF, "sFigure4_ni_heatmap.pdf"),
       p_sfig4, width = 10, height = 8, device = cairo_pdf)
cat("  >> sFigure 4 saved\n")

# =============================================================================
# sFIGURE 5: CLOSURE-AF Endpoint Cascade
# =============================================================================

cat("\n--- sFigure 5: CLOSURE-AF Endpoint Cascade ---\n")

cascade_data <- closure_cascade %>%
  mutate(
    label = factor(model,
      levels = rev(c("Published NACE", "Reconstructed NACE",
                     "Model 1: MACE", "Model 2: Thromboembolic"))),
    # Color by whether HR crosses 1
    favours = ifelse(hr_mc > 1, "Favours OAC", "Neutral"),
    hr_text = sprintf("HR %.2f | dRMST %.3f yr", hr_mc, drmst),
    events_text = sprintf("Device: %d  |  Control: %d", ev_device, ev_ctrl)
  )

p_sfig5 <- ggplot(cascade_data, aes(x = hr_mc, y = label)) +
  # Reference line at 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  # Bars
  geom_col(aes(fill = favours), width = 0.55, alpha = 0.85) +
  # HR text
  geom_text(aes(x = max(hr_mc) + 0.08, label = hr_text),
            hjust = 0, size = 4.5, fontface = "bold", color = "black") +
  # Event counts
  geom_text(aes(x = 0.8, label = events_text),
            hjust = 1, size = 4, color = "black") +
  # Scales
  scale_fill_manual(values = c("Favours OAC" = "#FC8D59", "Neutral" = "#91BFDB"),
                    guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.85), breaks = seq(0.6, 1.6, 0.2)) +
  labs(
    title = "sFigure 5: CLOSURE-AF Endpoint Cascade",
    subtitle = "Progressive attenuation of harm signal as bleeding components are removed",
    x = "Hazard Ratio (Reconstructed)",
    y = NULL
  ) +
  # Annotations
  annotate("segment", x = 1.05, xend = 1.4, y = 0.3, yend = 0.3,
           arrow = arrow(length = unit(0.15, "cm")),
           color = COL_CONTROL, linewidth = 0.6) +
  annotate("text", x = 1.42, y = 0.3, label = "Favours OAC",
           size = 4, color = COL_CONTROL, hjust = 0, fontface = "bold") +
  annotate("segment", x = 0.95, xend = 0.7, y = 0.3, yend = 0.3,
           arrow = arrow(length = unit(0.15, "cm")),
           color = COL_DEVICE, linewidth = 0.6) +
  annotate("text", x = 0.68, y = 0.3, label = "Favours Device",
           size = 4, color = COL_DEVICE, hjust = 1, fontface = "bold") +
  # Arrow showing attenuation
  annotate("segment", x = 1.5, y = 4, xend = 1.1, yend = 1,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           color = "black", linewidth = 0.7, linetype = "dotted") +
  annotate("text", x = 1.52, y = 2.5,
           label = "Progressive\nattenuation",
           size = 4, color = "black", fontface = "bold.italic", angle = -60) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title  = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin = margin(10, 120, 10, 10)
  )

ggsave(file.path(DIR_SF, "sFigure5_closure_cascade.pdf"),
       p_sfig5, width = 12, height = 6, device = cairo_pdf)
cat("  >> sFigure 5 saved\n")

# =============================================================================
# sFIGURE 6: Subgroup Interaction Forest
# =============================================================================

cat("\n--- sFigure 6: Subgroup Interaction Forest ---\n")

# Combine CHAMPION-AF and OPTION subgroup data
subgroup_all <- bind_rows(
  subgroup_champion_safety %>% mutate(trial = "CHAMPION-AF"),
  subgroup_option_safety   %>% mutate(trial = "OPTION")
) %>%
  mutate(
    # Y-axis label: "Trial — Category"
    row_label = paste0(trial, ":  ", category),
    hr_text   = sprintf("%.2f (%.2f\u2013%.2f)", hr, hr_lo, hr_hi),
    # Subgroup as facet (clean names)
    subgroup  = factor(subgroup,
                       levels = c("Age", "Sex", "HAS-BLED",
                                  "CHA2DS2-VASc", "Prior embolic"))
  )

# Within each facet, order by category then trial
subgroup_all <- subgroup_all %>%
  arrange(subgroup, category, trial) %>%
  mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

p_sfig6 <- ggplot(subgroup_all, aes(x = hr, y = row_label, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = hr_lo, xmax = hr_hi), linewidth = 1) +
  geom_point(size = 4, shape = 15) +
  # HR text on the right
  geom_text(aes(label = hr_text),
            x = log10(1.6), hjust = 0, size = 3.8,
            color = "black", show.legend = FALSE) +
  # Facet by subgroup — free y so each panel has only its rows
  facet_wrap(~ subgroup, ncol = 1, scales = "free_y",
             strip.position = "left") +
  scale_color_manual(
    values = c("CHAMPION-AF" = unname(COL_POOL["Contemporary"]),
               "OPTION"      = unname(COL_POOL["Expanded"])),
    name = "Trial"
  ) +
  scale_x_log10(breaks = c(0.25, 0.5, 1),
                limits = c(0.15, 2.2)) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "sFigure 6: Subgroup Analysis — Non-Procedural Bleeding",
    subtitle = "HR for LAAC vs OAC by pre-specified subgroups",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.title   = element_text(size = 12, color = "black"),
    strip.text.y.left = element_text(size = 12, face = "bold", color = "black",
                                     angle = 0, hjust = 1),
    strip.placement = "outside",
    legend.position = "top",
    legend.text  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    plot.title    = element_text(size = 15, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    plot.margin   = margin(10, 60, 10, 10),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(file.path(DIR_SF, "sFigure6_subgroup.pdf"),
       p_sfig6, width = 12, height = 10, device = cairo_pdf)
cat("  >> sFigure 6 saved\n")

# =============================================================================
# DONE
# =============================================================================

cat("\n>> 07_suppl_figures.R complete.\n")
cat(sprintf(">> All supplementary figures saved to: %s\n", DIR_SF))
