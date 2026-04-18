# ============================================================
# Forest plot: Stroke — device vs AC across time bands
#   Panel 1: Contemporary (CHAMPION + OPTION + PRAGUE-17)
#   Panel 2: Expanded    (+ CLOSURE-AF)
# Bands:
#   A: 0–7 days (peri-procedural)
#   B: 7–180 days (early)
#   C: > 180 days (late)
# Pooled HR: 1-stage meta-analysis via stratified Cox
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

ipd_dir <- "~/Desktop/manuel/json_csv/ipd_data"

# ===== 1) Load & harmonize =====
load_trial <- function(path, trial, ac_lab, dev_lab) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$trial <- trial
  df$arm_h <- ifelse(df$arm == ac_lab, "AC",
              ifelse(df$arm == dev_lab, "device", NA))
  df <- df[!is.na(df$arm_h), ]
  df$arm_h <- factor(df$arm_h, levels = c("AC", "device"))
  df[, c("trial", "arm_h", "time", "status")]
}

ch <- load_trial(file.path(ipd_dir, "CHAMPION_all_stroke.csv"),
                 "CHAMPION", "AC", "device")
op <- load_trial(file.path(ipd_dir, "OPTION_stroke_se.csv"),
                 "OPTION", "AC", "device")
pr <- load_trial(file.path(ipd_dir, "PRAGUE17_all_stroke_TIA.csv"),
                 "PRAGUE-17", "AC", "device")
cl <- load_trial(file.path(ipd_dir, "CLOSURE_stroke_se.csv"),
                 "CLOSURE-AF", "Medical Therapy", "LAA Closure")

cat("Data loaded:\n")
for (d in list(ch, op, pr, cl)) {
  cat(sprintf("  %-12s N=%d (AC %d, device %d), events=%d\n",
              d$trial[1], nrow(d),
              sum(d$arm_h == "AC"), sum(d$arm_h == "device"),
              sum(d$status == 1)))
}

# ===== 2) Band filter =====
d180 <- 180 / 365.25

apply_band <- function(df, band) {
  if (band == "A") {
    df2 <- df
    df2$status <- ifelse(df$time <= d180, df$status, 0)
    df2$time   <- pmin(df$time, d180)
  } else {
    df2 <- df[df$time > d180, ]
    df2$time   <- df2$time - d180
  }
  df2
}

# ===== 3) Cox helpers =====
run_cox_single <- function(df, label) {
  cx <- tryCatch(coxph(Surv(time, status) ~ arm_h, data = df),
                 error = function(e) NULL)
  ev_ac  <- sum(df$arm_h == "AC"     & df$status == 1)
  ev_dev <- sum(df$arm_h == "device" & df$status == 1)
  n_ac   <- sum(df$arm_h == "AC")
  n_dev  <- sum(df$arm_h == "device")
  if (is.null(cx) || ev_ac + ev_dev < 2) {
    return(data.frame(label = label, hr = NA, lo = NA, hi = NA,
                      n_ac = n_ac, n_dev = n_dev,
                      ev_ac = ev_ac, ev_dev = ev_dev))
  }
  data.frame(
    label = label,
    hr = as.numeric(exp(coef(cx))[1]),
    lo = as.numeric(exp(confint(cx))[1, 1]),
    hi = as.numeric(exp(confint(cx))[1, 2]),
    n_ac = n_ac, n_dev = n_dev,
    ev_ac = ev_ac, ev_dev = ev_dev
  )
}

run_cox_pooled <- function(df_list, label) {
  df_all <- do.call(rbind, df_list)
  cx <- tryCatch(coxph(Surv(time, status) ~ arm_h + strata(trial),
                        data = df_all),
                 error = function(e) NULL)
  ev_ac  <- sum(df_all$arm_h == "AC"     & df_all$status == 1)
  ev_dev <- sum(df_all$arm_h == "device" & df_all$status == 1)
  if (is.null(cx)) {
    return(data.frame(label = label, hr = NA, lo = NA, hi = NA,
                      n_ac = NA, n_dev = NA,
                      ev_ac = ev_ac, ev_dev = ev_dev, pooled = TRUE))
  }
  data.frame(
    label = label,
    hr = as.numeric(exp(coef(cx))[1]),
    lo = as.numeric(exp(confint(cx))[1, 1]),
    hi = as.numeric(exp(confint(cx))[1, 2]),
    n_ac = sum(df_all$arm_h == "AC"),
    n_dev = sum(df_all$arm_h == "device"),
    ev_ac = ev_ac, ev_dev = ev_dev, pooled = TRUE
  )
}

# ===== 4) Build results per panel × band =====
bands <- c("A" = "A. ≤180 days",
           "B" = "B. >180 days")

build_panel <- function(trial_list, panel_label) {
  rows <- list()
  for (b in names(bands)) {
    trial_dfs <- lapply(trial_list, apply_band, band = b)
    for (i in seq_along(trial_dfs)) {
      r <- run_cox_single(trial_dfs[[i]], label = trial_dfs[[i]]$trial[1])
      r$band <- bands[[b]]
      r$pooled <- FALSE
      rows[[length(rows) + 1]] <- r
    }
    rp <- run_cox_pooled(trial_dfs, label = "Pooled (stratified Cox)")
    rp$band <- bands[[b]]
    rows[[length(rows) + 1]] <- rp
  }
  do.call(rbind, rows) %>% mutate(panel = panel_label)
}

contemp  <- build_panel(list(ch, op, pr),         "Contemporary trials")
expanded <- build_panel(list(ch, op, pr, cl),     "Expanded contemporary trials")

plot_df <- rbind(contemp, expanded)
plot_df$band <- factor(plot_df$band, levels = bands)

cat("\n=== Results ===\n")
print(plot_df %>% select(panel, band, label, hr, lo, hi, ev_ac, ev_dev),
      row.names = FALSE)

# ===== 5) Plot =====
# Build ordered labels per panel for y-axis stacking
make_plot <- function(df, panel_title) {
  df <- df %>% filter(panel == panel_title)
  # Order within each band: trials alphabetical, pooled last
  df <- df %>% arrange(band, pooled, label)
  # Build unique row order
  df$row_id <- seq_len(nrow(df))
  df$y_label <- ifelse(df$pooled, "Pooled", df$label)
  df$shape_v <- ifelse(df$pooled, "diamond", "square")

  # Add band header rows by faceting
  # Use ggplot with facet by band (free y)
  # Cap HRs and CIs to plotting window, flag extremes with arrows
  xlo <- 0.1; xhi <- 10
  df$plot_hr <- pmin(pmax(df$hr, xlo), xhi)
  df$plot_lo <- pmin(pmax(df$lo, xlo), xhi)
  df$plot_hi <- pmin(pmax(df$hi, xlo), xhi)
  df$hr_label <- ifelse(
    is.na(df$hr) | df$ev_ac + df$ev_dev < 2,
    sprintf("— (insufficient events)  [%d vs %d]", df$ev_dev, df$ev_ac),
    ifelse(df$hi > xhi | df$lo < xlo | df$hr > xhi | df$hr < xlo,
           sprintf("%.2f (%.2f–%.2f)*  [%d vs %d]",
                   df$hr, df$lo, df$hi, df$ev_dev, df$ev_ac),
           sprintf("%.2f (%.2f–%.2f)  [%d vs %d]",
                   df$hr, df$lo, df$hi, df$ev_dev, df$ev_ac))
  )

  ggplot(df, aes(x = plot_hr, y = reorder(y_label, -row_id))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = plot_lo, xmax = plot_hi),
                  orientation = "y", width = 0.25, linewidth = 0.5) +
    geom_point(aes(shape = shape_v, size = shape_v,
                   fill = shape_v, color = shape_v)) +
    geom_text(aes(x = xhi * 1.4, label = hr_label),
              hjust = 0, size = 2.9, family = "mono") +
    scale_shape_manual(values = c("diamond" = 23, "square" = 15), guide = "none") +
    scale_size_manual(values = c("diamond" = 3.5, "square" = 2.5), guide = "none") +
    scale_fill_manual(values = c("diamond" = "black", "square" = "black"), guide = "none") +
    scale_color_manual(values = c("diamond" = "black", "square" = "black"), guide = "none") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                       limits = c(xlo, xhi * 12)) +
    coord_cartesian(clip = "off") +
    facet_wrap(~ band, ncol = 1, scales = "free_y", strip.position = "top") +
    labs(x = "Hazard ratio (device vs AC)  — log scale",
         y = NULL, title = panel_title) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", hjust = 0),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.y = element_text(size = 10),
      panel.spacing.y = unit(0.6, "lines"),
      plot.margin = margin(5, 130, 5, 5)
    )
}

p1 <- make_plot(plot_df, "Contemporary trials")
p2 <- make_plot(plot_df, "Expanded contemporary trials")

combined <- p1 | p2

ggsave("~/Desktop/manuel/JACC_figures/figure_forest_stroke.png",
       combined, width = 14, height = 7, dpi = 300)
ggsave("~/Desktop/manuel/JACC_figures/figure_forest_stroke.tif",
       combined, width = 14, height = 7, dpi = 300,
       device = "tiff", compression = "lzw")

cat("\nSaved:\n  ~/Desktop/manuel/forest_stroke.pdf\n  ~/Desktop/manuel/forest_stroke.png\n")

# Save tabular results too
write.csv(plot_df, "~/Desktop/manuel/forest_stroke_results.csv", row.names = FALSE)
cat("Table:  ~/Desktop/manuel/forest_stroke_results.csv\n")
