# ============================================================
# Modeled bleeding incidence rate over time — POOLED ONLY
# Poisson GAM: status ~ arm + s(time, by=arm) + s(trial, bs="re")
#
# Endpoint mapping (because definitions are heterogeneous):
#   CHAMPION-AF : major bleeding
#   OPTION      : major bleeding
#   PRAGUE-17   : clinically relevant bleeding (no major-only curve published)
#   CLOSURE-AF  : major bleeding (BARC 3–5, reconstructed)
#
# Primary analysis : Contemporary / Expanded, including PRAGUE-17 (CRB)
# Sensitivity       : same, excluding PRAGUE-17 (major-only trials)
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(mgcv)
  library(ggplot2)
  library(patchwork)
})

ipd_dir <- "~/Desktop/manuel/json_csv/ipd_data"

load_trial <- function(path, trial, ac_lab, dev_lab) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$trial <- trial
  df$arm_h <- ifelse(df$arm == ac_lab, "AC",
              ifelse(df$arm == dev_lab, "device", NA))
  df <- df[!is.na(df$arm_h), ]
  df$arm_h <- factor(df$arm_h, levels = c("AC", "device"))
  df[, c("trial", "arm_h", "time", "status")]
}

ch <- load_trial(file.path(ipd_dir, "CHAMPION_major_bleeding.csv"),
                 "CHAMPION", "AC", "device")
op <- load_trial(file.path(ipd_dir, "OPTION_major_bleeding.csv"),
                 "OPTION", "AC", "device")
pr <- load_trial(file.path(ipd_dir, "PRAGUE17_clinically_relevant_bleeding.csv"),
                 "PRAGUE-17", "AC", "device")
cl <- load_trial(file.path(ipd_dir, "CLOSURE_majorbleeding.csv"),
                 "CLOSURE-AF", "Medical Therapy", "LAA Closure")

cat("=== Loaded bleeding IPD ===\n")
for (d in list(ch, op, pr, cl)) {
  cat(sprintf("  %-12s N=%d (AC %d, dev %d), events AC=%d, dev=%d\n",
              d$trial[1], nrow(d),
              sum(d$arm_h == "AC"), sum(d$arm_h == "device"),
              sum(d$arm_h == "AC" & d$status == 1),
              sum(d$arm_h == "device" & d$status == 1)))
}

split_to_episodes <- function(df, cuts) {
  df$id <- seq_len(nrow(df))
  sp <- survSplit(Surv(time, status) ~ ., data = df, cut = cuts,
                  episode = "interval", start = "tstart",
                  end = "time", event = "status")
  sp$pt       <- sp$time - sp$tstart
  sp$mid_time <- (sp$tstart + sp$time) / 2
  sp[sp$pt > 0, ]
}

cuts <- c(seq(0.02, 0.5, by = 0.02),
          seq(0.55, 2,   by = 0.08),
          seq(2.15, 6,   by = 0.15))

ch_sp <- split_to_episodes(ch, cuts)
op_sp <- split_to_episodes(op, cuts)
pr_sp <- split_to_episodes(pr, cuts)
cl_sp <- split_to_episodes(cl, cuts)

# ===== Fit =====
fit_pooled <- function(sp, k_time = 8) {
  sp$trial <- factor(sp$trial)
  # If only 2 trials, s(trial, bs="re") is effectively a fixed 2-level
  # intercept; use plain trial fixed effect instead to avoid singularity.
  if (length(unique(sp$trial)) >= 3) {
    gam(status ~ arm_h + s(mid_time, by = arm_h, k = k_time) +
          s(trial, bs = "re") + offset(log(pt)),
        data = sp, family = poisson, method = "REML")
  } else {
    gam(status ~ arm_h + s(mid_time, by = arm_h, k = k_time) +
          trial + offset(log(pt)),
        data = sp, family = poisson, method = "REML")
  }
}

m_primary_contemp  <- fit_pooled(rbind(ch_sp, op_sp, pr_sp))
m_primary_expand   <- fit_pooled(rbind(ch_sp, op_sp, pr_sp, cl_sp))
m_sens_contemp     <- fit_pooled(rbind(ch_sp, op_sp))          # no PRAGUE
m_sens_expand      <- fit_pooled(rbind(ch_sp, op_sp, cl_sp))   # no PRAGUE

# ===== Predict, excluding trial RE =====
tgrid <- seq(0.02, 4, by = 0.02)
predict_pooled <- function(mod) {
  first_trial <- levels(mod$model$trial)[1]
  nd <- expand.grid(mid_time = tgrid, arm_h = c("AC", "device"),
                    pt = 1, trial = first_trial)
  p <- predict(mod, newdata = nd, se.fit = TRUE,
               exclude = "s(trial)")
  nd$rate <- exp(p$fit) * 100
  nd$lo   <- exp(p$fit - 1.96 * p$se.fit) * 100
  nd$hi   <- exp(p$fit + 1.96 * p$se.fit) * 100
  nd$trial <- NULL
  nd
}

pred_primary_c <- predict_pooled(m_primary_contemp) %>% mutate(panel = "Contemporary trials")
pred_primary_e <- predict_pooled(m_primary_expand)  %>% mutate(panel = "Expanded contemporary trials")
pred_sens_c    <- predict_pooled(m_sens_contemp)    %>% mutate(panel = "Contemporary trials\n(sensitivity: without PRAGUE-17)")
pred_sens_e    <- predict_pooled(m_sens_expand)     %>% mutate(panel = "Expanded contemporary trials\n(sensitivity: without PRAGUE-17)")

d180 <- 180 / 365.25

# ===== Plot =====
make_plot <- function(df, title) {
  df$hi <- pmin(df$hi, 50)
  ggplot(df, aes(x = mid_time, y = rate, color = arm_h, fill = arm_h)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.4) +
    geom_vline(xintercept = d180, linetype = "dotted",
               color = "grey40", linewidth = 0.5) +
    annotate("text", x = d180 + 0.05, y = Inf, label = "180 days",
             hjust = 0, vjust = 1.5, size = 3, color = "grey30") +
    scale_color_manual(values = c("AC" = "#1F77B4", "device" = "#D62728"),
                       labels = c("AC (anticoagulation)",
                                  "Device (LAA closure)"), name = NULL) +
    scale_fill_manual(values = c("AC" = "#1F77B4", "device" = "#D62728"),
                      guide = "none") +
    scale_x_continuous(breaks = c(0, d180, 1, 2, 3, 4),
                       labels = c("0", "0.5", "1", "2", "3", "4"),
                       limits = c(0, 4)) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = title,
         x = "Years since randomization",
         y = "Modeled bleeding rate (per 100 patient-years)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "top",
          legend.key.width = unit(1.2, "cm"))
}

# ===== Primary figure =====
p1 <- make_plot(pred_primary_c, "Contemporary trials")
p2 <- make_plot(pred_primary_e, "Expanded contemporary trials")
primary_fig <- p1 | p2

ggsave("~/Desktop/manuel/JACC_figures/figure_bleeding_incidence.png",
       primary_fig, width = 12, height = 5.5, dpi = 300)
ggsave("~/Desktop/manuel/JACC_figures/figure_bleeding_incidence.tif",
       primary_fig, width = 12, height = 5.5, dpi = 300,
       device = "tiff", compression = "lzw")

# ===== Sensitivity figure =====
p3 <- make_plot(pred_sens_c, "Contemporary trials\n(sensitivity: without PRAGUE-17)")
p4 <- make_plot(pred_sens_e, "Expanded contemporary trials\n(sensitivity: without PRAGUE-17)")
sens_fig <- p3 | p4

ggsave("~/Desktop/manuel/JACC_figures/figure_bleeding_incidence_sensitivity.png",
       sens_fig, width = 12, height = 5.5, dpi = 300)
ggsave("~/Desktop/manuel/JACC_figures/figure_bleeding_incidence_sensitivity.tif",
       sens_fig, width = 12, height = 5.5, dpi = 300,
       device = "tiff", compression = "lzw")

cat("\nSaved:\n")
cat("  ~/Desktop/manuel/model_incidence_bleeding.pdf\n")
cat("  ~/Desktop/manuel/model_incidence_bleeding_sensitivity.pdf\n")

# ===== Random-effect variance =====
re_sd <- function(mod, name) {
  if ("s(trial)" %in% sapply(mod$smooth, function(s) s$label)) {
    vc <- gam.vcomp(mod, rescale = FALSE)
    sd_trial <- vc[grepl("trial", rownames(vc)), "std.dev"]
    cat(sprintf("  %-40s SD(trial) = %.3f\n", name, sd_trial))
  } else {
    cat(sprintf("  %-40s trial as fixed effect (n_trials=2)\n", name))
  }
}
cat("\n=== Between-trial random effect (log-rate SD) ===\n")
re_sd(m_primary_contemp, "Primary – Contemporary")
re_sd(m_primary_expand,  "Primary – Expanded")
re_sd(m_sens_contemp,    "Sensitivity – Contemporary (no PRAGUE)")
re_sd(m_sens_expand,     "Sensitivity – Expanded (no PRAGUE)")

# ===== Landmark table =====
landmarks <- c("180d" = d180, "1y" = 1, "2y" = 2, "3y" = 3)
extract <- function(pred_df, at_t, panel) {
  pred_df %>%
    group_by(arm_h) %>%
    slice_min(abs(mid_time - at_t), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(at = names(landmarks)[landmarks == at_t],
           panel = panel) %>%
    select(panel, at, arm_h, rate, lo, hi)
}

smry <- do.call(rbind, lapply(landmarks, function(t) {
  rbind(
    extract(pred_primary_c, t, "Primary — Contemporary"),
    extract(pred_primary_e, t, "Primary — Expanded"),
    extract(pred_sens_c,    t, "Sensitivity (no PRAGUE) — Contemporary"),
    extract(pred_sens_e,    t, "Sensitivity (no PRAGUE) — Expanded")
  )
}))
cat("\n=== Modeled pooled bleeding incidence (per 100 PY), 95% CI ===\n")
print(smry, row.names = FALSE, n = 60)
write.csv(smry, "~/Desktop/manuel/pooled_incidence_modeled_bleeding.csv",
          row.names = FALSE)
