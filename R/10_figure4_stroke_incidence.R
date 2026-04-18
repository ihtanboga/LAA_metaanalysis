# ============================================================
# Modeled stroke incidence rate over time — POOLED ONLY
# Poisson GAM: status ~ arm + s(time, by = arm) + s(trial, bs = "re")
#   trial enters as a random intercept
# Predictions exclude the random effect (population-average shape)
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

ch <- load_trial(file.path(ipd_dir, "CHAMPION_all_stroke.csv"),
                 "CHAMPION", "AC", "device")
op <- load_trial(file.path(ipd_dir, "OPTION_stroke_se.csv"),
                 "OPTION", "AC", "device")
pr <- load_trial(file.path(ipd_dir, "PRAGUE17_all_stroke_TIA.csv"),
                 "PRAGUE-17", "AC", "device")
cl <- load_trial(file.path(ipd_dir, "CLOSURE_stroke_se.csv"),
                 "CLOSURE-AF", "Medical Therapy", "LAA Closure")

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

pool_contemp <- rbind(ch_sp, op_sp, pr_sp)
pool_expand  <- rbind(ch_sp, op_sp, pr_sp, cl_sp)
pool_contemp$trial <- factor(pool_contemp$trial)
pool_expand$trial  <- factor(pool_expand$trial)

# ===== Fit Poisson GAM with trial random intercept =====
fit_pooled <- function(sp, k_time = 8) {
  gam(status ~ arm_h + s(mid_time, by = arm_h, k = k_time) +
        s(trial, bs = "re") + offset(log(pt)),
      data = sp, family = poisson, method = "REML")
}

m_c <- fit_pooled(pool_contemp)
m_e <- fit_pooled(pool_expand)

cat("=== Contemporary pooled model ===\n")
print(summary(m_c))
cat("\n=== Expanded pooled model ===\n")
print(summary(m_e))

# Variance component of the trial random effect
var_trial <- function(mod) {
  vc <- gam.vcomp(mod, rescale = FALSE)
  if (is.null(dim(vc))) vc else vc[grepl("trial", rownames(vc)), ]
}
cat("\nRandom-effect variance components:\n")
cat("Contemporary:\n"); print(var_trial(m_c))
cat("Expanded:\n");     print(var_trial(m_e))

# ===== Predict, excluding the trial random effect =====
tgrid <- seq(0.02, 4, by = 0.02)

predict_pooled <- function(mod, label) {
  nd <- expand.grid(mid_time = tgrid, arm_h = c("AC", "device"),
                    pt = 1, trial = levels(mod$model$trial)[1])
  p <- predict(mod, newdata = nd, se.fit = TRUE,
               exclude = "s(trial)")
  nd$rate  <- exp(p$fit) * 100
  nd$lo    <- exp(p$fit - 1.96 * p$se.fit) * 100
  nd$hi    <- exp(p$fit + 1.96 * p$se.fit) * 100
  nd$label <- label
  nd$trial <- NULL
  nd
}

pred_c <- predict_pooled(m_c, "Contemporary trials")
pred_e <- predict_pooled(m_e, "Expanded contemporary trials")
pred_all <- rbind(pred_c, pred_e)

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
                                  "Device (LAA closure)"),
                       name = NULL) +
    scale_fill_manual(values = c("AC" = "#1F77B4", "device" = "#D62728"),
                      guide = "none") +
    scale_x_continuous(breaks = c(0, d180, 1, 2, 3, 4),
                       labels = c("0", "0.5", "1", "2", "3", "4"),
                       limits = c(0, 4)) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = title,
         x = "Years since randomization",
         y = "Modeled stroke rate (per 100 patient-years)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "top",
          legend.key.width = unit(1.2, "cm"))
}

p1 <- make_plot(pred_c, "Contemporary trials")
p2 <- make_plot(pred_e, "Expanded contemporary trials")
combined <- p1 | p2

# JACC spec: ≥2100 px wide at 300 DPI → double-panel 12" = 3600 px
ggsave("~/Desktop/manuel/JACC_figures/figure_stroke_incidence.png",
       combined, width = 12, height = 5.5, dpi = 300)
ggsave("~/Desktop/manuel/JACC_figures/figure_stroke_incidence.tif",
       combined, width = 12, height = 5.5, dpi = 300,
       device = "tiff", compression = "lzw")

cat("\nSaved:\n")
cat("  ~/Desktop/manuel/model_incidence_stroke.pdf\n")
cat("  ~/Desktop/manuel/model_incidence_stroke.png\n")

# ===== Landmarks =====
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
  rbind(extract(pred_c, t, "Contemporary"),
        extract(pred_e, t, "Expanded"))
}))
cat("\n=== Modeled pooled incidence (per 100 PY), 95% CI ===\n")
print(smry, row.names = FALSE)
write.csv(smry, "~/Desktop/manuel/pooled_incidence_modeled_stroke.csv",
          row.names = FALSE)
