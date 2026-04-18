# ============================================================
# Modeled NON-PROCEDURAL bleeding incidence rate over time — POOLED
# Poisson GAM: status ~ arm + s(time, by=arm) + s(trial, bs="re")
#
# Trials & endpoints (non-procedural flavor):
#   CHAMPION-AF : non-procedure bleeding
#   OPTION      : non–procedure-related major or CRNM bleeding
#   PRAGUE-17   : non-procedural clinically relevant bleeding
# CLOSURE-AF:   NOT included — the source publication does not report a
#               non-procedural bleeding endpoint (only composite BARC 3–5).
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(mgcv)
  library(ggplot2)
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

ch <- load_trial(file.path(ipd_dir, "CHAMPION_nonprocedure_bleeding.csv"),
                 "CHAMPION", "AC", "device")
op <- load_trial(file.path(ipd_dir, "OPTION_nonprocedure_major_or_CRNM_bleeding.csv"),
                 "OPTION", "AC", "device")
pr <- load_trial(file.path(ipd_dir, "PRAGUE17_nonprocedural_clinically_relevant_bleeding.csv"),
                 "PRAGUE-17", "AC", "device")

cat("=== Loaded non-procedural bleeding IPD ===\n")
for (d in list(ch, op, pr)) {
  cat(sprintf("  %-10s N=%d (AC %d, dev %d), events AC=%d, dev=%d\n",
              d$trial[1], nrow(d),
              sum(d$arm_h == "AC"), sum(d$arm_h == "device"),
              sum(d$arm_h == "AC"     & d$status == 1),
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

pool_sp <- rbind(split_to_episodes(ch, cuts),
                 split_to_episodes(op, cuts),
                 split_to_episodes(pr, cuts))
pool_sp$trial <- factor(pool_sp$trial)

# ===== Fit =====
m <- gam(status ~ arm_h + s(mid_time, by = arm_h, k = 8) +
           s(trial, bs = "re") + offset(log(pt)),
         data = pool_sp, family = poisson, method = "REML")

cat("\n=== Model summary ===\n")
print(summary(m))

vc <- gam.vcomp(m, rescale = FALSE)
sd_trial <- vc[grepl("trial", rownames(vc)), "std.dev"]
cat(sprintf("\nBetween-trial SD on log-rate scale = %.3f\n", sd_trial))

# ===== Predict, excluding RE =====
tgrid <- seq(0.02, 4, by = 0.02)
nd <- expand.grid(mid_time = tgrid, arm_h = c("AC", "device"),
                  pt = 1, trial = levels(pool_sp$trial)[1])
p <- predict(m, newdata = nd, se.fit = TRUE, exclude = "s(trial)")
nd$rate <- exp(p$fit) * 100
nd$lo   <- exp(p$fit - 1.96 * p$se.fit) * 100
nd$hi   <- exp(p$fit + 1.96 * p$se.fit) * 100
nd$trial <- NULL

d180 <- 180 / 365.25

# ===== Plot (single panel) =====
pl <- ggplot(nd, aes(x = mid_time, y = rate, color = arm_h, fill = arm_h)) +
  geom_ribbon(aes(ymin = lo, ymax = pmin(hi, 50)), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
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
  labs(title = "Modeled non-procedural bleeding incidence\nCHAMPION-AF + OPTION + PRAGUE-17 (CLOSURE-AF not available)",
       x = "Years since randomization",
       y = "Modeled bleeding rate (per 100 patient-years)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        legend.position = "top",
        legend.key.width = unit(1.2, "cm"))

# JACC spec: ≥2100 px wide at 300 DPI → single panel 7.5" = 2250 px
ggsave("~/Desktop/manuel/JACC_figures/figure_nonproc_bleeding_incidence.png",
       pl, width = 7.5, height = 5.5, dpi = 300)
ggsave("~/Desktop/manuel/JACC_figures/figure_nonproc_bleeding_incidence.tif",
       pl, width = 7.5, height = 5.5, dpi = 300,
       device = "tiff", compression = "lzw")

cat("\nSaved:\n")
cat("  ~/Desktop/manuel/model_incidence_nonproc_bleeding.pdf\n")
cat("  ~/Desktop/manuel/model_incidence_nonproc_bleeding.png\n")

# ===== Landmark table =====
landmarks <- c("180d" = d180, "1y" = 1, "2y" = 2, "3y" = 3)
smry <- do.call(rbind, lapply(landmarks, function(t) {
  nd %>%
    group_by(arm_h) %>%
    slice_min(abs(mid_time - t), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(at = names(landmarks)[landmarks == t]) %>%
    select(at, arm_h, rate, lo, hi)
}))
cat("\n=== Non-procedural bleeding rate (per 100 PY), 95% CI ===\n")
print(smry, row.names = FALSE, n = 20)
write.csv(smry, "~/Desktop/manuel/pooled_incidence_modeled_nonproc_bleeding.csv",
          row.names = FALSE)
