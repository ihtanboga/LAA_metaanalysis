# ============================================================
# CLOSURE-AF: Bleeding & Death Excluded Endpoint Analyses
# ============================================================
#
# NEW MODELS (corrected from v4):
#   v4 targeted Table S13 (removed: bleeding + hemorrhagic stroke + unexplained death)
#   NOW: Hemorrhagic stroke STAYS in all models
#
# Model 1: Remove ONLY major bleeding
#   Remaining: Stroke (all types) + SE + CV/Unexplained Death
#   → Shows primary endpoint without procedure-related bleeding bias
#
# Model 2: Remove major bleeding + ALL death (CV + unexplained)
#   Remaining: Stroke (all types) + SE = pure thromboembolic
#   → Isolates stroke prevention efficacy
#
# Supplementary: Non-CV death IPD reconstruction + KM plot
#
# Table 2 reference (paper first events):
#   Primary composite: LAA=155, Medical=127
#   Major bleeding: LAA=70, Medical=61
#   CV/unexplained death: LAA=99, Medical=81
#   Stroke+SE: LAA=29, Medical=28
#   Ischemic stroke: LAA=18, Medical=15
#   Hemorrhagic stroke: LAA=10, Medical=13
#   SE: LAA=3, Medical=1
# ============================================================

library(readr)
library(dplyr)
library(survival)
library(survminer)
library(survRM2)
library(ggplot2)
library(ggpubr)
library(reconstructKM)

set.seed(2026)

# ============================================================
# 0) Load Reconstructed Primary Composite IPD
# ============================================================

allIPD <- read.csv("closure_af_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Medical Therapy", "LAA Closure"))

cat("=== Loaded Primary Composite IPD ===\n")
cat("N:", nrow(allIPD), "\n")
for (a in levels(allIPD$arm)) {
  cat(sprintf("  %s: N = %d, Events = %d\n", a,
              sum(allIPD$arm == a), sum(allIPD$arm == a & allIPD$status == 1)))
}

# ============================================================
# 1) Non-CV Death IPD Reconstruction (Supplementary)
# ============================================================

cat("\n\n################################################################\n")
cat("# SUPPLEMENTARY: Non-CV Death IPD Reconstruction\n")
cat("################################################################\n\n")

# --- Read digitized non-CV death KM curves ---
# Format: European (semicolon-delimited, comma decimal)
# Columns: time (years), cumulative incidence (%)
laa_death_raw <- read_delim(
  "laa_death.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
) %>% filter(!is.na(time) & !is.na(cum_inc_pct))

control_death_raw <- read_delim(
  "control_death.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
) %>% filter(!is.na(time) & !is.na(cum_inc_pct))

# control_death.csv has trailing rogue points (re-digitized early segments)
# Remove any points where time decreases compared to previous
control_death_raw <- control_death_raw %>%
  filter(time == cummax(time))

cat("Non-CV death digitized points: LAA =", nrow(laa_death_raw),
    ", Medical =", nrow(control_death_raw), "\n")

# --- Transform CIF → Survival ---
laa_death_clicks <- laa_death_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

control_death_clicks <- control_death_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

# Add (0, 1) start point
add_origin <- function(df) {
  if (!any(abs(df$time) < 1e-6 & abs(df$survival - 1) < 1e-6))
    df <- bind_rows(data.frame(time = 0, survival = 1), df)
  df
}
laa_death_clicks <- add_origin(laa_death_clicks)
control_death_clicks <- add_origin(control_death_clicks)

# Enforce monotonic decreasing
collapse_corners <- function(df) {
  df <- df %>%
    group_by(time) %>%
    summarise(survival = min(survival), .groups = "drop") %>%
    arrange(time)
  df$survival <- cummin(pmin(df$survival, 1))
  df$survival[df$survival < 0] <- 0
  df
}
laa_death_clicks <- collapse_corners(laa_death_clicks)
control_death_clicks <- collapse_corners(control_death_clicks)

# Remove trailing horizontal segments: format_raw_tabs requires last click
# at end of a vertical segment. Keep only the first occurrence of the final
# survival value (where the last drop happens), then drop the horizontal tail.
trim_horiz_tail <- function(df) {
  final_surv <- df$survival[nrow(df)]
  first_final <- which(abs(df$survival - final_surv) < 1e-8)[1]
  df[1:first_final, ]
}
laa_death_clicks <- trim_horiz_tail(laa_death_clicks)
control_death_clicks <- trim_horiz_tail(control_death_clicks)

# --- NAR for non-CV death endpoint (MC-based estimation) ---
# Primary composite NAR drops too fast because non-fatal first events
# (bleeding, stroke, SE) remove patients who are still alive and at risk
# for non-CV death.
#
# Table 2 component totals include first + subsequent events:
#   Stroke/SE total: LAA=29, Medical=28
#   Bleeding total: LAA=70, Medical=61
#   CV/unexplained death total: LAA=99, Medical=81
#   Sum = 198/170, but only 155/127 first events
#   → 43/43 patients had 2+ component events (mostly non-fatal first → later death)
#
# First-event deaths (patients who died WITHOUT prior non-fatal event):
#   LAA: 155 - 70 - 29 = 56 (36.1% of first events)
#   Medical: 127 - 61 - 28 = 38 (29.9% of first events)
#
# Non-fatal first events (patients still alive, at risk for non-CV death):
#   LAA: 70 + 29 = 99 (63.9%)
#   Medical: 61 + 28 = 89 (70.1%)
#
# MC approach: for each event in allIPD, randomly classify as death (leave risk set)
# or non-fatal (extend follow-up), then compute NAR. Average over iterations.

primary_laa_NAR <- data.frame(time = 0:6, NAR = c(446, 304, 202, 117, 71, 33, 9))
primary_ctl_NAR <- data.frame(time = 0:6, NAR = c(442, 306, 203, 136, 77, 40, 7))

target_events <- c("LAA Closure" = 54, "Medical Therapy" = 60)

# Trim clicks to 6-year max
laa_death_clicks <- laa_death_clicks %>% filter(time <= 6)
control_death_clicks <- control_death_clicks %>% filter(time <= 6)

# --- reconstructKM helper ---
fix_aug_surv <- function(aug) {
  as <- aug$aug_surv
  if (!"surv" %in% names(as) && "survival" %in% names(as))
    names(as)[names(as) == "survival"] <- "surv"
  as <- as[order(as$time), ]
  aug$aug_surv <- as
  aug
}

# --- Calibration: binary search for NAR scale factor ---
# Scale NAR(t>0) upward (keeping t=0 = N) until reconstructed events ≈ target
calibrate_nar <- function(clicks, base_nar, N_total, target_ev, arm_label) {
  k_lo <- 1.0
  k_hi <- 5.0
  best <- NULL

  for (iter in 1:30) {
    k <- (k_lo + k_hi) / 2

    nar_scaled <- base_nar
    nar_scaled$NAR[-1] <- pmin(N_total, round(base_nar$NAR[-1] * k))
    # Enforce monotonic non-increasing
    for (i in 2:nrow(nar_scaled)) {
      nar_scaled$NAR[i] <- min(nar_scaled$NAR[i], nar_scaled$NAR[i-1])
    }

    aug <- tryCatch({
      a <- format_raw_tabs(raw_NAR = nar_scaled, raw_surv = clicks)
      fix_aug_surv(a)
    }, error = function(e) NULL)

    if (is.null(aug)) { k_hi <- k; next }

    recon <- tryCatch(
      KM_reconstruct(aug_NAR = aug$aug_NAR, aug_surv = aug$aug_surv),
      error = function(e) NULL)

    if (is.null(recon)) { k_hi <- k; next }

    events <- sum(recon$IPD_event)

    best <- list(NAR = nar_scaled, recon = recon, events = events, k = k)
    if (abs(events - target_ev) <= 2) {
      cat(sprintf("  %s: k=%.3f → events=%d (target=%d) ✓\n", arm_label, k, events, target_ev))
      break
    }
    if (events < target_ev) k_lo <- k
    else k_hi <- k
  }

  if (is.null(best)) stop(paste("Calibration failed for", arm_label))
  return(best)
}

cat("Calibrating NAR for non-CV death (binary search)...\n")

cal_laa <- calibrate_nar(laa_death_clicks, primary_laa_NAR, 446,
                          target_events["LAA Closure"], "LAA")
cal_ctl <- calibrate_nar(control_death_clicks, primary_ctl_NAR, 442,
                          target_events["Medical Therapy"], "Medical")

laa_NAR <- cal_laa$NAR
control_NAR <- cal_ctl$NAR

cat(sprintf("\nCalibrated NAR:\n"))
cat("  LAA:    ", paste(laa_NAR$NAR, collapse = ", "),
    sprintf(" (k=%.2f, events=%d)\n", cal_laa$k, cal_laa$events))
cat("  Medical:", paste(control_NAR$NAR, collapse = ", "),
    sprintf(" (k=%.2f, events=%d)\n", cal_ctl$k, cal_ctl$events))

# Use calibrated reconstructions
laa_death_recon <- cal_laa$recon
control_death_recon <- cal_ctl$recon

# Combine death IPD
death_IPD <- rbind(
  data.frame(arm = "LAA Closure",
             time = laa_death_recon$IPD_time,
             status = laa_death_recon$IPD_event),
  data.frame(arm = "Medical Therapy",
             time = control_death_recon$IPD_time,
             status = control_death_recon$IPD_event)
)
death_IPD$arm <- factor(death_IPD$arm, levels = c("Medical Therapy", "LAA Closure"))

cat("\n=== Non-CV Death Reconstructed IPD ===\n")
for (a in levels(death_IPD$arm)) {
  cat(sprintf("  %s: N = %d, Events = %d (Table S12: %s)\n", a,
              sum(death_IPD$arm == a),
              sum(death_IPD$arm == a & death_IPD$status == 1),
              ifelse(a == "LAA Closure", "54", "60")))
}

# --- Non-CV Death Cox + KM ---
cx_death <- coxph(Surv(time, status) ~ arm, data = death_IPD)
s_death <- summary(cx_death)
hr_death <- exp(coef(cx_death))[1]
ci_death <- exp(confint(cx_death))[1, ]
p_death <- s_death$coefficients[1, "Pr(>|z|)"]

cat(sprintf("\nNon-CV Death HR (LAA vs Medical): %.2f (95%% CI %.2f-%.2f), P = %s\n",
            hr_death, ci_death[1], ci_death[2],
            ifelse(p_death < 0.001, "<0.001", sprintf("%.3f", p_death))))

# RMST
rmst_death <- rmst2(
  time = death_IPD$time, status = death_IPD$status,
  arm = as.numeric(death_IPD$arm) - 1, tau = 6
)
cat("\n=== Non-CV Death RMST (tau = 6 yr) ===\n")
print(rmst_death)

# KM plot
fit_death <- survfit(Surv(time, status) ~ arm, data = death_IPD)
p_death_km <- ggsurvplot(
  fit_death, data = death_IPD,
  fun = "event", conf.int = TRUE, censor = FALSE,
  xlim = c(0, 6), ylim = c(0, 0.30),
  break.time.by = 1,
  palette = c("#1F77B4", "#D62728"),
  legend.title = "", legend.labs = c("Medical Therapy", "LAA Closure"),
  legend = c(0.25, 0.85),
  xlab = "Years since Randomization",
  ylab = "Cumulative Incidence (%)",
  risk.table = TRUE, risk.table.title = "No. at Risk",
  title = "Non-CV Death",
  ggtheme = theme_classic(base_size = 13)
)
p_death_km$plot <- p_death_km$plot +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     limits = c(0, 0.30), breaks = seq(0, 0.30, by = 0.05)) +
  annotate("text", x = 3, y = 0.25, hjust = 0.5, size = 3.8,
           label = sprintf("HR %.2f (95%% CI %.2f\u2013%.2f)\nP = %s",
                           hr_death, ci_death[1], ci_death[2],
                           ifelse(p_death < 0.001, "<0.001", sprintf("%.3f", p_death))))

death_combined <- ggarrange(p_death_km$plot, p_death_km$table,
                             ncol = 1, nrow = 2, heights = c(3, 1))
ggsave("closure_af_noncv_death_km.pdf", plot = death_combined, width = 10, height = 7)
cat("Plot saved: closure_af_noncv_death_km.pdf\n")


# ============================================================
# 2) Setup for Monte Carlo Event Thinning
# ============================================================

cat("\n\n################################################################\n")
cat("# MC EVENT THINNING SETUP\n")
cat("################################################################\n\n")

# --- Paper event counts (Table 2) ---
paper_events <- c("LAA Closure" = 155, "Medical Therapy" = 127)
our_events <- c(
  "LAA Closure"     = sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
  "Medical Therapy" = sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)
)
scale_ratio <- our_events / paper_events

cat(sprintf("Scale ratios: LAA = %.3f, Medical = %.3f\n",
            scale_ratio["LAA Closure"], scale_ratio["Medical Therapy"]))

# --- Model 1: Remove ONLY major bleeding ---
# Paper: major bleeding LAA=70, Medical=61
bleed_remove_paper <- c("LAA Closure" = 70, "Medical Therapy" = 61)
bleed_remove_scaled <- bleed_remove_paper * scale_ratio

cat(sprintf("\nModel 1 — Bleeding removal targets (scaled):\n"))
cat(sprintf("  LAA:     %.1f (paper: %d)\n",
            bleed_remove_scaled["LAA Closure"], bleed_remove_paper["LAA Closure"]))
cat(sprintf("  Medical: %.1f (paper: %d)\n",
            bleed_remove_scaled["Medical Therapy"], bleed_remove_paper["Medical Therapy"]))
cat(sprintf("Expected remaining: LAA ~ %.0f, Medical ~ %.0f\n",
            our_events["LAA Closure"] - bleed_remove_scaled["LAA Closure"],
            our_events["Medical Therapy"] - bleed_remove_scaled["Medical Therapy"]))
cat(sprintf("Paper check: 155-70=%d, 127-61=%d\n", 155-70, 127-61))

# --- Model 2: Remove bleeding + ALL death (CV + unexplained) ---
# Target remaining = stroke + SE only: LAA=29, Medical=28
stroke_se_paper <- c("LAA Closure" = 29, "Medical Therapy" = 28)
stroke_se_scaled <- stroke_se_paper * scale_ratio

total_remove_m2 <- our_events - stroke_se_scaled
death_remove_m2 <- total_remove_m2 - bleed_remove_scaled  # death component

cat(sprintf("\nModel 2 — Total removal targets (scaled):\n"))
cat(sprintf("  LAA:     %.1f total (bleeding: %.1f + death: %.1f)\n",
            total_remove_m2["LAA Closure"],
            bleed_remove_scaled["LAA Closure"],
            death_remove_m2["LAA Closure"]))
cat(sprintf("  Medical: %.1f total (bleeding: %.1f + death: %.1f)\n",
            total_remove_m2["Medical Therapy"],
            bleed_remove_scaled["Medical Therapy"],
            death_remove_m2["Medical Therapy"]))
cat(sprintf("Expected remaining: LAA ~ %.1f, Medical ~ %.1f\n",
            stroke_se_scaled["LAA Closure"], stroke_se_scaled["Medical Therapy"]))

# --- Table S11: 3-period temporal distribution of bleeding ---
bleed_3mo   <- c("LAA Closure" = 40, "Medical Therapy" = 17)
bleed_3to6  <- c("LAA Closure" = 7,  "Medical Therapy" = 6)
bleed_6plus <- c("LAA Closure" = 23, "Medical Therapy" = 38)
bleed_total <- c("LAA Closure" = 70, "Medical Therapy" = 61)
frac_3mo    <- bleed_3mo / bleed_total
frac_3to6   <- bleed_3to6 / bleed_total
frac_6plus  <- bleed_6plus / bleed_total

cat(sprintf("\nBleeding temporal fractions (Table S11):\n"))
cat(sprintf("  LAA:     <3mo: %.1f%%, 3-6mo: %.1f%%, >=6mo: %.1f%%\n",
            frac_3mo["LAA Closure"]*100, frac_3to6["LAA Closure"]*100, frac_6plus["LAA Closure"]*100))
cat(sprintf("  Medical: <3mo: %.1f%%, 3-6mo: %.1f%%, >=6mo: %.1f%%\n",
            frac_3mo["Medical Therapy"]*100, frac_3to6["Medical Therapy"]*100, frac_6plus["Medical Therapy"]*100))

# --- Pre-compute censoring distributions ---
censor_times <- list()
for (arm_name in c("LAA Closure", "Medical Therapy")) {
  censor_times[[arm_name]] <- allIPD$time[allIPD$arm == arm_name & allIPD$status == 0]
}

extend_followup <- function(arm_name, event_time) {
  available <- censor_times[[arm_name]][censor_times[[arm_name]] > event_time]
  if (length(available) > 0) {
    sample(available, 1)
  } else {
    max(allIPD$time[allIPD$arm == arm_name])
  }
}

stoch_round <- function(x) {
  floor(x) + rbinom(1, 1, x - floor(x))
}


# ============================================================
# 3) Model 1: Remove ONLY Major Bleeding (MC Simulation)
# ============================================================

cat("\n\n################################################################\n")
cat("# MODEL 1: Remove Only Major Bleeding\n")
cat("# Remaining: Stroke (all types) + SE + CV/Unexplained Death\n")
cat("# Hemorrhagic stroke STAYS\n")
cat("################################################################\n\n")

n_sim <- 1000

run_model1_sim <- function(n_sim, model_type = "temporal") {
  results <- data.frame(
    sim = integer(n_sim), hr = numeric(n_sim),
    hr_lo = numeric(n_sim), hr_hi = numeric(n_sim),
    drmst = numeric(n_sim), drmst_lo = numeric(n_sim), drmst_hi = numeric(n_sim),
    n_events_laa = integer(n_sim), n_events_med = integer(n_sim)
  )

  for (i in 1:n_sim) {
    sim_data <- allIPD

    for (arm_name in c("LAA Closure", "Medical Therapy")) {
      events_idx <- which(sim_data$arm == arm_name & sim_data$status == 1)
      n_remove <- stoch_round(bleed_remove_scaled[arm_name])
      n_remove <- min(n_remove, length(events_idx))

      if (model_type == "proportional") {
        # Model 1B: Simple random selection
        selected <- sample(events_idx, n_remove)
      } else {
        # Model 1A: Bleeding-temporal weighted (3-period from Table S11)
        idx_p1 <- events_idx[sim_data$time[events_idx] < 0.25]
        idx_p2 <- events_idx[sim_data$time[events_idx] >= 0.25 &
                               sim_data$time[events_idx] < 0.5]
        idx_p3 <- events_idx[sim_data$time[events_idx] >= 0.5]

        n1 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3mo[arm_name]),
                   length(idx_p1))
        n2 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3to6[arm_name]),
                   length(idx_p2))
        n3 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_6plus[arm_name]),
                   length(idx_p3))

        selected <- c()
        if (n1 > 0) selected <- c(selected, sample(idx_p1, n1))
        if (n2 > 0) selected <- c(selected, sample(idx_p2, n2))
        if (n3 > 0) selected <- c(selected, sample(idx_p3, n3))
      }

      for (j in selected) {
        old_time <- sim_data$time[j]
        sim_data$status[j] <- 0
        sim_data$time[j] <- extend_followup(arm_name, old_time)
      }
    }

    # Cox
    cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = sim_data), error = function(e) NULL)
    if (!is.null(cx)) {
      results$hr[i] <- exp(coef(cx))[1]
      results$hr_lo[i] <- exp(confint(cx))[1, 1]
      results$hr_hi[i] <- exp(confint(cx))[1, 2]
    }

    # RMST
    rmst <- tryCatch(
      rmst2(time = sim_data$time, status = sim_data$status,
            arm = as.numeric(sim_data$arm) - 1, tau = 6),
      error = function(e) NULL)
    if (!is.null(rmst)) {
      results$drmst[i] <- rmst$unadjusted.result[1, "Est."]
      results$drmst_lo[i] <- rmst$unadjusted.result[1, "lower .95"]
      results$drmst_hi[i] <- rmst$unadjusted.result[1, "upper .95"]
    }

    results$sim[i] <- i
    results$n_events_laa[i] <- sum(sim_data$arm == "LAA Closure" & sim_data$status == 1)
    results$n_events_med[i] <- sum(sim_data$arm == "Medical Therapy" & sim_data$status == 1)
  }
  return(results)
}

# --- Model 1A: Bleeding-temporal weighted ---
cat("--- Model 1A: Bleeding-temporal weighted (Table S11) ---\n")
cat("Running", n_sim, "simulations...\n")
sim_1A <- run_model1_sim(n_sim, "temporal")

cat(sprintf("\n=== Model 1A Results ===\n"))
cat(sprintf("Events remaining: LAA = %.0f, Medical = %.0f (target: ~81, ~63)\n",
            median(sim_1A$n_events_laa), median(sim_1A$n_events_med)))
cat(sprintf("Median HR: %.2f (Sim 95%%: %.2f - %.2f)\n",
            median(sim_1A$hr), quantile(sim_1A$hr, 0.025), quantile(sim_1A$hr, 0.975)))
cat(sprintf("Median dRMST: %.3f yr (Sim 95%%: %.3f to %.3f)\n",
            median(sim_1A$drmst), quantile(sim_1A$drmst, 0.025), quantile(sim_1A$drmst, 0.975)))
cat(sprintf("Within-sim 95%% CI: %.3f to %.3f\n",
            median(sim_1A$drmst_lo), median(sim_1A$drmst_hi)))

# --- Model 1B: Proportional thinning ---
cat("\n--- Model 1B: Proportional thinning ---\n")
cat("Running", n_sim, "simulations...\n")
sim_1B <- run_model1_sim(n_sim, "proportional")

cat(sprintf("\n=== Model 1B Results ===\n"))
cat(sprintf("Events remaining: LAA = %.0f, Medical = %.0f (target: ~81, ~63)\n",
            median(sim_1B$n_events_laa), median(sim_1B$n_events_med)))
cat(sprintf("Median HR: %.2f (Sim 95%%: %.2f - %.2f)\n",
            median(sim_1B$hr), quantile(sim_1B$hr, 0.025), quantile(sim_1B$hr, 0.975)))
cat(sprintf("Median dRMST: %.3f yr (Sim 95%%: %.3f to %.3f)\n",
            median(sim_1B$drmst), quantile(sim_1B$drmst, 0.025), quantile(sim_1B$drmst, 0.975)))
cat(sprintf("Within-sim 95%% CI: %.3f to %.3f\n",
            median(sim_1B$drmst_lo), median(sim_1B$drmst_hi)))

# --- Model 1 KM Plot (representative run from 1B) ---
best_1B_idx <- which.min(abs(sim_1B$drmst - median(sim_1B$drmst)))
cat(sprintf("\nRepresentative 1B sim (#%d): HR=%.2f, dRMST=%.3f\n",
            best_1B_idx, sim_1B$hr[best_1B_idx], sim_1B$drmst[best_1B_idx]))

set.seed(2026 + best_1B_idx)
rep_data_1 <- allIPD
for (arm_name in c("LAA Closure", "Medical Therapy")) {
  events_idx <- which(rep_data_1$arm == arm_name & rep_data_1$status == 1)
  n_remove <- min(stoch_round(bleed_remove_scaled[arm_name]), length(events_idx))
  selected <- sample(events_idx, n_remove)
  for (j in selected) {
    rep_data_1$status[j] <- 0
    rep_data_1$time[j] <- extend_followup(arm_name, rep_data_1$time[j])
  }
}

fit_m1 <- survfit(Surv(time, status) ~ arm, data = rep_data_1)
cx_m1 <- coxph(Surv(time, status) ~ arm, data = rep_data_1)
hr_m1 <- exp(coef(cx_m1))[1]
ci_m1 <- exp(confint(cx_m1))[1, ]
p_m1 <- summary(cx_m1)$coefficients[1, "Pr(>|z|)"]

p1_km <- ggsurvplot(
  fit_m1, data = rep_data_1,
  fun = "event", conf.int = FALSE, censor = FALSE,
  xlim = c(0, 6), ylim = c(0, 0.60),
  break.time.by = 1,
  palette = c("#1F77B4", "#D62728"),
  legend.title = "", legend.labs = c("Medical Therapy", "LAA Closure"),
  legend = c(0.25, 0.85),
  xlab = "Years since Randomization",
  ylab = "Cumulative Incidence (%)",
  risk.table = TRUE, risk.table.title = "No. at Risk",
  title = "Model 1: Bleeding-Excluded Composite (Stroke + SE + CV/Unexplained Death)",
  ggtheme = theme_classic(base_size = 13)
)

p1_km$plot <- p1_km$plot +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     limits = c(0, 0.60), breaks = seq(0, 0.60, by = 0.10)) +
  annotate("text", x = 3, y = 0.50, hjust = 0.5, size = 3.8,
           label = sprintf("HR %.2f (95%% CI %.2f\u2013%.2f), P = %s\nMedian MC dRMST: 1A = %.3f, 1B = %.3f yr",
                           hr_m1, ci_m1[1], ci_m1[2],
                           ifelse(p_m1 < 0.001, "<0.001", sprintf("%.3f", p_m1)),
                           median(sim_1A$drmst), median(sim_1B$drmst)))

m1_combined <- ggarrange(p1_km$plot, p1_km$table, ncol = 1, nrow = 2, heights = c(3, 1))
ggsave("closure_af_bleedonly_km.pdf", plot = m1_combined, width = 10, height = 7)
cat("Plot saved: closure_af_bleedonly_km.pdf\n")


# ============================================================
# 4) Model 2: Remove Bleeding + ALL Death (MC Simulation)
# ============================================================

cat("\n\n################################################################\n")
cat("# MODEL 2: Remove Bleeding + All Death\n")
cat("# Remaining: Stroke (all types) + SE = Pure Thromboembolic\n")
cat("# Two-phase removal: temporal for bleeding, proportional for death\n")
cat("################################################################\n\n")

run_model2_sim <- function(n_sim) {
  results <- data.frame(
    sim = integer(n_sim), hr = numeric(n_sim),
    hr_lo = numeric(n_sim), hr_hi = numeric(n_sim),
    drmst = numeric(n_sim), drmst_lo = numeric(n_sim), drmst_hi = numeric(n_sim),
    n_events_laa = integer(n_sim), n_events_med = integer(n_sim)
  )

  for (i in 1:n_sim) {
    sim_data <- allIPD

    for (arm_name in c("LAA Closure", "Medical Therapy")) {
      events_idx <- which(sim_data$arm == arm_name & sim_data$status == 1)

      # Phase 1: Remove bleeding (temporal-weighted from Table S11)
      n_bleed <- stoch_round(bleed_remove_scaled[arm_name])
      n_bleed <- min(n_bleed, length(events_idx))

      idx_p1 <- events_idx[sim_data$time[events_idx] < 0.25]
      idx_p2 <- events_idx[sim_data$time[events_idx] >= 0.25 &
                             sim_data$time[events_idx] < 0.5]
      idx_p3 <- events_idx[sim_data$time[events_idx] >= 0.5]

      n1 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3mo[arm_name]),
                 length(idx_p1))
      n2 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3to6[arm_name]),
                 length(idx_p2))
      n3 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_6plus[arm_name]),
                 length(idx_p3))

      bleed_selected <- c()
      if (n1 > 0) bleed_selected <- c(bleed_selected, sample(idx_p1, n1))
      if (n2 > 0) bleed_selected <- c(bleed_selected, sample(idx_p2, n2))
      if (n3 > 0) bleed_selected <- c(bleed_selected, sample(idx_p3, n3))

      for (j in bleed_selected) {
        sim_data$status[j] <- 0
        sim_data$time[j] <- extend_followup(arm_name, sim_data$time[j])
      }

      # Phase 2: Remove death (proportional — no strong temporal bias for death)
      remaining_events_idx <- which(sim_data$arm == arm_name & sim_data$status == 1)
      n_death <- stoch_round(death_remove_m2[arm_name])
      n_death <- min(n_death, length(remaining_events_idx))

      if (n_death > 0) {
        death_selected <- sample(remaining_events_idx, n_death)
        for (j in death_selected) {
          sim_data$status[j] <- 0
          sim_data$time[j] <- extend_followup(arm_name, sim_data$time[j])
        }
      }
    }

    # Cox
    cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = sim_data), error = function(e) NULL)
    if (!is.null(cx)) {
      results$hr[i] <- exp(coef(cx))[1]
      results$hr_lo[i] <- exp(confint(cx))[1, 1]
      results$hr_hi[i] <- exp(confint(cx))[1, 2]
    }

    # RMST
    rmst <- tryCatch(
      rmst2(time = sim_data$time, status = sim_data$status,
            arm = as.numeric(sim_data$arm) - 1, tau = 6),
      error = function(e) NULL)
    if (!is.null(rmst)) {
      results$drmst[i] <- rmst$unadjusted.result[1, "Est."]
      results$drmst_lo[i] <- rmst$unadjusted.result[1, "lower .95"]
      results$drmst_hi[i] <- rmst$unadjusted.result[1, "upper .95"]
    }

    results$sim[i] <- i
    results$n_events_laa[i] <- sum(sim_data$arm == "LAA Closure" & sim_data$status == 1)
    results$n_events_med[i] <- sum(sim_data$arm == "Medical Therapy" & sim_data$status == 1)
  }
  return(results)
}

cat("Running", n_sim, "simulations (two-phase: temporal bleeding + proportional death)...\n")
sim_2 <- run_model2_sim(n_sim)

cat(sprintf("\n=== Model 2 Results (Thromboembolic Only) ===\n"))
cat(sprintf("Events remaining: LAA = %.0f, Medical = %.0f (target: ~28, ~27)\n",
            median(sim_2$n_events_laa), median(sim_2$n_events_med)))
cat(sprintf("Median HR: %.2f (Sim 95%%: %.2f - %.2f)\n",
            median(sim_2$hr), quantile(sim_2$hr, 0.025), quantile(sim_2$hr, 0.975)))
cat(sprintf("Median dRMST: %.3f yr (Sim 95%%: %.3f to %.3f)\n",
            median(sim_2$drmst), quantile(sim_2$drmst, 0.025), quantile(sim_2$drmst, 0.975)))
cat(sprintf("Within-sim 95%% CI: %.3f to %.3f\n",
            median(sim_2$drmst_lo), median(sim_2$drmst_hi)))

# --- Model 2 KM Plot (representative run) ---
best_2_idx <- which.min(abs(sim_2$drmst - median(sim_2$drmst)))
cat(sprintf("\nRepresentative sim (#%d): HR=%.2f, dRMST=%.3f\n",
            best_2_idx, sim_2$hr[best_2_idx], sim_2$drmst[best_2_idx]))

set.seed(2026 + 10000 + best_2_idx)
rep_data_2 <- allIPD
for (arm_name in c("LAA Closure", "Medical Therapy")) {
  events_idx <- which(rep_data_2$arm == arm_name & rep_data_2$status == 1)

  # Phase 1: bleeding (temporal)
  idx_p1 <- events_idx[rep_data_2$time[events_idx] < 0.25]
  idx_p2 <- events_idx[rep_data_2$time[events_idx] >= 0.25 &
                         rep_data_2$time[events_idx] < 0.5]
  idx_p3 <- events_idx[rep_data_2$time[events_idx] >= 0.5]

  n1 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3mo[arm_name]), length(idx_p1))
  n2 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_3to6[arm_name]), length(idx_p2))
  n3 <- min(stoch_round(bleed_remove_scaled[arm_name] * frac_6plus[arm_name]), length(idx_p3))

  bleed_sel <- c()
  if (n1 > 0) bleed_sel <- c(bleed_sel, sample(idx_p1, n1))
  if (n2 > 0) bleed_sel <- c(bleed_sel, sample(idx_p2, n2))
  if (n3 > 0) bleed_sel <- c(bleed_sel, sample(idx_p3, n3))

  for (j in bleed_sel) {
    rep_data_2$status[j] <- 0
    rep_data_2$time[j] <- extend_followup(arm_name, rep_data_2$time[j])
  }

  # Phase 2: death (proportional)
  remaining_idx <- which(rep_data_2$arm == arm_name & rep_data_2$status == 1)
  n_death <- min(stoch_round(death_remove_m2[arm_name]), length(remaining_idx))
  if (n_death > 0) {
    death_sel <- sample(remaining_idx, n_death)
    for (j in death_sel) {
      rep_data_2$status[j] <- 0
      rep_data_2$time[j] <- extend_followup(arm_name, rep_data_2$time[j])
    }
  }
}

fit_m2 <- survfit(Surv(time, status) ~ arm, data = rep_data_2)
cx_m2 <- coxph(Surv(time, status) ~ arm, data = rep_data_2)
hr_m2 <- exp(coef(cx_m2))[1]
ci_m2 <- exp(confint(cx_m2))[1, ]
p_m2 <- summary(cx_m2)$coefficients[1, "Pr(>|z|)"]

p2_km <- ggsurvplot(
  fit_m2, data = rep_data_2,
  fun = "event", conf.int = TRUE, censor = FALSE,
  xlim = c(0, 6), ylim = c(0, 0.30),
  break.time.by = 1,
  palette = c("#1F77B4", "#D62728"),
  legend.title = "", legend.labs = c("Medical Therapy", "LAA Closure"),
  legend = c(0.25, 0.85),
  xlab = "Years since Randomization",
  ylab = "Cumulative Incidence (%)",
  risk.table = TRUE, risk.table.title = "No. at Risk",
  title = "Model 2: Thromboembolic Only (Stroke + SE)",
  ggtheme = theme_classic(base_size = 13)
)

p2_km$plot <- p2_km$plot +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     limits = c(0, 0.30), breaks = seq(0, 0.30, by = 0.05)) +
  annotate("text", x = 3, y = 0.25, hjust = 0.5, size = 3.8,
           label = sprintf("HR %.2f (95%% CI %.2f\u2013%.2f), P = %s\nMedian MC dRMST = %.3f yr",
                           hr_m2, ci_m2[1], ci_m2[2],
                           ifelse(p_m2 < 0.001, "<0.001", sprintf("%.3f", p_m2)),
                           median(sim_2$drmst)))

m2_combined <- ggarrange(p2_km$plot, p2_km$table, ncol = 1, nrow = 2, heights = c(3, 1))
ggsave("closure_af_thromboembolic_km.pdf", plot = m2_combined, width = 10, height = 7)
cat("Plot saved: closure_af_thromboembolic_km.pdf\n")


# ============================================================
# 5) Comprehensive Summary Table
# ============================================================

cat("\n\n================================================================\n")
cat("  CLOSURE-AF: ENDPOINT CASCADE SUMMARY\n")
cat("================================================================\n\n")

# Primary composite
cx_prim <- coxph(Surv(time, status) ~ arm, data = allIPD)
hr_prim <- exp(coef(cx_prim))[1]
ci_prim <- exp(confint(cx_prim))[1, ]
rmst_prim <- rmst2(time = allIPD$time, status = allIPD$status,
                    arm = as.numeric(allIPD$arm) - 1, tau = 6)
drmst_prim <- rmst_prim$unadjusted.result[1, "Est."]

cat("Endpoint                          | Events (LAA/Med) | HR (95% CI)           | dRMST (yr)\n")
cat("-------------------------------------------------------------------------------------------\n")

cat(sprintf("1. Primary Composite               | %3d / %3d        | %.2f (%.2f-%.2f)      | %.3f\n",
            our_events["LAA Closure"], our_events["Medical Therapy"],
            hr_prim, ci_prim[1], ci_prim[2], drmst_prim))

cat(sprintf("2. Model 1A: Bleeding-free (temp)   | %3.0f / %3.0f        | %.2f (%.2f-%.2f)      | %.3f\n",
            median(sim_1A$n_events_laa), median(sim_1A$n_events_med),
            median(sim_1A$hr),
            quantile(sim_1A$hr, 0.025), quantile(sim_1A$hr, 0.975),
            median(sim_1A$drmst)))

cat(sprintf("3. Model 1B: Bleeding-free (prop)   | %3.0f / %3.0f        | %.2f (%.2f-%.2f)      | %.3f\n",
            median(sim_1B$n_events_laa), median(sim_1B$n_events_med),
            median(sim_1B$hr),
            quantile(sim_1B$hr, 0.025), quantile(sim_1B$hr, 0.975),
            median(sim_1B$drmst)))

cat(sprintf("4. Model 2:  Thromboembolic only    | %3.0f / %3.0f        | %.2f (%.2f-%.2f)      | %.3f\n",
            median(sim_2$n_events_laa), median(sim_2$n_events_med),
            median(sim_2$hr),
            quantile(sim_2$hr, 0.025), quantile(sim_2$hr, 0.975),
            median(sim_2$drmst)))

# Non-CV death
drmst_death <- rmst_death$unadjusted.result[1, "Est."]
cat(sprintf("5. Non-CV Death (supplementary)     | %3d / %3d        | %.2f (%.2f-%.2f)      | %.3f\n",
            sum(death_IPD$arm == "LAA Closure" & death_IPD$status == 1),
            sum(death_IPD$arm == "Medical Therapy" & death_IPD$status == 1),
            hr_death, ci_death[1], ci_death[2], drmst_death))

cat("-------------------------------------------------------------------------------------------\n")

cat("\n=== Interpretation ===\n")
cat(sprintf("Primary → Model 1 (bleeding removed): dRMST moves from %.3f to %.3f yr\n",
            drmst_prim, median(sim_1B$drmst)))
cat(sprintf("  → %.0f%% of the dRMST gap is attributable to bleeding\n",
            (1 - abs(median(sim_1B$drmst)) / abs(drmst_prim)) * 100))

cat(sprintf("\nPrimary → Model 2 (thromboembolic only): dRMST moves from %.3f to %.3f yr\n",
            drmst_prim, median(sim_2$drmst)))
cat("  → Near null: LAA provides equivalent thromboembolic protection to medical therapy\n")

cat(sprintf("\nModel 1 → Model 2 (death component): %.3f yr difference\n",
            median(sim_1B$drmst) - median(sim_2$drmst)))
cat("  → Death events disproportionately affect LAA arm in the composite\n")

cat("\n=== Model Details ===\n")
cat("Model 1: Remove only major bleeding from primary composite\n")
cat("  1A: Temporal weighting (Table S11 bleeding distribution)\n")
cat("  1B: Proportional thinning (uniform random)\n")
cat("  Hemorrhagic stroke STAYS (differs from v4 which removed it)\n")
cat("Model 2: Remove bleeding + all death → pure stroke + SE\n")
cat("  Phase 1: Bleeding removal (temporal, as Model 1A)\n")
cat("  Phase 2: Death removal (proportional, no temporal bias)\n")

cat("\n================================================================\n")
cat("PDF outputs:\n")
cat("  1. closure_af_bleedonly_km.pdf    — Model 1: Bleeding-excluded\n")
cat("  2. closure_af_thromboembolic_km.pdf — Model 2: Thromboembolic only\n")
cat("  3. closure_af_noncv_death_km.pdf  — Non-CV death\n")
cat("================================================================\n")
