# ============================================================
# CLOSURE-AF: Stroke + Systemic Embolism Endpoint IPD (stroke_se.csv)
# ============================================================
# Logic (Model 2 from closure_af_death_free.R):
#   1. Load primary composite IPD (closure_af_allIPD.csv)
#   2. Phase 1: Remove bleeding events (temporal-weighted, Table S11)
#         LAA 70 × 0.948 = 66.4, Med 61 × 0.953 = 58.1
#   3. Phase 2: Remove CV/unexplained death events (proportional)
#         death_remove = total_remove − bleeding_remove
#         LAA: 119.5 − 66.4 = 53.1
#         Med:  94.3 − 58.1 = 36.2
#   4. Remaining events = stroke (all types) + SE
#         Target: LAA 29 × 0.948 = 27.5, Med 28 × 0.953 = 26.7
#   5. Output: stroke_se.csv (arm, time, status) — status=1 means stroke or SE
#
# Paper reference (Table 2):
#   Stroke or SE: LAA 29 (18 ischemic + 10 hemorrhagic + 3 SE), Med 28 (15 + 13 + 1)
# ============================================================

library(dplyr)
library(survival)
library(survRM2)

set.seed(2026)

# ===== Load primary composite IPD =====
allIPD <- read.csv("closure_af_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Medical Therapy", "LAA Closure"))

cat("=== Loaded primary composite IPD ===\n")
cat(sprintf("N = %d (LAA %d, Med %d)\n", nrow(allIPD),
            sum(allIPD$arm == "LAA Closure"),
            sum(allIPD$arm == "Medical Therapy")))
cat(sprintf("Primary events: LAA %d, Med %d (paper 155/127)\n",
            sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
            sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)))

# ===== Scaling =====
paper_events <- c("LAA Closure" = 155, "Medical Therapy" = 127)
our_events   <- c(
  "LAA Closure"    = sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
  "Medical Therapy" = sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)
)
scale_ratio <- our_events / paper_events

bleed_paper   <- c("LAA Closure" = 70, "Medical Therapy" = 61)
bleed_scaled  <- bleed_paper * scale_ratio

stroke_paper  <- c("LAA Closure" = 29, "Medical Therapy" = 28)
stroke_scaled <- stroke_paper * scale_ratio

total_remove <- our_events - stroke_scaled
death_remove <- total_remove - bleed_scaled

cat(sprintf("\nTargets (scaled):\n"))
cat(sprintf("  Bleeding remove: LAA %.1f, Med %.1f\n",
            bleed_scaled["LAA Closure"], bleed_scaled["Medical Therapy"]))
cat(sprintf("  Death remove:    LAA %.1f, Med %.1f\n",
            death_remove["LAA Closure"], death_remove["Medical Therapy"]))
cat(sprintf("  Stroke+SE keep:  LAA %.1f, Med %.1f\n",
            stroke_scaled["LAA Closure"], stroke_scaled["Medical Therapy"]))

# Table S11 bleeding temporal
bleed_tot   <- c("LAA Closure" = 70, "Medical Therapy" = 61)
frac_3mo   <- c("LAA Closure" = 40, "Medical Therapy" = 17) / bleed_tot
frac_3to6  <- c("LAA Closure" = 7,  "Medical Therapy" = 6)  / bleed_tot
frac_6plus <- c("LAA Closure" = 23, "Medical Therapy" = 38) / bleed_tot

censor_times <- list()
for (a in c("LAA Closure", "Medical Therapy")) {
  censor_times[[a]] <- allIPD$time[allIPD$arm == a & allIPD$status == 0]
}
extend_followup <- function(a, t) {
  avail <- censor_times[[a]][censor_times[[a]] > t]
  if (length(avail) > 0) sample(avail, 1) else max(allIPD$time[allIPD$arm == a])
}
stoch_round <- function(x) floor(x) + rbinom(1, 1, x - floor(x))

# ===== MC: build stroke+SE dataset (remove bleeding + death) =====
run_sim <- function() {
  sim <- allIPD
  for (a in c("LAA Closure", "Medical Therapy")) {

    # Phase 1: bleeding removal (temporal, Table S11)
    ev_idx <- which(sim$arm == a & sim$status == 1)
    idx_p1 <- ev_idx[sim$time[ev_idx] <  0.25]
    idx_p2 <- ev_idx[sim$time[ev_idx] >= 0.25 & sim$time[ev_idx] < 0.5]
    idx_p3 <- ev_idx[sim$time[ev_idx] >= 0.5]

    n1 <- min(stoch_round(bleed_scaled[a] * frac_3mo[a]),   length(idx_p1))
    n2 <- min(stoch_round(bleed_scaled[a] * frac_3to6[a]),  length(idx_p2))
    n3 <- min(stoch_round(bleed_scaled[a] * frac_6plus[a]), length(idx_p3))

    bleed_sel <- c()
    if (n1 > 0) bleed_sel <- c(bleed_sel, sample(idx_p1, n1))
    if (n2 > 0) bleed_sel <- c(bleed_sel, sample(idx_p2, n2))
    if (n3 > 0) bleed_sel <- c(bleed_sel, sample(idx_p3, n3))

    for (j in bleed_sel) {
      sim$status[j] <- 0
      sim$time[j]   <- extend_followup(a, sim$time[j])
    }

    # Phase 2: death removal (proportional)
    rem_idx <- which(sim$arm == a & sim$status == 1)
    n_death <- min(stoch_round(death_remove[a]), length(rem_idx))
    if (n_death > 0) {
      death_sel <- sample(rem_idx, n_death)
      for (j in death_sel) {
        sim$status[j] <- 0
        sim$time[j]   <- extend_followup(a, sim$time[j])
      }
    }
  }
  sim
}

# ===== Run MC & pick representative =====
n_sim <- 1000
hrs <- numeric(n_sim); drmsts <- numeric(n_sim)
n_ev_laa <- integer(n_sim); n_ev_med <- integer(n_sim)

for (i in seq_len(n_sim)) {
  s <- run_sim()
  n_ev_laa[i] <- sum(s$arm == "LAA Closure" & s$status == 1)
  n_ev_med[i] <- sum(s$arm == "Medical Therapy" & s$status == 1)
  cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = s), error = function(e) NULL)
  if (!is.null(cx)) hrs[i] <- exp(coef(cx))[1]
  r <- tryCatch(rmst2(s$time, s$status, as.numeric(s$arm) - 1, tau = 6),
                error = function(e) NULL)
  if (!is.null(r)) drmsts[i] <- r$unadjusted.result[1, "Est."]
}

cat(sprintf("\n=== MC summary (n=%d) ===\n", n_sim))
cat(sprintf("Events per sim: LAA %.0f (target %.1f), Med %.0f (target %.1f)\n",
            median(n_ev_laa), stroke_scaled["LAA Closure"],
            median(n_ev_med), stroke_scaled["Medical Therapy"]))
cat(sprintf("Median HR: %.2f (sim 95%%: %.2f-%.2f)\n",
            median(hrs), quantile(hrs, 0.025), quantile(hrs, 0.975)))
cat(sprintf("Median dRMST: %.3f yr (sim 95%%: %.3f to %.3f)\n",
            median(drmsts), quantile(drmsts, 0.025), quantile(drmsts, 0.975)))

best_idx <- which.min(abs(drmsts - median(drmsts)))
set.seed(2026 + 50000 + best_idx)
stroke_IPD <- run_sim()

cat(sprintf("\nRepresentative sim #%d:\n", best_idx))
cat(sprintf("  Events: LAA %d, Med %d\n",
            sum(stroke_IPD$arm == "LAA Closure" & stroke_IPD$status == 1),
            sum(stroke_IPD$arm == "Medical Therapy" & stroke_IPD$status == 1)))

cx <- coxph(Surv(time, status) ~ arm, data = stroke_IPD)
cat(sprintf("  HR = %.2f (95%% CI %.2f-%.2f)\n",
            exp(coef(cx))[1],
            exp(confint(cx))[1, 1], exp(confint(cx))[1, 2]))

write.csv(stroke_IPD, "stroke_se.csv", row.names = FALSE)
cat("\nSaved: stroke_se.csv (stroke + SE endpoint)\n")
