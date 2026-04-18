# ============================================================
# CLOSURE-AF: Pure Bleeding Endpoint IPD (bleeding.csv)
# ============================================================
# Inverse of closure_af_bleeding_free.R:
#   bleeding_free.R: REMOVES bleeding from primary composite → bleeding-free
#   THIS:           KEEPS only bleeding from primary composite → bleeding-only
#
# Logic:
#   1. Load primary composite IPD (closure_af_allIPD.csv)
#   2. Among events (status=1), probabilistically flag ~70 LAA + ~61 Medical
#      as bleeding, using Table S11 temporal distribution:
#         LAA:    <3mo 57.1%, 3-6mo 10.0%, ≥6mo 32.9%
#         Med:    <3mo 27.9%, 3-6mo 9.8%,  ≥6mo 62.3%
#   3. Non-bleeding events → status=0, extend follow-up to admin-censor time
#   4. Output: bleeding.csv (arm, time, status)  — status=1 means bleeding
#
# Targets (scaled to reconstructed IPD):
#   Paper:        LAA 70 / Medical 61
#   Our ratios:   LAA 147/155=0.948, Med 121/127=0.953
#   Scaled:       LAA 66.4 / Medical 58.1
# ============================================================

library(dplyr)
library(survival)
library(survRM2)

set.seed(2026)

# ===== Load primary composite IPD =====
allIPD <- read.csv("closure_af_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Medical Therapy", "LAA Closure"))

cat("=== Loaded primary composite IPD ===\n")
cat(sprintf("N = %d (LAA %d, Med %d)\n",
            nrow(allIPD),
            sum(allIPD$arm == "LAA Closure"),
            sum(allIPD$arm == "Medical Therapy")))
cat(sprintf("Primary events: LAA %d, Med %d (paper 155/127)\n",
            sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
            sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)))

# ===== Scaling & temporal weights (Table S11) =====
paper_events <- c("LAA Closure" = 155, "Medical Therapy" = 127)
our_events   <- c(
  "LAA Closure"    = sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
  "Medical Therapy" = sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)
)
scale_ratio <- our_events / paper_events

bleed_paper  <- c("LAA Closure" = 70, "Medical Therapy" = 61)
bleed_scaled <- bleed_paper * scale_ratio

bleed_3mo   <- c("LAA Closure" = 40, "Medical Therapy" = 17)
bleed_3to6  <- c("LAA Closure" = 7,  "Medical Therapy" = 6)
bleed_6plus <- c("LAA Closure" = 23, "Medical Therapy" = 38)
bleed_tot   <- c("LAA Closure" = 70, "Medical Therapy" = 61)
frac_3mo   <- bleed_3mo / bleed_tot
frac_3to6  <- bleed_3to6 / bleed_tot
frac_6plus <- bleed_6plus / bleed_tot

cat(sprintf("\nBleeding targets (scaled): LAA %.1f, Med %.1f\n",
            bleed_scaled["LAA Closure"], bleed_scaled["Medical Therapy"]))

# ===== Extend follow-up for non-bleeding events =====
censor_times <- list()
for (a in c("LAA Closure", "Medical Therapy")) {
  censor_times[[a]] <- allIPD$time[allIPD$arm == a & allIPD$status == 0]
}
extend_followup <- function(a, t) {
  avail <- censor_times[[a]][censor_times[[a]] > t]
  if (length(avail) > 0) sample(avail, 1) else max(allIPD$time[allIPD$arm == a])
}
stoch_round <- function(x) floor(x) + rbinom(1, 1, x - floor(x))

# ===== MC: find representative bleeding-only dataset =====
n_sim <- 1000

run_sim <- function() {
  sim <- allIPD
  for (a in c("LAA Closure", "Medical Therapy")) {
    ev_idx <- which(sim$arm == a & sim$status == 1)

    # Temporal-weighted selection of bleeding events
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

    # Non-selected events are NON-bleeding → censor + extend follow-up
    non_bleed <- setdiff(ev_idx, bleed_sel)
    for (j in non_bleed) {
      sim$status[j] <- 0
      sim$time[j] <- extend_followup(a, sim$time[j])
    }
  }
  sim
}

# Summary across simulations
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
            median(n_ev_laa), bleed_scaled["LAA Closure"],
            median(n_ev_med), bleed_scaled["Medical Therapy"]))
cat(sprintf("Median HR: %.2f (sim 95%%: %.2f-%.2f)\n",
            median(hrs), quantile(hrs, 0.025), quantile(hrs, 0.975)))
cat(sprintf("Median dRMST: %.3f yr (sim 95%%: %.3f to %.3f)\n",
            median(drmsts), quantile(drmsts, 0.025), quantile(drmsts, 0.975)))

# Pick representative sim (median dRMST) and regenerate
best_idx <- which.min(abs(drmsts - median(drmsts)))
set.seed(2026 + best_idx)
bleeding_IPD <- run_sim()

cat(sprintf("\nRepresentative sim #%d:\n", best_idx))
cat(sprintf("  Events: LAA %d, Med %d\n",
            sum(bleeding_IPD$arm == "LAA Closure" & bleeding_IPD$status == 1),
            sum(bleeding_IPD$arm == "Medical Therapy" & bleeding_IPD$status == 1)))

cx <- coxph(Surv(time, status) ~ arm, data = bleeding_IPD)
cat(sprintf("  HR = %.2f (95%% CI %.2f-%.2f)\n",
            exp(coef(cx))[1],
            exp(confint(cx))[1, 1], exp(confint(cx))[1, 2]))

write.csv(bleeding_IPD, "bleeding.csv", row.names = FALSE)
cat("\nSaved: bleeding.csv (bleeding-only endpoint)\n")
