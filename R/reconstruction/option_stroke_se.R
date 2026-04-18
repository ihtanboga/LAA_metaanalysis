# ============================================================
# OPTION Trial: Stroke + SE IPD (remove death from composite)
# ============================================================
# Input:  OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism (AC + device)
#         Composite endpoint: Death from any cause + Stroke + SE
# Method: Reconstruct composite IPD via Guyot, then remove death events
#         uniformly (proportional thinning) using paper targets:
#           Device: 29 deaths (3.8% KM at 36 mo)
#           AC:     34 deaths (4.5% KM at 36 mo)
# Output: OPTION_stroke_se.csv (arm, time, status) — status=1 means stroke or SE
# ============================================================

library(readr)
library(dplyr)
library(reconstructKM)
library(survival)

set.seed(2026)

json_dir <- "~/Desktop/manuel/json_csv"
out_dir  <- "~/Desktop/manuel/json_csv/ipd_data"

read_option_arm <- function(csv_file) {
  df <- read_csv(file.path(json_dir, csv_file), comment = "#",
                 show_col_types = FALSE)
  df
}

# ===== 1) Read digitized composite CIF curves =====
ac_raw  <- read_option_arm("OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism_AC.csv")
dev_raw <- read_option_arm("OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism_device.csv")

cat("AC points:", nrow(ac_raw), "| device points:", nrow(dev_raw), "\n")

# Convert cumulative incidence (%) -> survival
prep <- function(df) {
  df %>%
    transmute(time = pmax(0, x), survival = pmax(0, 1 - y / 100)) %>%
    bind_rows(data.frame(time = 0, survival = 1), .) %>%
    group_by(time) %>%
    summarise(survival = min(survival), .groups = "drop") %>%
    arrange(time) %>%
    mutate(survival = cummin(pmin(survival, 1)))
}

ac_clicks  <- prep(ac_raw)
dev_clicks <- prep(dev_raw)

# Trim trailing horizontal tail (required by format_raw_tabs)
trim_tail <- function(df) {
  final <- df$survival[nrow(df)]
  first_final <- which(abs(df$survival - final) < 1e-8)[1]
  df[1:first_final, ]
}
ac_clicks  <- trim_tail(ac_clicks)
dev_clicks <- trim_tail(dev_clicks)

# ===== 2) NAR tables (from JSON) =====
ac_NAR  <- data.frame(time = 0:3, NAR = c(797, 754, 740, 701))
dev_NAR <- data.frame(time = 0:3, NAR = c(803, 764, 748, 717))  # placeholder; verify

# Verify device NAR from json
dev_json_nar <- jsonlite::fromJSON(
  file.path(json_dir, "OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism_device.json")
)
dev_NAR <- data.frame(
  time = dev_json_nar$arms$nar[[1]]$time,
  NAR  = dev_json_nar$arms$nar[[1]]$n
)
cat("Device NAR:\n"); print(dev_NAR)
cat("AC NAR:\n"); print(ac_NAR)

# ===== 3) Reconstruct IPD =====
fix_aug <- function(aug) {
  as <- aug$aug_surv
  if (!"surv" %in% names(as) && "survival" %in% names(as))
    names(as)[names(as) == "survival"] <- "surv"
  aug$aug_surv <- as[order(as$time), ]
  aug
}

ac_aug  <- fix_aug(format_raw_tabs(raw_NAR = ac_NAR,  raw_surv = ac_clicks))
dev_aug <- fix_aug(format_raw_tabs(raw_NAR = dev_NAR, raw_surv = dev_clicks))

ac_recon  <- KM_reconstruct(aug_NAR = ac_aug$aug_NAR,  aug_surv = ac_aug$aug_surv)
dev_recon <- KM_reconstruct(aug_NAR = dev_aug$aug_NAR, aug_surv = dev_aug$aug_surv)

composite_IPD <- rbind(
  data.frame(arm = "AC",     time = ac_recon$IPD_time,  status = ac_recon$IPD_event),
  data.frame(arm = "device", time = dev_recon$IPD_time, status = dev_recon$IPD_event)
)
composite_IPD$arm <- factor(composite_IPD$arm, levels = c("AC", "device"))

n_composite <- c(
  AC     = sum(composite_IPD$arm == "AC"     & composite_IPD$status == 1),
  device = sum(composite_IPD$arm == "device" & composite_IPD$status == 1)
)

cat("\n=== Composite (death + stroke + SE) reconstruction ===\n")
cat(sprintf("Events: AC = %d, device = %d\n", n_composite["AC"], n_composite["device"]))
cat(sprintf("N:      AC = %d, device = %d\n",
            sum(composite_IPD$arm == "AC"),
            sum(composite_IPD$arm == "device")))

# ===== 4) Remove deaths uniformly =====
# Paper targets (Kaplan-Meier at 36 months):
#   Device: 29 deaths (3.8%)
#   AC:     34 deaths (4.5%)
death_target <- c(AC = 34, device = 29)

# Scale to reconstructed event count
# Paper composite events are unknown exactly — estimate from KM at 3y
# AC:     5.9% × 797 = 47 → scale = n_composite["AC"] / 47
# device: 5.3% × 803 = 42 → scale = n_composite["device"] / 42
paper_composite_est <- c(AC = round(0.059 * 797), device = round(0.053 * 803))
scale_ratio <- n_composite / paper_composite_est
death_scaled <- death_target * scale_ratio

cat(sprintf("\nEstimated paper composite: AC = %d, device = %d\n",
            paper_composite_est["AC"], paper_composite_est["device"]))
cat(sprintf("Scale ratio: AC = %.3f, device = %.3f\n",
            scale_ratio["AC"], scale_ratio["device"]))
cat(sprintf("Death removal target (scaled): AC = %.1f, device = %.1f\n",
            death_scaled["AC"], death_scaled["device"]))

# Extend follow-up helper
censor_times <- list()
for (a in c("AC", "device")) {
  censor_times[[a]] <- composite_IPD$time[composite_IPD$arm == a & composite_IPD$status == 0]
}
extend_fu <- function(a, t) {
  avail <- censor_times[[a]][censor_times[[a]] > t]
  if (length(avail) > 0) sample(avail, 1) else max(composite_IPD$time[composite_IPD$arm == a])
}
stoch_round <- function(x) floor(x) + rbinom(1, 1, x - floor(x))

run_sim <- function() {
  sim <- composite_IPD
  for (a in c("AC", "device")) {
    ev_idx <- which(sim$arm == a & sim$status == 1)
    n_rm   <- min(stoch_round(death_scaled[a]), length(ev_idx))
    if (n_rm > 0) {
      sel <- sample(ev_idx, n_rm)
      for (j in sel) {
        sim$status[j] <- 0
        sim$time[j]   <- extend_fu(a, sim$time[j])
      }
    }
  }
  sim
}

# ===== 5) MC + representative pick =====
n_sim <- 1000
hrs <- numeric(n_sim); drmsts <- numeric(n_sim)
n_ev_ac <- integer(n_sim); n_ev_dev <- integer(n_sim)

library(survRM2)
for (i in seq_len(n_sim)) {
  s <- run_sim()
  n_ev_ac[i]  <- sum(s$arm == "AC"     & s$status == 1)
  n_ev_dev[i] <- sum(s$arm == "device" & s$status == 1)
  cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = s), error = function(e) NULL)
  if (!is.null(cx)) hrs[i] <- exp(coef(cx))[1]
  r <- tryCatch(rmst2(s$time, s$status, as.numeric(s$arm) - 1, tau = 3),
                error = function(e) NULL)
  if (!is.null(r)) drmsts[i] <- r$unadjusted.result[1, "Est."]
}

cat(sprintf("\n=== MC (n=%d) ===\n", n_sim))
cat(sprintf("Stroke+SE events per sim: AC %.0f, device %.0f\n",
            median(n_ev_ac), median(n_ev_dev)))
cat(sprintf("Median HR (device vs AC): %.2f (sim 95%%: %.2f-%.2f)\n",
            median(hrs), quantile(hrs, 0.025), quantile(hrs, 0.975)))
cat(sprintf("Median dRMST (tau=3): %.3f yr\n", median(drmsts)))

best_idx <- which.min(abs(drmsts - median(drmsts)))
set.seed(2026 + best_idx)
stroke_IPD <- run_sim()

cat(sprintf("\nRepresentative sim #%d: AC %d events, device %d events\n",
            best_idx,
            sum(stroke_IPD$arm == "AC"     & stroke_IPD$status == 1),
            sum(stroke_IPD$arm == "device" & stroke_IPD$status == 1)))

out_path <- file.path(out_dir, "OPTION_stroke_se.csv")
write.csv(stroke_IPD, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
