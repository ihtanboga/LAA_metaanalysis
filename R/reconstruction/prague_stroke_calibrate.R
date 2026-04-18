# ============================================================
# PRAGUE-17 All Stroke/TIA — NAR-calibrated reconstruction
# Arms were mislabeled during digitization; swap and recalibrate.
# Targets: AC = 13 events, device (LAA closure) = 14 events
# Source files (with swap applied):
#   "device.csv" (final CIF 9.56%)  → true AC  (target 13)
#   "AC.csv"     (final CIF 14.93%) → true device (target 14)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
  library(reconstructKM)
})

set.seed(2026)

json_dir <- "~/Desktop/manuel/json_csv"
out_dir  <- "~/Desktop/manuel/json_csv/ipd_data"

prep_clicks <- function(df) {
  df %>%
    transmute(time = pmax(0, x), survival = pmax(0, 1 - y / 100)) %>%
    filter(time > 1e-8) %>%
    bind_rows(data.frame(time = 0, survival = 1), .) %>%
    arrange(time) %>%
    group_by(time) %>%
    summarise(survival = min(survival), .groups = "drop") %>%
    arrange(time) %>%
    mutate(survival = cummin(pmin(survival, 1)))
}

trim_tail <- function(df) {
  final <- df$survival[nrow(df)]
  first_final <- which(abs(df$survival - final) < 1e-8)[1]
  df[1:first_final, ]
}

fix_aug <- function(aug) {
  as <- aug$aug_surv
  if (!"surv" %in% names(as) && "survival" %in% names(as))
    names(as)[names(as) == "survival"] <- "surv"
  aug$aug_surv <- as[order(as$time), ]
  aug
}

interpolate_clicks <- function(df, n_points = 80) {
  tg <- seq(0, max(df$time), length.out = n_points)
  sg <- approx(df$time, df$survival, xout = tg, rule = 2)$y
  data.frame(time = tg, survival = cummin(pmin(sg, 1)))
}

calibrate <- function(clicks, base_nar, N_total, target_ev, arm_label) {
  k_lo <- 0.1; k_hi <- 10.0
  best <- NULL
  for (iter in 1:40) {
    k <- (k_lo + k_hi) / 2
    nar_scaled <- base_nar
    nar_scaled$NAR[-1] <- pmin(N_total, round(base_nar$NAR[-1] * k))
    for (i in 2:nrow(nar_scaled)) {
      nar_scaled$NAR[i] <- min(nar_scaled$NAR[i], nar_scaled$NAR[i-1])
    }
    aug <- tryCatch({
      a <- format_raw_tabs(raw_NAR = nar_scaled, raw_surv = clicks)
      fix_aug(a)
    }, error = function(e) NULL)
    if (is.null(aug)) { k_hi <- k; next }
    recon <- tryCatch(
      KM_reconstruct(aug_NAR = aug$aug_NAR, aug_surv = aug$aug_surv),
      error = function(e) NULL)
    if (is.null(recon)) { k_hi <- k; next }
    events <- sum(recon$IPD_event)
    best <- list(NAR = nar_scaled, recon = recon, events = events, k = k)
    if (abs(events - target_ev) <= 1) {
      cat(sprintf("  %s: k=%.3f → events=%d (target=%d) ✓\n",
                  arm_label, k, events, target_ev))
      break
    }
    if (events < target_ev) k_lo <- k else k_hi <- k
  }
  best
}

# ===== SWAP: file "device.csv" → true AC (target 13) =====
cat("=== True AC (from file tagged 'device') ===\n")
ac_csv  <- read_csv(file.path(json_dir,
    "PRAGUE17_ITT_allstroketransientischemicattack_device.csv"),
    comment = "#", show_col_types = FALSE)
ac_json <- fromJSON(file.path(json_dir,
    "PRAGUE17_ITT_allstroketransientischemicattack_device.json"))

ac_clicks <- interpolate_clicks(trim_tail(prep_clicks(ac_csv)), 80)
ac_nar <- data.frame(time = ac_json$arms$nar[[1]]$time,
                      NAR = ac_json$arms$nar[[1]]$n)
ac_clicks <- ac_clicks %>% filter(time <= max(ac_nar$time))
cat("Final CIF:", round((1 - min(ac_clicks$survival)) * 100, 2), "%\n")
cat("Original NAR:\n"); print(ac_nar)

ac_cal <- calibrate(ac_clicks, ac_nar, N_total = max(ac_nar$NAR),
                    target_ev = 13, "AC")
cat("Calibrated NAR:\n"); print(ac_cal$NAR)

# ===== SWAP: file "AC.csv" → true device (target 14) =====
cat("\n=== True device (from file tagged 'AC') ===\n")
dev_csv  <- read_csv(file.path(json_dir,
    "PRAGUE17_ITT_allstroketransientischemicattack_AC.csv"),
    comment = "#", show_col_types = FALSE)
dev_json <- fromJSON(file.path(json_dir,
    "PRAGUE17_ITT_allstroketransientischemicattack_AC.json"))

dev_clicks <- interpolate_clicks(trim_tail(prep_clicks(dev_csv)), 80)
dev_nar <- data.frame(time = dev_json$arms$nar[[1]]$time,
                       NAR = dev_json$arms$nar[[1]]$n)
dev_clicks <- dev_clicks %>% filter(time <= max(dev_nar$time))
cat("Final CIF:", round((1 - min(dev_clicks$survival)) * 100, 2), "%\n")
cat("Original NAR:\n"); print(dev_nar)

dev_cal <- calibrate(dev_clicks, dev_nar, N_total = max(dev_nar$NAR),
                     target_ev = 14, "device")
cat("Calibrated NAR:\n"); print(dev_cal$NAR)

# ===== Combine =====
combined <- rbind(
  data.frame(arm = "AC",
             time = ac_cal$recon$IPD_time,
             status = ac_cal$recon$IPD_event),
  data.frame(arm = "device",
             time = dev_cal$recon$IPD_time,
             status = dev_cal$recon$IPD_event)
)

cat(sprintf("\n=== Final PRAGUE-17 all_stroke_TIA IPD (swap applied) ===\n"))
cat(sprintf("AC:     N=%d, events=%d (target 13)\n",
            sum(combined$arm == "AC"),
            sum(combined$arm == "AC" & combined$status == 1)))
cat(sprintf("device: N=%d, events=%d (target 14)\n",
            sum(combined$arm == "device"),
            sum(combined$arm == "device" & combined$status == 1)))

out_file <- file.path(out_dir, "PRAGUE17_all_stroke_TIA.csv")
write.csv(combined, out_file, row.names = FALSE)
cat("\nOverwritten:", out_file, "\n")
