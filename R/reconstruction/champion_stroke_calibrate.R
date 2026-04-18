# ============================================================
# CHAMPION All Stroke — NAR-calibrated reconstruction
# Target: device = 50 events, AC = 50 events
# Method: binary search on NAR scale factor until reconstructed
#         event count matches target (same approach as CLOSURE non-CV death)
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

# Linear interpolation to densify sparse digitized curves
interpolate_clicks <- function(df, n_points = 80) {
  tg <- seq(0, max(df$time), length.out = n_points)
  sg <- approx(df$time, df$survival, xout = tg, rule = 2)$y
  data.frame(time = tg, survival = cummin(pmin(sg, 1)))
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

# ===== NAR calibration (binary search, scale NAR values upward) =====
calibrate <- function(clicks, base_nar, N_total, target_ev, arm_label) {
  # Scale NAR values: higher k → NAR stays closer to N_total → fewer
  # drops → Guyot allocates MORE events to match the KM curve shape.
  k_lo <- 0.3
  k_hi <- 10.0
  best <- NULL

  for (iter in 1:40) {
    k <- (k_lo + k_hi) / 2

    nar_scaled <- base_nar
    nar_scaled$NAR[-1] <- pmin(N_total, round(base_nar$NAR[-1] * k))
    # Monotonic non-increasing
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

  if (is.null(best)) stop(paste("Calibration failed for", arm_label))
  best
}

# ===== Device arm =====
cat("=== CHAMPION device arm ===\n")
dev_csv  <- read_csv(file.path(json_dir, "champion_ITT_AllStroke_device.csv"),
                     comment = "#", show_col_types = FALSE)
dev_json <- fromJSON(file.path(json_dir, "champion_ITT_AllStroke_device.json"))

dev_clicks <- interpolate_clicks(trim_tail(prep_clicks(dev_csv)), 80)
dev_nar_raw <- data.frame(time = dev_json$arms$nar[[1]]$time,
                           NAR = dev_json$arms$nar[[1]]$n)
dev_clicks <- dev_clicks %>% filter(time <= max(dev_nar_raw$time))

cat("Original NAR: "); print(dev_nar_raw)

dev_cal <- calibrate(dev_clicks, dev_nar_raw, N_total = 1499,
                     target_ev = 50, "device")

cat("Calibrated NAR:\n"); print(dev_cal$NAR)

# ===== AC arm =====
cat("\n=== CHAMPION AC arm ===\n")
ac_csv  <- read_csv(file.path(json_dir, "champion_ITT_AllStroke_AC.csv"),
                    comment = "#", show_col_types = FALSE)
ac_json <- fromJSON(file.path(json_dir, "champion_ITT_AllStroke_AC.json"))

ac_clicks <- interpolate_clicks(trim_tail(prep_clicks(ac_csv)), 80)
ac_nar_raw <- data.frame(time = ac_json$arms$nar[[1]]$time,
                          NAR = ac_json$arms$nar[[1]]$n)
ac_clicks <- ac_clicks %>% filter(time <= max(ac_nar_raw$time))

cat("Original NAR: "); print(ac_nar_raw)

ac_cal <- calibrate(ac_clicks, ac_nar_raw, N_total = 1501,
                    target_ev = 33, "AC")

cat("Calibrated NAR:\n"); print(ac_cal$NAR)

# ===== Combine =====
combined <- rbind(
  data.frame(arm = "AC",
             time = ac_cal$recon$IPD_time,
             status = ac_cal$recon$IPD_event),
  data.frame(arm = "device",
             time = dev_cal$recon$IPD_time,
             status = dev_cal$recon$IPD_event)
)

cat(sprintf("\n=== Final CHAMPION all_stroke IPD ===\n"))
cat(sprintf("AC:     N=%d, events=%d (target 50)\n",
            sum(combined$arm == "AC"),
            sum(combined$arm == "AC" & combined$status == 1)))
cat(sprintf("device: N=%d, events=%d (target 50)\n",
            sum(combined$arm == "device"),
            sum(combined$arm == "device" & combined$status == 1)))

out_file <- file.path(out_dir, "CHAMPION_all_stroke.csv")
write.csv(combined, out_file, row.names = FALSE)
cat("\nOverwritten:", out_file, "\n")
