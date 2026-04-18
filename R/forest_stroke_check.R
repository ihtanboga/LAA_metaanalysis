# Sanity check: overall HR per trial + pooled (no time bands)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
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

run_cox <- function(df, label) {
  cx <- coxph(Surv(time, status) ~ arm_h, data = df)
  data.frame(
    label = label,
    hr = as.numeric(exp(coef(cx))[1]),
    lo = as.numeric(exp(confint(cx))[1, 1]),
    hi = as.numeric(exp(confint(cx))[1, 2]),
    ev_ac  = sum(df$arm_h == "AC"     & df$status == 1),
    ev_dev = sum(df$arm_h == "device" & df$status == 1),
    n_ac   = sum(df$arm_h == "AC"),
    n_dev  = sum(df$arm_h == "device")
  )
}

run_pooled <- function(df_list, label) {
  df_all <- do.call(rbind, df_list)
  cx <- coxph(Surv(time, status) ~ arm_h + strata(trial), data = df_all)
  data.frame(
    label = label,
    hr = as.numeric(exp(coef(cx))[1]),
    lo = as.numeric(exp(confint(cx))[1, 1]),
    hi = as.numeric(exp(confint(cx))[1, 2]),
    ev_ac  = sum(df_all$arm_h == "AC"     & df_all$status == 1),
    ev_dev = sum(df_all$arm_h == "device" & df_all$status == 1),
    n_ac   = sum(df_all$arm_h == "AC"),
    n_dev  = sum(df_all$arm_h == "device")
  )
}

overall <- rbind(
  run_cox(ch, "CHAMPION"),
  run_cox(op, "OPTION"),
  run_cox(pr, "PRAGUE-17"),
  run_cox(cl, "CLOSURE-AF"),
  run_pooled(list(ch, op, pr), "Pooled — Contemporary (CHAMPION+OPTION+PRAGUE)"),
  run_pooled(list(ch, op, pr, cl), "Pooled — Expanded (+ CLOSURE-AF)")
)

cat("=== Overall HR (no time bands) — device vs AC ===\n\n")
overall$out <- sprintf("HR %.2f (95%% CI %.2f–%.2f)  | events: %d vs %d  | N: %d vs %d",
                       overall$hr, overall$lo, overall$hi,
                       overall$ev_dev, overall$ev_ac,
                       overall$n_dev, overall$n_ac)
for (i in seq_len(nrow(overall))) {
  cat(sprintf("%-50s  %s\n", overall$label[i], overall$out[i]))
}

# ===== Sum of banded pooled events (should match overall pooled) =====
cat("\n\n=== Sanity check: pooled events reconcile with banded sum ===\n")
cat("Contemporary: overall pooled events = ",
    sum(overall[overall$label == "Pooled — Contemporary (CHAMPION+OPTION+PRAGUE)", c("ev_ac","ev_dev")]),
    "\n")
cat("  Banded sum: A (11+14) + B (58+57) = ", 11+14+58+57, "\n")

cat("Expanded: overall pooled events = ",
    sum(overall[overall$label == "Pooled — Expanded (+ CLOSURE-AF)", c("ev_ac","ev_dev")]),
    "\n")
cat("  Banded sum: A (22+22) + B (72+76) = ", 22+22+72+76, "\n")

write.csv(overall[, 1:8], "~/Desktop/manuel/forest_stroke_overall.csv", row.names = FALSE)
cat("\nSaved: ~/Desktop/manuel/forest_stroke_overall.csv\n")
