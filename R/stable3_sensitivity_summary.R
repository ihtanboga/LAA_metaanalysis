###############################################################################
# Supplemental Table 3 — Sensitivity Analyses Summary (REGENERATED, harmonized)
#
# This is the ACTUAL "Supplemental Table 3" in 03_OnlineAppendix.docx
# (title: "Sensitivity Analyses Summary"). It was previously typed from a stale
# pre-/intermediate-harmonization input, so its IS+SE rows disagreed with main
# Table 3. This script rebuilds it from the FINAL harmonized IS+SE counts (the
# same ones that feed main Table 3 = data/endpoint_events.xlsx).
#
# Arm sizes: CHAMPION-AF uses n_device=1499, n_control=1501 EXACTLY as in
# endpoint_events.xlsx, so the regenerated rows reproduce main Table 3 to the
# decimal (Primary Contemporary RR 1.29; All Trials RR 1.30). (trial_design.csv
# lists 1501/1499; the 2-patient swap is immaterial except at the rounding edge.)
#
# Scope per author decision: ONLY this table is regenerated. endpoint_events.csv
# and all other outputs (Figure 2, BRR/Figure 5) are left untouched.
#
# IS+SE  -> RR, Mantel-Haenszel, event-count based (THIS is what changed)
# NonProcBleed -> HR, inverse-variance (unchanged; counts identical to csv)
#
# Output: supplement/stable3_sensitivity_summary.csv  (+ console, with old vs new)
###############################################################################

suppressMessages(library(meta))

# ---- FINAL harmonized IS+SE counts (from endpoint_events.xlsx) + correct N ----
# study, ev_device, n_device, ev_control, n_control, pool, single-study effect
isse <- data.frame(
  study      = c("CHAMPION-AF","OPTION","PRAGUE-17","CLOSURE-AF","PROTECT-AF/PREVAIL"),
  ev_device  = c(50, 11, 15, 27, 45),
  n_device   = c(1499, 803, 201, 446, 732),   # CHAMPION 1499 per xlsx (matches main Table 3)
  ev_control = c(35, 11, 13, 25, 14),
  n_control  = c(1501, 797, 201, 442, 382),   # CHAMPION 1501 per xlsx
  pool       = c("contemporary","contemporary","contemporary","expanded","historical"),
  eff = c(NA, 0.97, 1.38, NA, 1.71), lo = c(NA,0.42,0.63,NA,0.94), hi = c(NA,2.25,3.03,NA,3.11),
  stringsAsFactors = FALSE)

# ---- Non-procedural bleeding (HR-based; identical in csv & xlsx) --------------
d <- read.csv("data/endpoint_events.csv", stringsAsFactors = FALSE)
bleed <- d[d$endpoint == "NonProcBleed",
           c("study","pool","effect","effect_lo","effect_hi")]

# ---- pooling helpers --------------------------------------------------------
rr_re <- function(s) {  # RR random-effects (MH, REML, HKSJ)
  m <- metabin(s$ev_device, s$n_device, s$ev_control, s$n_control, studlab = s$study,
               sm = "RR", method = "MH", method.tau = "REML",
               method.random.ci = "HK", incr = 0.5)
  c(est = exp(m$TE.random), lo = exp(m$lower.random), hi = exp(m$upper.random),
    p = m$pval.random, i2 = m$I2 * 100)
}
rr_fe <- function(s) {  # RR fixed/common-effect (MH)
  m <- metabin(s$ev_device, s$n_device, s$ev_control, s$n_control, studlab = s$study,
               sm = "RR", method = "MH", incr = 0.5)
  c(est = exp(m$TE.common), lo = exp(m$lower.common), hi = exp(m$upper.common),
    p = m$pval.common, i2 = m$I2 * 100)
}
hr_re <- function(s) {  # HR random-effects (IV, REML, HKSJ)
  s <- s[!is.na(s$effect), ]
  se <- (log(s$effect_hi) - log(s$effect_lo)) / (2 * qnorm(0.975))
  m <- metagen(log(s$effect), se, studlab = s$study, sm = "HR",
               method.tau = "REML", method.random.ci = "HK")
  c(est = exp(m$TE.random), lo = exp(m$lower.random), hi = exp(m$upper.random),
    p = m$pval.random, i2 = m$I2 * 100)
}
hr_fe <- function(s) {  # HR common-effect (IV)
  s <- s[!is.na(s$effect), ]
  se <- (log(s$effect_hi) - log(s$effect_lo)) / (2 * qnorm(0.975))
  m <- metagen(log(s$effect), se, studlab = s$study, sm = "HR")
  c(est = exp(m$TE.common), lo = exp(m$lower.common), hi = exp(m$upper.common),
    p = m$pval.common, i2 = m$I2 * 100)
}
fmt <- function(v, sm) sprintf("%s %.2f (%.2f-%.2f)", sm, v["est"], v["lo"], v["hi"])
row <- function(ep, lab, k, v, sm) data.frame(
  Endpoint = ep, Sensitivity = lab, k = k, Estimate = fmt(v, sm),
  P = ifelse(is.na(v["p"]), "—", sprintf("%.4f", v["p"])),
  I2 = ifelse(is.na(v["i2"]), "—", sprintf("%.0f", v["i2"])), stringsAsFactors = FALSE)
single <- function(ep, lab, eff, lo, hi, sm) data.frame(
  Endpoint = ep, Sensitivity = lab, k = 1,
  Estimate = sprintf("%s %.2f (%.2f-%.2f)", sm, eff, lo, hi),
  P = "—", I2 = "—", stringsAsFactors = FALSE)

# ---- IS+SE block (RR) -------------------------------------------------------
con <- isse[isse$pool == "contemporary", ]
exp4 <- isse[isse$pool %in% c("contemporary","expanded"), ]
out <- rbind(
  row("IS+SE", "Primary (RE, REML, HKSJ), Contemporary", 3, rr_re(con), "RR"),
  row("IS+SE", "LOO Contemporary: -CHAMPION-AF", 2, rr_re(con[con$study != "CHAMPION-AF", ]), "RR"),
  row("IS+SE", "LOO Contemporary: -OPTION",      2, rr_re(con[con$study != "OPTION", ]), "RR"),
  row("IS+SE", "LOO Contemporary: -PRAGUE-17",   2, rr_re(con[con$study != "PRAGUE-17", ]), "RR"),
  row("IS+SE", "Fixed-effect (Mantel-Haenszel), Contemporary", 3, rr_fe(con), "RR"),
  row("IS+SE", "FLX-only (CHAMPION-AF + OPTION)", 2, rr_re(con[con$study %in% c("CHAMPION-AF","OPTION"), ]), "RR"),
  single("IS+SE", "Historical Warfarin Benchmark (separate)", 1.71, 0.94, 3.11, "RR"),
  row("IS+SE", "Expanded Contemporary (RE, REML, HKSJ)", 4, rr_re(exp4), "RR"),
  row("IS+SE", "All Trials (RE, REML, HKSJ)", 5, rr_re(isse), "RR")
)

# ---- Non-proc bleeding block (HR) -------------------------------------------
cb  <- bleed[bleed$pool == "contemporary", ]
eb4 <- bleed[bleed$pool %in% c("contemporary","expanded"), ]
out <- rbind(out,
  row("NonProcBleed", "Primary (RE, REML, HKSJ), Contemporary", 3, hr_re(cb), "HR"),
  row("NonProcBleed", "LOO Contemporary: -CHAMPION-AF", 2, hr_re(cb[cb$study != "CHAMPION-AF", ]), "HR"),
  row("NonProcBleed", "LOO Contemporary: -OPTION",      2, hr_re(cb[cb$study != "OPTION", ]), "HR"),
  row("NonProcBleed", "LOO Contemporary: -PRAGUE-17",   2, hr_re(cb[cb$study != "PRAGUE-17", ]), "HR"),
  row("NonProcBleed", "Fixed-effect (IV), Contemporary", 3, hr_fe(cb), "HR"),
  row("NonProcBleed", "FLX-only (CHAMPION-AF + OPTION)", 2, hr_re(cb[cb$study %in% c("CHAMPION-AF","OPTION"), ]), "HR"),
  single("NonProcBleed", "Historical Warfarin Benchmark (separate)", 0.48, 0.32, 0.71, "HR"),
  row("NonProcBleed", "Expanded Contemporary (RE, REML, HKSJ)", 4, hr_re(eb4), "HR"),
  row("NonProcBleed", "All Trials (RE, REML, HKSJ)", 5, hr_re(bleed), "HR")
)

write.csv(out, "supplement/stable3_sensitivity_summary.csv", row.names = FALSE)

cat("Supplemental Table 3 (Sensitivity Analyses Summary) — REGENERATED (harmonized)\n")
cat(strrep("=", 92), "\n")
print(out, row.names = FALSE)
cat("\nKey check — IS+SE Primary Contemporary should now match main Table 3 = RR 1.29 (0.83-1.99);\n")
cat("All Trials should be RR 1.30 (1.01-1.67), P~0.046.\n")
