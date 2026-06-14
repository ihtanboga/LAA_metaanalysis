###############################################################################
# RBA Analysis — Severity-weighted benefit-risk ratio (BRR), Kaul framework
#
# Computes per-1,000 patient-year (PY) absolute rate differences DIRECTLY from
# the harmonized event data (no pre-baked rba_inputs.csv), then applies three
# pre-specified severity-weighting scenarios.
#
# Method (verified to reproduce the manuscript estimates):
#   * Patient-years per arm/pool  = Sum(n_arm * follow_up_yr)   [trial_design.csv]
#   * Rate (per 1,000 PY)         = 1000 * events / PY
#   * Signed rate difference      = LAAC - OAC  (all endpoints "lower is better":
#                                   negative = benefit, positive = harm)
#   * Procedural complications    = total_complications_pct * implants_attempted,
#                                   expressed over device PY (control arm = 0)
#   * BRR  = weighted benefit / weighted harm
#   * INB  = weighted benefit - weighted harm   (per 1,000 PY)
#
# Two follow-up (PY) assumptions are reported side by side:
#   (1) "real"  : actual mean follow-up per trial (PRIMARY, consistent w/ text)
#   (2) "cap3y" : follow-up capped at 3.0 years for every trial (sensitivity)
#
# Inputs : data/endpoint_events.csv, data/trial_design.csv,
#          data/procedural_complications_corrected.csv
# Outputs: rba/rba_rate_differences.csv      (PRIMARY = real follow-up)
#          rba/rba_brr_results.csv           (PRIMARY = real follow-up)
#          rba/rba_brr_results_BOTH_PY.csv   (real vs cap3y comparison)
###############################################################################

ev <- read.csv("data/endpoint_events.csv",                  stringsAsFactors = FALSE)
td <- read.csv("data/trial_design.csv",                     stringsAsFactors = FALSE)
pc <- read.csv("data/procedural_complications_corrected.csv", stringsAsFactors = FALSE)

# --- Harmonized IS+SE overlay -------------------------------------------------
# endpoint_events.csv is stale for IS+SE; the FINAL harmonized counts (the ones
# that feed main Table 3, RR 1.29/1.30) live in data/endpoint_events.xlsx.
# Overlay them here so the BRR's IS+SE harm matches the rest of the manuscript
# (IS+SE excess +2.2/1,000 PY; Contemporary harm 8.5; equal BRR 3.19).
# (csv NOT modified, per scope.) Counts: device / control.
isse_h <- data.frame(
  study      = c("CHAMPION-AF","OPTION","PRAGUE-17","CLOSURE-AF","PROTECT-AF/PREVAIL"),
  ev_device  = c(50, 11, 15, 27, 45),
  ev_control = c(35, 11, 13, 25, 14), stringsAsFactors = FALSE)
for (i in seq_len(nrow(isse_h))) {
  r <- ev$endpoint == "IS_SE" & ev$study == isse_h$study[i]
  ev$ev_device[r]  <- isse_h$ev_device[i]
  ev$ev_control[r] <- isse_h$ev_control[i]
}

# ---- Pool definitions (event-data study labels) -----------------------------
pool_studies <- list(
  "Contemporary" = c("CHAMPION-AF", "OPTION", "PRAGUE-17"),
  "All Trials"   = c("CHAMPION-AF", "OPTION", "PRAGUE-17",
                     "CLOSURE-AF", "PROTECT-AF/PREVAIL")
)
study_to_design <- list(
  "CHAMPION-AF"        = "CHAMPION-AF",
  "OPTION"             = "OPTION",
  "PRAGUE-17"          = "PRAGUE-17",
  "CLOSURE-AF"         = "CLOSURE-AF",
  "PROTECT-AF/PREVAIL" = c("PROTECT-AF", "PREVAIL")
)

# ---- Endpoints + per-scenario severity weights (vs non-proc bleeding = 1) ----
endpoints <- c("NonProcBleed", "HemStroke", "AllCauseDeath", "IS_SE", "Procedural")
endpoint_label <- c(
  NonProcBleed  = "Non-procedural bleeding", HemStroke = "Haemorrhagic stroke",
  AllCauseDeath = "All-cause death", IS_SE = "Ischaemic stroke + SE",
  Procedural    = "Procedural complications")
weights <- list(
  Equal             = c(NonProcBleed=1, HemStroke=1, AllCauseDeath=1,  IS_SE=1, Procedural=1),
  Moderate          = c(NonProcBleed=1, HemStroke=3, AllCauseDeath=5,  IS_SE=3, Procedural=1.5),
  "Stroke-priority" = c(NonProcBleed=1, HemStroke=5, AllCauseDeath=10, IS_SE=5, Procedural=2))

# Procedural complication COUNT (device arm) = total_complications_pct * attempts
attempts_n     <- as.integer(sub(".*/", "", pc$implant_success))
pc$proc_events <- round(pc$total_complications_pct / 100 * attempts_n)

# ---- Core computation, parameterized by follow-up cap -----------------------
compute <- function(fu_cap = Inf) {
  fu <- pmin(td$follow_up_yr, fu_cap); names(fu) <- td$study
  py_for_pool <- function(studies) {
    d <- unlist(study_to_design[studies])
    list(device  = sum(td$n_device [match(d, td$study)] * fu[d]),
         control = sum(td$n_control[match(d, td$study)] * fu[d]))
  }
  proc_events_for_pool <- function(studies)
    sum(pc$proc_events[pc$study %in% unlist(study_to_design[studies])])

  rd <- data.frame()
  for (pool in names(pool_studies)) {
    studies <- pool_studies[[pool]]; py <- py_for_pool(studies)
    for (e in endpoints) {
      if (e == "Procedural") { ev_dev <- proc_events_for_pool(studies); ev_ctrl <- 0
      } else {
        sub <- ev[ev$endpoint == e & ev$study %in% studies, ]
        ev_dev <- sum(sub$ev_device); ev_ctrl <- sum(sub$ev_control)
      }
      rdif <- 1000*ev_dev/py$device - 1000*ev_ctrl/py$control
      rd <- rbind(rd, data.frame(pool=pool, endpoint=e, label=endpoint_label[e],
        ev_device=ev_dev, ev_control=ev_ctrl,
        py_device=round(py$device), py_control=round(py$control),
        rate_device=round(1000*ev_dev/py$device,2),
        rate_control=round(1000*ev_ctrl/py$control,2),
        rate_diff=round(rdif,3), side=ifelse(rdif<0,"benefit","harm"),
        stringsAsFactors=FALSE))
    }
  }
  brr <- data.frame()
  for (pool in names(pool_studies)) {
    sub <- rd[rd$pool==pool,]; diffs <- setNames(sub$rate_diff, sub$endpoint)
    for (sc in names(weights)) {
      contrib <- weights[[sc]][names(diffs)] * diffs
      benefit <- -sum(contrib[contrib<0]); harm <- sum(contrib[contrib>0])
      isse <- diffs["IS_SE"]; cni <- contrib[names(contrib)!="IS_SE"]
      be_isse <- (-sum(cni[cni<0]) - sum(cni[cni>0])) / isse
      brr <- rbind(brr, data.frame(pool=pool, scenario=sc,
        weighted_benefit=round(benefit,1), weighted_harm=round(harm,1),
        BRR=round(benefit/harm,2), INB=round(benefit-harm,1),
        breakeven_overall=round(benefit/harm,2), breakeven_ISSE=round(be_isse,1),
        verdict=ifelse(benefit/harm>1,"Favours LAAC","Favours OAC"),
        stringsAsFactors=FALSE))
    }
  }
  list(rd=rd, brr=brr)
}

real  <- compute(Inf)    # PRIMARY
cap3y <- compute(3.0)    # sensitivity

# ---- Write PRIMARY (real follow-up) outputs ---------------------------------
write.csv(real$rd,  "rba/rba_rate_differences.csv", row.names = FALSE)
write.csv(real$brr, "rba/rba_brr_results.csv",      row.names = FALSE)

# ---- Combined comparison ----------------------------------------------------
cmp <- merge(
  transform(real$brr [, c("pool","scenario","BRR","INB","weighted_benefit","weighted_harm")],
            PY="real_followup"),
  transform(cap3y$brr[, c("pool","scenario","BRR","INB","weighted_benefit","weighted_harm")],
            PY="cap_3y"), all = TRUE)
write.csv(rbind(transform(real$brr,  PY="real_followup"),
                transform(cap3y$brr, PY="cap_3y")),
          "rba/rba_brr_results_BOTH_PY.csv", row.names = FALSE)

# ---- Report -----------------------------------------------------------------
show_py <- function(res, cap) { cat(sprintf("\n--- PY assumption: %s ---\n", cap))
  for (pool in names(pool_studies)) {
    p <- res$rd[res$rd$pool==pool,][1,]
    cat(sprintf("  %-13s device=%.0f PY  control=%.0f PY\n", pool, p$py_device, p$py_control)) } }

cat("================= PER-1,000-PY RATE DIFFERENCES (real follow-up) =================\n")
print(real$rd[, c("pool","label","rate_device","rate_control","rate_diff","side")], row.names=FALSE)

cat("\n================= SEVERITY-WEIGHTED BRR — REAL FOLLOW-UP (PRIMARY) =================\n")
print(real$brr, row.names = FALSE)
cat("\n================= SEVERITY-WEIGHTED BRR — FOLLOW-UP CAPPED AT 3y (sensitivity) =====\n")
print(cap3y$brr, row.names = FALSE)
show_py(real, "real follow-up"); show_py(cap3y, "capped at 3.0y")

cat("\n----- Side-by-side BRR -----\n")
sbs <- data.frame(pool=real$brr$pool, scenario=real$brr$scenario,
                  BRR_real=real$brr$BRR, BRR_cap3y=cap3y$brr$BRR,
                  INB_real=real$brr$INB, INB_cap3y=cap3y$brr$INB)
print(sbs, row.names = FALSE)
