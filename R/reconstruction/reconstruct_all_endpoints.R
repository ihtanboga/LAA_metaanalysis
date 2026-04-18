# ============================================================
# Reconstruct IPD for each trial × endpoint, combining AC + device
# ============================================================
# For every JSON pair in json_csv/ (AC + device), reconstruct IPD via
# Guyot algorithm and save a single combined CSV per endpoint under
# json_csv/ipd_data/<TRIAL>_<endpoint_slug>.csv
# ============================================================

library(readr)
library(dplyr)
library(jsonlite)
library(reconstructKM)

set.seed(2026)

json_dir <- "~/Desktop/manuel/json_csv"
out_dir  <- "~/Desktop/manuel/json_csv/ipd_data"

# Endpoint slugs (trial, endpoint_slug, AC_basename, device_basename)
endpoints <- list(
  # ---- CHAMPION ----
  list("CHAMPION", "major_bleeding",
       "champion_ITT_majorbleeding_AC",
       "champion_ITT_majorbleeding_device"),
  list("CHAMPION", "nonprocedural_ISTH_major_bleeding",
       "CHAMPION_ITT_NonproceduralISTHMajorBleeding_AC",
       "CHAMPION_ITT_NonproceduralISTHMajorBleeding_device"),
  list("CHAMPION", "nonprocedure_bleeding",
       "champion_ITT_nonprocedurebleeding_AC",
       "champion_ITT_nonprocedurebleeding_device"),
  list("CHAMPION", "nonprocedural_modified_ISTH_CRNM_bleeding",
       "champion_ITT_NonproceduralModifiedISTHClinicallyRelevantNonmajorBleeding_AC",
       "champion_ITT_NonproceduralModifiedISTHClinicallyRelevantNonmajorBleeding_device"),
  list("CHAMPION", "procedural_nonprocedural_ISTH_major_and_CRNM_bleeding",
       "champion_ITT_ProceduralandNonproceduralISTHMajorandModifiedISTHClinicallyRelevantNonmajorBleeding_AC",
       "champion_ITT_ProceduralandNonproceduralISTHMajorandModifiedISTHClinicallyRelevantNonmajorBleeding_device"),
  list("CHAMPION", "all_stroke",
       "champion_ITT_AllStroke_AC",
       "champion_ITT_AllStroke_device"),

  # ---- OPTION ----
  list("OPTION", "major_bleeding",
       "OPTION_ITT_majorbleeding_AC",
       "OPTION_ITT_majorbleeding_device"),
  list("OPTION", "ISTH_major_bleeding_including_procedure_related",
       "OPTION_ITT_ISTHmajorbleedingincludingprocedurerelatedbleeding_AC",
       "OPTION_ITT_ISTHmajorbleedingincludingprocedurerelatedbleeding_device"),
  list("OPTION", "nonprocedure_major_or_CRNM_bleeding",
       "OPTION_ITT_NonProcedureRelatedMajorBleedingorClinicallyRelevantNonmajorBleeding_AC",
       "OPTION_ITT_NonProcedureRelatedMajorBleedingorClinicallyRelevantNonmajorBleeding_device"),
  list("OPTION", "death_stroke_SE_composite",
       "OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism_AC",
       "OPTION_ITT_DeathfromAnyCauseStrokeorSystemicEmbolism_device"),

  # ---- PRAGUE-17 ----
  list("PRAGUE17", "all_stroke_TIA",
       "PRAGUE17_ITT_allstroketransientischemicattack_AC",
       "PRAGUE17_ITT_allstroketransientischemicattack_device"),
  list("PRAGUE17", "clinically_relevant_bleeding",
       "PRAGUE17_ITT_clinicallyrelevantbleeding_AC",
       "PRAGUE17_ITT_clinicallyrelevantbleeding_device"),
  list("PRAGUE17", "nonprocedural_clinically_relevant_bleeding",
       "PRAGUE17_ITT_nonproceduralclinicallyrelevantbleeding_AC",
       "PRAGUE17_ITT_nonproceduralclinicallyrelevantbleeding_device")
)

# ===== Helpers =====
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

reconstruct_arm <- function(base_name, arm_label) {
  csv_path  <- file.path(json_dir, paste0(base_name, ".csv"))
  json_path <- file.path(json_dir, paste0(base_name, ".json"))

  df   <- read_csv(csv_path, comment = "#", show_col_types = FALSE)
  meta <- fromJSON(json_path)

  nar_raw <- meta$arms$nar[[1]]
  nar <- data.frame(time = nar_raw$time, NAR = nar_raw$n)

  clicks <- trim_tail(prep_clicks(df))
  # Trim clicks to NAR window
  clicks <- clicks %>% filter(time <= max(nar$time))

  aug   <- fix_aug(format_raw_tabs(raw_NAR = nar, raw_surv = clicks))
  recon <- KM_reconstruct(aug_NAR = aug$aug_NAR, aug_surv = aug$aug_surv)

  data.frame(arm = arm_label,
             time = recon$IPD_time,
             status = recon$IPD_event)
}

# ===== Main loop =====
summary_rows <- list()

for (ep in endpoints) {
  trial       <- ep[[1]]
  slug        <- ep[[2]]
  ac_base     <- ep[[3]]
  dev_base    <- ep[[4]]
  out_file    <- file.path(out_dir, sprintf("%s_%s.csv", trial, slug))

  cat(sprintf("\n=== %s | %s ===\n", trial, slug))

  ac  <- tryCatch(reconstruct_arm(ac_base,  "AC"),     error = function(e) e)
  dev <- tryCatch(reconstruct_arm(dev_base, "device"), error = function(e) e)

  if (inherits(ac, "error") || inherits(dev, "error")) {
    cat("  ERROR:\n")
    if (inherits(ac, "error"))  cat("    AC: ", conditionMessage(ac), "\n")
    if (inherits(dev, "error")) cat("    device: ", conditionMessage(dev), "\n")
    next
  }

  combined <- rbind(ac, dev)
  write.csv(combined, out_file, row.names = FALSE)

  n_ac     <- sum(combined$arm == "AC")
  n_dev    <- sum(combined$arm == "device")
  ev_ac    <- sum(combined$arm == "AC"     & combined$status == 1)
  ev_dev   <- sum(combined$arm == "device" & combined$status == 1)

  cat(sprintf("  AC:     N=%d, events=%d\n", n_ac, ev_ac))
  cat(sprintf("  device: N=%d, events=%d\n", n_dev, ev_dev))
  cat(sprintf("  Saved:  %s\n", basename(out_file)))

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trial = trial, endpoint = slug,
    N_AC = n_ac, events_AC = ev_ac,
    N_device = n_dev, events_device = ev_dev,
    file = basename(out_file)
  )
}

cat("\n\n===== SUMMARY =====\n")
summary_df <- do.call(rbind, summary_rows)
print(summary_df, row.names = FALSE)
write.csv(summary_df,
          file.path(out_dir, "_reconstruction_summary.csv"),
          row.names = FALSE)
cat("\nSummary saved: _reconstruction_summary.csv\n")
