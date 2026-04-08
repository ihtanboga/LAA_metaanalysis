###############################################################################
# run_all.R
# LAA Closure vs OAC Meta-Analysis — Master Script
#
# Runs the full analysis pipeline:
#   1. Load packages and data
#   2. Primary analysis (IS+SE, non-proc bleeding, all major bleeding)
#   3. Secondary analysis (ACM, CV death, HS, etc.)
#   4. Sensitivity analyses & benefit-risk
#   5. Generate tables (main + supplementary)
#   6. Generate figures (main + supplementary)
#   7. Render manuscript & supplement to DOCX
###############################################################################

cat("================================================================\n")
cat("  LAA CLOSURE vs OAC — META-ANALYSIS PIPELINE\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

BASE <- normalizePath("~/Desktop/LAA/LAA_meta")

run_script <- function(script_name) {
  path <- file.path(BASE, "R", script_name)
  cat(sprintf("\n>> Running %s ...\n", script_name))
  t0 <- Sys.time()
  tryCatch(
    source(path, local = new.env(parent = globalenv())),
    error = function(e) {
      cat(sprintf("!! ERROR in %s: %s\n", script_name, conditionMessage(e)))
    }
  )
  dt <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf(">> %s completed in %.1f seconds\n", script_name, as.numeric(dt)))
}

# ---- STEP 1: Data ----
cat("=== STEP 1: Loading data ===\n")
run_script("01_data.R")

# ---- STEP 2: Primary Analysis ----
cat("\n=== STEP 2: Primary Analysis ===\n")
run_script("02_primary_analysis.R")

# ---- STEP 3: Secondary Analysis ----
cat("\n=== STEP 3: Secondary Analysis ===\n")
run_script("03_secondary_analysis.R")

# ---- STEP 4: Sensitivity Analysis ----
cat("\n=== STEP 4: Sensitivity Analysis & Benefit-Risk ===\n")
run_script("04_sensitivity.R")

# ---- STEP 5: Tables ----
cat("\n=== STEP 5: Generating Tables ===\n")
run_script("08_tables.R")

# ---- STEP 6: Figures ----
cat("\n=== STEP 6: Generating Main Figures ===\n")
run_script("06_main_figures.R")

cat("\n=== STEP 7: Generating Supplementary Figures ===\n")
run_script("07_suppl_figures.R")

# ---- STEP 8: Render Manuscripts to DOCX ----
cat("\n=== STEP 8: Rendering Manuscripts ===\n")

render_md_to_docx <- function(md_file, docx_file) {
  md_path   <- file.path(BASE, "manuscript", md_file)
  docx_path <- file.path(BASE, "manuscript", docx_file)

  if (!file.exists(md_path)) {
    cat(sprintf("!! %s not found, skipping DOCX render\n", md_file))
    return(invisible(NULL))
  }

  # Try pandoc first, then rmarkdown
  pandoc_cmd <- sprintf('pandoc "%s" -o "%s" --reference-doc="" 2>/dev/null || true',
                        md_path, docx_path)

  tryCatch({
    # Method 1: pandoc directly
    if (Sys.which("pandoc") != "") {
      system2("pandoc", args = c(md_path, "-o", docx_path,
                                  "--wrap=none",
                                  "--toc",
                                  "--toc-depth=2"),
              stdout = TRUE, stderr = TRUE)
      if (file.exists(docx_path)) {
        cat(sprintf(">> %s rendered via pandoc\n", docx_file))
        return(invisible(TRUE))
      }
    }

    # Method 2: rmarkdown::render
    if (requireNamespace("rmarkdown", quietly = TRUE)) {
      tmpRmd <- tempfile(fileext = ".Rmd")
      writeLines(c("---",
                    sprintf("title: ''"),
                    "output: word_document",
                    "---", "",
                    readLines(md_path)),
                 tmpRmd)
      rmarkdown::render(tmpRmd, output_file = docx_path, quiet = TRUE)
      file.remove(tmpRmd)
      cat(sprintf(">> %s rendered via rmarkdown\n", docx_file))
      return(invisible(TRUE))
    }

    cat(sprintf("!! Cannot render %s — neither pandoc nor rmarkdown available\n", md_file))
  }, error = function(e) {
    cat(sprintf("!! DOCX render error for %s: %s\n", md_file, conditionMessage(e)))
  })
}

render_md_to_docx("manuscript.md", "manuscript.docx")
render_md_to_docx("supplement.md", "supplement.docx")

# ---- SUMMARY ----
cat("\n================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("Output files:\n")
cat("  Analysis results:  output/*.RData\n")
cat("  Main tables:       tables/*.png, tables/*.docx\n")
cat("  Suppl tables:      suppl_tables/*.png, suppl_tables/*.docx\n")
cat("  Main figures:      figures/*.pdf\n")
cat("  Suppl figures:     suppl_figures/*.pdf\n")
cat("  Manuscript:        manuscript/manuscript.md + .docx\n")
cat("  Supplement:        manuscript/supplement.md + .docx\n")

# List generated files
cat("\n--- Generated files ---\n")
for (dir_name in c("output", "tables", "suppl_tables", "figures", "suppl_figures", "manuscript")) {
  dir_path <- file.path(BASE, dir_name)
  files <- list.files(dir_path, recursive = FALSE)
  if (length(files) > 0) {
    cat(sprintf("\n  %s/\n", dir_name))
    for (f in files) cat(sprintf("    %s\n", f))
  }
}
