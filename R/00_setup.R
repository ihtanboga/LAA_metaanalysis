###############################################################################
# 00_setup.R
# LAA Closure vs OAC Meta-Analysis — Package Loading & Global Settings
###############################################################################

# ---- CRAN packages ----
pkgs <- c(
 "meta",           # generic meta-analysis (metagen, metabin, forest)
 "metafor",        # advanced meta-analysis (rma)
 "survival",       # Cox, KM
 "survRM2",        # RMST
 "tidyverse",      # dplyr, ggplot2, tidyr, readr, stringr, purrr, forcats
 "patchwork",      # combine ggplot panels
 "scales",         # axis formatting
 "gridExtra",      # table grobs
 "gt",             # publication tables
 "flextable",      # DOCX-compatible tables
 "officer",        # DOCX manipulation
 "rmarkdown",      # markdown rendering
 "knitr",          # table output
 "ggtext",         # rich text in ggplot
 "RColorBrewer",   # colour palettes
 "forestploter",   # publication-quality forest plots (CRAN fallback)
 "checkmate"       # input validation
)

for (p in pkgs) {
 if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
 library(p, character.only = TRUE)
}

# ---- GitHub packages (optional) ----
HAS_FORESTER <- FALSE
tryCatch({
  if (!requireNamespace("forester", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("rdboyes/forester", quiet = TRUE)
  }
  library(forester)
  HAS_FORESTER <- TRUE
  cat(">> forester package loaded\n")
}, error = function(e) {
  cat(">> forester not available, using forestploter fallback\n")
})

# ---- Global paths ----
BASE   <- normalizePath("~/Desktop/LAA/LAA_meta")
DIR_R  <- file.path(BASE, "R")
DIR_D  <- file.path(BASE, "data")
DIR_O  <- file.path(BASE, "output")
DIR_F  <- file.path(BASE, "figures")
DIR_SF <- file.path(BASE, "suppl_figures")
DIR_T  <- file.path(BASE, "tables")
DIR_ST <- file.path(BASE, "suppl_tables")
DIR_M  <- file.path(BASE, "manuscript")

for (d in c(DIR_D, DIR_O, DIR_F, DIR_SF, DIR_T, DIR_ST, DIR_M))
 dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Global options ----
options(digits = 4)
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank()))

# Colour palette for device vs control
COL_DEVICE  <- "#2166AC"
COL_CONTROL <- "#B2182B"
COL_POOL    <- c("Contemporary" = "#1B9E77",
                 "Expanded"     = "#D95F02",
                 "Historical"   = "#7570B3")

cat(">> 00_setup.R loaded successfully\n")
