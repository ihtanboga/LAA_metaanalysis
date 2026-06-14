###############################################################################
# Figure 5 — Severity-weighted benefit-risk ratio (BRR) and incremental net
# benefit (INB) under three pre-specified weighting scenarios.
#
# Two pools: Contemporary DOAC-era (primary) and All Trials.
# Source numbers: rba/rba_brr_results.csv  (run rba/rba_analysis.R first)
#
# Outputs:
#   main/fig5_benefit_risk.pdf                          (editable, cairo)
#   final_heart_rhythm/Figures_300dpi/Figure_5_BenefitRisk.png  (300 dpi)
#   ... .tif / .jpg are produced from the PNG by run_figure5.sh (ImageMagick)
###############################################################################

br <- read.csv("rba/rba_brr_results.csv", stringsAsFactors = FALSE)

scen_levels <- c("Equal", "Moderate", "Stroke-priority")
scen_short  <- c(Equal = "Equal", Moderate = "Moderate", "Stroke-priority" = "Stroke-\npriority")
pool_levels <- c("Contemporary", "All Trials")
pool_col    <- c("Contemporary" = "#1F6FB2", "All Trials" = "#E08214")

# matrices: rows = pool, cols = scenario
get_mat <- function(field) {
  m <- matrix(NA, nrow = length(pool_levels), ncol = length(scen_levels),
              dimnames = list(pool_levels, scen_levels))
  for (i in seq_along(pool_levels))
    for (j in seq_along(scen_levels))
      m[i, j] <- br[[field]][br$pool == pool_levels[i] & br$scenario == scen_levels[j]]
  m
}
brr_m <- get_mat("BRR")
inb_m <- get_mat("INB")

draw_figure <- function() {
  par(mfrow = c(1, 2), mar = c(5.5, 4.8, 4, 1.2), mgp = c(3, 0.7, 0),
      cex.axis = 0.95, cex.lab = 1.05)

  ## ---- Panel A: BRR ----
  ymaxA <- max(brr_m) * 1.18
  bpA <- barplot(brr_m, beside = TRUE, col = pool_col[pool_levels], border = "white",
                 ylim = c(0, ymaxA), names.arg = scen_short[scen_levels],
                 ylab = "Benefit-risk ratio (BRR)", axes = TRUE)
  abline(h = 1, lty = 2, lwd = 1.6, col = "#B22222")
  text(bpA[1, 1], 1, "BRR = 1 (equipoise)", pos = 3, offset = 0.25, adj = 0,
       col = "#B22222", cex = 0.78, font = 3, xpd = NA)
  text(bpA, brr_m, sprintf("%.2f", brr_m), pos = 3, offset = 0.3,
       cex = 0.9, font = 2, col = pool_col[pool_levels])
  mtext("A", side = 3, line = 1.6, at = par("usr")[1], font = 2, cex = 1.4, xpd = NA)
  title(main = "Severity-weighted benefit-risk ratio", cex.main = 1.05, line = 0.6)

  ## ---- Panel B: INB ----
  ymaxB <- max(inb_m) * 1.20; yminB <- min(0, min(inb_m) * 1.2)
  bpB <- barplot(inb_m, beside = TRUE, col = pool_col[pool_levels], border = "white",
                 ylim = c(yminB, ymaxB), names.arg = scen_short[scen_levels],
                 ylab = "Incremental net benefit (events / 1,000 PY)", axes = TRUE)
  abline(h = 0, lwd = 1.2)
  text(bpB, inb_m, sprintf("%+.1f", inb_m),
       pos = ifelse(inb_m >= 0, 3, 1), offset = 0.3,
       cex = 0.9, font = 2, col = pool_col[pool_levels])
  mtext("B", side = 3, line = 1.6, at = par("usr")[1], font = 2, cex = 1.4, xpd = NA)
  title(main = "Incremental net benefit", cex.main = 1.05, line = 0.6)

  ## ---- shared legend + footnote ----
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend(x = 0.5, y = 0.095, xjust = 0.5, horiz = TRUE, bty = "n",
         legend = pool_levels, fill = pool_col[pool_levels], border = "white",
         cex = 0.95, x.intersp = 0.6)
  text(0.5, 0.022, paste0(
       "Weights vs non-procedural bleeding (=1): Moderate  IS+SE 3x, death 5x, ",
       "haemorrhagic stroke 3x, procedural 1.5x;  Stroke-priority  IS+SE 5x, ",
       "death 10x, haemorrhagic stroke 5x, procedural 2x."),
       cex = 0.62, col = "grey30", xpd = NA)
}

# ---- render: editable PDF + 300-dpi PNG (cairo) -----------------------------
W <- 11; H <- 6.0   # inches
cairo_pdf("main/fig5_benefit_risk.pdf", width = W, height = H)
draw_figure(); dev.off()

png("final_heart_rhythm/Figures_300dpi/Figure_5_BenefitRisk.png",
    width = W, height = H, units = "in", res = 300, type = "cairo")
draw_figure(); dev.off()

cat("Figure 5 written:\n  main/fig5_benefit_risk.pdf\n",
    " final_heart_rhythm/Figures_300dpi/Figure_5_BenefitRisk.png\n", sep = "")
cat("\nBRR matrix:\n"); print(round(brr_m, 2))
cat("\nINB matrix:\n"); print(round(inb_m, 1))
