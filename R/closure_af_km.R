# ============================================================
# CLOSURE-AF: Reconstructed Kaplan-Meier
# LAA Closure vs Physician-Directed Best Medical Care
# Primary Endpoint: Stroke + Systemic Embolism + Major Bleeding + CV/Unexplained Death
# Source: Landmesser et al., NEJM 2026
# ============================================================

library(readr)
library(dplyr)
library(reconstructKM)
library(survival)
library(survminer)
library(ggpubr)

# ===== 1) Digitized CSV dosyalarını oku =====
laa_raw <- read_delim(
  "laa.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
)

control_raw <- read_delim(
  "control.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
)

# Son satır boş olabilir
laa_raw <- laa_raw %>% filter(!is.na(time) & !is.na(cum_inc_pct))
control_raw <- control_raw %>% filter(!is.na(time) & !is.na(cum_inc_pct))

cat("LAA raw noktalar:", nrow(laa_raw), "\n")
cat("Control raw noktalar:", nrow(control_raw), "\n")

# ===== 2) Kümülatif insidans (%) → Survival (0-1) =====
laa_clicks <- laa_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

control_clicks <- control_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

# (0, 1) başlangıç noktasını ekle
if (!any(abs(laa_clicks$time) < 1e-6 & abs(laa_clicks$survival - 1) < 1e-6)) {
  laa_clicks <- bind_rows(data.frame(time = 0, survival = 1), laa_clicks)
}
if (!any(abs(control_clicks$time) < 1e-6 & abs(control_clicks$survival - 1) < 1e-6)) {
  control_clicks <- bind_rows(data.frame(time = 0, survival = 1), control_clicks)
}

# Monoton azalan garanti
collapse_corners <- function(df) {
  df <- df %>%
    group_by(time) %>%
    summarise(survival = min(survival), .groups = "drop") %>%
    arrange(time)
  df$survival <- cummin(pmin(df$survival, 1))
  df$survival[df$survival < 0] <- 0
  df
}

laa_clicks <- collapse_corners(laa_clicks)
control_clicks <- collapse_corners(control_clicks)

cat("LAA clicks (processed):", nrow(laa_clicks), "\n")
cat("Control clicks (processed):", nrow(control_clicks), "\n")

# ===== 3) Number at Risk (NAR) tabloları =====
# Figure 2'den: yıl 0, 1, 2, 3, 4, 5, 6
laa_NAR <- data.frame(
  time = 0:6,
  NAR  = c(446, 304, 202, 117, 71, 33, 9)
)

control_NAR <- data.frame(
  time = 0:6,
  NAR  = c(442, 306, 203, 136, 77, 40, 7)
)

# Clicks'i NAR max zamanına kırp
max_time <- max(laa_NAR$time)
laa_clicks <- laa_clicks %>% filter(time <= max_time)
control_clicks <- control_clicks %>% filter(time <= max_time)

# ===== 4) reconstructKM: format_raw_tabs =====
laa_aug <- format_raw_tabs(raw_NAR = laa_NAR, raw_surv = laa_clicks)
control_aug <- format_raw_tabs(raw_NAR = control_NAR, raw_surv = control_clicks)

# Sütun adı düzeltme
fix_aug_surv <- function(aug) {
  as <- aug$aug_surv
  if (!"surv" %in% names(as) && "survival" %in% names(as))
    names(as)[names(as) == "survival"] <- "surv"
  as <- as[order(as$time), ]
  aug$aug_surv <- as
  aug
}
laa_aug <- fix_aug_surv(laa_aug)
control_aug <- fix_aug_surv(control_aug)

# ===== 5) IPD Rekonstrüksiyonu =====
laa_recon <- KM_reconstruct(aug_NAR = laa_aug$aug_NAR, aug_surv = laa_aug$aug_surv)
control_recon <- KM_reconstruct(aug_NAR = control_aug$aug_NAR, aug_surv = control_aug$aug_surv)

# ===== 6) Arm'ları birleştir =====
laa_IPD <- data.frame(
  arm    = "LAA Closure",
  time   = laa_recon$IPD_time,
  status = laa_recon$IPD_event
)
control_IPD <- data.frame(
  arm    = "Medical Therapy",
  time   = control_recon$IPD_time,
  status = control_recon$IPD_event
)

allIPD <- rbind(laa_IPD, control_IPD)
allIPD$arm <- factor(allIPD$arm, levels = c("Medical Therapy", "LAA Closure"))

# Kontrol
cat("\n=== Reconstructed IPD ===\n")
cat("Toplam IPD satır:", nrow(allIPD), "\n")
cat("  LAA Closure    :", sum(allIPD$arm == "LAA Closure"), "\n")
cat("  Medical Therapy:", sum(allIPD$arm == "Medical Therapy"), "\n")
cat("Events: LAA =", sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
    ", Medical =", sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1), "\n")

# ===== 7) IPD'yi kaydet =====
write.csv(allIPD, "closure_af_allIPD.csv", row.names = FALSE)
cat("IPD kaydedildi: closure_af_allIPD.csv\n")

# ===== 8) Kaplan-Meier Çizimi =====
col_laa     <- "#D62728"   # kırmızı - device
col_control <- "#1F77B4"   # mavi - medical therapy

fit <- survfit(Surv(time, status) ~ arm, data = allIPD)

p <- ggsurvplot(
  fit,
  data         = allIPD,
  fun          = "event",
  conf.int     = FALSE,
  censor       = FALSE,
  xlim         = c(0, 6),
  ylim         = c(0, 1.0),
  break.time.by = 1,
  palette      = c(col_control, col_laa),
  legend.title = "",
  legend.labs  = c("Medical Therapy", "LAA Closure"),
  legend       = c(0.25, 0.85),
  xlab         = "Years since Randomization",
  ylab         = "Cumulative Incidence",
  risk.table       = TRUE,
  risk.table.title = "No. at Risk",
  ggtheme      = theme_classic(base_size = 13)
)

p$plot <- p$plot +
  scale_y_continuous(
    labels = function(x) paste0(x * 100, "%"),
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.10)
  )

# ===== 9) Cox Regresyonu =====
coxfit <- coxph(Surv(time, status) ~ arm, data = allIPD)
sum_cx <- summary(coxfit)
hr  <- exp(coef(coxfit))[1]
ci  <- exp(confint(coxfit))[1, ]
pvl <- sum_cx$coefficients[1, "Pr(>|z|)"]
fmt_p <- ifelse(pvl < 0.001, "<0.001", sprintf("%.3f", pvl))

cat("\n=== Cox Regression ===\n")
cat(sprintf("HR (LAA Closure vs Medical Therapy): %.2f (95%% CI %.2f-%.2f), P = %s\n",
            hr, ci[1], ci[2], fmt_p))

# Grafiğe HR anotasyonu
hr_lab <- sprintf("HR %.2f (95%% CI %.2f\u2013%.2f)\nP = %s", hr, ci[1], ci[2], fmt_p)
p$plot <- p$plot +
  annotate("text", x = 3, y = 0.60, hjust = 0.5, label = hr_lab, size = 4)

print(p)

# ===== 10) PH Varsayımı Testi =====
ph_test <- cox.zph(coxfit)
cat("\n=== Proportional Hazards Test ===\n")
print(ph_test)

# ===== 11) PDF kaydet =====
combined <- ggarrange(
  p$plot, p$table,
  ncol = 1, nrow = 2,
  heights = c(3, 1)
)
ggsave("closure_af_km_plot.pdf", plot = combined, width = 9, height = 7)
cat("\nPDF kaydedildi: closure_af_km_plot.pdf\n")

# ===== 12) RMST Analizi =====
# Makaledeki RMST farkını doğrulamak için
library(survRM2)
rmst_res <- rmst2(
  time   = allIPD$time,
  status = allIPD$status,
  arm    = as.numeric(allIPD$arm) - 1,  # 0 = Medical, 1 = LAA
  tau    = 6  # restrict to 6 years
)
cat("\n=== RMST Analysis (tau = 6 years) ===\n")
print(rmst_res)

# ===== 13) Özet Tablo =====
cat("\n=== CLOSURE-AF Reconstructed Analysis Summary ===\n")
cat("Trial: CLOSURE-AF (Landmesser et al., NEJM 2026)\n")
cat("Primary Endpoint: Stroke + SE + Major Bleeding + CV/Unexplained Death\n")
cat(sprintf("N: LAA Closure = %d, Medical Therapy = %d\n",
            sum(allIPD$arm == "LAA Closure"),
            sum(allIPD$arm == "Medical Therapy")))
cat(sprintf("Events: LAA = %d, Medical = %d\n",
            sum(allIPD$arm == "LAA Closure" & allIPD$status == 1),
            sum(allIPD$arm == "Medical Therapy" & allIPD$status == 1)))
cat(sprintf("HR (LAA vs Medical): %.2f (95%% CI %.2f-%.2f), P = %s\n",
            hr, ci[1], ci[2], fmt_p))
cat(sprintf("PH assumption p-value: %.3f\n", ph_test$table[1, "p"]))
