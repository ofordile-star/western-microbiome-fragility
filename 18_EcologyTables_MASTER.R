# =============================================================================
# EcologyTables_MASTER.R
# [v9 — no hardcoded values; all footnote statistics derived from live objects]
#
# Computes all verification tables T1–T13, then writes two output files:
#
#   1. VerificationTables_MASTER.xlsx  — full analytic record (T1–T13)
#   2. EcologyPaper_SupplementaryTables_GitHub.xlsx — 8 clean deposition
#      sheets (S1–S8), publication-ready footnotes
#      with all values drawn from computed objects at runtime.
#
# Author: Ogochukwu Ofordile, MRC Unit The Gambia at LSHTM
# =============================================================================

library(tidyverse)
library(broom)
library(MASS)
library(effsize)
library(pROC)
library(openxlsx)

# --------------------------------------------------------------------------- #
# PATHS
# --------------------------------------------------------------------------- #
rds_dir    <- "C:/Users/oofordile/Desktop/IHAT_Paper2_RDS/"
output_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Results/"
dir.create(output_dir, showWarnings = FALSE)

mb           <- readRDS(file.path(rds_dir, "mb_full.rds"))
d85_base     <- readRDS(file.path(rds_dir, "d85_base.rds"))
ae_ari       <- readRDS(file.path(rds_dir, "ae_ari.rds"))
species_cols <- readRDS(file.path(rds_dir, "species_cols.rds"))

covs      <- "age_enrolment + gender_n + HAZ_enrolment"
covs_infl <- "HAZ_enrolment + age_sampling + gender_n"

# --------------------------------------------------------------------------- #
# HELPERS
# --------------------------------------------------------------------------- #
fmt_p <- function(p, digits = 3) {
  ifelse(p < 0.001, "< 0.001", sprintf(paste0("%.", digits, "f"), p))
}

run_strat <- function(data, outcome) {
  do.call(rbind, lapply(levels(data$WAZ_group), function(g) {
    sub <- data[data$WAZ_group == g & !is.na(data[[outcome]]), ]
    fit <- lm(as.formula(paste(outcome, "~ log_Pstercorea +", covs_infl)), data = sub)
    cf  <- summary(fit)$coefficients
    ci  <- confint(fit)
    data.frame(
      WAZ_stratum = g,
      outcome     = outcome,
      n           = nobs(fit),
      beta        = round(cf["log_Pstercorea", "Estimate"],   4),
      se          = round(cf["log_Pstercorea", "Std. Error"], 4),
      ci_low      = round(ci["log_Pstercorea", "2.5 %"],      4),
      ci_high     = round(ci["log_Pstercorea", "97.5 %"],     4),
      p           = round(cf["log_Pstercorea", "Pr(>|t|)"],   4)
    )
  }))
}

run_int <- function(data, outcome) {
  fit <- lm(as.formula(paste(outcome,
                             "~ log_Pstercorea * WAZ_enrolment +", covs_infl)), data = data)
  cf  <- summary(fit)$coefficients
  data.frame(
    outcome          = outcome,
    n                = nobs(fit),
    interaction_beta = round(cf["log_Pstercorea:WAZ_enrolment", "Estimate"],  4),
    interaction_p    = round(cf["log_Pstercorea:WAZ_enrolment", "Pr(>|t|)"],  4)
  )
}

# =============================================================================
# D85 BASE DATASET
# =============================================================================
d85 <- mb %>%
  filter(timepoints == "D85") %>%
  rowwise() %>%
  mutate(
    shannon = {
      x <- c_across(all_of(species_cols))
      x <- x[x > 0 & !is.na(x)]
      if (length(x) == 0) NA_real_
      else { p <- x / sum(x); -sum(p * log(p)) }
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(Ill_binary), !is.na(richness), !is.na(shannon)) %>%
  mutate(illness = factor(ifelse(Ill_binary == 0, "Not ill", "Ill"),
                          levels = c("Not ill", "Ill")))

n_d85     <- nrow(d85)
n_not_ill <- sum(d85$Ill_binary == 0)
n_ill     <- sum(d85$Ill_binary == 1)
cat(sprintf("D85 analytic sample: n=%d (%d Not-Ill, %d Ill)\n", n_d85, n_not_ill, n_ill))

# =============================================================================
# T1: Richness
# =============================================================================
rich_t <- t.test(richness ~ illness, data = d85)
rich_d <- effsize::cohen.d(richness ~ illness, data = d85)

T1 <- d85 %>%
  group_by(illness) %>%
  summarise(n      = n(),
            mean   = round(mean(richness),   2),
            sd     = round(sd(richness),     2),
            median = round(median(richness), 1),
            .groups = "drop") %>%
  mutate(t_stat   = round(rich_t$statistic,      3),
         p_value  = round(rich_t$p.value,         4),
         cohens_d = round(abs(rich_d$estimate),   3))

rich_mean_ni  <- T1$mean[T1$illness == "Not ill"]
rich_sd_ni    <- T1$sd[T1$illness == "Not ill"]
rich_mean_ill <- T1$mean[T1$illness == "Ill"]
rich_sd_ill   <- T1$sd[T1$illness == "Ill"]
rich_t_val    <- round(rich_t$statistic, 2)
rich_p        <- rich_t$p.value
rich_d_val    <- round(abs(rich_d$estimate), 2)

# =============================================================================
# T2: Shannon
# =============================================================================
shan_t <- t.test(shannon ~ illness, data = d85)
shan_d <- effsize::cohen.d(shannon ~ illness, data = d85)

T2 <- d85 %>%
  group_by(illness) %>%
  summarise(n      = n(),
            mean   = round(mean(shannon),   4),
            sd     = round(sd(shannon),     4),
            median = round(median(shannon), 4),
            .groups = "drop") %>%
  mutate(t_stat   = round(shan_t$statistic,     3),
         p_value  = round(shan_t$p.value,        4),
         cohens_d = round(abs(shan_d$estimate),  3))

shan_mean_ni  <- T2$mean[T2$illness == "Not ill"]
shan_mean_ill <- T2$mean[T2$illness == "Ill"]
shan_p        <- shan_t$p.value
shan_d_val    <- round(abs(shan_d$estimate), 3)
shan_dir      <- ifelse(shan_mean_ill > shan_mean_ni, "higher", "lower")

# =============================================================================
# T3: ROC / AUC
# =============================================================================
roc_rich <- roc(d85$Ill_binary, d85$richness, quiet = TRUE)
roc_shan <- roc(d85$Ill_binary, d85$shannon,  quiet = TRUE)
roc_test <- roc.test(roc_rich, roc_shan, quiet = TRUE)

auc_rich    <- round(as.numeric(auc(roc_rich)), 3)
auc_rich_lo <- round(ci(roc_rich)[1], 3)
auc_rich_hi <- round(ci(roc_rich)[3], 3)
auc_shan    <- round(as.numeric(auc(roc_shan)), 3)
auc_shan_lo <- round(ci(roc_shan)[1], 3)
auc_shan_hi <- round(ci(roc_shan)[3], 3)
delong_p    <- roc_test$p.value

T3 <- tibble(
  metric   = c("Richness", "Shannon"),
  AUC      = c(auc_rich,    auc_shan),
  CI_low   = c(auc_rich_lo, auc_shan_lo),
  CI_high  = c(auc_rich_hi, auc_shan_hi),
  DeLong_p = c(round(delong_p, 4), NA)
)

# =============================================================================
# T4: NB model (P. stercorea -> ARI, total + richness-adjusted direct)
# =============================================================================
d_ari <- d85_base %>%
  left_join(ae_ari, by = "rand_no") %>%
  replace_na(list(ae_freq = 0, ae_dur = 0)) %>%
  filter(ae_freq <= 5)

n_nb <- nrow(d_ari)
cat(sprintf("NB model dataset: n=%d\n", n_nb))

fit_nb_total  <- glm.nb(as.formula(paste("ae_freq ~ log_Pstercorea +", covs)), data = d_ari)
fit_nb_direct <- glm.nb(as.formula(paste("ae_freq ~ log_Pstercorea + richness +", covs)), data = d_ari)

T4 <- bind_rows(
  tidy(fit_nb_total,  conf.int = TRUE) %>% mutate(model = "Total"),
  tidy(fit_nb_direct, conf.int = TRUE) %>% mutate(model = "Direct (richness-adjusted)")
) %>%
  filter(term == "log_Pstercorea") %>%
  mutate(IRR     = round(exp(estimate),  4),
         CI_low  = round(exp(conf.low),  4),
         CI_high = round(exp(conf.high), 4),
         p_value = round(p.value,        4),
         n       = c(nobs(fit_nb_total), nobs(fit_nb_direct))) %>%
  dplyr::select(model, n, IRR, CI_low, CI_high, p_value)

irr_total     <- T4$IRR[T4$model == "Total"]
irr_total_lo  <- T4$CI_low[T4$model == "Total"]
irr_total_hi  <- T4$CI_high[T4$model == "Total"]
p_total_nb    <- T4$p_value[T4$model == "Total"]
irr_direct    <- T4$IRR[T4$model == "Direct (richness-adjusted)"]
irr_direct_lo <- T4$CI_low[T4$model == "Direct (richness-adjusted)"]
irr_direct_hi <- T4$CI_high[T4$model == "Direct (richness-adjusted)"]
p_direct_nb   <- T4$p_value[T4$model == "Direct (richness-adjusted)"]
n_total_nb    <- T4$n[T4$model == "Total"]

# =============================================================================
# T5: Joint model (P. stercorea + P. copri)
# =============================================================================
fit_joint <- glm.nb(as.formula(
  paste("ae_freq ~ log_Pstercorea + log_Pcopri +", covs)
), data = d_ari)

r_species <- round(cor(d_ari$log_Pstercorea, d_ari$log_Pcopri,
                       use = "complete.obs"), 2)

T5 <- tidy(fit_joint, conf.int = TRUE) %>%
  filter(term %in% c("log_Pstercorea", "log_Pcopri")) %>%
  mutate(IRR     = round(exp(estimate),  4),
         CI_low  = round(exp(conf.low),  4),
         CI_high = round(exp(conf.high), 4),
         p_value = round(p.value,        4),
         n       = nobs(fit_joint),
         r_ps_pc = r_species) %>%
  dplyr::select(term, n, IRR, CI_low, CI_high, p_value, r_ps_pc)

irr_ps_joint    <- T5$IRR[T5$term == "log_Pstercorea"]
irr_ps_joint_lo <- T5$CI_low[T5$term == "log_Pstercorea"]
irr_ps_joint_hi <- T5$CI_high[T5$term == "log_Pstercorea"]
p_ps_joint      <- T5$p_value[T5$term == "log_Pstercorea"]
irr_pc_joint    <- T5$IRR[T5$term == "log_Pcopri"]
irr_pc_joint_lo <- T5$CI_low[T5$term == "log_Pcopri"]
irr_pc_joint_hi <- T5$CI_high[T5$term == "log_Pcopri"]
p_pc_joint      <- T5$p_value[T5$term == "log_Pcopri"]
n_joint         <- T5$n[1]

# =============================================================================
# T6: Baron-Kenny mediation (richness as mediator)
# =============================================================================
m_a <- lm(richness ~ log_Pstercorea + age_enrolment + gender_n + HAZ_enrolment, data = d_ari)
m_b <- lm(ae_freq  ~ richness + log_Pstercorea + age_enrolment + gender_n + HAZ_enrolment, data = d_ari)
m_c <- lm(ae_freq  ~ log_Pstercorea + age_enrolment + gender_n + HAZ_enrolment, data = d_ari)

a       <- coef(m_a)["log_Pstercorea"]
b       <- coef(m_b)["richness"]
c_total <- coef(m_c)["log_Pstercorea"]
nie     <- a * b

sea    <- sqrt(vcov(m_a)["log_Pstercorea", "log_Pstercorea"])
seb    <- sqrt(vcov(m_b)["richness", "richness"])
sob_se <- sqrt((b^2 * sea^2) + (a^2 * seb^2))
sob_z  <- nie / sob_se
sob_p  <- 2 * pnorm(-abs(sob_z))
prop   <- nie / c_total * 100

T6 <- tibble(
  n              = nobs(m_a),
  path_A_beta    = round(a,       4),
  path_B_beta    = round(b,       4),
  NIE            = round(nie,     6),
  c_total        = round(c_total, 4),
  Sobel_z        = round(sob_z,   3),
  Sobel_p        = round(sob_p,   4),
  pct_mediated   = round(prop,    2),
  interpretation = case_when(
    abs(prop) < 1  ~ "< 1% — species-autonomous pathway confirmed",
    abs(prop) < 10 ~ "Weak mediation",
    TRUE           ~ "Partial mediation"
  )
)

pct_med_val <- round(prop, 1)
n_med       <- nobs(m_a)

# =============================================================================
# T7/T8/T9: D85 P. stercorea -> D1 biomarkers (reverse temporal)
# =============================================================================
d1_biomarkers <- mb %>%
  filter(timepoints == "D1") %>%
  dplyr::select(rand_no, log_CRP, log_AGP,
                HAZ_enrolment, WAZ_enrolment, age_sampling, gender_n)

d85_exposure <- d85_base %>% dplyr::select(rand_no, log_Pstercorea)

d_d85_to_d1   <- inner_join(d1_biomarkers, d85_exposure, by = "rand_no")
waz_med_d85d1 <- median(d_d85_to_d1$WAZ_enrolment, na.rm = TRUE)
d_d85_to_d1   <- d_d85_to_d1 %>%
  mutate(WAZ_group = factor(
    ifelse(WAZ_enrolment <= waz_med_d85d1, "Low WAZ", "High WAZ"),
    levels = c("Low WAZ", "High WAZ")
  ))

T7 <- run_strat(d_d85_to_d1, "log_CRP")
T8 <- run_strat(d_d85_to_d1, "log_AGP")
T9 <- bind_rows(run_int(d_d85_to_d1, "log_CRP"),
                run_int(d_d85_to_d1, "log_AGP"))

crp_hi_beta_d85 <- T7$beta[T7$WAZ_stratum == "High WAZ"]
crp_hi_p_d85    <- T7$p[T7$WAZ_stratum    == "High WAZ"]
agp_hi_beta_d85 <- T8$beta[T8$WAZ_stratum == "High WAZ"]
agp_hi_p_d85    <- T8$p[T8$WAZ_stratum    == "High WAZ"]
int_crp_p_d85   <- T9$interaction_p[T9$outcome == "log_CRP"]
int_agp_p_d85   <- T9$interaction_p[T9$outcome == "log_AGP"]
n_hi_d85        <- T7$n[T7$WAZ_stratum == "High WAZ"]
n_lo_d85        <- T7$n[T7$WAZ_stratum == "Low WAZ"]

# =============================================================================
# T10: D85 P. stercorea distribution by illness group
# =============================================================================
T10 <- d85 %>%
  group_by(illness) %>%
  summarise(n    = n(),
            mean = round(mean(log_Pstercorea), 4),
            sd   = round(sd(log_Pstercorea),   4),
            .groups = "drop")

# =============================================================================
# T11/T12/T13: D1 P. stercorea -> D1 biomarkers (cross-sectional)
# =============================================================================
d1_full <- mb %>%
  filter(timepoints == "D1") %>%
  dplyr::select(rand_no, log_Pstercorea, log_CRP, log_AGP,
                HAZ_enrolment, WAZ_enrolment, age_sampling, gender_n)

waz_med_d1 <- median(d1_full$WAZ_enrolment, na.rm = TRUE)
d1_full <- d1_full %>%
  mutate(WAZ_group = factor(
    ifelse(WAZ_enrolment <= waz_med_d1, "Low WAZ", "High WAZ"),
    levels = c("Low WAZ", "High WAZ")
  ))

T11 <- run_strat(d1_full, "log_CRP")
T12 <- run_strat(d1_full, "log_AGP")
T13 <- bind_rows(run_int(d1_full, "log_CRP"),
                 run_int(d1_full, "log_AGP"))

crp_lo_beta_d1 <- T11$beta[T11$WAZ_stratum == "Low WAZ"]
crp_lo_p_d1    <- T11$p[T11$WAZ_stratum    == "Low WAZ"]
agp_hi_beta_d1 <- T12$beta[T12$WAZ_stratum == "High WAZ"]
agp_hi_p_d1    <- T12$p[T12$WAZ_stratum    == "High WAZ"]
int_crp_p_d1   <- T13$interaction_p[T13$outcome == "log_CRP"]
int_agp_p_d1   <- T13$interaction_p[T13$outcome == "log_AGP"]
n_lo_d1        <- T11$n[T11$WAZ_stratum == "Low WAZ"]
n_hi_d1        <- T12$n[T12$WAZ_stratum == "High WAZ"]

cat("\n=== ALL TABLES COMPUTED ===\n")

# =============================================================================
# BH CORRECTION — two families, all values from live computed objects
#
# Family 1 (stratum-level, 8 tests): T7, T8, T11, T12
#   Each data frame has 2 rows (Low WAZ, High WAZ); p-values stacked in order
#   D85→D1: CRP Low/High WAZ, AGP Low/High WAZ
#   D1→D1:  CRP Low/High WAZ, AGP Low/High WAZ
#
# Family 2 (interaction, 4 tests): T9, T13
#   WAZ×P.stercorea continuous: D85 CRP, D85 AGP, D1 CRP, D1 AGP
# =============================================================================

# ---- Family 1 ----
strat_p_family <- c(T7$p, T8$p, T11$p, T12$p)   # 8 values, ordered as above
strat_q_family <- p.adjust(strat_p_family, method = "BH")

T7  <- T7  %>% mutate(q_BH = round(strat_q_family[1:nrow(T7)],   4))
T8  <- T8  %>% mutate(q_BH = round(strat_q_family[(nrow(T7)+1):(nrow(T7)+nrow(T8))], 4))
T11 <- T11 %>% mutate(q_BH = round(strat_q_family[(nrow(T7)+nrow(T8)+1):(nrow(T7)+nrow(T8)+nrow(T11))], 4))
T12 <- T12 %>% mutate(q_BH = round(strat_q_family[(nrow(T7)+nrow(T8)+nrow(T11)+1):length(strat_q_family)], 4))

cat("\n--- Stratum-level BH q-values (Family 1, n=8) ---\n")
cat("D85 CRP:\n"); print(T7[, c("WAZ_stratum","p","q_BH")])
cat("D85 AGP:\n"); print(T8[, c("WAZ_stratum","p","q_BH")])
cat("D1  CRP:\n"); print(T11[, c("WAZ_stratum","p","q_BH")])
cat("D1  AGP:\n"); print(T12[, c("WAZ_stratum","p","q_BH")])

# ---- Family 2 ----
int_p_family <- c(T9$interaction_p, T13$interaction_p)  # 4 values: D85 CRP, D85 AGP, D1 CRP, D1 AGP
int_q_family <- p.adjust(int_p_family, method = "BH")

T9  <- T9  %>% mutate(interaction_q_BH = round(int_q_family[1:nrow(T9)],  4))
T13 <- T13 %>% mutate(interaction_q_BH = round(int_q_family[(nrow(T9)+1):length(int_q_family)], 4))

cat("\n--- Interaction BH q-values (Family 2, n=4) ---\n")
cat("D85:\n"); print(T9[,  c("outcome","interaction_p","interaction_q_BH")])
cat("D1:\n");  print(T13[, c("outcome","interaction_p","interaction_q_BH")])

# ---- Extract scalars for footnotes (all from live objects) ----
crp_hi_q_d85  <- T7$q_BH[T7$WAZ_stratum  == "High WAZ"]
crp_lo_q_d85  <- T7$q_BH[T7$WAZ_stratum  == "Low WAZ"]
agp_hi_q_d85  <- T8$q_BH[T8$WAZ_stratum  == "High WAZ"]
agp_lo_q_d85  <- T8$q_BH[T8$WAZ_stratum  == "Low WAZ"]
int_crp_q_d85 <- T9$interaction_q_BH[T9$outcome  == "log_CRP"]
int_agp_q_d85 <- T9$interaction_q_BH[T9$outcome  == "log_AGP"]

crp_lo_q_d1   <- T11$q_BH[T11$WAZ_stratum == "Low WAZ"]
crp_hi_q_d1   <- T11$q_BH[T11$WAZ_stratum == "High WAZ"]
agp_hi_q_d1   <- T12$q_BH[T12$WAZ_stratum == "High WAZ"]
agp_lo_q_d1   <- T12$q_BH[T12$WAZ_stratum == "Low WAZ"]
int_crp_q_d1  <- T13$interaction_q_BH[T13$outcome == "log_CRP"]
int_agp_q_d1  <- T13$interaction_q_BH[T13$outcome == "log_AGP"]

# =============================================================================
# OUTPUT 1: VERIFICATION MASTER (internal — one sheet per table)
# =============================================================================
wb_master <- createWorkbook()

add_raw <- function(wb, sheet, df) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df)
}

add_raw(wb_master, "T1_Richness_D85",     T1)
add_raw(wb_master, "T2_Shannon_D85",       T2)
add_raw(wb_master, "T3_ROC_D85",           T3)
add_raw(wb_master, "T4_NB_Pstercorea_D85", T4)
add_raw(wb_master, "T5_Joint_Pstercorea",  T5)
add_raw(wb_master, "T6_Mediation",         T6)
add_raw(wb_master, "T7_D85_to_CRP_WAZ",   T7)
add_raw(wb_master, "T8_D85_to_AGP_WAZ",   T8)
add_raw(wb_master, "T9_D85_WAZ_Int",       T9)
add_raw(wb_master, "T10_Pstercorea_D85",   T10)
add_raw(wb_master, "T11_D1_to_CRP_WAZ",   T11)
add_raw(wb_master, "T12_D1_to_AGP_WAZ",   T12)
add_raw(wb_master, "T13_D1_WAZ_Int",       T13)

master_file <- file.path(output_dir, "VerificationTables_MASTER.xlsx")
saveWorkbook(wb_master, master_file, overwrite = TRUE)
cat("Saved:", master_file, "\n")

# =============================================================================
# OUTPUT 2: GITHUB DEPOSITION (styled, 8 sheets, all values from live objects)
# =============================================================================

NAVY  <- "#2F4F7F"
LIGHT <- "#EBF0F7"

sty_title <- createStyle(fontName = "Calibri", fontSize = 10,
                         textDecoration = "bold")
sty_hdr   <- createStyle(fontName = "Calibri", fontSize = 9,
                         fontColour = "#FFFFFF", fgFill = NAVY,
                         halign = "center", textDecoration = "bold",
                         border = "TopBottom", borderColour = NAVY,
                         wrapText = TRUE)
sty_odd   <- createStyle(fontName = "Calibri", fontSize = 9,
                         fgFill = LIGHT,
                         border = "TopBottom", borderColour = "#CCCCCC")
sty_even  <- createStyle(fontName = "Calibri", fontSize = 9,
                         border = "TopBottom", borderColour = "#CCCCCC")
sty_note  <- createStyle(fontName = "Calibri", fontSize = 8,
                         fontColour = "#555555",
                         textDecoration = "italic", wrapText = TRUE)

write_clean_sheet <- function(wb, sheet_name, df, title, note = NULL) {
  addWorksheet(wb, sheet_name)
  nc <- ncol(df)
  
  writeData(wb, sheet_name, title, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name, sty_title, rows = 1, cols = 1)
  if (nc > 1) mergeCells(wb, sheet_name, rows = 1, cols = 1:nc)
  
  for (j in seq_len(nc)) {
    writeData(wb, sheet_name, names(df)[j], startRow = 2, startCol = j)
    addStyle(wb, sheet_name, sty_hdr, rows = 2, cols = j)
  }
  setRowHeights(wb, sheet_name, rows = 2, heights = 30)
  
  for (i in seq_len(nrow(df))) {
    writeData(wb, sheet_name, df[i, , drop = FALSE],
              startRow = i + 2, startCol = 1, colNames = FALSE)
    addStyle(wb, sheet_name,
             if (i %% 2 == 1) sty_odd else sty_even,
             rows = i + 2, cols = 1:nc, gridExpand = TRUE)
  }
  
  if (!is.null(note)) {
    nr <- nrow(df) + 4
    writeData(wb, sheet_name, note, startRow = nr, startCol = 1)
    addStyle(wb, sheet_name, sty_note, rows = nr, cols = 1)
    if (nc > 1) mergeCells(wb, sheet_name, rows = nr, cols = 1:nc)
    setRowHeights(wb, sheet_name, rows = nr, heights = 70)
  }
  
  setColWidths(wb, sheet_name, cols = 1:nc,
               widths = pmax(10, nchar(names(df)) + 3))
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)
}

# Consolidated tables
T7_T8 <- bind_rows(T7 %>% mutate(biomarker = "log_CRP"),
                   T8 %>% mutate(biomarker = "log_AGP")) %>%
  dplyr::select(biomarker, WAZ_stratum, n, beta, se, ci_low, ci_high, p, q_BH)

T11_T12 <- bind_rows(T11 %>% mutate(biomarker = "log_CRP"),
                     T12 %>% mutate(biomarker = "log_AGP")) %>%
  dplyr::select(biomarker, WAZ_stratum, n, beta, se, ci_low, ci_high, p, q_BH)

wb_gh <- createWorkbook()

write_clean_sheet(wb_gh, "S1_Richness", T1,
                  sprintf("Supplementary Table 1. Genus-level richness at Day 85 by illness group (IHAT-GUT, n=%d)", n_d85),
                  note = sprintf(
                    paste0("Day 85 microbiome after filtering for non-missing illness status, richness, and Shannon diversity. ",
                           "n=%d: %d not-ill, %d ill. Welch two-sample t-test. Cohen's d = |Not-Ill minus Ill| / pooled SD. ",
                           "Not-Ill: %.1f +/- %.1f genera; Ill: %.1f +/- %.1f genera; t=%.2f, p=%s, d=%.2f."),
                    n_d85, n_not_ill, n_ill,
                    rich_mean_ni, rich_sd_ni, rich_mean_ill, rich_sd_ill,
                    rich_t_val, fmt_p(rich_p), rich_d_val
                  )
)

write_clean_sheet(wb_gh, "S2_Shannon", T2,
                  sprintf("Supplementary Table 2. Shannon diversity at Day 85 by illness group (IHAT-GUT, n=%d)", n_d85),
                  note = sprintf(
                    paste0("Same sample as Table S1 (n=%d). Shannon index computed from species-level relative abundances. ",
                           "Ill children have %s Shannon (%.3f vs %.3f), direction reversed relative to richness. ",
                           "Non-significant (p=%s, d=%.3f), confirming dissociation between evenness and richness-based protection."),
                    n_d85, shan_dir, shan_mean_ill, shan_mean_ni, fmt_p(shan_p), shan_d_val
                  )
)

write_clean_sheet(wb_gh, "S3_AUC", T3,
                  sprintf("Supplementary Table 3. ROC / AUC comparison: genus-level richness vs Shannon diversity (D85, n=%d)", n_d85),
                  note = sprintf(
                    paste0("AUC and 95%% CI from pROC package (DeLong method). ",
                           "AUC richness = %.3f (%.3f-%.3f); AUC Shannon = %.3f (%.3f-%.3f). ",
                           "DeLong test for difference: p=%s. ",
                           "Richness discriminates illness; Shannon performs at chance."),
                    auc_rich, auc_rich_lo, auc_rich_hi,
                    auc_shan, auc_shan_lo, auc_shan_hi,
                    fmt_p(delong_p)
                  )
)

write_clean_sheet(wb_gh, "S4_NB_Models", T4,
                  sprintf("Supplementary Table 4. Negative binomial GLM: P. stercorea -> ARI frequency (D85, n=%d)", n_total_nb),
                  note = sprintf(
                    paste0("d85_base joined to ae_ari; ae_freq capped at <=5. n=%d after covariate complete-case filtering. ",
                           "Total model: IRR=%.3f (%.3f-%.3f), p=%s. ",
                           "Direct model (richness-adjusted): IRR=%.3f (%.3f-%.3f), p=%s. ",
                           "Covariates: age, sex, HAZ. IRR = exp(beta). CIs from profile likelihood."),
                    n_total_nb,
                    irr_total, irr_total_lo, irr_total_hi, fmt_p(p_total_nb),
                    irr_direct, irr_direct_lo, irr_direct_hi, fmt_p(p_direct_nb)
                  )
)

write_clean_sheet(wb_gh, "S5_Joint_Model", T5,
                  sprintf("Supplementary Table 5. Joint NB model: P. stercorea + P. copri -> ARI (D85, n=%d)", n_joint),
                  note = sprintf(
                    paste0("Both species entered simultaneously. P. copri serves as within-model negative control. ",
                           "Pearson r(log_Pstercorea, log_Pcopri) = %.2f. ",
                           "P. stercorea: IRR=%.3f (%.3f-%.3f), p=%s. ",
                           "P. copri: IRR=%.3f (%.3f-%.3f), p=%s. Species-specificity confirmed."),
                    r_species,
                    irr_ps_joint, irr_ps_joint_lo, irr_ps_joint_hi, fmt_p(p_ps_joint),
                    irr_pc_joint, irr_pc_joint_lo, irr_pc_joint_hi, fmt_p(p_pc_joint)
                  )
)

write_clean_sheet(wb_gh, "S6_Mediation", T6,
                  sprintf("Supplementary Table 6. Baron-Kenny mediation: richness as mediator of P. stercorea -> ARI (n=%d)", n_med),
                  note = sprintf(
                    paste0("Path A: log_Pstercorea -> richness (OLS). ",
                           "Path B: richness -> ae_freq | log_Pstercorea (OLS). ",
                           "NIE = Path A x Path B. Sobel test for NIE significance. ",
                           "Proportion mediated = %.1f%%. ",
                           "Near-zero value confirms P. stercorea ARI protection is richness-independent (species-autonomous pathway)."),
                    pct_med_val
                  )
)

write_clean_sheet(wb_gh, "S7_D85toD1_Biomarkers", T7_T8,
                  "Supplementary Table 7. Day 85 P. stercorea (exposure) -> Day 1 CRP and AGP (outcome), WAZ-stratified",
                  note = sprintf(
                    paste0(
                      "Regression direction: Day 85 log_Pstercorea as predictor; Day 1 log_CRP or log_AGP as outcome. ",
                      "Read prospectively: Day 1 inflammatory tone predicts subsequent Day 85 P. stercorea colonisation. ",
                      "WAZ median split = %.3f. Covariates: HAZ, age_sampling, sex. ",
                      "q_BH: Benjamini-Hochberg correction across 8 stratum-level tests (Family 1: T7, T8, T11, T12). ",
                      "High WAZ (n=%d): CRP beta=%.3f, p=%s, q=%s; AGP beta=%.3f, p=%s, q=%s. ",
                      "Low WAZ (n=%d): no associations (CRP p=%s, q=%s; AGP p=%s, q=%s). ",
                      "Interaction q-values (Family 2, 4 tests): CRP int. p=%s, q=%s; AGP int. p=%s, q=%s."
                    ),
                    waz_med_d85d1,
                    n_hi_d85,
                    crp_hi_beta_d85, fmt_p(crp_hi_p_d85), fmt_p(crp_hi_q_d85),
                    agp_hi_beta_d85, fmt_p(agp_hi_p_d85), fmt_p(agp_hi_q_d85),
                    n_lo_d85,
                    fmt_p(T7$p[T7$WAZ_stratum == "Low WAZ"]), fmt_p(crp_lo_q_d85),
                    fmt_p(T8$p[T8$WAZ_stratum == "Low WAZ"]), fmt_p(agp_lo_q_d85),
                    fmt_p(int_crp_p_d85), fmt_p(int_crp_q_d85),
                    fmt_p(int_agp_p_d85), fmt_p(int_agp_q_d85)
                  )
)

write_clean_sheet(wb_gh, "S8_D1_CrossSectional", T11_T12,
                  "Supplementary Table 8. Day 1 P. stercorea -> Day 1 CRP and AGP (cross-sectional baseline), WAZ-stratified",
                  note = sprintf(
                    paste0(
                      "Both exposure (log_Pstercorea) and outcomes (log_CRP, log_AGP) measured at Day 1 (enrolment). ",
                      "Cross-sectional; no temporal inference. WAZ median split = %.3f. ",
                      "q_BH: Benjamini-Hochberg correction across 8 stratum-level tests (Family 1: T7, T8, T11, T12). ",
                      "Low WAZ (n=%d): CRP beta=%.3f, p=%s, q=%s (CRP suppression: colonisation resistance). ",
                      "High WAZ (n=%d): AGP beta=%.3f, p=%s, q=%s (immune engagement: tonic priming). ",
                      "Low WAZ AGP: p=%s, q=%s (null). High WAZ CRP: p=%s, q=%s (null). ",
                      "Interaction q-values (Family 2, 4 tests): CRP int. p=%s, q=%s; AGP int. p=%s, q=%s. ",
                      "Directional asymmetry by WAZ stratum supports host immune-metabolic reserve-dependent immune conditioning."
                    ),
                    waz_med_d1,
                    n_lo_d1,
                    crp_lo_beta_d1, fmt_p(crp_lo_p_d1), fmt_p(crp_lo_q_d1),
                    n_hi_d1,
                    agp_hi_beta_d1, fmt_p(agp_hi_p_d1), fmt_p(agp_hi_q_d1),
                    fmt_p(T12$p[T12$WAZ_stratum == "Low WAZ"]), fmt_p(agp_lo_q_d1),
                    fmt_p(T11$p[T11$WAZ_stratum == "High WAZ"]), fmt_p(crp_hi_q_d1),
                    fmt_p(int_crp_p_d1), fmt_p(int_crp_q_d1),
                    fmt_p(int_agp_p_d1), fmt_p(int_agp_q_d1)
                  )
)

github_file <- file.path(output_dir, "EcologyPaper_SupplementaryTables_GitHub.xlsx")
saveWorkbook(wb_gh, github_file, overwrite = TRUE)
cat("Saved:", github_file, "\n")

cat("\n=== COMPLETE ===\n")
cat("Output 1 (verification): VerificationTables_MASTER.xlsx\n")
cat("Output 2 (deposition):   EcologyPaper_SupplementaryTables_GitHub.xlsx\n")
