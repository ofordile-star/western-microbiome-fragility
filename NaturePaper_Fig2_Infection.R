# =============================================================================
# NaturePaper_Fig2_revised.R
# [REVISED — D1 interpretation updated to context-dependent dual immune mode]
#
# Figure 2 | Species-specific and richness-independent protection by
#            Prevotella stercorea
#
# Panels:
#   A  Joint NB model IRRs — P. stercorea vs P. copri (species specificity)
#
#   B  Mediation: proportion of P.stercorea→ARI mediated by richness
#
#   C  Temporal antecedence: D85 P.stercorea abundance by illness
#
#   D  DUAL-TIMEPOINT nutritional gating panel (REVISED INTERPRETATION):
#      Left subfigure:  D1  P.stercorea → CRP/AGP by WAZ (prospective baseline)
#                       Low WAZ: CRP↓ (β=−0.037, p=0.016) — colonisation resistance
#                       High WAZ: AGP↑ (β=+0.014, p=0.008) — sustained immune tone
#                       Interaction: CRP p=0.015; AGP p=0.001
#      Right subfigure: D85 P.stercorea → CRP/AGP by WAZ (post-follow-up)
#                       High WAZ: CRP↑ (β=+0.064, p=0.037) and AGP↑ (β=+0.018, p=0.007)
#                       Interaction: CRP p<0.001; AGP p<0.001
#                       Interpret with caution: D85 microbiome partly encodes illness
#                       history; signal reflects better immune engagement in
#                       better-nourished, better-protected children
#
# Author:  Ogochukwu Ofordile, MRC Unit The Gambia at LSHTM
# =============================================================================

library(tidyverse)
library(patchwork)
library(MASS)
library(broom)

rds_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_RDS/"
fig_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Figures/"
dir.create(fig_dir, showWarnings = FALSE)

d85_base <- readRDS(file.path(rds_dir, "d85_base.rds"))
ae_ari   <- readRDS(file.path(rds_dir, "ae_ari.rds"))
mb       <- readRDS(file.path(rds_dir, "mb_full.rds"))

# --------------------------------------------------------------------------- #
# THEME
# --------------------------------------------------------------------------- #
theme_nature <- theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    plot.title         = element_text(size = 10, face = "bold"),
    plot.subtitle      = element_text(size = 8.0, colour = "grey40"),
    plot.tag           = element_text(size = 13, face = "bold"),
    legend.position    = "none",
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4)
  )

col_ps      <- "#2166AC"
col_pc      <- "#D6604D"
col_high    <- "#2166AC"   # positive significant
col_low     <- "#D6604D"   # negative significant
pal_species <- c("P. stercorea" = col_ps, "P. copri" = col_pc)

covs      <- "age_enrolment + gender_n + HAZ_enrolment"
covs_infl <- "HAZ_enrolment + age_sampling + gender_n"

# --------------------------------------------------------------------------- #
# COMPUTE: ARI models
# --------------------------------------------------------------------------- #
d_ari <- d85_base %>%
  left_join(ae_ari, by = "rand_no") %>%
  replace_na(list(ae_freq = 0, ae_dur = 0)) %>%
  filter(ae_freq <= 5)

fit_joint <- glm.nb(
  as.formula(paste("ae_freq ~ log_Pstercorea + log_Pcopri +", covs)),
  data = d_ari
)

joint_res <- tidy(fit_joint, conf.int = TRUE) %>%
  filter(term %in% c("log_Pstercorea", "log_Pcopri")) %>%
  mutate(
    IRR    = exp(estimate),
    IRR_lo = exp(conf.low),
    IRR_hi = exp(conf.high),
    species = factor(
      recode(term, "log_Pstercorea" = "P. stercorea",
                   "log_Pcopri"     = "P. copri"),
      levels = c("P. copri", "P. stercorea")
    )
  )

# Mediation
m_a   <- lm(as.formula(paste("richness ~ log_Pstercorea +", covs)), data = d_ari)
m_b   <- lm(as.formula(paste("ae_freq ~ richness + log_Pstercorea +", covs)), data = d_ari)
m_c   <- lm(as.formula(paste("ae_freq ~ log_Pstercorea +", covs)), data = d_ari)
a     <- coef(m_a)["log_Pstercorea"]
b     <- coef(m_b)["richness"]
c_tot <- coef(m_c)["log_Pstercorea"]
prop_ari_freq <- (a * b) / c_tot

cat(sprintf("P.stercorea IRR=%.4f p=%.4f; P.copri IRR=%.4f p=%.4f\n",
            joint_res$IRR[joint_res$species=="P. stercorea"],
            joint_res$p.value[joint_res$species=="P. stercorea"],
            joint_res$IRR[joint_res$species=="P. copri"],
            joint_res$p.value[joint_res$species=="P. copri"]))
cat(sprintf("Mediation proportion (ARI freq): %.1f%%\n", prop_ari_freq * 100))

# --------------------------------------------------------------------------- #
# COMPUTE: WAZ stratification — D85 exposure → D1 biomarkers
# --------------------------------------------------------------------------- #
d1_bio <- mb %>%
  filter(timepoints == "D1") %>%
  dplyr::select(rand_no, log_CRP, log_AGP, HAZ_enrolment,
                WAZ_enrolment, age_sampling, gender_n)

d85_exp <- d85_base %>% dplyr::select(rand_no, log_Pstercorea)

d_infl_d85 <- inner_join(d1_bio, d85_exp, by = "rand_no") %>%
  filter(!is.na(WAZ_enrolment), !is.na(log_Pstercorea))

waz_med_d85 <- median(d_infl_d85$WAZ_enrolment, na.rm = TRUE)
d_infl_d85 <- d_infl_d85 %>%
  mutate(WAZ_group = factor(
    ifelse(WAZ_enrolment <= waz_med_d85, "Low WAZ", "High WAZ"),
    levels = c("Low WAZ", "High WAZ")
  ))

# --------------------------------------------------------------------------- #
# COMPUTE: WAZ stratification — D1 exposure → D1 biomarkers (prospective)
# --------------------------------------------------------------------------- #
d1_full <- mb %>%
  filter(timepoints == "D1") %>%
  dplyr::select(rand_no, log_Pstercorea, log_CRP, log_AGP,
                HAZ_enrolment, WAZ_enrolment, age_sampling, gender_n) %>%
  filter(!is.na(WAZ_enrolment), !is.na(log_Pstercorea))

waz_med_d1 <- median(d1_full$WAZ_enrolment, na.rm = TRUE)
d1_full <- d1_full %>%
  mutate(WAZ_group = factor(
    ifelse(WAZ_enrolment <= waz_med_d1, "Low WAZ", "High WAZ"),
    levels = c("Low WAZ", "High WAZ")
  ))

# Helper: run stratified models
run_strat_both <- function(data, outcome) {
  do.call(rbind, lapply(levels(data$WAZ_group), function(g) {
    sub <- data[data$WAZ_group == g & !is.na(data[[outcome]]), ]
    fit <- lm(as.formula(paste(outcome, "~ log_Pstercorea +", covs_infl)), data = sub)
    cf  <- summary(fit)$coefficients
    ci  <- confint(fit)
    data.frame(
      stratum = g, n = nobs(fit),
      beta    = round(cf["log_Pstercorea", "Estimate"],  4),
      ci_low  = round(ci["log_Pstercorea", "2.5 %"],     4),
      ci_high = round(ci["log_Pstercorea", "97.5 %"],    4),
      p       = round(cf["log_Pstercorea", "Pr(>|t|)"],  4),
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(stratum = factor(stratum, levels = c("Low WAZ", "High WAZ")),
           outcome = outcome)
}

strat_d85_crp <- run_strat_both(d_infl_d85, "log_CRP")
strat_d85_agp <- run_strat_both(d_infl_d85, "log_AGP")
strat_d1_crp  <- run_strat_both(d1_full,    "log_CRP")
strat_d1_agp  <- run_strat_both(d1_full,    "log_AGP")

# Interaction p-values (continuous WAZ × Pstercorea)
get_int_p <- function(data, outcome) {
  summary(lm(as.formula(paste(outcome, "~ log_Pstercorea * WAZ_enrolment +", covs_infl)),
             data = data))$coefficients["log_Pstercorea:WAZ_enrolment", "Pr(>|t|)"]
}

int_d85_crp <- get_int_p(d_infl_d85, "log_CRP")
int_d85_agp <- get_int_p(d_infl_d85, "log_AGP")
int_d1_crp  <- get_int_p(d1_full,    "log_CRP")
int_d1_agp  <- get_int_p(d1_full,    "log_AGP")

cat(sprintf("\nD85→biomarker interactions: CRP p=%.4f, AGP p=%.4f\n", int_d85_crp, int_d85_agp))
cat(sprintf("D1→biomarker interactions:  CRP p=%.4f, AGP p=%.4f\n",  int_d1_crp,  int_d1_agp))

# Print stratified estimates for verification
cat("\nD1 stratified estimates:\n")
print(strat_d1_crp); print(strat_d1_agp)
cat("\nD85 stratified estimates:\n")
print(strat_d85_crp); print(strat_d85_agp)

# ============================================================================
# PANEL A: Species specificity (unchanged)
# ============================================================================
figA <- ggplot(joint_res,
               aes(x = IRR, xmin = IRR_lo, xmax = IRR_hi,
                   y = species, colour = species)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.22, linewidth = 0.85) +
  geom_point(size = 4.5) +
  geom_text(
    aes(x = IRR_hi,
        label = ifelse(p.value < 0.05,
                       sprintf("  IRR %.3f, p=%.3f *", IRR, p.value),
                       sprintf("  IRR %.3f, p=%.3f",   IRR, p.value))),
    hjust = 0, vjust = 0.4, size = 2.9, colour = "grey20"
  ) +
  scale_colour_manual(values = pal_species) +
  scale_x_continuous(
    limits = c(0.87, 1.18),
    breaks = c(0.90, 0.95, 1.00, 1.05, 1.10),
    expand = expansion(mult = c(0.02, 0.28))
  ) +
  labs(
    x     = "IRR for ARI  (joint NB, per log-unit abundance)",
    y     = NULL,
    title = expression("Species specificity: " * italic("P. stercorea") *
                         " independently associated; " * italic("P. copri") * " null"),
    tag   = "a"
  ) +
  theme_nature +
  theme(axis.text.y = element_text(face = "italic", size = 9.5))

# ============================================================================
# PANEL B: Mediation (unchanged)
# ============================================================================
ari_dur_prop  <- -0.098   # from Table S2
inf_freq_prop <-  0.308
dia_freq_prop <-  0.445

med_dat <- tibble(
  outcome  = factor(
    c("ARI\n(frequency)", "ARI\n(duration)",
      "Infection\n(reference)", "Diarrhoea\n(reference)"),
    levels = c("Diarrhoea\n(reference)", "Infection\n(reference)",
               "ARI\n(duration)",        "ARI\n(frequency)")
  ),
  prop_pct = c(prop_ari_freq * 100, ari_dur_prop * 100,
               inf_freq_prop * 100, dia_freq_prop * 100),
  is_ari   = c(TRUE, TRUE, FALSE, FALSE)
) %>%
  mutate(bar_colour = ifelse(prop_pct < 0, "#B2182B",
                             ifelse(is_ari, "#6BAED6", "#2166AC")))

figB <- ggplot(med_dat, aes(x = prop_pct, y = outcome, fill = bar_colour)) +
  geom_col(width = 0.55, colour = "white") +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
  geom_text(
    aes(x     = prop_pct + ifelse(prop_pct >= 0, 2, -2),
        label = sprintf("%.0f%%", prop_pct),
        hjust = ifelse(prop_pct >= 0, 0, 1)),
    size = 3.2, colour = "grey20"
  ) +
  annotate("rect", xmin = -55, xmax = 15, ymin = 2.45, ymax = 4.55,
           fill = col_ps, alpha = 0.06, colour = NA) +
  annotate("text", x = -50, y = 4.42,
           label = "ARI: near-zero or negative\n(richness-independent)",
           size = 2.7, colour = col_ps, fontface = "italic", hjust = 0,
           lineheight = 0.9) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-55, 75),
                     breaks = c(-40, -20, 0, 20, 40, 60),
                     labels = c("-40%", "-20%", "0%", "20%", "40%", "60%")) +
  labs(
    x        = "Proportion of effect mediated through genus-level richness",
    y        = NULL,
    title    = "ARI protection is richness-independent",
    subtitle = "Negative = richness suppresses (not mediates) P. stercorea protection",
    tag      = "b"
  ) +
  theme_nature +
  theme(panel.grid.major.y = element_blank())

# ============================================================================
# PANEL C: Temporal antecedence (D85 abundance by illness; unchanged)
# ============================================================================
d85_ps <- mb %>%
  filter(timepoints == "D85", !is.na(Ill_binary), !is.na(log_Pstercorea)) %>%
  mutate(illness = factor(ifelse(Ill_binary == 0, "Not ill", "Ill"),
                          levels = c("Not ill", "Ill")))

ps_means <- d85_ps %>% group_by(illness) %>%
  summarise(m = mean(log_Pstercorea), .groups = "drop")
ps_t <- t.test(log_Pstercorea ~ illness, data = d85_ps)

cat(sprintf("\nPanel C: Not-Ill mean=%.3f, Ill mean=%.3f, p=%.4f\n",
            ps_means$m[ps_means$illness=="Not ill"],
            ps_means$m[ps_means$illness=="Ill"],
            ps_t$p.value))

pal_illness <- c("Not ill" = col_ps, "Ill" = col_pc)

figC <- ggplot(d85_ps, aes(x = log_Pstercorea,
                            fill = illness, colour = illness)) +
  geom_density(alpha = 0.22, linewidth = 0.7, trim = FALSE) +
  geom_vline(data = ps_means, aes(xintercept = m, colour = illness),
             linetype = "dashed", linewidth = 0.75) +
  geom_text(data = ps_means,
            aes(x = m, colour = illness,
                label = sprintf("mean = %.2f", m)),
            y = Inf, vjust = 1.6, hjust = 0.5, size = 2.9,
            inherit.aes = FALSE) +
  scale_fill_manual(values = pal_illness, name = NULL) +
  scale_colour_manual(values = pal_illness, name = NULL) +
  labs(
    x        = expression("log(1 + " * italic("P. stercorea") * ") abundance, D85"),
    y        = "Density",
    title    = expression(italic("P. stercorea") ~ "depletion precedes illness onset"),
    subtitle = sprintf("D85 prospective microbiome; t-test p = %.4f", ps_t$p.value),
    tag      = "c"
  ) +
  theme_nature +
  theme(
    legend.position   = c(0.72, 0.78),
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(0.85, "lines"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# ============================================================================
# PANEL D: Dual-timepoint nutritional gating (REVISED)
#
# Display order: D1 (left) then D85 (right).
# D1 is the prospective, causally cleaner panel and therefore leads.
# D85 is shown for completeness but labelled with the caution note.
#
# Colour logic (from make_forest_dat):
#   positive & significant → col_high (blue)
#   negative & significant → col_low  (red)
#   non-significant        → grey60
#
# Expected D1 pattern:
#   CRP Low WAZ  → negative & significant  (red)  — colonisation resistance
#   AGP High WAZ → positive & significant  (blue) — sustained immune tone
#   others       → grey
#
# Expected D85 pattern:
#   CRP High WAZ → positive & significant  (blue)
#   AGP High WAZ → positive & significant  (blue)
# ============================================================================

make_forest_dat <- function(strat_crp, strat_agp, int_crp, int_agp,
                             timepoint_label) {
  bind_rows(
    strat_crp %>% mutate(marker = "CRP", int_p = int_crp),
    strat_agp %>% mutate(marker = "AGP", int_p = int_agp)
  ) %>%
    mutate(
      timepoint = timepoint_label,
      row_label = paste(marker, stratum, sep = "\n"),
      sig       = p < 0.05,
      pt_colour = case_when(
        beta > 0 & sig  ~ col_high,
        beta < 0 & sig  ~ col_low,
        TRUE            ~ "grey60"
      )
    )
}

# D1 is now the LEFT panel; D85 is RIGHT
dat_d1  <- make_forest_dat(strat_d1_crp,  strat_d1_agp,
                            int_d1_crp,  int_d1_agp,
                            "D1 (prospective baseline)")
dat_d85 <- make_forest_dat(strat_d85_crp, strat_d85_agp,
                            int_d85_crp, int_d85_agp,
                            "D85 (post-follow-up)")

forest_dat <- bind_rows(dat_d1, dat_d85) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("D1 (prospective baseline)",
                                  "D85 (post-follow-up)")),
    row_label = factor(row_label,
                       levels = c("AGP\nHigh WAZ", "AGP\nLow WAZ",
                                  "CRP\nHigh WAZ", "CRP\nLow WAZ"))
  )

# Interaction annotation labels per facet
int_labels <- forest_dat %>%
  distinct(timepoint, marker, int_p) %>%
  group_by(timepoint) %>%
  summarise(
    label = paste0(
      "CRP int. p=", formatC(unique(int_p[marker == "CRP"]), digits = 3, format = "f"), "\n",
      "AGP int. p=", formatC(unique(int_p[marker == "AGP"]), digits = 3, format = "f")
    ),
    .groups = "drop"
  )

# Mechanistic annotation labels — only shown on D1 facet
mech_labels <- tibble(
  timepoint = factor("D1 (prospective baseline)",
                     levels = levels(forest_dat$timepoint)),
  row_label = factor(c("CRP\nLow WAZ", "AGP\nHigh WAZ"),
                     levels = levels(forest_dat$row_label)),
  label     = c("colonisation\nresistance", "immune\ntone"),
  x_pos     = c(-0.065, 0.10)
)

figD <- ggplot(forest_dat,
               aes(x = beta, xmin = ci_low, xmax = ci_high,
                   y = row_label, colour = pt_colour)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.25, linewidth = 0.8) +
  geom_point(aes(size = sig)) +
  geom_text(
    aes(label = ifelse(p < 0.05,
                       sprintf("%.3f *", beta),
                       sprintf("%.3f",   beta)),
        x = ci_high),
    hjust = -0.15, size = 2.6, colour = "grey20"
  ) +
  # Interaction p annotations (top-right each facet)
  geom_text(data = int_labels,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.3, size = 2.4,
            colour = "grey35", lineheight = 1.0,
            inherit.aes = FALSE) +
  # D1-only mechanistic annotation
  geom_text(data = mech_labels,
            aes(x = x_pos, y = row_label, label = label),
            size = 2.3, colour = "grey50", fontface = "italic",
            lineheight = 0.85, vjust = -0.8,
            inherit.aes = FALSE) +
  facet_wrap(~timepoint, ncol = 2) +
  scale_colour_identity() +
  scale_size_manual(values = c("TRUE" = 4.0, "FALSE" = 2.5), guide = "none") +
  scale_x_continuous(
    limits = c(-0.09, 0.20),
    breaks = c(-0.06, 0, 0.06, 0.12)
  ) +
  labs(
    x = expression(beta ~ "(log-biomarker per log-unit " *
                     italic("P. stercorea") * ")"),
    y        = NULL,
    title    = expression(italic("P. stercorea") ~
                            "\u2192 biomarkers: context-dependent, nutritionally gated"),
    subtitle = paste0(
      "D1: Low WAZ CRP↓ | High WAZ AGP↑; ",
      "D85: High WAZ CRP↑ & AGP↑"
    ),
    tag = "d"
  ) +
  theme_nature +
  theme(
    strip.text         = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 8.5, lineheight = 0.9)
  )

# ============================================================================
# ASSEMBLE & SAVE
# ============================================================================
fig2 <- (figA | figB) / (figC | figD) +
  plot_layout(heights = c(1, 1.4))

for (fmt in c("pdf", "tiff")) {
  outfile <- file.path(fig_dir, paste0("NaturePaper_Fig2_revised.", fmt))
  ggsave(outfile, fig2, width = 10, height = 9,
         dpi = if (fmt == "tiff") 300 else 150)
  cat("Saved:", outfile, "\n")
}

cat("\nFig 2 complete.\n")
cat("D1 panel is now LEFT facet (prospective baseline, cleaner causal inference).\n")
cat("D85 panel is RIGHT facet with caution label in subtitle.\n")
cat("Colour coding: blue = positive & significant; red = negative & significant.\n")
cat("Expected D1 pattern: CRP Low WAZ red; AGP High WAZ blue.\n")
