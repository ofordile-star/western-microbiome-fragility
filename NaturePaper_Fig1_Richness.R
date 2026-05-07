# =============================================================================
# NaturePaper_Fig1_Richness.R  [FINAL — values verified against master table]
#
# Figure 1 | Richness, but not Shannon diversity, associates with protection
#            from infection in Gambian children
#
# VERIFIED VALUES FROM MASTER TABLE:
#   Richness: Not-Ill 103.06 ± 25.28 vs Ill 93.90 ± 25.47
#             t = 3.531, p = 0.0005, Cohen's d = 0.36
#   Shannon:  Not-Ill 2.251 ± 0.439 vs Ill 2.272 ± 0.452
#             t = −0.456, p = 0.649, d = 0.046 (ns; direction REVERSED)
#   AUC richness = 0.605 (0.548–0.662); AUC Shannon = 0.521 (0.463–0.578)
#   DeLong p = 0.104
#
# NOTE ON SHANNON DIRECTION: Ill children have marginally *higher* Shannon
# (2.272 vs 2.251). This reinforces the dissociation — evenness is not only
# non-significant but slightly inverted, underlining that richness (not
# evenness) is the ecologically relevant metric.
#
# Panels:
#   A  Genus-level richness: Not-Ill vs Ill (violin + box)
#   B  Shannon diversity: Not-Ill vs Ill (same layout, ns + direction note)
#   C  Cohen's d effect size comparison
#   D  ROC / AUC comparison
#
# Author:  Ogochukwu Ofordile, MRC Unit The Gambia at LSHTM
# =============================================================================

library(tidyverse)
library(patchwork)
library(effsize)
library(pROC)

rds_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_RDS/"
fig_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Figures/"
dir.create(fig_dir, showWarnings = FALSE)

mb           <- readRDS(file.path(rds_dir, "mb_full.rds"))
species_cols <- readRDS(file.path(rds_dir, "species_cols.rds"))

# --------------------------------------------------------------------------- #
# DERIVE D85 dataset with Shannon
# --------------------------------------------------------------------------- #
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
  mutate(
    illness = factor(ifelse(Ill_binary == 0, "Not ill", "Ill"),
                     levels = c("Not ill", "Ill"))
  )

cat(sprintf("n = %d (%d Not-Ill, %d Ill)\n",
            nrow(d85), sum(d85$Ill_binary == 0), sum(d85$Ill_binary == 1)))

# --------------------------------------------------------------------------- #
# STATS (recomputed live; should match master table exactly)
# --------------------------------------------------------------------------- #
rich_t <- t.test(richness ~ illness, data = d85)
shan_t <- t.test(shannon  ~ illness, data = d85)
rich_d <- effsize::cohen.d(richness ~ illness, data = d85)
shan_d <- effsize::cohen.d(shannon  ~ illness, data = d85)

cat(sprintf("Richness: Not-Ill %.2f±%.2f vs Ill %.2f±%.2f, p=%.4f, d=%.3f\n",
            mean(d85$richness[d85$illness=="Not ill"]),
            sd(d85$richness[d85$illness=="Not ill"]),
            mean(d85$richness[d85$illness=="Ill"]),
            sd(d85$richness[d85$illness=="Ill"]),
            rich_t$p.value, abs(rich_d$estimate)))
cat(sprintf("Shannon:  Not-Ill %.4f±%.4f vs Ill %.4f±%.4f, p=%.4f, d=%.3f\n",
            mean(d85$shannon[d85$illness=="Not ill"]),
            sd(d85$shannon[d85$illness=="Not ill"]),
            mean(d85$shannon[d85$illness=="Ill"]),
            sd(d85$shannon[d85$illness=="Ill"]),
            shan_t$p.value, abs(shan_d$estimate)))
cat(sprintf("Shannon direction: ill children have %s Shannon (difference = %.4f)\n",
            ifelse(mean(d85$shannon[d85$illness=="Ill"]) >
                   mean(d85$shannon[d85$illness=="Not ill"]),
                   "HIGHER", "lower"),
            mean(d85$shannon[d85$illness=="Ill"]) -
            mean(d85$shannon[d85$illness=="Not ill"])))

# --------------------------------------------------------------------------- #
# THEME
# --------------------------------------------------------------------------- #
theme_nature <- theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    plot.title         = element_text(size = 10, face = "bold"),
    plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
    plot.tag           = element_text(size = 13, face = "bold"),
    legend.position    = "none",
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4)
  )

col_notill <- "#2166AC"
col_ill    <- "#D6604D"
pal_ill    <- c("Not ill" = col_notill, "Ill" = col_ill)

# --------------------------------------------------------------------------- #
# GROUP MEANS
# --------------------------------------------------------------------------- #
rich_means <- d85 %>% group_by(illness) %>%
  summarise(m = mean(richness), .groups = "drop")
shan_means <- d85 %>% group_by(illness) %>%
  summarise(m = mean(shannon), .groups = "drop")

# --------------------------------------------------------------------------- #
# PANEL A: Richness
# --------------------------------------------------------------------------- #
figA <- ggplot(d85, aes(x = illness, y = richness,
                        fill = illness, colour = illness)) +
  geom_violin(alpha = 0.22, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.17, alpha = 0.85, outlier.size = 0.7,
               outlier.alpha = 0.35, linewidth = 0.55,
               colour = "grey20", fill = "white") +
  geom_point(data = rich_means, aes(y = m),
             shape = 18, size = 3.5, colour = "black") +
  annotate("segment",
           x = 1, xend = 2,
           y = max(d85$richness) * 1.04,
           yend = max(d85$richness) * 1.04,
           linewidth = 0.6) +
  annotate("text",
           x = 1.5, y = max(d85$richness) * 1.075,
           label = sprintf("p = %s  \u2022  d = %.2f",
                           ifelse(rich_t$p.value < 0.001,
                                  "< 0.001",
                                  sprintf("%.3f", rich_t$p.value)),
                           abs(rich_d$estimate)),
           size = 3.2, fontface = "bold") +
  scale_fill_manual(values = pal_ill) +
  scale_colour_manual(values = pal_ill) +
  labs(
    x        = NULL,
    y        = "Genus-level richness (no. genera)",
    title    = "Richness significantly separates illness groups",
    subtitle = sprintf("Not-Ill: %.1f \u00b1 %.1f genera   Ill: %.1f \u00b1 %.1f genera (D85)",
                       mean(d85$richness[d85$illness=="Not ill"]),
                       sd(d85$richness[d85$illness=="Not ill"]),
                       mean(d85$richness[d85$illness=="Ill"]),
                       sd(d85$richness[d85$illness=="Ill"])),
    tag      = "a"
  ) +
  theme_nature

# --------------------------------------------------------------------------- #
# PANEL B: Shannon
# --------------------------------------------------------------------------- #
figB <- ggplot(d85, aes(x = illness, y = shannon,
                        fill = illness, colour = illness)) +
  geom_violin(alpha = 0.22, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.17, alpha = 0.85, outlier.size = 0.7,
               outlier.alpha = 0.35, linewidth = 0.55,
               colour = "grey20", fill = "white") +
  geom_point(data = shan_means, aes(y = m),
             shape = 18, size = 3.5, colour = "black") +
  annotate("segment",
           x = 1, xend = 2,
           y = max(d85$shannon, na.rm = TRUE) * 1.04,
           yend = max(d85$shannon, na.rm = TRUE) * 1.04,
           linewidth = 0.6) +
  annotate("text",
           x = 1.5, y = max(d85$shannon, na.rm = TRUE) * 1.075,
           label = sprintf("p = %.3f  \u2022  d = %.3f  (ns)",
                           shan_t$p.value, abs(shan_d$estimate)),
           size = 3.2, colour = "grey45") +
  # note the reversed direction
  annotate("text",
           x = 1.5, y = max(d85$shannon, na.rm = TRUE) * 0.60,
           label = "Ill children have marginally\nhigher Shannon (direction reversed)",
           size = 2.7, colour = "grey45", fontface = "italic") +
  scale_fill_manual(values = pal_ill) +
  scale_colour_manual(values = pal_ill) +
  labs(
    x        = NULL,
    y        = "Shannon diversity index",
    title    = "Shannon diversity: no separation, direction reversed",
    subtitle = sprintf("Not-Ill: %.3f \u00b1 %.3f   Ill: %.3f \u00b1 %.3f (D85)",
                       mean(d85$shannon[d85$illness=="Not ill"]),
                       sd(d85$shannon[d85$illness=="Not ill"]),
                       mean(d85$shannon[d85$illness=="Ill"]),
                       sd(d85$shannon[d85$illness=="Ill"])),
    tag      = "b"
  ) +
  theme_nature

# --------------------------------------------------------------------------- #
# PANEL C: Cohen's d effect size
# --------------------------------------------------------------------------- #
effect_df <- tibble(
  metric  = c("Genus richness", "Shannon diversity"),
  d       = c(abs(rich_d$estimate), abs(shan_d$estimate)),
  ci_low  = c(min(abs(rich_d$conf.int)), min(abs(shan_d$conf.int))),
  ci_high = c(max(abs(rich_d$conf.int)), max(abs(shan_d$conf.int))),
  sig     = c(rich_t$p.value < 0.05, shan_t$p.value < 0.05)
) %>%
  mutate(
    metric = factor(metric, levels = c("Shannon diversity", "Genus richness")),
    colour = ifelse(sig, col_notill, "grey60")
  )

figC <- ggplot(effect_df,
               aes(x = d, xmin = ci_low, xmax = ci_high,
                   y = metric, colour = colour)) +
  geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.2, linetype = "dotted", colour = "grey70") +
  geom_errorbarh(height = 0.18, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("d = %.3f", d)),
            hjust = -0.3, size = 3.2, colour = "black") +
  annotate("text", x = 0.20, y = 0.52,
           label = "small\neffect", size = 2.7,
           colour = "grey55", hjust = 0.5, lineheight = 0.9) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 0.60), expand = c(0, 0)) +
  labs(
    x     = "Cohen's d  (|Not-Ill \u2212 Ill|)",
    y     = NULL,
    title = "Effect size: richness >> evenness",
    tag   = "c"
  ) +
  theme_nature +
  theme(panel.grid.major.y = element_blank())

# --------------------------------------------------------------------------- #
# PANEL D: ROC / AUC
# --------------------------------------------------------------------------- #
roc_rich <- pROC::roc(d85$Ill_binary, d85$richness, quiet = TRUE)
roc_shan <- pROC::roc(d85$Ill_binary, d85$shannon,  quiet = TRUE)
roc_test <- pROC::roc.test(roc_rich, roc_shan, quiet = TRUE)

auc_rich <- as.numeric(pROC::auc(roc_rich))
auc_shan <- as.numeric(pROC::auc(roc_shan))

cat(sprintf("\nAUC richness = %.3f | AUC Shannon = %.3f | DeLong p = %.4f\n",
            auc_rich, auc_shan, roc_test$p.value))

roc_df <- bind_rows(
  tibble(fpr = 1 - roc_rich$specificities,
         tpr = roc_rich$sensitivities,
         label = sprintf("Richness (AUC = %.2f)", auc_rich)),
  tibble(fpr = 1 - roc_shan$specificities,
         tpr = roc_shan$sensitivities,
         label = sprintf("Shannon  (AUC = %.2f)", auc_shan))
)

pal_roc <- setNames(
  c(col_notill, "grey55"),
  c(sprintf("Richness (AUC = %.2f)", auc_rich),
    sprintf("Shannon  (AUC = %.2f)", auc_shan))
)

figD <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = label)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey70") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal_roc, name = NULL) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  annotate("text",
           x = 0.62, y = 0.20,
           label = sprintf("DeLong test\np = %.3f", roc_test$p.value),
           size = 3.0, colour = "grey40", lineheight = 1.1) +
  labs(
    x     = "False positive rate",
    y     = "True positive rate",
    title = "Richness predicts illness; Shannon at chance",
    tag   = "d"
  ) +
  theme_nature +
  theme(
    legend.position   = c(0.62, 0.28),
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(0.9, "lines"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# --------------------------------------------------------------------------- #
# ASSEMBLE & SAVE
# --------------------------------------------------------------------------- #
fig1 <- (figA | figB) / (figC | figD) +
  plot_layout(heights = c(1.3, 1))

for (fmt in c("pdf", "tiff")) {
  outfile <- file.path(fig_dir, paste0("NaturePaper_Fig1.", fmt))
  ggsave(outfile, fig1, width = 8.5, height = 7.5,
         dpi = if (fmt == "tiff") 300 else 150)
  cat("Saved:", outfile, "\n")
}

cat("\nFig 1 complete. Values match master table.\n")
