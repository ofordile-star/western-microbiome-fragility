# =============================================================================
# NaturePaper_Fig4_ResilienceLandscape.R  [REVISED v3 for ISME Communications]
#
# Figure 4 | Conceptual, illustrative model of microbiome resilience and
#             functional redundancy. No panel is fitted to data (see Methods).
#
# v3 CHANGE (per editorial review of the resubmission):
# - Panel B rebuilt entirely. It previously plotted real IHAT-GUT/Western
#   richness values against an axis literally labelled "functional
#   redundancy" -- even flagged as illustrative, this invited the exact
#   R3 objection it was trying to defuse ("you haven't measured redundancy,
#   so how can you position populations on a redundancy axis?"). The new
#   panel keeps real richness values on the x-axis (observed) but plots NO
#   data points against redundancy; the y-axis carries only a wide,
#   explicitly labelled uncertainty band around a hypothesised, untested
#   trend line. No population is assigned a redundancy value anywhere.
#
# v2 CHANGES (post local test render):
# - Dropped ggnewscale entirely -- it was never actually needed (fill and
#   colour are independent ggplot2 aesthetics; panel b only ever used one
#   fill scale). Removes an install dependency for no visual cost.
# - On-panel titles/subtitles/captions shortened; long explanatory text
#   moved to a single figure-level caption via plot_annotation(), which is
#   also what was causing per-panel captions to be clipped at panel edges.
# - Panel d: legend simplified from two long wrapped strings ("Gastro-
#   intestinal (higher redundancy)" / "Skin / nares (lower redundancy)") to
#   short labels ("GI" / "Skin/nares"), with the redundancy explanation
#   moved to the figure caption. This also removes the point-label text
#   that was overlapping the legend glyphs near the top of the panel.
# - Panel tags (a-d) switched to patchwork's automatic tagging instead of
#   manual ggplot `tag =` labels for consistent placement.
# =============================================================================

library(tidyverse)
library(patchwork)

fig_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Figures/"
dir.create(fig_dir, showWarnings = FALSE)

theme_concept <- theme_classic(base_size = 11) +
  theme(
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9),
    plot.title      = element_text(size = 9.5, face = "bold"),
    plot.subtitle   = element_text(size = 7.5, colour = "grey40"),
    legend.position = "none"
  )

col_western <- "#D6604D"
col_nonwestern <- "#2166AC"
col_gi   <- "#2166AC"
col_skin <- "#D6604D"
col_tip  <- "#B2182B"

# ============================================================================
# PANEL A -- stability landscape (illustrative)
# ============================================================================
# Tilted double-well potential: U = k*((x-xm)^2 - a^2)^2 + s*(x-xm).
# This is the standard construction for a two-minima landscape with
# unequal well depth (a symmetric quartic double-well, tilted linearly).
# The prior version (two subtracted quartic bumps) did not reliably
# produce a genuine second minimum -- the "non-Western" point ended up
# sitting on a rising slope rather than in a basin. Parameters below were
# chosen (see tuning check) to give two clearly separated local minima
# with a visible depth contrast between them.
richness_seq <- seq(55, 130, length.out = 1000)
x1 <- 75; x2 <- 108              # nominal reference positions (labels only)
xm <- (x1 + x2) / 2
a  <- (x2 - x1) / 2
k  <- 8e-6
s  <- -0.04

potential <- k * ((richness_seq - xm)^2 - a^2)^2 + s * (richness_seq - xm)
potential <- (potential - min(potential)) / (max(potential) - min(potential))
landscape_df <- tibble(richness = richness_seq, U = potential)

find_local_minima <- function(x, y) {
  n <- length(y)
  idx <- which(y[2:(n-1)] < y[1:(n-2)] & y[2:(n-1)] < y[3:n]) + 1
  tibble(x = x[idx], y = y[idx])
}
minima <- find_local_minima(landscape_df$richness, landscape_df$U)
ball_w_x  <- minima$x[which.min(minima$x)]
ball_nw_x <- minima$x[which.max(minima$x)]
ball_w_U  <- minima$y[which.min(minima$x)] + 0.02
ball_nw_U <- minima$y[which.max(minima$x)] + 0.02

figA <- ggplot(landscape_df, aes(x = richness, y = U)) +
  geom_ribbon(aes(ymin = 0, ymax = U), fill = "grey88", alpha = 0.55) +
  geom_line(linewidth = 1.1, colour = "grey45") +
  annotate("point", x = ball_w_x,  y = ball_w_U,  size = 5, colour = col_western) +
  annotate("point", x = ball_nw_x, y = ball_nw_U, size = 5, colour = col_nonwestern) +
  annotate("text", x = ball_w_x,  y = ball_w_U  + 0.10,
           label = sprintf("Western\n(~%.0f genera)", ball_w_x), size = 2.7,
           colour = col_western, lineheight = 0.9) +
  annotate("text", x = ball_nw_x - 1, y = ball_nw_U + 0.10,
           label = sprintf("Non-Western\n(~%.0f genera)", ball_nw_x), size = 2.7,
           colour = col_nonwestern, lineheight = 0.9) +
  scale_x_continuous(breaks = round(c(60, ball_w_x, ball_nw_x, 125)),
                     labels = as.character(round(c(60, ball_w_x, ball_nw_x, 125)))) +
  scale_y_continuous(limits = c(0, 0.78)) +
  labs(
    x = "Genus-level richness (no. genera)",
    y = "Ecological potential\n(lower = more stable)",
    title = "Richer systems occupy a deeper\nstability basin (illustrative)"
  ) +
  theme_concept

# ============================================================================
# PANEL B -- REBUILT: fully conceptual, no data points on a redundancy axis.
# Per manuscript Methods: "no panel plots empirical richness data against a
# redundancy axis as though that relationship were known." Richness ranges
# on the x-axis are real (literature/IHAT-GUT derived, shown as a rug), but
# the y-axis (redundancy) carries no data points -- only a wide, explicitly
# labelled uncertainty band around a hypothesised (untested) trend.
# ============================================================================
richness_rug <- tribble(
  ~population,            ~richness,
  "Western\n(typical)",     79,
  "IHAT-GUT\n(ill)",         94,
  "IHAT-GUT\n(not ill)",    103,
  "Non-Western\nagrarian",  111
)

rich_range <- seq(40, 140, length.out = 200)
hyp_mid   <- 20 + 0.55 * (rich_range - 40)
hyp_lo    <- pmax(0, hyp_mid - 28)
hyp_hi    <- pmin(100, hyp_mid + 28)
hyp_df <- tibble(richness = rich_range, mid = hyp_mid, lo = hyp_lo, hi = hyp_hi)

figB <- ggplot(hyp_df, aes(x = richness)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = col_tip, alpha = 0.12) +
  geom_line(aes(y = mid), colour = col_tip, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 105, y = 85, label = "hypothesised relationship\n(direction untested)",
           size = 2.6, colour = col_tip, fontface = "italic", lineheight = 0.9) +
  geom_rug(data = richness_rug, aes(x = richness), sides = "b",
           length = unit(0.05, "npc"), linewidth = 0.9, colour = "grey30",
           inherit.aes = FALSE) +
  scale_x_continuous(limits = c(40, 140), breaks = c(50, 75, 100, 125)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  labs(
    x = "Genus-level richness (observed; literature/IHAT-GUT-derived)",
    y = "Functional redundancy\n(not measured -- axis is conceptual)",
    title = "Richness is observed; redundancy is\nhypothesised, not assigned",
    caption = "Ticks: Western, IHAT-GUT (ill/not ill), non-Western reference richness (Table 1)"
  ) +
  theme_concept +
  theme(plot.caption = element_text(size = 6.5, colour = "grey45", hjust = 0))

# ============================================================================

# PANEL C -- body-site recovery: speed, not success
# ============================================================================
t_vec <- seq(0, 14, length.out = 400)
perturb_on <- 3.0; perturb_off <- 6.5

gi_traj <- case_when(
  t_vec < perturb_on  ~ 100,
  t_vec < perturb_off ~ 100 - 12 * (t_vec - perturb_on) / (perturb_off - perturb_on),
  TRUE ~ pmin(100, 88 + 12 * (1 - exp(-0.9 * (t_vec - perturb_off))))
)
skin_traj <- case_when(
  t_vec < perturb_on  ~ 100,
  t_vec < perturb_off ~ 100 - 18 * (t_vec - perturb_on) / (perturb_off - perturb_on),
  TRUE ~ pmin(100, 82 + 18 * (1 - exp(-0.28 * (t_vec - perturb_off))))
)

phase_df <- bind_rows(
  tibble(t = t_vec, composition = gi_traj,   system = "GI"),
  tibble(t = t_vec, composition = skin_traj, system = "Skin/nares")
) %>% mutate(system = factor(system, levels = c("GI", "Skin/nares")))

pal_phase <- c("GI" = col_gi, "Skin/nares" = col_skin)

figC <- ggplot(phase_df, aes(x = t, y = composition, colour = system)) +
  annotate("rect", xmin = perturb_on, xmax = perturb_off, ymin = 70, ymax = 104,
           fill = "grey50", alpha = 0.10) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pal_phase, name = NULL) +
  scale_y_continuous(limits = c(70, 104), breaks = c(75, 85, 95),
                     labels = c("75%", "85%", "95%")) +
  scale_x_continuous(breaks = c(0, perturb_on, perturb_off, 14),
                     labels = c("Baseline", "Onset", "Offset", "Long-term")) +
  labs(
    x = "Time", y = "Composition vs pre-flight",
    title = "Both body sites recover,\nbut at different rates"
  ) +
  theme_concept +
  theme(
    legend.position   = c(0.78, 0.22),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.8, "lines"),
    legend.background = element_rect(fill = "white", colour = NA),
    axis.text.x       = element_text(size = 8)
  )

# ============================================================================
# PANEL D -- recovery time vs perturbation magnitude, GI vs skin/nares
# ============================================================================
pmag <- seq(0, 100, by = 0.5)
skin_rec <- 2 + 0.10 * pmag + 0.012 * pmag^2
gi_rec   <- 1.5 + 0.04 * pmag + 0.003 * pmag^2

basin_df <- bind_rows(
  tibble(perturb = pmag, rec_time = skin_rec, system = "Skin/nares"),
  tibble(perturb = pmag, rec_time = gi_rec,   system = "GI")
) %>% mutate(system = factor(system, levels = c("GI", "Skin/nares")))

emp_pts <- tibble(
  perturb  = c(30, 30, 80),
  rec_time = c(8,  24, 55),
  system   = c("GI", "Skin/nares", "Skin/nares"),
  label    = c("Spaceflight\n(GI)", "Spaceflight\n(skin/nares)", "Antibiotics*")
) %>% mutate(system = factor(system, levels = c("GI", "Skin/nares")))

pal_basin <- c("GI" = col_gi, "Skin/nares" = col_skin)

figD <- ggplot(basin_df, aes(x = perturb, y = rec_time, colour = system)) +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(data = emp_pts, aes(x = perturb, y = rec_time, colour = system),
             shape = 17, size = 3.2, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    data = emp_pts, aes(x = perturb, y = rec_time, label = label, colour = system),
    size = 2.4, lineheight = 0.85, inherit.aes = FALSE,
    seed = 1, segment.size = 0.3, min.segment.length = 0, show.legend = FALSE
  ) +
  scale_colour_manual(values = pal_basin, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%","25%","50%","75%","100%")) +
  scale_y_continuous(limits = c(0, 65), breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(
    x = "Perturbation magnitude (illustrative)",
    y = "Recovery time (a.u.)",
    title = "Lower-redundancy sites take longer\nto recover, at a given magnitude"
  ) +
  theme_concept +
  theme(
    legend.position   = c(0.22, 0.85),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.8, "lines"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# ============================================================================
# ASSEMBLE
# ============================================================================
fig4 <- ((figA | figB) / (figC | figD)) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "All panels are conceptual illustrations of resilience/redundancy theory; none are fitted to data (see Methods).",
      "Panel b: x-axis richness values are observed (literature/IHAT-GUT); no population is assigned a redundancy value -- the y-axis relationship shown is hypothesised, not measured.",
      "Panel c/d: GI vs skin/nares is a body-site comparison within one cohort, not a Western/non-Western population comparison (see Limitations).",
      "*Antibiotics point (Wipperman et al. 2017) reflects a larger perturbation magnitude than spaceflight, not equal severity.",
      sep = "\n"
    )
  ) &
  theme(
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.caption = element_text(size = 7, colour = "grey45", hjust = 0, lineheight = 1.2)
  )

for (fmt in c("pdf", "tiff")) {
  outfile <- file.path(fig_dir, paste0("NaturePaper_Fig4_REVISED_v2.", fmt))
  ggsave(outfile, fig4, width = 10, height = 8.5,
         dpi = if (fmt == "tiff") 300 else 150)
  cat("Saved:", outfile, "\n")
}
cat("\nFig 4 v2: dropped ggnewscale dependency, shortened on-panel text,\n")
cat("simplified panel d legend/labels (needs 'ggrepel' package), auto-tags.\n")
