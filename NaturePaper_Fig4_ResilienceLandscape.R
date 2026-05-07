# =============================================================================
# NaturePaper_Fig4_ResilienceLandscape.R  [FIXED]
# Ecological Collapse Perspective — Nature
#
# FIX: approx() was called with x and U vectors of the same name as the
#      outer variables x and U, causing a length mismatch when the tibble
#      had already been created. Fixed by using distinct variable names for
#      the potential function computation, and by computing ball positions
#      from the data frame rather than the raw vectors.
#
# Figure 4 | Conceptual model: microbiome resilience as a function of
#             ecological redundancy
#
# All panels are conceptual/theoretical. No primary data required.
#
# Author:  Ogochukwu Ofordile
# =============================================================================

library(tidyverse)
library(patchwork)

fig_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Figures/"
dir.create(fig_dir, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# THEME
# --------------------------------------------------------------------------- #
theme_concept <- theme_classic(base_size = 11) +
  theme(
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9),
    plot.title      = element_text(size = 10, face = "bold"),
    plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
    plot.tag        = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )

col_western    <- "#D6604D"
col_nonwestern <- "#2166AC"
col_tip        <- "#B2182B"

# ============================================================================
# PANEL A: Stability landscape (double-well potential)
#
# FIX: renamed inner vectors to richness_seq / potential to avoid collision
#      with the landscape_df tibble columns; ball positions computed from
#      the data frame using which.min() rather than approx() on raw vectors.
# ============================================================================

richness_seq <- seq(55, 130, length.out = 500)

x1 <- 75    # Western attractor (genera)
x2 <- 108   # Non-Western attractor
xt <- 89    # Tipping point

# Asymmetric double-well: Western basin deliberately shallower
potential <- 0.0012 * (richness_seq - x1)^2 * (richness_seq - xt)^2 / 200 -
             0.0022 * 0.3 * (richness_seq - x2)^2 * (richness_seq - xt)^2 / 900 +
             0.00005 * (richness_seq - 90)^2

# Rescale to [0, 1]
potential <- (potential - min(potential)) / (max(potential) - min(potential))

landscape_df <- tibble(richness = richness_seq, U = potential)

# Ball and tipping-point positions: use which.min / which.min on landscape_df
# to avoid any vector-length mismatch
find_U <- function(target_x, df) {
  idx <- which.min(abs(df$richness - target_x))
  df$U[idx]
}

ball_w_x  <- x1 + 0.8
ball_nw_x <- x2 - 0.5
ball_w_U  <- find_U(ball_w_x,  landscape_df) + 0.025
ball_nw_U <- find_U(ball_nw_x, landscape_df) + 0.025
tip_x     <- xt
tip_U     <- find_U(xt, landscape_df)

figA <- ggplot(landscape_df, aes(x = richness, y = U)) +
  geom_ribbon(aes(ymin = 0, ymax = U), fill = "grey88", alpha = 0.55) +
  geom_line(linewidth = 1.1, colour = "grey45") +
  # balls
  annotate("point", x = ball_w_x,  y = ball_w_U,
           size = 5, colour = col_western,    shape = 19) +
  annotate("point", x = ball_nw_x, y = ball_nw_U,
           size = 5, colour = col_nonwestern, shape = 19) +
  # tipping point
  annotate("point", x = tip_x, y = tip_U,
           size = 3, shape = 25, fill = col_tip, colour = col_tip) +
  annotate("text",  x = tip_x, y = tip_U + 0.07,
           label = "tipping\npoint", size = 2.8, colour = col_tip,
           fontface = "italic", hjust = 0.5, lineheight = 0.9) +
  # ball labels
  annotate("text", x = ball_w_x,  y = ball_w_U  + 0.12,
           label = "Western\n(~75 genera)", size = 2.8,
           colour = col_western, hjust = 0.5, lineheight = 0.9) +
  annotate("text", x = ball_nw_x - 1, y = ball_nw_U + 0.12,
           label = "Non-Western\n(~108 genera)", size = 2.8,
           colour = col_nonwestern, hjust = 0.5, lineheight = 0.9) +
  # basin depth arrows
  annotate("segment",
           x = ball_w_x, xend = ball_w_x, y = 0.01, yend = ball_w_U - 0.025,
           colour = col_western, linewidth = 0.5,
           arrow = arrow(ends = "both", length = unit(0.07, "inches"))) +
  annotate("text", x = ball_w_x + 4, y = 0.05,
           label = "shallow basin", size = 2.5, colour = col_western, hjust = 0) +
  annotate("segment",
           x = ball_nw_x, xend = ball_nw_x, y = 0.01, yend = ball_nw_U - 0.025,
           colour = col_nonwestern, linewidth = 0.5,
           arrow = arrow(ends = "both", length = unit(0.07, "inches"))) +
  annotate("text", x = ball_nw_x - 10, y = 0.07,
           label = "deep basin", size = 2.5, colour = col_nonwestern, hjust = 0) +
  scale_x_continuous(
    breaks = c(60, x1, xt, x2, 125),
    labels = c("60", sprintf("%d\n(Western)", x1),
               sprintf("%d\n(threshold)", xt),
               sprintf("%d\n(non-Western)", x2), "125")
  ) +
  scale_y_continuous(limits = c(0, 0.85)) +
  labs(
    x        = "Genus-level richness (no. genera)",
    y        = "Ecological potential (lower = more stable)",
    title    = "Alternative stable states: Western system sits in shallower basin",
    subtitle = "Asymmetric double-well potential; constructed to illustrate resilience theory",
    tag      = "a"
  ) +
  theme_concept +
  theme(axis.text.x = element_text(size = 8.5, lineheight = 0.9))

# ============================================================================
# PANEL B: Redundancy gradient — population positions
# ============================================================================
redundancy_df <- tribble(
  ~population,             ~redundancy, ~richness, ~rich_sd, ~colour,
  "Western\ndysbiotic",     10,          52,        15,       "#B2182B",
  "Western\n(typical)",     28,          79,        12,       col_western,
  "IHAT-GUT\n(ill)",        55,          94,        25,       "#969696",
  "IHAT-GUT\n(not ill)",    70,         103,        25,       col_nonwestern,
  "Non-Western\nagrarian",  88,         111,        18,       "#08519C"
) %>%
  mutate(
    population = factor(population,
                        levels = c("Western\ndysbiotic", "Western\n(typical)",
                                   "IHAT-GUT\n(ill)", "IHAT-GUT\n(not ill)",
                                   "Non-Western\nagrarian"))
  )

tip_redundancy <- 42

figB <- ggplot(redundancy_df,
               aes(x = redundancy, y = richness, colour = colour)) +
  annotate("rect",
           xmin = tip_redundancy - 8, xmax = tip_redundancy + 8,
           ymin = 35, ymax = 135,
           fill = col_tip, alpha = 0.07) +
  annotate("text", x = tip_redundancy, y = 133,
           label = "vulnerability\nzone", size = 2.7,
           colour = col_tip, hjust = 0.5, fontface = "italic", lineheight = 0.9) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              colour = "grey70", linewidth = 0.7, linetype = "dashed") +
  geom_errorbar(aes(ymin = richness - rich_sd, ymax = richness + rich_sd),
                width = 3, linewidth = 0.7) +
  geom_point(size = 5) +
  geom_text(aes(label = population),
            vjust = -1.1, size = 2.7, lineheight = 0.9, colour = "grey20") +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100")
  ) +
  scale_y_continuous(limits = c(30, 145)) +
  labs(
    x        = "Functional redundancy (arbitrary units)",
    y        = "Mean genus-level richness \u00b1 SD",
    title    = "Western microbiomes cluster in low-redundancy region",
    subtitle = "IHAT-GUT not-ill children: 103 genera (non-Western range)\nRichness values grounded in IHAT-GUT data; redundancy axis conceptual",
    tag      = "b"
  ) +
  theme_concept

# ============================================================================
# PANEL C: Perturbation trajectory — low vs high redundancy
# ============================================================================
t_vec        <- seq(0, 14, length.out = 400)
perturb_on   <- 3.0
perturb_off  <- 6.5
tip_threshold <- 65

w_traj <- case_when(
  t_vec < perturb_on  ~ 79,
  t_vec < perturb_off ~ 79 - 16 * (t_vec - perturb_on) / (perturb_off - perturb_on),
  TRUE                ~ pmax(63 - 4.5 * (t_vec - perturb_off), 44)
)

nw_traj <- case_when(
  t_vec < perturb_on  ~ 100,
  t_vec < perturb_off ~ 100 - 12 * (t_vec - perturb_on) / (perturb_off - perturb_on),
  TRUE                ~ 88 + 12 * (1 - exp(-0.6 * (t_vec - perturb_off)))
)

phase_df <- bind_rows(
  tibble(t = t_vec, richness = w_traj,  system = "Western (low redundancy)"),
  tibble(t = t_vec, richness = nw_traj, system = "Non-Western (high redundancy)")
) %>%
  mutate(system = factor(system,
                         levels = c("Non-Western (high redundancy)",
                                    "Western (low redundancy)")))

pal_phase <- c("Western (low redundancy)"      = col_western,
               "Non-Western (high redundancy)" = col_nonwestern)

# ISS empirical anchors (approximate; verify against primary sources)
iss_pts <- tibble(
  t        = c(3.0, 4.5, 6.5, 8.0, 11.0),
  richness = c(79,   67,  57,  62,   68),
  label    = c("Pre-\nflight", "Wk 2", "Wk 10", "R+1\nwk", "R+6\nmo")
)

figC <- ggplot(phase_df, aes(x = t, y = richness, colour = system)) +
  annotate("rect",
           xmin = perturb_on, xmax = perturb_off, ymin = 30, ymax = 108,
           fill = "grey50", alpha = 0.10) +
  annotate("text",
           x = (perturb_on + perturb_off) / 2, y = 107,
           label = "Perturbation", size = 2.9, colour = "grey40", hjust = 0.5) +
  geom_hline(yintercept = tip_threshold, linetype = "dashed",
             colour = col_tip, linewidth = 0.55) +
  annotate("text", x = 13.9, y = tip_threshold + 2.5,
           label = "tipping\nthreshold", size = 2.7, colour = col_tip,
           hjust = 1, lineheight = 0.9) +
  geom_line(linewidth = 1.1) +
  geom_point(data = iss_pts, aes(x = t, y = richness),
             shape = 17, size = 3, colour = col_western, inherit.aes = FALSE) +
  geom_text(data = iss_pts, aes(x = t, y = richness, label = label),
            vjust = -0.6, hjust = 0.5, size = 2.5, colour = col_western,
            lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(values = pal_phase, name = NULL) +
  scale_y_continuous(
    limits = c(30, 112),
    breaks = c(40, 60, 80, 100),
    labels = c("40%", "60%", "80%", "100%")
  ) +
  scale_x_continuous(
    breaks = c(0, perturb_on, perturb_off, 14),
    labels = c("Baseline", "Onset", "Offset", "Long-term")
  ) +
  labs(
    x        = "Time",
    y        = "Microbiome integrity (% baseline)",
    title    = "Low-redundancy system crosses tipping point; incomplete recovery",
    subtitle = "Triangles = approximate ISS data; lines = model prediction",
    tag      = "c"
  ) +
  theme_concept +
  theme(
    legend.position   = c(0.78, 0.82),
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(0.85, "lines"),
    legend.background = element_rect(fill = "white", colour = NA),
    axis.text.x       = element_text(size = 8.5)
  )

# ============================================================================
# PANEL D: Recovery basin — perturbation magnitude vs recovery time
# ============================================================================
pmag <- seq(0, 100, by = 0.5)

w_rec  <- ifelse(pmag < 55,
                 2 + 0.08 * pmag + 0.01 * pmag^2,
                 NA_real_)
nw_rec <- 1.5 + 0.04 * pmag + 0.003 * pmag^2

basin_df <- bind_rows(
  tibble(perturb = pmag, rec_time = w_rec,
         system = "Western (low redundancy)"),
  tibble(perturb = pmag, rec_time = nw_rec,
         system = "Non-Western (high redundancy)")
) %>%
  mutate(system = factor(system,
                         levels = c("Non-Western (high redundancy)",
                                    "Western (low redundancy)")))

# Empirical anchors (approximate; see source notes)
emp_pts <- tibble(
  perturb  = c(30,  30,  20),
  rec_time = c(24,   8,   6),
  system   = c("Western (low redundancy)",
               "Non-Western (high redundancy)",
               "Western (low redundancy)"),
  label    = c("Spaceflight\n(Western)",
               "Spaceflight\n(non-Western est.)",
               "Antibiotics\n(Western)")
) %>%
  mutate(system = factor(system,
                         levels = c("Non-Western (high redundancy)",
                                    "Western (low redundancy)")))

pal_basin <- c("Western (low redundancy)"      = col_western,
               "Non-Western (high redundancy)" = col_nonwestern)

figD <- ggplot(basin_df, aes(x = perturb, y = rec_time, colour = system)) +
  annotate("rect",
           xmin = 53, xmax = 100, ymin = 0, ymax = 65,
           fill = col_tip, alpha = 0.06) +
  annotate("text", x = 76, y = 61,
           label = "Western: past\ntipping point", size = 2.7,
           colour = col_tip, hjust = 0.5, fontface = "italic", lineheight = 0.9) +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(data = emp_pts,
             aes(x = perturb, y = rec_time, colour = system),
             shape = 17, size = 3.5, inherit.aes = FALSE) +
  geom_text(data = emp_pts,
            aes(x = perturb, y = rec_time, label = label, colour = system),
            vjust = -0.55, size = 2.5, lineheight = 0.9, hjust = 0.5,
            inherit.aes = FALSE) +
  scale_colour_manual(values = pal_basin, name = NULL) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_y_continuous(limits = c(0, 65), breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(
    x        = "Perturbation magnitude (% diversity removed)",
    y        = "Expected recovery time (arbitrary units)",
    title    = "Shallower basin = slower recovery at any perturbation magnitude",
    subtitle = "Conceptual; triangles = approximate empirical anchors",
    tag      = "d"
  ) +
  theme_concept +
  theme(
    legend.position   = c(0.28, 0.78),
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(0.85, "lines"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# ============================================================================
# ASSEMBLE & SAVE
# ============================================================================
fig4 <- (figA | figB) / (figC | figD) +
  plot_layout(heights = c(1.1, 1))

for (fmt in c("pdf", "tiff")) {
  outfile <- file.path(fig_dir, paste0("NaturePaper_Fig4.", fmt))
  ggsave(outfile, fig4,
         width  = 10,
         height = 8,
         dpi    = if (fmt == "tiff") 300 else 150)
  cat("Saved:", outfile, "\n")
}

cat("\nFig 4 complete. All panels are conceptual.\n")
cat("Panel B richness values grounded in IHAT-GUT data.\n")
cat("ISS anchor points in panel C are approximate; verify before submission.\n")
