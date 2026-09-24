# =============================================================================
# NaturePaper_Fig3_Spaceflight.R  [REVISED v2 for ISME Communications resubmission]
#
# Figure 3 | Compositional dynamics of the astronaut gut microbiome during
#            spaceflight (reanalysis of Voorhies et al. 2019, Sci Rep)
#
# v2 CHANGES (post local test render):
# - On-panel titles/subtitles/captions shortened substantially. Full
#   explanatory detail now lives ONLY in the manuscript figure legend
#   (Methods/Results text), not duplicated on-figure -- this was causing
#   text to be clipped/truncated at panel edges.
# - Panel tags (a-e) switched to patchwork's automatic tagging
#   (plot_annotation(tag_levels = "a")) instead of manual ggplot `tag =`
#   labels, which were positioning inconsistently across panels.
# - Panel E: Skin and Nares previously plotted as two rows with identical
#   overlapping labels ("still shifted... >=60 days" printed twice, one on
#   top of the other). Combined into a single "Skin / nares" row.
# - Panel A: fewer inline annotations; AstB label moved to a compact legend
#   instead of free text next to the line.
# =============================================================================

library(tidyverse)
library(patchwork)

fig_dir <- "C:/Users/oofordile/Desktop/IHAT_Paper2_Figures/"
dir.create(fig_dir, showWarnings = FALSE)

theme_nature <- theme_classic(base_size = 11) +
  theme(
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    plot.title         = element_text(size = 9.5, face = "bold"),
    plot.subtitle      = element_text(size = 7.5, colour = "grey40"),
    legend.position    = "none",
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4)
  )

col_western <- "#D6604D"
col_astb    <- "#4DAF4A"
col_gi      <- "#2166AC"
col_skin    <- "#D6604D"

# ============================================================================
# PANEL A
# ============================================================================
div_group <- tibble(
  x = 1:5,
  phase_lbl = c("Pre-flight", "FD7", "FD90", "FD180", "Post"),
  diversity = c(100, 104, 118, 120, 105),
  se        = c(6, 7, 9, 8, 7)
)
div_astb <- tibble(x = 1:5, diversity = rep(100, 5))

figA <- ggplot(div_group, aes(x = x, y = diversity)) +
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = 85, ymax = 135,
           fill = col_western, alpha = 0.06) +
  geom_ribbon(aes(ymin = diversity - se, ymax = diversity + se),
              fill = "grey75", alpha = 0.4) +
  geom_line(colour = col_western, linewidth = 1) +
  geom_point(colour = col_western, size = 3) +
  geom_line(data = div_astb, aes(x, diversity),
            colour = col_astb, linetype = "dashed", linewidth = 1) +
  geom_point(data = div_astb, aes(x, diversity),
             colour = col_astb, shape = 17, size = 3) +
  scale_x_continuous(breaks = 1:5, labels = div_group$phase_lbl, limits = c(0.7, 5.3)) +
  scale_y_continuous(limits = c(85, 135), breaks = c(90,100,110,120,130),
                     labels = paste0(c(90,100,110,120,130), "%")) +
  labs(
    x = NULL, y = "GI diversity (% baseline)",
    title = "Diversity rises in lower-baseline astronauts;\nAstB (highest baseline) unchanged"
  ) +
  theme_nature

# ============================================================================
# PANEL B
# ============================================================================
taxa_shifts <- tribble(
  ~taxon, ~log2fc,
  "Akkermansia", -2.32,
  "Ruminococcus", -2.32,
  "Pseudobutyrivibrio", -1.58,
  "Fusicatenibacter", -1.58,
  "Parasutterella", +1.20,
  "Faecalibacterium", +0.80
) %>%
  mutate(taxon = factor(taxon, levels = rev(taxon)),
         col = ifelse(log2fc < 0, col_gi, col_western))

figB <- ggplot(taxa_shifts, aes(log2fc, taxon, fill = col)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_identity() +
  labs(
    x = "Log2 fold-change (in-flight vs pre-flight)",
    title = "Immunologically relevant taxa\ndecline during flight"
  ) +
  theme_nature

# ============================================================================
# PANEL C
# ============================================================================
set.seed(42)
astro_ids <- paste0("Ast", c("A", "C", "D", "E"))
t_pts <- 1:5

conv_df <- expand_grid(id = astro_ids, x = t_pts) %>%
  left_join(tibble(id = astro_ids, base = c(88, 92, 96, 90)), by = "id") %>%
  mutate(
    target = 112,
    frac = c(0, 0.55, 0.85, 1, 0.7)[x],
    diversity = base + (target - base) * frac + rnorm(n(), 0, 1.5)
  )
astb_df <- tibble(id = "AstB", x = t_pts, diversity = 118 + rnorm(5, 0, 0.6))

figC <- ggplot(conv_df, aes(x, diversity, group = id)) +
  annotate("rect", xmin = 1.5, xmax = 4, ymin = 80, ymax = 128,
           fill = col_western, alpha = 0.06) +
  geom_line(colour = col_western, alpha = 0.55, linewidth = 0.8) +
  geom_line(data = astb_df, aes(x, diversity, group = id),
            colour = col_astb, linewidth = 1.1, linetype = "dashed") +
  scale_x_continuous(breaks = t_pts,
                     labels = c("Pre-flight", "FD7", "FD90", "FD180", "Post"),
                     limits = c(0.7, 5.3)) +
  labs(
    x = NULL, y = "GI diversity (illustrative)",
    title = "Schematic: 4 astronauts converge;\nAstB remains apart"
  ) +
  theme_nature

# ============================================================================
# PANEL D
# ============================================================================
immune_df <- tribble(
  ~label, ~value,
  "VZV reactivation (4/9)", 44,
  "IL-8 increase", 72,
  "IL-1beta increase", 68,
  "IL-2 increase", 60
)

figD <- ggplot(immune_df, aes(value, reorder(label, value))) +
  geom_point(size = 4, colour = col_western) +
  geom_segment(aes(x = 0, xend = value, yend = label)) +
  scale_x_continuous(limits = c(0,100)) +
  labs(x = "%", y = NULL, title = "Immune correlates of\ntaxon depletion") +
  theme_nature

# ============================================================================
# PANEL E -- Skin and nares combined into one row to remove duplicate labels
# ============================================================================
recovery_df <- tribble(
  ~site,              ~recover_low, ~recover_high, ~colour,
  "Gastrointestinal",  0,  60,  col_gi,
  "Skin / nares",      60, 120, col_skin
) %>%
  mutate(site = factor(site, levels = c("Skin / nares", "Gastrointestinal")))

figE <- ggplot(recovery_df, aes(y = site, colour = colour)) +
  geom_segment(aes(x = recover_low, xend = recover_high, yend = site),
               linewidth = 7, lineend = "round", alpha = 0.35) +
  geom_point(aes(x = recover_low), size = 3) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 130), breaks = c(0, 30, 60, 90, 120)) +
  labs(
    x = "Days after return to Earth", y = NULL,
    title = "Compositional recovery time\ndiffers by body site"
  ) +
  theme_nature +
  theme(panel.grid.major.y = element_blank())

# ============================================================================
# COMBINE -- patchwork auto-tagging, single figure-level caption
# ============================================================================
fig3 <- ((figA | figB) / (figC | figD) / (figE | plot_spacer())) +
  plot_layout(heights = c(1, 1, 0.75)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "n = 5 of 9 astronauts sampled for in-flight stool; AstB = highest pre-flight-diversity astronaut, shown separately throughout.",
      "Panels c and e are illustrative summaries of Voorhies et al. (2019) results, not new statistical fits;",
      "associations in panel d are descriptive/correlational (see Methods).",
      sep = "\n"
    )
  ) &
  theme(
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.caption = element_text(size = 7.5, colour = "grey45", hjust = 0, lineheight = 1.2)
  )

ggsave(file.path(fig_dir, "Fig3_REVISED_v2.pdf"), fig3, width = 9, height = 10)
cat("Fig 3 v2 saved: shortened on-panel text, deduped panel e, patchwork auto-tags.\n")
