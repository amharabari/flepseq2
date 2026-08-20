## ==================================================================
## assemble_panels.R
## Loads the individually saved ggplot objects and combines them into
## the four composite report figures using patchwork.
## Run AFTER all four .Rmd scripts have been knitted at least once
## (so every r_plots/objects/*.rds file exists).
## ==================================================================

library(patchwork)
library(ggplot2)

dir.create("r_plots/panels", showWarnings = FALSE, recursive = TRUE)

## Small helper to strip a plot down for use as a sub-panel:
## drops the title (panel letter + external caption will identify it)
## and optionally the legend, so it can be shared once across the figure
strip_for_panel <- function(p, keep_legend = FALSE, keep_title = FALSE) {
  if (!keep_title)  p <- p + labs(title = NULL)
  if (!keep_legend) p <- p + theme(legend.position = "none")
  p
}

## ------------------------------------------------------------------
## Helper 2: force a uniform y-axis (0-100%) and a uniform gray strip
## style, overriding whatever each original .Rmd set individually
## ------------------------------------------------------------------
unify_axis_and_strip <- function(p, ylim = c(0, 100), breaks = seq(0, 100, 20)) {
  p +
    scale_y_continuous(limits = ylim, breaks = breaks) +   # replaces old y scale
    coord_cartesian(ylim = ylim) +                          # replaces old coord zoom
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text       = element_text(face = "bold", color = "black")
    )
}

## ------------------------------------------------------------------
## FIGURE 1 - Global overview across tail categories
## (Medium acts as the natural, built-in negative control)
## ------------------------------------------------------------------
short_global  <- readRDS("r_plots/objects/boxshort_global.rds")
medium_global <- readRDS("r_plots/objects/box_Medium_global.rds")
long_global   <- readRDS("r_plots/objects/boxLong_global.rds")
uri_global    <- readRDS("r_plots/objects/boxplot_U_per_gene.rds")

## apply the y-axis / strip fix to each panel FIRST
short_global  <- unify_axis_and_strip(short_global)
medium_global <- unify_axis_and_strip(medium_global)
long_global   <- unify_axis_and_strip(long_global)
uri_global    <- unify_axis_and_strip(uri_global)

fig1 <- (strip_for_panel(short_global)  + labs(title = "Short poly(A) tails")) +
        (strip_for_panel(medium_global) + labs(title = "Medium poly(A) tails")) +
        (strip_for_panel(long_global)   + labs(title = "Long poly(A) tails")) +
  #  (strip_for_panel(uri_global, keep_legend = TRUE) + labs(title = "Uridylation")) +
  plot_layout(nrow = 1, ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        plot.tag = element_text(face = "bold", size = 16))


ggsave("r_plots/panels/Figure1_global_overview.pdf", fig1, width = 11, height = 5, dpi = 30)
ggsave("r_plots/panels/Figure1_global_overview.png", fig1, width = 11, height = 5, dpi = 300)
ggsave("r_plots/panels/Figure1_global_overview.svg", fig1, width = 12, height = 5, dpi = 300, scale = 1.3)


## ------------------------------------------------------------------
## FIGURE 2 - Short poly(A) tail differential analysis (A/B/C)
## ------------------------------------------------------------------
short_up   <- readRDS("r_plots/objects/short_diffcount_up.rds")
short_down <- readRDS("r_plots/objects/short_diffcount_down.rds")
short_diff <- readRDS("r_plots/objects/short_diff_box.rds")
short_comb <- readRDS("r_plots/objects/short_combined_box.rds")

short_diff   <- unify_axis_and_strip(short_diff)
short_comb   <- unify_axis_and_strip(short_comb)


panelA_short <- strip_for_panel(short_up)   / strip_for_panel(short_down)   # stacked up/down bars
panelB_short <- strip_for_panel(short_diff, keep_legend = TRUE)
panelC_short <- strip_for_panel(short_comb)

fig2 <- (panelA_short | panelB_short) / panelC_short
fig2 <- fig2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))


ggsave("r_plots/panels/Figure2_short_differential.pdf", fig2, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure2_short_differential.png", fig2, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure2_short_differential.svg", fig2, width = 16, height = 6, dpi = 300)

## Explicit grid layout: A and B occupy the top row side by side,
## C spans the FULL width of the row underneath both of them.
## Each letter's number of cells controls its relative width/height.
design <- "
AB
AB
CC
"

fig2 <- wrap_plots(A = panelA_short, B = panelB_short, C = panelC_short, design = design)

fig2 <- fig2 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))


## ------------------------------------------------------------------
## FIGURE 3 - Long poly(A) tail differential analysis (A/B/C)
## ------------------------------------------------------------------
long_up   <- readRDS("r_plots/objects/long_diffcount_up.rds")
long_down <- readRDS("r_plots/objects/long_diffcount_down.rds")
long_diff <- readRDS("r_plots/objects/long_diff_box.rds")
long_comb <- readRDS("r_plots/objects/long_combined_box.rds")

long_diff   <- unify_axis_and_strip(long_diff)
long_comb   <- unify_axis_and_strip(long_comb)


panelA_long <- strip_for_panel(long_up) / strip_for_panel(long_down)
panelB_long <- strip_for_panel(long_diff, keep_legend = TRUE)
panelC_long <- strip_for_panel(long_comb)

fig3 <- panelA_long | panelB_long | panelC_long
fig3 <- fig3 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("r_plots/panels/Figure3_long_differential.pdf", fig3, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure3_long_differential.png", fig3, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure3_long_differential.svg", fig3, width = 16, height = 6, dpi = 300)


## ------------------------------------------------------------------
## FIGURE 4 - Uridylation differential analysis (A/B/C)
## ------------------------------------------------------------------
uri_up   <- readRDS("r_plots/objects/uridylation_diffcount_up.rds")
uri_down <- readRDS("r_plots/objects/uridylation_diffcount_down.rds")
uri_diff <- readRDS("r_plots/objects/uridylation_diff_box.rds")
uri_comb <- readRDS("r_plots/objects/uridylation_combined_box.rds")

uri_diff   <- unify_axis_and_strip(uri_diff)
uri_comb   <- unify_axis_and_strip(uri_comb)


panelA_uri <- strip_for_panel(uri_up) / strip_for_panel(uri_down)
panelB_uri <- strip_for_panel(uri_diff, keep_legend = TRUE)
panelC_uri <- strip_for_panel(uri_comb)

fig4 <- panelA_uri | panelB_uri | panelC_uri
fig4 <- fig4 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("r_plots/panels/Figure4_uridylation_differential.pdf", fig4, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure4_uridylation_differential.png", fig4, width = 16, height = 6, dpi = 300)
ggsave("r_plots/panels/Figure4_uridylation_differential.svg", fig4, width = 16, height = 6, dpi = 300)


cat("All four composite figures written to r_plots/panels/\n")
