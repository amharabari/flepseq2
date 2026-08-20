## ------------------------------------------------------------------
## FIGURE 2 - Differential analysis 
## ------------------------------------------------------------------

short_up   <- readRDS("r_plots/objects/short_diffcount_up.rds")
short_down <- readRDS("r_plots/objects/short_diffcount_down.rds")

long_up   <- readRDS("r_plots/objects/long_diffcount_up.rds")
long_down <- readRDS("r_plots/objects/long_diffcount_down.rds")

uri_up   <- readRDS("r_plots/objects/uridylation_diffcount_up.rds")
uri_down <- readRDS("r_plots/objects/uridylation_diffcount_down.rds")

short_diff <- readRDS("r_plots/objects/short_diffcount_ALL.rds")
long_diff <- readRDS("r_plots/objects/long_diffcount_ALL.rds")
uri_diff <- readRDS("r_plots/objects/uridylation_diffcount_ALL.rds")


panelA_diff <-  wrap_elements(short_diff)   # stacked up/down bars
panelB_diff <- wrap_elements(long_diff)   
panelC_diff <- wrap_elements(uri_diff)

library(patchwork) 

fig2 <- panelA_diff | panelB_diff | panelC_diff
fig2 <- fig2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))
ggsave("r_plots/panels/test.svg", fig2,width = 10, height = 5, dpi = 300, scale = 1.3)
ggsave("r_plots/panels/test.png", fig2,width = 10, height = 5, dpi = 300, scale = 1.3)



panelA_diff <- short_up / short_down
panelB_diff <- long_up  / long_down
panelC_diff <- uri_up   / uri_down

fig2 <- panelA_diff | panelB_diff | panelC_diff

fig2 <- fig2 +
  plot_layout(guides = "collect") +          # merge the 6 duplicate legends into 1
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag        = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text       = element_text(size = 12),
    axis.title      = element_text(size = 12),
    strip.text      = element_text(size = 12)
  )

ggsave("r_plots/panels/differential_counts_overview.svg", fig2,
       width = 9, height = 5, dpi = 300,scale = 1.3 )
ggsave("r_plots/panels/differential_counts_overview.png", fig2,
       width = 14, height = 6, dpi = 300,c)