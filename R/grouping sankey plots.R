library(ggpubr)
library(png)
library(grid)

# 1. Combine plots with panel tags and a common legend if applicable
multi_plot <- ggarrange(
    p1, p2, p3, p4,
    labels = "AUTO", 
    ncol = 2, 
    nrow = 2,
    common.legend = TRUE,
    legend = "bottom"
)

# 2. Add global title or annotations
multi_plot_annotated <- annotate_figure(
    multi_plot,
    top = text_grob("Figure 1: Taxa-Pathways links in camda cohort", face = "bold", size = 14)
)

# 3. Export to publication-quality vector PDF
ggsave(
    filename = "figure1_sankey_grid1.pdf",
    plot = multi_plot_annotated,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300,
    device = cairo_pdf # Ensures embedded fonts and clean rendering
)