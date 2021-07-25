load("Figures/resubmission/onset_correlogram_plot.Rdata")
load("Figures/resubmission/offset_correlogram_plot.Rdata")
load("Figures/resubmission/duration_correlogram_plot.Rdata")

cp <- cowplot::plot_grid(onset_correlogram_plot, offset_correlogram_plot,
                         duration_correlogram_plot,
                         ncol = 1, labels = c("Emergence",
                                              "Termination",
                                              "Duration"),
                         label_x = 0.099)
cp

library(ggplot2)
ggsave("Figures/resubmission/MoransI.png", width = 8, height = 6)
