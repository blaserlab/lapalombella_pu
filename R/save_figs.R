# source("R/Figure_6.R")
# source("R/Figure_7.R")

save_plot("F6A.pdf", plot = F6A, base_height = 5.2, base_width = 8, path = figs_out)
save_plot("F6E.pdf", plot = F6E, base_height = 7, base_width = 5, path = figs_out)

save_plot("F7_microenv_map.pdf", plot = F7_microenv_map, path = figs_out)
save_plot("F7_microenv_bp.pdf", plot = F7_microenv_bp, path = figs_out)
save_plot("gb_plot2.pdf", plot = gb_plot2, path = figs_out)
save_plot("gb_plot.pdf", plot = gb_plot2, path = figs_out)

save_plot("SF9A.pdf", plot = leiden_clust, base_height = 6, base_width = 6, path = figs_out)
