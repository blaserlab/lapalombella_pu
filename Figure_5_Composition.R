#source("R/figs/Ethan_Figs_081622.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Figures"

# F1top <- plot_grid(NULL, labels = "A")

F5row1 <- plot_grid(
  F5A,
  NULL,
  ncol = 2,
  rel_widths = c(1.75, 1),
  labels = c("A", "B")
)

F5row2 <- plot_grid(
  aml_gexp_umap, labels = c("C")
)
F5row3 <- plot_grid(F5hm, NULL, ncol=2, rel_widths = c(1,1), labels = "D")

F5 <- plot_grid(F5row1,
                       F5row2,
                       F5row3,
                       nrow = 3,
                       rel_heights = c(0.75, 1.2, 1))

save_plot(F5,
          filename = "tempF5.png",
          base_width = 7.5,
          base_height = 9.75)
#ggsave("F5.pdf", path = T_Figs, width = 7.5, height = 9.75)

#Supplemental Fig1

S1top <- plot_grid(
  S1A,
  S1B,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("A", "B")
)

S1rowtwo <- plot_grid(
  S1C,
  S1D,
  ncol = 2,
  rel_widths = c(1, 2.8),
  labels = c("C", "D")
)
S1rowthree <- plot_grid(S1E, labels = "E")


S1 <- plot_grid(S1top,
                S1rowtwo,
                S1rowthree,
                nrow = 3,
                rel_heights = c(1, 0.7, 1.2))

save_plot(S1,
          filename = "tempS1.png",
          base_width = 7.5,
          base_height = 9.75)

ggsave("S1_leiden.pdf", path = T_Figs, width = 7.5, height = 9.75)
