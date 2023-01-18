#source("R/figs/Ethan_Figs_081622.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Figures"

# F1top <- plot_grid(NULL, labels = "A")

F6row1 <- plot_grid(
  NULL,
  NULL,
  ncol = 2,
  rel_widths = c(1.5, 1),
  labels = c("A", "B")
)

F6row2 <- plot_grid(
  NULL, labels = c("C")
)
F6row3 <- plot_grid(F6heatmap, NULL, ncol=2, rel_widths = c(0.9,1))#, labels = c("D","E"))

F6hm <- plot_grid(F6row1,
                       F6row2,
                       F6row3,
                       nrow = 3,
                       rel_heights = c(0.75, 1, 1.5)) #0.75,1.2,1

save_plot(F6hm,
          filename = "tempF6.png",
          base_width = 7.5,
          base_height = 9.75)
#ggsave("F6hm.pdf", path = T_Figs, width = 7.5, height = 9.75)

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
