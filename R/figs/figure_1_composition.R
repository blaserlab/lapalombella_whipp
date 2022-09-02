# source("R/Ethan_Figs_081622.R")

F1top <- plot_grid(NULL, labels = "A")

F1rowtwo <- plot_grid(
  F1A,
  NULL,
  ncol = 2,
  rel_widths = c(1.3, 1),
  labels = c("B", "C")
)

F1rowthree <- plot_grid(
  F1D,
  NULL,
  #heatmap goes here
  ncol = 2,
  rel_widths = c(1.3, 1),
  labels = c("D", "E")
)

F1 <- plot_grid(F1top,
                F1rowtwo,
                F1rowthree,
                nrow = 3,
                rel_heights = c(0.1, 1, 1.5))

save_plot(F1,
          filename = "temp1.png",
          base_width = 7.5,
          base_height = 9.75)
