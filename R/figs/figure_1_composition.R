#source("R/figs/Ethan_Figs_081622.R")
#T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"
#T_Tables <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables"

F1top <- plot_grid(NULL, labels = "A")

F1rowtwo <- plot_grid(
  F1A,
  NULL,
  ncol = 2,
  rel_widths = c(1.4, 1),
  labels = c("B", "C")
)

F1rowthree <- plot_grid(
  F1D,
  F1E,
  ncol = 2,
  rel_widths = c(1.4, 1),
  labels = c("D", "E")
)

F1_leiden <- plot_grid(F1top,
                F1rowtwo,
                F1rowthree,
                nrow = 3,
                rel_heights = c(0.15, 1, 1.25))

save_plot(F1_leiden,
          filename = "tempF1.png",
          base_width = 7.5,
          base_height = 9.75)
ggsave("F1_leiden.pdf", path = T_Figs, width = 7.5, height = 9.75)

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

