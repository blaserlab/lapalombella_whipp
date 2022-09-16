#source("R/figs/Ethan_Figs_081622.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"


#F5top <- plot_grid(NULL, labels = "A")

F5top <- plot_grid(
  NULL,
  NULL,
  NULL,
  ncol = 3,
  rel_widths = c(1.4, 1, 1),
  labels = c("A", "B", "C")
)

F5rowtwo <- plot_grid(
  F5D,
  NULL,
  ncol = 2,
  rel_widths = c(1.2, 1),
  labels = c("D", "E")
)
F5rowthree <- plot_grid(
  F5F,
  NULL,
  ncol = 2,
  rel_widths = c(1.2, 1),
  labels = c("F", NULL)
)

F5D_F5F <- plot_grid(F5top,
                F5rowtwo,
                F5rowthree,
                nrow = 3,
                rel_heights = c(0.1, 1, 2.25, 3))

save_plot(F5D_F5F,
          filename = "tempF5.png",
          base_width = 7.5,
          base_height = 9.75)
#ggsave("F5D_F5F.pdf", path = T_Figs, width = 7.5, height = 9.75)

F5Etop <- plot_grid(
  F5E_1,
 labels = "E"
)
F5Etwo <- plot_grid(
  F5E_2,
  labels = NULL
)
F5Ethree <- plot_grid(
  F5E_3,
  labels = NULL
)
F5Efour <- plot_grid(
  F5E_2,
  labels = NULL
)
F5Espacer1 <- plot_grid(
  NULL)
F5E <- plot_grid(F5Espacer1,F5Etop,
                     F5Etwo, F5Espacer1,
                     F5Ethree,
                     F5Efour,
                     nrow = 6,
                     rel_heights = c(0.25,1, 1,0.25, 1, 1))

save_plot(F5E,
          filename = "tempF5E.png",
          base_width = 7.5,
          base_height = 9.75)
ggsave("F5E.pdf", path = T_Figs, width = 7.5, height = 9.75)

#Supplemental Fig4 E, F, & dotplot

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
                rel_heights = c(0.5, 1, 1))

save_plot(S1,
          filename = "tempS1.png",
          base_width = 7.5,
          base_height = 9.75)

#ggsave("S1_leiden.pdf", path = T_Figs, width = 7.5, height = 9.75)