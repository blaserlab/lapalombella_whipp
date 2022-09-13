#source("R/figs/Ethan_Figs_081622.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"
T_tables <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables"



F3top <- plot_grid(
  F3A,
  F3B,
  ncol = 2,
  rel_widths = c(1.4,1),
  labels = c("A", "B"))

F3rowtwo <- plot_grid(
  F3C,
  F3D,
  ncol = 2,
  rel_widths = c(1, 1.3),
  labels = c("C", "D" )
)

F3 <- plot_grid(F3top,
                F3rowtwo,
                #F3rowthree,
                nrow = 2,
                rel_heights = c(0.4, 1))

save_plot(F3,
          filename = "tempF3.png",
          base_width = 7.5,
          base_height = 9.75)
ggsave("F3.pdf", path = T_Figs, width = 7.5, height = 9.75)

#Supplemental Fig1

S3top <- plot_grid(
  NULL,
  NULL,
  ncol = 2,
  rel_widths = c(0.2, 1),
  labels = c("A", "B")
)

S3rowtwo <- plot_grid(
  NULL,
  NULL,
  S3Ea,
  ncol = 3,
  rel_widths = c(0.2, 0.2, 1),
  labels = c("C", "D", "E")
)
S3rowthree <- plot_grid(
  NULL,
  S3Eb,
  ncol = 2,
  rel_widths = c(0.2, 1),
  labels = c("F", NULL ))
  
S3rowfour <- plot_grid(
  NULL,
  NULL,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c(NULL, "G")
)


S3 <- plot_grid(S3top,
                S3rowtwo,
                S3rowthree,
                S3rowfour,
                nrow = 4,
                rel_heights = c(1, 1, 1, 1))

save_plot(S3,
          filename = "tempS3.png",
          base_width = 7.5,
          base_height = 9.75)

#ggsave("S3.pdf", path = T_Figs, width = 7.5, height = 9.75)


