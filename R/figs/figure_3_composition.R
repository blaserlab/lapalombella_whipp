#source("R/figs/Ethan_Figs_081622.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"
#T_tables <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables"



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
#ggsave("F3.pdf", path = T_Figs, width = 7.5, height = 9.75)

#Supplemental Fig2
S2E_1 <- (S3E1 | S3E2 | S3E3 | S3E4)
# save_plot(S2E_1,
#           filename = "tempS2E_1.png",
#           base_width = 7.5,
#           base_height = 3)
#ggsave("S2E_1.pdf", path = T_Figs, width = 7.5, height = 3)

S2E_2 <- (S3E5 | S3E6 | S3E7 | S3E8 | S3E9)
# save_plot(S2E_2,
#           filename = "tempS2E_1.png",
#           base_width = 7.5,
#           base_height = 3)
#ggsave("S2E_2.pdf", path = T_Figs, width = 7.5, height = 3)




