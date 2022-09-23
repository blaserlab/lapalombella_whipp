# #source("R/figs/Supp_fig5_staging.R")
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Supp5"

S5top <- plot_grid(
  S5A,
  NULL,
  ncol = 2,
  rel_widths = c(1.5,1),
  labels = c("A", NULL))

S5rowtwo <- plot_grid(S5B_1, labels = "B")
S5rowthree <- plot_grid(S5B_2, labels = NULL)


S5 <- plot_grid(S5top,
                S5rowtwo,
                S5rowthree,
                nrow = 3,
                rel_heights = c(1.75, 1, 1))

save_plot(S5,
          filename = "tempS5.png",
          base_width = 7.5,
          base_height = 9.75)

ggsave("S5.pdf", path = T_Figs, width = 7.5, height = 9.75)

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

