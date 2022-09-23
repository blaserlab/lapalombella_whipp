WalkerAccess <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs"
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Supp5"
unique(colData(mouse_cds_list[[4]])$tissue)
unique(colData(mouse_cds_list[[4]])$genotype)

#Lymph Node - PRMT5xTCL1 v TCL1
colData(mouse_cds_list[[4]])$genotype <- recode(colData(mouse_cds_list[[4]])$genotype,
                                                "PRMT5" = "PRMT5/TCL1",
                                                "TCL1" = "TCL1",
                                                "P/T" = "PRMT5/TCL1")
unique(colData(mouse_cds_list[[4]])$genotype)

# Number of ea genotype
# length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "TCL1"])) #n=3
# length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "PRMT5/TCL1"])) #n=4

#Harmoize clusters
colData(mouse_cds_list[[4]])$kmeans_10_harmonized <- recode(colData(mouse_cds_list[[4]])$kmeans_10_clusters, 
                                                            "1" = "LN.1",
                                                            "2" = "LN.2",
                                                            "3" = "LN.3", 
                                                            "4" = "LN.4",
                                                            "5" = "LN.5",
                                                            "6" = "LN.6", 
                                                            "7" = "LN.7",
                                                            "8" = "LN.8",
                                                            "9" = "LN.9", 
                                                            "10" = "LN.10")

colData(mouse_cds_list[[4]])$kmeans_10_harmonized <- factor(colData(mouse_cds_list[[4]])$kmeans_10_harmonized, 
                                                            levels = paste0("LN.", 1:10))
#Cell type assignment
# bb_genebubbles(
#   obj = filter_cds(mouse_cds_list[[4]], cells = bb_cellmeta(mouse_cds_list[[4]])),
#   genes = c("Cd19", "Ms4a1", "Cd79a", "Cd3d","Cd4","Cd14","Itgam","Cd8a","Cd177","Pdcd1","Foxp3"), cell_grouping = "kmeans_10_clusters")
# 
# 
# LN.8_tm <- monocle3::top_markers(filter_cds(mouse_cds_list[[4]], cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
#                                    filter(kmeans_10_harmonized == "LN.8")), group_cells_by = "genotype", genes_to_test_per_group = 50, cores = 12)

colData(mouse_cds_list[[4]])$k10_assignment <- recode(colData(mouse_cds_list[[4]])$kmeans_10_clusters, "1" = "Cd8+ T", "2" = "B", "3" = "Cd8+ T",
                                             "4" = "B", "5" = "B", "6" = "Cd4+ T","7" = "B", "8" = "Other",
                                             "9" = "Neutrophils", "10" = "Monocytes")
# make logical values for CD19+CD5+ cells
mat <- monocle3::exprs(mouse_cds_list[[4]])

cd19_tbl <- colnames(mat[ ,mat["ENSMUSG00000030724", ] > 0]) |> as_tibble() |> mutate(Cd19_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[4]] <- bb_tbl_to_coldata(mouse_cds_list[[4]], min_tbl = cd19_tbl)

cd5_tbl <- colnames(mat[ ,mat["ENSMUSG00000024669", ] > 0]) |> as_tibble() |> mutate(Cd5_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[4]] <- bb_tbl_to_coldata(mouse_cds_list[[4]], min_tbl = cd5_tbl)


colData(mouse_cds_list[[4]])$cd19_cd5_pos <- colData(mouse_cds_list[[4]])$Cd19_pos & colData(mouse_cds_list[[4]])$Cd5_pos  

mouse_cds_list[[4]] <- bb_cellmeta(mouse_cds_list[[4]]) |> 
  filter(cd19_cd5_pos) |> 
  select(cell_id) |> 
  mutate(cd19_cd5_label = "CD19+/CD5+ cells") |> 
  bb_tbl_to_coldata(mouse_cds_list[[4]], min_tbl = _)

###########################Supp 5 Figs
S5A1<-bb_var_umap(mouse_cds_list[[4]], "k10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F, facet_by = "genotype")+
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank())

S5A2 <- bb_var_umap(mouse_cds_list[[4]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2") + 
  theme(strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.key.size = unit(3, 'mm'))

S5A <- (S5A1/S5A2) 
S5A <- grid.arrange(patchworkGrob(S5A1/S5A2), left = textGrob("UMAP 2", rot = 90, vjust = 1.5), bottom = textGrob("UMAP 1", vjust = -1))
######S5B
S5B_plotlist <- map(.x = c("Myc","Ly6a","Mki67","Ccr7","Cxcr5","Il10","Ctla4","Egr1","Cd274","Cd93"),
                    .f = \(x, dat = mouse_cds_list[[4]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) +
                        facet_wrap( ~ genotype) + scale_x_continuous(breaks = c(-10,0, 10)) + 
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white")) +
                        theme(axis.line = element_blank()) + 
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text.x = element_text()) +
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(legend.position = "none") + theme(strip.text = element_blank())
                      p 
                    })

# S5B_plotlist[[1]]+ theme(panel.spacing = unit(1, "lines"))
#S5B_plotlist[[1]]/S5B_plotlist[[2]]/S5B_plotlist[[3]]/S5B_plotlist[[4]]/S5B_plotlist[[5]]
#S5B_plotlist[[6]]/S5B_plotlist[[7]]/S5B_plotlist[[8]]/S5B_plotlist[[9]]/S5B_plotlist[[10]]

S5B1 <-
  ((S5B_plotlist[[1]]+theme(legend.position = "left",
                           legend.key.size = unit(4, "mm"),
                           legend.margin = margin(c(0, -9, 0, -6)),
                           legend.title = element_blank()))/plot_spacer()/ bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression")) + plot_layout(heights = c(1, -0.08 ,1))

S5B2 <- S5B_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                        cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                          filter(k10_assignment == "B")), variable = "genotype",
                                             genes_to_plot = "Ly6a", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B3 <- S5B_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                        cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                          filter(k10_assignment == "B")), variable = "genotype",
                                             genes_to_plot = "Mki67", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B4 <-S5B_plotlist[[4]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Ccr7", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B5 <-S5B_plotlist[[5]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cxcr5", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))

S5B6 <-
  S5B_plotlist[[6]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.08 , 1))
S5B7 <-S5B_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Ctla4", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B8 <-S5B_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Egr1", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B9 <-S5B_plotlist[[9]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cd274", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))
S5B10 <-(S5B_plotlist[[10]])/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cd93", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.08 ,1))

S5B_1 <- 
  S5B1 |
  S5B2 |
  S5B3 |
  S5B4 |
  S5B5

S5B_2 <- 
  S5B6 |
  S5B7 |
  S5B8 |
  S5B9 |
  S5B10
S5B 

S5B<- S5B_1/plot_spacer()/S5B_2 + plot_layout(heights = c(1, -0.15 ,1))

F3A <- (F3A1/F3A2) 
S5Ba <- grid.arrange(patchworkGrob(S5B))                                          
ggsave("S5B.pdf", path = T_Figs, width = 9.5, height = 5)

#####S5C -Heatmap -RT Clust2&7 and TCL1 Clust4
S5C1<-bb_var_umap(mouse_cds_list[[4]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")
S5A1/S5C1/S5A2  

#Alt CD19/CD5
S5_alt <- bb_var_umap(
  mouse_cds_list[[4]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  labs(y = "UMAP 2", x = "UMAP 1") +
  theme(legend.justification = "center")

bb_gene_umap(mouse_cds_list[[4]], gene_or_genes = "Mki67", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)
bb_gene_umap(mouse_cds_list[[4]], gene_or_genes = "Myc", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)
bb_gene_umap(filter_cds(mouse_cds_list[[4]], cells = bb_cellmeta(mouse_cds_list[[4]]) |> filter(k10_assignment =="B")),  gene_or_genes = "Cd93", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)

#Top Markers of all B cells by genotype (top 50)
S5A1/S5C1  
S5_LN_Bclust_Top50markers<- monocle3::top_markers(
  filter_cds(
    mouse_cds_list[[4]],
    cells = bb_cellmeta(mouse_cds_list[[4]]) |> filter(k10_assignment == "B")
  ),
  group_cells_by = "kmeans_10_harmonized",
  genes_to_test_per_group = 50,
  cores = 12
)

#write_csv(S5_LN_Bclust_Top50markers, file = file.path(WalkerTables, "S5_LN_Bclust_Top50markers.csv"))
S5_LN_Bclust_Top50markers <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables//Supp5_PRMT5xTCL1_vs_TCL1_LN/S5_LN_Bclust_Top50markers.csv")

