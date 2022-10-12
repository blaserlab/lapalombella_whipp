
#WalkerAccess <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"
#WalkerTables <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs//Tables/Heatmap Tables"

colData(mouse_cds_list[[2]])$genotype <- recode(colData(mouse_cds_list[[2]])$genotype,
                                                "PRMT5" = "PRMT5/TCL1",
                                                "TCL1" = "TCL1",
                                                "P/T" = "PRMT5/TCL1")
#unique(mouse_cds_list[[2]]$tissue)

#Fig5 Figs
#recode and harmoize clusters
colData(mouse_cds_list[[2]])$kmeans_10_harmonized <- recode(colData(mouse_cds_list[[2]])$kmeans_10_clusters, 
                                                            "1" = "5.1",
                                                            "2" = "5.2",
                                                            "3" = "5.3", 
                                                            "4" = "5.4",
                                                            "5" = "5.5",
                                                            "6" = "5.6", 
                                                            "7" = "5.7",
                                                            "8" = "5.8",
                                                            "9" = "5.9", 
                                                            "10" = "5.10")

colData(mouse_cds_list[[2]])$kmeans_10_harmonized <- factor(colData(mouse_cds_list[[2]])$kmeans_10_harmonized, 
                                                            levels = paste0("5.", 1:10))
#cell type calling
# bb_genebubbles(
#   obj = mouse_cds_list[[2]],
#   genes = c(
#     "Cd19",
#     "Ms4a1",
#     "Cd79a",
#     "Cd3d",
#     "Cd4",
#     "Cd14",
#     "Itgam",
#     "Cd8a",
#     "Cd177",
#     "Pdcd1",
#     "Foxp3"
#   ),
#   cell_grouping = "kmeans_10_harmonized"
# ) 

#cluster/cell type assignment
colData(mouse_cds_list[[2]])$k_10_assignment <- recode(colData(mouse_cds_list[[2]])$kmeans_10_harmonized, "5.1" = "B", "5.2" = "Cd8+ T", "5.3" = "B", "5.4" = "Cd4+ T", "5.5" = "B", "5.6" = "B", "5.7" = "Neutro", "5.8" = "B", "5.9" = "B", "5.10" = "Mono")

# make logical values for CD19+CD5+ cells
mat <- monocle3::exprs(mouse_cds_list[[2]])

cd19_tbl <- colnames(mat[ ,mat["ENSMUSG00000030724", ] > 0]) |> as_tibble() |> mutate(Cd19_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[2]] <- bb_tbl_to_coldata(mouse_cds_list[[2]], min_tbl = cd19_tbl)

cd5_tbl <- colnames(mat[ ,mat["ENSMUSG00000024669", ] > 0]) |> as_tibble() |> mutate(Cd5_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[2]] <- bb_tbl_to_coldata(mouse_cds_list[[2]], min_tbl = cd5_tbl)


colData(mouse_cds_list[[2]])$cd19_cd5_pos <- colData(mouse_cds_list[[2]])$Cd19_pos & colData(mouse_cds_list[[2]])$Cd5_pos  

mouse_cds_list[[2]] <- bb_cellmeta(mouse_cds_list[[2]]) |> 
  filter(cd19_cd5_pos) |> 
  select(cell_id) |> 
  mutate(cd19_cd5_label = "CD19+/CD5+ cells") |> 
  bb_tbl_to_coldata(mouse_cds_list[[2]], min_tbl = _)
#colData(mouse_cds_list[[2]])$cd19_cd5_label

######################### Figure 5: main figs
F5D1 <- bb_var_umap(mouse_cds_list[[2]], "k_10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")

F5D2 <- bb_var_umap(mouse_cds_list[[2]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
F5D <- F5D1 / F5D2 + plot_layout(heights = c(1.5,1))

############F5E
F5E_plotlist <- map(.x = c("Myc", "Mki67", "Egr1", "Cxcr5", "Ccr7", "Il10", "Ctla4", "Cd274"),
                    .f = \(x, dat = mouse_cds_list[[2]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25 #adjusted cell size - default is 0.5
                      ) +
                        scale_color_distiller(palette = "Oranges", #?scale_color_distiller(...)
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) + #fix scale limit 
                        facet_wrap( ~ genotype) + scale_y_continuous(breaks = c(-10,0, 10)) + 
                        scale_x_continuous(breaks = c(-5,5, 15))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white")) +
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        #labs(x = NULL, y = x) +
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(strip.text = element_blank()) + theme(legend.position = "none") #+
                        #theme(legend.title = element_blank()) 
                      #if (x != "Myc") p <- p + theme(strip.text = element_blank())
                      p 
                    })

#F5E_plotlist[[1]] | F5E_plotlist[[2]] | F5E_plotlist[[3]] | F5E_plotlist[[4]] + theme(legend.position = "top")

F5E1<-
  F5E_plotlist[[1]]/ plot_spacer() / bb_gene_violinplot(
    filter_cds(mouse_cds_list[[2]],
               cells = bb_cellmeta(mouse_cds_list[[2]]) |>
                 filter(k_10_assignment == "B")),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1 , 1))

F5E2 <- F5E_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Mki67",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1, -0.1 ,1))
  
F5E3 <-F5E_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Egr1",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1, -0.1 ,1))
  
F5E4 <- (F5E_plotlist[[4]]+theme(legend.position = "right",
                               legend.key.size = unit(4, "mm"),
                               legend.margin = margin(c(0, -7, 0, -6.5)),
                               legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cxcr5",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

F5E_1 <- 
  F5E1 |
  F5E2 |
  F5E3 |
  F5E4
#ggsave("F5E_1.pdf", path = WalkerAccess, width = 11.6, height = 4.4)

F5E5 <-
  F5E_plotlist[[5]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Ccr7",
    pseudocount = 0
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1, 1))
F5E6 <- F5E_plotlist[[6]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
F5E7 <- F5E_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Ctla4",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1, 1))
F5E8 <- F5E_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cd274",
    pseudocount = 0
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
F5E_2 <- 
  F5E5 |
  F5E6 |
  F5E7 |
  F5E8
ggsave("F5E_2.pdf", path = WalkerAccess, width = 11.6, height = 4.4)

######Fig 5F
#kmeans10 UMAP
F5F1<- 
  bb_var_umap(mouse_cds_list[[2]], var = "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) +
  facet_wrap(~genotype)
#ggsave("fig5_k10_harmonized.pdf", path = figures_out, width = 7.8, height = 4.7)

#Heatmap/topmarkers
#top markers
F5_k10_Top50_tm <-monocle3::top_markers(filter_cds(mouse_cds_list[[2]], 
                                                   cells = bb_cellmeta(mouse_cds_list[[2]]) |> 
                                                     filter(kmeans_10_harmonized %in% c("5.1", "5.3", "5.5", "5.6", "5.9"))), group_cells_by = "kmeans_10_harmonized", genes_to_test_per_group = 50, cores = 10)
# F5_k10_Top50_clust5.1_bygenotype <- monocle3::top_markers(filter_cds(mouse_cds_list[[2]],
#                                                      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
#                                                        filter(kmeans_10_harmonized == "5.1")), group_cells_by = "genotype", genes_to_test_per_group = 50, cores = 10)
#write_csv(F5_k10_Top50_tm, file = file.path(WalkerTables, "F5_k10_Top50_tm.csv"))
F5_k10_Top50_tm <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Fig5_PRMT5xTCL1_vs_TCL1/F5_k10_Top50_tm.csv")

markers <- F5_k10_Top50_tm |> 
  filter(cell_group %in% c("5.1","5.6")) |> 
  pull(gene_short_name)

fig5_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[2]], 
                                          cells = bb_cellmeta(mouse_cds_list[[2]]) |> 
                                            filter(kmeans_10_harmonized %in% c("5.1", "5.3", "5.5", "5.6", "5.8","5.9")),
                                          genes = bb_rowmeta(mouse_cds_list[[2]]) |> 
                                            filter(gene_short_name %in% markers)), 
                         cell_group_df = bb_cellmeta(mouse_cds_list[[2]]) |> 
                           select(cell_id, kmeans_10_harmonized)) |> 
  t() |> 
  scale() |> 
  t()

rownames(fig5_mat) <- tibble(feature_id = rownames(fig5_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[2]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

# fig5_heatmap_mat <- as.data.frame(fig5_mat)|> mutate(GeneID = rownames(fig5_mat)) |> select(GeneID, everything())
# write_csv(fig5_heatmap_mat, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Fig5_PRMT5xTCL1_vs_TCL1", "fig5_heatmap_mat.csv"))
# write_csv(fig5_heatmap_mat, file = file.path("~/network/T/Labs/EHL/Rosa/EXPERIMENTS/PRMT5/Manuscripts/Drafts/Drafts_2022/Source Data", "fig5_heatmap_mat.csv"))



fig5_mat <- tibble(feature_id = rownames(fig5_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[2]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

# fig5_heatmap_mat <- fig5_mat |> as.data.frame()
# write_csv(fig5_heatmap_mat, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Fig5_PRMT5xTCL1_vs_TCL1", "fig5_heatmap_mat.csv"))
fig5_colfun = circlize::colorRamp2(breaks = c(min(fig5_mat),
                                              0,
                                              max(fig5_mat)),
                                   colors = heatmap_3_colors)

#highlights <- c("Cd83", "Tubb5", "Tuba1b","Wnt10a", "Ccr7", "Cdk1", "Birc5","Junb","Atf4","Cd69","Nfkbia","Nfkbid","Rel","Ubc","Ccna2", "Blnk")
#highlights <-c("Ccr7", "Cdk4", "Cxcr5", "Birc5", "Il4","Npm1","Jun","Junb", "Fos","Fosb","Atf3","Atf4","Myc","Cd69","Il10", "Top2a", "Hmgb1", "Hmgb2", "Cd83","Ube2a", "Tubb5", "Tuba1b")
manuscript_highlights <- c("Tubb5", "Tuba1b", "Npm1") #manuscript highlights in dataset
#Additional lymphoma associated genes: DisGeNet.org - Lymphoma; CUI: C0024299
lymphoma_genes<- readxl::read_excel("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/disgenet.org_LymphomaGenes_C0024299_disease_gda_summary.xlsx")[3]
lymphoma_genes$Gene <- str_to_title(lymphoma_genes$Gene)
filt<- F5_k10_Top50_tm |> filter(F5_k10_Top50_tm$gene_short_name %in% lymphoma_genes$Gene)
lymphoma_gois <- filter(filt, cell_group %in% c('5.6'))[["gene_short_name"]]
#Human top markers
RTLN_tm <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Fig1 human RT data/leiden clustering/LN_B_leiden_Top50_tm.csv")
RTLN_tm$gene_short_name <- str_to_title(RTLN_tm$gene_short_name)
RTLN_tm <- filter(RTLN_tm, cell_group == c('3','11'))
human_overlap<- F5_k10_Top50_tm |> filter(F5_k10_Top50_tm$gene_short_name %in% RTLN_tm$gene_short_name)
human_overlap <- filter(human_overlap, cell_group %in% c('5.6'))[["gene_short_name"]]
highlights <- unique(c(manuscript_highlights, lymphoma_gois, human_overlap))


fig5_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig5_mat) %in% highlights),
  labels = rownames(fig5_mat)[rownames(fig5_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2.5
))

F5F2 <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(fig5_mat, 
                        col = fig5_colfun, 
                        name = "Expression", 
                        show_row_names = F, 
                        right_annotation = fig5_anno,
                        row_dend_width = unit(4, "mm"),
                        column_dend_height = unit(4, "mm"),
                        heatmap_legend_param = list(legend_direction = "vertical",
                                                    #legend_width = unit(1, "mm"),
                                                    title_position = "lefttop-rot", 
                                                    title_gp = gpar(fontsize = 6)
                        ))))

F5F <- F5F1 / F5F2 + plot_layout(heights = c(1, 1.75))
F5F

# F5_k10_Top50_tm |> filter(cell_group == "5.6") |> arrange(desc(marker_score))
# fig5_kmeans_10_tm_Top100 |> filter(cell_group == "5.6") |> arrange(desc(mean_expression))

##########################Supp Fig 4
S4E <- bb_var_umap(
  mouse_cds_list[[2]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  labs(x = "UMAP 1", y= "UMAP 2") +
  theme(legend.justification = "center")
#ggsave("S4E.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs", height = 3.8, width = 5.4)

S4F_plotlist <- map(.x = c("Cd93", "Il4"),
                    .f = \(x, dat = mouse_cds_list[[2]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25 #adjusted cell size - default is 0.5
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,1)) + #set fixed scale (limits)
                        facet_wrap( ~ genotype) + scale_y_continuous(breaks = c(-10,0, 10)) + 
                        scale_x_continuous(breaks = c(-5,5, 15))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white")) +
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        #labs(x = NULL, y = x) +
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(strip.text = element_blank()) + theme(legend.position = "none")
                      p 
                    })

#Pseudobulk
#B cell CDS

colData(mouse_cds_list[[2]])$genotype <-
  recode(
    colData(mouse_cds_list[[2]])$genotype,
    "PRMT5/TCL1" = "PT",
    "TCL1" = "TCL1")

colData(mouse_cds_list[[2]])$geno_k10assign <-
  paste0(colData(mouse_cds_list[[2]])$genotype, "_", colData(mouse_cds_list[[2]])$k_10_assignment)

f5_k10_B_cds<- filter_cds(mouse_cds_list[[2]],
                          cells = bb_cellmeta(mouse_cds_list[[2]]) |>
                            filter(k_10_assignment == "B"))
unique(colData(f5_k10_B_cds)$kmeans_10_harmonized)

exp_design <- 
  bb_cellmeta(f5_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()
exp_design

exp_design <- 
  bb_cellmeta(f5_k10_B_cds) |>
  group_by(genotype, k_10_assignment) |>
  summarise()
#exp_design
f5_k10_B_cds_copy <- f5_k10_B_cds 
rowData(f5_k10_B_cds_copy)$id <- rownames(rowData(f5_k10_B_cds_copy))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = f5_k10_B_cds_copy,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PT", "TCL1"))

#less conservative approach
#bb_monocle_regression(cds = f5_k10_B_cds, gene_or_genes = "Cd93", form = "~genotype")

pseudobulk_res$Header 

pseudobulk_res$Result |> filter(gene_short_name == "Cd93")
pseudobulk_res$Result |> filter(gene_short_name == "Il4")
Fig5_pseudobulk<- pseudobulk_res$Result
write.csv(Fig5_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Fig5_pseudobulk.csv")
#symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

S4F1 <- S4F_plotlist[[1]]/plot_spacer()/bb_gene_violinplot(
  filter_cds(
    mouse_cds_list[[2]],
    cells = bb_cellmeta(mouse_cds_list[[2]]) |>
      filter(k_10_assignment == "B")),
  variable = "genotype",
  genes_to_plot = "Cd93",
  pseudocount = 0
) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
  t = 0,
  r = -1.5,
  b = 0,
  l = -7
))) + stat_compare_means(method = "t.test") + labs(y = "B cell Expression") +plot_layout(heights = c(1,-0.1,1))
S4F2 <- (S4F_plotlist[[2]]+theme(legend.position = "right",
                                legend.key.size = unit(4, "mm"),
                                legend.margin = margin(c(0, -7, 0, -6.5)),
                                legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(
  filter_cds(
    mouse_cds_list[[2]],
    cells = bb_cellmeta(mouse_cds_list[[2]]) |>
      filter(k_10_assignment == "B")
  ),
  variable = "genotype",
  genes_to_plot = "Il4",
  pseudocount = 0) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) +plot_layout(heights = c(1,-0.1,1))

S4F<- S4F1|S4F2
S4F
#ggsave("S4F.pdf", path = WalkerAccess, height = 2.93, width = 4.59)

S4_genebubble <- 
  bb_genebubbles(
    mouse_cds_list[[2]],
    genes = c(
      "Ms4a1",
      "Cd79a",
      "Cd19",
      "Cd3d",
      "Cd4",
      "Cd14",
      "Itgam",
      "Cd8a",
      "Cd177",
      "Pdcd1",
      "Foxp3"
    ),
    cell_grouping = c("kmeans_10_harmonized", "k_10_assignment"),
    return_value = "data"
  ) |> 
  ggplot(mapping = aes(x = kmeans_10_harmonized, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~k_10_assignment, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
#ggsave("S4_genebubble.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs", height = 6.1, width = 5.8)

#TODO Supp GO terms Supp Fig 4G

####################################################################################################################################

#Graph based Analysis

#fig5_graphbased_UMAP<- bb_var_umap(mouse_cds_list[[2]], var = "graphclust", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) +
#  facet_wrap(~genotype)
#ggsave("fig5_graphbased_UMAP.pdf", path = WalkerAccess, width = 7.8, height = 4.7)
# cds_sub <- mouse_cds_list[[2]]
# fig5_graphbased_genebubbles<- bb_genebubbles(
#   obj = filter_cds(cds_sub, cells = bb_cellmeta(cds_sub)),
#   genes = c("Cd19", "Ms4a1", "Cd79a", "Cd3d","Cd4","Cd14","Itgam","Cd8a","Cd177","Pdcd1","Foxp3"),
#   cell_grouping = "graphclust")
# #ggsave("fig5_graphbased_genebubbles.pdf", path = WalkerAccess, width = 8, height = 5.5)
# 
# 
# colData(cds_sub)$gbased_assignment <- recode(colData(cds_sub)$graphclust, "1" = "Cd8+ T", "2" = "B", "3" = "B",
#                                              "4" = "B", "5" = "B", "6" = "B","7" = "B", "8" = "Neutrophils",
#                                              "9" = "Cd8+ T", "10" = "Cd4+ T", "11" = "B", "12" = "Monocytes",
#                                              "13" = "B", "14" = "Cd8+ T", "15" = "B", "16" = "Cd8+ T", "17" = "B",
#                                              "18" = "B", "19" = "B", "20" = "Cd4+ T", "21" = "B", "22" = "B", 
#                                              "23" = "Monocytes", "24" = "Cd4+ T", "25" = "B", "26" = "Neutrophils", 
#                                              "27" = "B", "28" = "Cd8+ T", "29" = "Monocytes", "30" = "B",
#                                              "31" = "B", "32" = "Monocytes", "33" = "B")
# fig5d <- bb_var_umap(cds_sub, "gbased_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")
# 
# fig5_gb_UMAP <- bb_var_umap(cds_sub, "graphclust", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")
# ggsave("fig5_gb_UMAP.pdf", path = WalkerAccess, width = 10, height = 3.9)
# 
# 
# #Heatmap/topmarkers
# #top markers
# unique(cds_sub$genotype)
# SplnB_gb<-cds_sub[,colData(cds_sub)$gbased_assignment == "B"]
#  SplnB_gb_tm<- monocle3::top_markers(SplnB_gb,
#                                         group_cells_by = "graphclust", genes_to_test_per_group = 100, cores = 10)
# SplnB_gb_300tm<- monocle3::top_markers(SplnB_gb,
#                       group_cells_by = "graphclust", genes_to_test_per_group = 300, cores = 10)
# unique(SplnB_gb$genotype)
# write_csv(SplnB_gb_300tm, file = file.path(WalkerTables, "SplnB_gb_300tm.csv"))
# # Fig5_Spln_AllB_gb_Top300_byGenotype<- monocle3::top_markers(SplnB_gb,
# #                                     group_cells_by = "genotype", genes_to_test_per_group = 300, cores = 10)
# # unique(SplnB_gb$genotype)
# # unique(Fig5_Spln_AllB_gb_Top100_byGenotype$genotype)
# # write_csv(Fig5_Spln_AllB_gb_Top300_byGenotype, file = file.path(WalkerTables, "Fig5_Spln_AllB_gb_Top300_byGenotype.csv"))
# 
# 
# #All B comparison btwn genotypes Top Markers
# Fig5_Spln_AllB_gb_Top100 <- monocle3::top_markers(SplnB_gb,
#                                             group_cells_by = "genotype", genes_to_test_per_group = 100, cores = 10)
# # 
# # 
# write_csv(Fig5_Spln_AllB_gb_Top100, file = file.path(WalkerTables, "Fig5_Spln_AllB_gb_Top100.csv"))
# 
# markers <- SplnB_gb_tm |> 
#   filter(cell_group %in% c("6", "21", "15", "5", "2", "17", "4", "3", "11", "7")) |>#c("2","3","4","5","6","7","11","13","15","17","18","19","21","22","25","27","30","31","33")) |> 
#   pull(gene_short_name)
# 
# fig5_gb_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[2]], 
#                                           cells = bb_cellmeta(mouse_cds_list[[2]]) |> 
#                                             filter(graphclust %in% c("6", "21", "15", "5", "2", "17", "4", "3", "11", "7")),#c("2","3","4","5","6","7","11","13","15","17","18","19","21","22","25","27","30","31","33")),
#                                           genes = bb_rowmeta(mouse_cds_list[[2]]) |> 
#                                             filter(gene_short_name %in% markers)), 
#                          cell_group_df = bb_cellmeta(mouse_cds_list[[2]]) |> 
#                            select(cell_id, graphclust)) |> 
#   t() |> 
#   scale() |> 
#   t()
# 
# rownames(fig5_gb_mat) <- tibble(feature_id = rownames(fig5_gb_mat)) |> 
#   left_join(bb_rowmeta(mouse_cds_list[[2]]) |> 
#               select(feature_id, gene_short_name)) |> 
#   pull(gene_short_name)
# 
# fig5_gb_colfun = circlize::colorRamp2(breaks = c(min(fig5_gb_mat),
#                                               0,
#                                               max(fig5_gb_mat)),
#                                    colors = heatmap_3_colors)
# 
# highlights <- c("Myc")
# 
# fig5_gb_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
#   at = which(rownames(fig5_gb_mat) %in% highlights),
#   labels = rownames(fig5_gb_mat)[rownames(fig5_gb_mat) %in% highlights],
#   labels_gp = gpar(fontsize = 8),
#   padding = 0
# ))
# 
# ComplexHeatmap::Heatmap(fig5_gb_mat, 
#                         col = fig5_gb_colfun, 
#                         name = "Expression", show_row_names = F, right_annotation = fig5_gb_anno)
# 
# SplnB_gb_tm |> filter(cell_group == "5.6") |> arrange(desc(marker_score))
# SplnB_gb_tm |> filter(cell_group == "21") |> arrange(desc(mean_expression))
# 
# bb_gene_umap(
#   mouse_cds_list[[2]],
#   gene_or_genes = "Mki67",
#   alt_dim_x = "aggr_UMAP_1",
#   alt_dim_y = "aggr_UMAP_2"
# ) + facet_wrap(~genotype)