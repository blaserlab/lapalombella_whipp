WalkerAccess <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs"
T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs"
#unique(colData(mouse_cds_list[[4]])$tissue)
#unique(colData(mouse_cds_list[[4]])$genotype)

#Lymph Node - PRMT5xTCL1 v TCL1
colData(mouse_cds_list[[4]])$genotype <- recode(colData(mouse_cds_list[[4]])$genotype,
                                                "PRMT5" = "PRMT5/TCL1",
                                                "TCL1" = "TCL1",
                                                "P/T" = "PRMT5/TCL1")
#unique(colData(mouse_cds_list[[4]])$genotype)

# Number of ea genotype
#length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "TCL1"])) #n=3
#length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "PRMT5/TCL1"])) #n=4

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
S5A1<-bb_var_umap(mouse_cds_list[[4]], "k10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")+
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
#ggsave("S5A.pdf", path = T_Figs)
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
                        #theme(axis.text.y = element_blank()) +       #Strip text
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(legend.position = "none") + theme(strip.text = element_blank())
                      p 
                    })

# S5B_plotlist[[1]]+ theme(panel.spacing = unit(1, "lines"))
#S5B_plotlist[[1]]/S5B_plotlist[[2]]/S5B_plotlist[[3]]/S5B_plotlist[[4]]/S5B_plotlist[[5]]
#S5B_plotlist[[6]]/S5B_plotlist[[7]]/S5B_plotlist[[8]]/S5B_plotlist[[9]]/S5B_plotlist[[10]]

#Pseudobulk
#B cell cds
s5_k10_B_cds<- filter_cds(mouse_cds_list[[4]],
                          cells = bb_cellmeta(mouse_cds_list[[4]]) |>
                            filter(k10_assignment == "B"))

exp_design <- 
  bb_cellmeta(s5_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()
exp_design

#exp_design
s5_k10_B_cds_copy <- s5_k10_B_cds 

rowData(s5_k10_B_cds_copy)$id <- rownames(rowData(s5_k10_B_cds_copy))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = s5_k10_B_cds_copy,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PRMT5/TCL1", "TCL1"))
#less conservative approach (pseudobulk is a very conservative approach)
#bb_monocle_regression(cds = f5_k10_B_cds, gene_or_genes = "Cd93", form = "~genotype")

pseudobulk_res$Header 

Supp5_PRMT5xTCL1vsTCL1_LN_B_pseudobulk<- pseudobulk_res$Result
#write.csv(Supp5_PRMT5xTCL1vsTCL1_LN_B_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Supp5_PRMT5xTCL1vsTCL1_LN_B_pseudobulk.csv")

S5B1 <-
  (S5B_plotlist[[1]]/plot_spacer()/ bb_gene_violinplot(        #+theme(axis.text.y = element_text())
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression")) + plot_layout(heights = c(1, -0.1 ,1))

S5B2 <- S5B_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                        cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                          filter(k10_assignment == "B")), variable = "genotype",
                                             genes_to_plot = "Ly6a", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B3 <- S5B_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                        cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                          filter(k10_assignment == "B")), variable = "genotype",
                                             genes_to_plot = "Mki67", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B4 <-S5B_plotlist[[4]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Ccr7", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B5 <-(S5B_plotlist[[5]]+theme(legend.position = "right",
                                legend.key.size = unit(4, "mm"),
                                legend.margin = margin(c(0, -7, 0, -6.5)),
                                legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cxcr5", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

S5B6 <-
  (S5B_plotlist[[6]]+theme(axis.text.y = element_text())) / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1 , 1))
S5B7 <-S5B_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Ctla4", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B8 <-S5B_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Egr1", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B9 <-S5B_plotlist[[9]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cd274", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5B10 <-(S5B_plotlist[[10]])/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                       cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                         filter(k10_assignment == "B")), variable = "genotype",
                                            genes_to_plot = "Cd93", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

# S5B_1 <- 
#   (S5B1 |plot_spacer()|
#   S5B2 |plot_spacer()|
#   S5B3 |plot_spacer()|
#   S5B4 |plot_spacer()|
#   S5B5) + plot_layout(widths = c(1,-0.125,1,-0.125,1,-0.125,1,-0.125,1))

S5B_1 <- S5B1|S5B2|S5B3|S5B4|S5B5
ggsave("S5B_1.pdf", path = T_Figs, width = 11.6, height = 3.3)

S5B_2 <- 
  (S5B6 |plot_spacer()|
  S5B7 |plot_spacer()|
  S5B8 |plot_spacer()|
  S5B9 |plot_spacer()|
  S5B10)+ plot_layout(widths = c(1,-0.125,1,-0.125,1,-0.125,1,-0.125,1))

S5B_2 <-S5B6|S5B7|S5B8|S5B9|S5B10
ggsave("S5B_2.pdf", path = T_Figs, width = 11.6, height = 3.3)

#S5B<- plot_spacer()/S5B_1/plot_spacer()/S5B_2/plot_spacer() + plot_layout(heights = c(-0.1,1, -0.175 ,1,-0.15))

#Alt CD19/CD5
S5_alt_Cd19Cd5 <- bb_var_umap(
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
ggsave("S5_alt_Cd19Cd5.pdf", path = T_Figs)

#bb_gene_umap(mouse_cds_list[[4]], gene_or_genes = "Mki67", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)
#bb_gene_umap(mouse_cds_list[[4]], gene_or_genes = "Myc", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)
#bb_gene_umap(filter_cds(mouse_cds_list[[4]], cells = bb_cellmeta(mouse_cds_list[[4]]) |> filter(k10_assignment =="B")),  gene_or_genes = "Cd93", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)

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
#S5_LN_Bclust_Top50markers <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables//Supp5_PRMT5xTCL1_vs_TCL1_LN/S5_LN_Bclust_Top50markers.csv")

#####S5C -Heatmap -RT Clust2&7 and TCL1 Clust4
S5C1<-bb_var_umap(mouse_cds_list[[4]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")+labs(x="UMAP 1", y="UMAP 2")
S5A1/S5C1  
markers <- S5_LN_Bclust_Top50markers |> 
  filter(cell_group %in% c("LN.2", "LN.7","LN.4")) |> 
  pull(gene_short_name)

S5_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[4]], 
                                          cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                            filter(kmeans_10_harmonized %in% c("LN.2","LN.4","LN.5","LN.7")),
                                          genes = bb_rowmeta(mouse_cds_list[[4]]) |> 
                                            filter(gene_short_name %in% markers)), 
                         cell_group_df = bb_cellmeta(mouse_cds_list[[4]]) |> 
                           select(cell_id, kmeans_10_harmonized)) |> 
  t() |> 
  scale() |> 
  t()

rownames(S5_mat) <- tibble(feature_id = rownames(S5_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[4]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

# Supp5_heatmap_mat <- as.data.frame(S5_mat)|> mutate(GeneID = rownames(S5_mat)) |> select(GeneID, everything())
#  write_csv(Supp5_heatmap_mat, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Supp5_PRMT5xTCL1_vs_TCL1_LN", "Supp5_heatmap_mat.csv"))
#  write_csv(Supp5_heatmap_mat, file = file.path("~/network/T/Labs/EHL/Rosa/EXPERIMENTS/PRMT5/Manuscripts/Drafts/Drafts_2022/Source Data", "Supp5_heatmap_mat.csv"))

S5_colfun = circlize::colorRamp2(breaks = c(min(S5_mat),
                                              0,
                                              max(S5_mat)),
                                   colors = heatmap_3_colors)

F3_F5_analysis_highlights <- str_to_title(c("CCR7", "BIRC5", "JUNB", "ATF4", "CD69", "TOP2A", "HMGB1", "HMGB2", "CD83", "TUBB5", "TUBA1B", "BLNK", "HSPA5", 
                               "CCDC34", "UBE2C", "IL2RG", "CKS1B", "SMC2", "STMN1", "TYMS", "UBC", "GAPDH", "APOE", "XRCC1",
                               "CD79A","ZFP36","PPP1R15A", "LDHA", "MKI67","CDK1", "GPX4", "BTG1", "ITM2B", "XBP1", "REL", "AURKB",
                               "EIF5A","C1QBP", "CD79B", "TK1","NFKBIA","HSP90AA1", "CENPM", "NR4A1","PFDN5","EIF4A2", "PIM1", "NPM1", 
                               "MARCKSL1", "ACTB", "NAPSA","CD37", "CD19", "CORO1A", "KLF2", "HNRNPK","RAC2", "MS4A1"))
#Additional lymphoma associated genes: DisGeNet.org - Lymphoma; CUI: C0024299
lymphoma_genes<- readxl::read_excel("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/disgenet.org_LymphomaGenes_C0024299_disease_gda_summary.xlsx")[3]
lymphoma_genes$Gene <- str_to_title(lymphoma_genes$Gene)
filt<- S5_LN_Bclust_Top50markers |> filter(S5_LN_Bclust_Top50markers$gene_short_name %in% lymphoma_genes$Gene)
lymphoma_gois <- filter(filt, cell_group %in% c("LN.4"))[["gene_short_name"]]
#Human top markers
RTLN_tm <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap_IPA_Tables/Fig1 human RT data/leiden clustering/LN_B_leiden_Top50_tm.csv")
RTLN_tm$gene_short_name <- str_to_title(RTLN_tm$gene_short_name)
RTLN_tm <- filter(RTLN_tm, cell_group == c('3','11'))
human_overlap<- S5_LN_Bclust_Top50markers |> filter(S5_LN_Bclust_Top50markers$gene_short_name %in% RTLN_tm$gene_short_name)
human_overlap <- filter(human_overlap, cell_group %in% c("LN.2","LN.7"))[["gene_short_name"]]
Supp5 <- c("Ms4a1","Coro1a","Cd52","Cdk1","Ccna2","Cd37","Rac2","Fau","Rhoa","Cd79b","Cd74","Cd83","Birc5","Hmgb1","Tubb5","Top2a","Hmgb2","Tuba1b")

highlights <- unique(c(F3_F5_analysis_highlights, human_overlap, Supp5))

S5_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(S5_mat) %in% highlights),
  labels = rownames(S5_mat)[rownames(S5_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2.5
))

S5C2 <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(S5_mat, 
                                                   col = S5_colfun, 
                                                   name = "Expression", 
                                                   show_row_names = F, 
                                                   right_annotation = S5_anno,
                                                   row_dend_width = unit(4, "mm"),
                                                   column_dend_height = unit(4, "mm"),
                                                   heatmap_legend_param = list(legend_direction = "vertical",
                                                                               #legend_width = unit(1, "mm"),
                                                                               title_position = "lefttop-rot", 
                                                                               title_gp = gpar(fontsize = 6)
                                                   ))))

S5C <- S5C1 / S5C2 + plot_layout(heights = c(1, 2))
S5C
ggsave("S5C.pdf", path = T_Figs)

# F5_k10_Top50_tm |> filter(cell_group == "5.6") |> arrange(desc(marker_score))
# fig5_kmeans_10_tm_Top100 |> filter(cell_group == "5.6") |> arrange(desc(mean_expression))

S5_genebubble <- 
  bb_genebubbles(
    mouse_cds_list[[4]],
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
    cell_grouping = c("kmeans_10_harmonized", "k10_assignment"),
    return_value = "data"
  ) |> 
  ggplot(mapping = aes(x = kmeans_10_harmonized, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~k10_assignment, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
ggsave("S5_genebubble.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Supp5", height = 4.23, width = 5.27)