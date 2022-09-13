#Create Disease_tissue column in cds
colData(cds_main)$disease_tissue <-
  paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
# Reordering group factor levels - For following graphing order
cds_main$disease_tissue <- factor(cds_main$disease_tissue,
                                  levels = c("CLL PBMC", "RT PBMC", "RT LN"))

#partition assignment
# bb_genebubbles(
#   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
#   genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8A", "CD4"), cell_grouping = "partition") + labs(x = "Partition Cluster", y = NULL) 
# colData(cds_main)$partition_assignment_1 <-
#   recode(
#     colData(cds_main)$partition,
#     "1" = "B",
#     "2" = "B",
#     "3" = "T",
#     "4" = "T",
#     "5" = "Mono",
#     "6" = "B",
#     "7" = "B"
#   )

#leiden assignment
# bb_genebubbles(
#   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
#   genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8A", "CD4"), cell_grouping = "leiden") + labs(x = "leiden Cluster", y = NULL)
colData(cds_main)$leiden_assignment_1 <-
  recode(
    colData(cds_main)$leiden,
    "1" = "B",
    "2" = "B",
    "3" = "B",
    "4" = "T",
    "5" = "B",
    "6" = "B",
    "7" = "T",
    "8" = "T",
    "9" = "B",
    "10" = "T",
    "11" = "B",
    "12" = "T",
    "13" = "Mono",
    "14" = "T",
    "15" = "B",
    "16" = "T",
    "17" = "B",
    "18" = "B"
  )

colData(cds_main)$f1d_lbl <- paste0(colData(cds_main)$disease_tissue, "-", recode(colData(cds_main)$patient, "pt_2712" = "Pt 1", "pt_1245" = "Pt 2"))
colData(cds_main)$patient <- recode(colData(cds_main)$patient, "pt_2712" = "Pt 1", "pt_1245" = "Pt 2")


#Compare to Lit data - Nadeu et al
nadeu_11b <-
  readxl::read_excel(
    "~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx",
    sheet = "Supplementary Table 11b",
    skip = 5,
    col_names = c(
      "feature_id",
      "gene_short_name",
      "mean",
      "l2fc",
      "se",
      "p",
      "padj",
      "direction"
    )
  )

cds_main <- nadeu_11b |> 
  filter(direction == "Up") |> 
  filter(padj < 0.05) |> 
  mutate(feature_id = str_remove(feature_id, "\\..*")) |> 
  mutate(nadeu_RT_gene = TRUE) |> 
  bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)

cds_main <- nadeu_11b |>
 filter(direction == "Down") |>
 filter(padj < 0.05) |>
 mutate(feature_id = str_remove(feature_id, "\\..*")) |>
 mutate(nadeu_CLL_gene = TRUE) |>
 bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)
#rowData(cds_main)

ln_cds<-cds_main[,colData(cds_main)$disease_tissue == "RT LN" &
                   colData(cds_main)$leiden_assignment_1 == "B"]


#supplemental figs - to show clustering pattern is due to clustering of both pts
#S1A <- bb_var_umap(cds_main, "patient", facet_by = "value", legend_pos = "none", foreground_alpha = 0.2)
#V(D)J

# possible alt S1A
 S1A <-bb_var_umap(cds_main, "patient", value_to_highlight = "Pt 1", legend_pos = "none", plot_title = "Pt 1", palette = "#EF8A62")+
 bb_var_umap(cds_main, "patient", value_to_highlight = "Pt 2", legend_pos = "none", plot_title = "Pt 2", palette = "#67A9CF")
 
 S1B <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(patient == "Pt 1")), "clonotype_id", facet_by = "disease_tissue")
 S1B<- S1B + theme(legend.position = "bottom")

#S1C <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "Pt 1")), "partition")
S1C <- bb_var_umap(cds_main, "leiden", overwrite_labels = T)
#S1C<- S1C + theme(legend.position = "bottom")

#possible alt S1B -----
#S1B <- bb_var_umap(cds_main, var = "partition", facet_by = "pt_recode")

S1D <- bb_genebubbles(
  cds_main,
  genes = c("CD14","CD4","CD8A", "CD3E", "CD79A", "MS4A1", "CD19"
 ),
  cell_grouping = c("leiden", "leiden_assignment_1"),
  return_value = "data"
) |> 
  ggplot(mapping = aes(x = leiden, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~leiden_assignment_1, scales = "free_x", ) +
  theme_minimal_grid(font_size = 8) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")

#Supp Fig - Nadeu et al RT UP aggregate gene expression mapping 
#Nadeu et al RT UP aggregate gene expression mapping 
# bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
# ) + 
#   facet_grid(row = vars(patient), col = (vars(disease_tissue)))

S1EP1 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
) + 
  facet_wrap(~disease_tissue)+
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*RT Gene*<br>Expression") +
  theme(legend.title = ggtext::element_markdown())
S1EP2 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_CLL_gene)
) +
  facet_wrap(~disease_tissue)+
  theme(strip.text = element_blank()) +
   theme(
  #   axis.line.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   axis.text.y = element_blank(),
     axis.title.x = element_blank(),
     axis.title.y = element_blank()
   ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*CLL Gene*<br>Expression") +
  theme(legend.title = ggtext::element_markdown())
#S1D <- S1DP1/S1DP2

S1E <- as_ggplot(grid.arrange(patchworkGrob(S1EP1/S1EP2), left = S1EP1$labels$y, bottom = textGrob(S1EP1$labels$x, hjust = 0.85, vjust = -1)))

#ggsave("suppfig_1A.pdf", path = figures_out, width = 4.95, height = 2.45)

# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "pt_2712")), "partition", facet_by = "disease_tissue") 

# LN_B_clst2 <-cds_main[,colData(cds_main)$partition == "2"]
# colData(LN_B_clst2)

#Top Markers Partition Cluster 2
# Fig1_tm <- monocle3::top_markers(cds_main,
#                      group_cells_by = "partition", genes_to_test_per_group = 300, cores = 10)
# write_csv(Fig1_tm, file = file.path(WalkerTables, "Fig1_300tm.csv"))


# main Fig 1

F1AP1 <-
  bb_var_umap(
    filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
    "leiden_assignment_1", overwrite_labels = T,
    facet_by = "disease_tissue"
  ) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1AP2 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
                     "log_local_n",
                     facet_by = "disease_tissue") +
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "Log<sub>10</sub><br>Cells") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.4,"cm"))+
  theme(legend.title = ggtext::element_markdown())
F1AP3 <- bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
                      gene_or_genes = "PRMT5") +
  facet_wrap( ~ disease_tissue) +
  theme(strip.text = element_blank()) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*PRMT5*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.4,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
  
F1A <- F1AP1/F1AP2/F1AP3

# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "partition_assignment_1"
# ) + facet_wrap(~disease_tissue)
# 
# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_1245")), "partition_assignment_1"
# ) + facet_wrap(~disease_tissue)
# 
#  bb_gene_umap(
#    filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), "PRMT5"
#  ) + facet_grid(row = vars(patient), col = (vars(disease_tissue)))

#Gene dot plot
# F1A_optional_dotplot <-bb_gene_dotplot(
#   cds_main[, colData(cds_main)$patient == "Pt 1" &
#              colData(cds_main)$clonotype_id %in% "clonotype1" &
#              colData(cds_main)$leiden_assignment_1 %in% "B"],
#   markers = c("PRMT5", "MYC", "MKI67"),
#   group_cells_by = "disease_tissue",
#   group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
#   colorscale_name = "Expression",
#   sizescale_name = "Proportion\nExpressing",
# ) + labs(x = NULL, y = NULL)
# ggsave("F1A_optional_dotplot.pdf", path = T_Figs, width = 3.57, height = 2.98)

#Fig1D
F1D0 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                    "leiden_assignment_1", facet_by = "f1d_lbl", overwrite_labels = T) +
    theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(panel.background = element_rect(color = "black"))+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1D1<- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                   "density", facet_by = "f1d_lbl") + 
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*Density*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D2<- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("PRMT5") 
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*PRMT5*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1D3<- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("MYC")
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*MYC*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1D4 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("MKI67")
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    # axis.line.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*MKI67*") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# F1D5<- bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
#   gene_or_genes = c("BIRC5")
# ) + facet_wrap(~f1d_lbl)+
#   theme(strip.text = element_blank()) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )+
#   theme(panel.background = element_rect(color = "black")) +
#   labs(color = "*BIRC5*") +
#   theme(legend.title = ggtext::element_markdown())

F1D <-
  as_ggplot(grid.arrange(
    patchworkGrob(F1D0 / F1D1 / F1D2 / F1D3 / F1D4),
    left = textGrob(F1D3$labels$y, rot=90, vjust = 1.5, hjust=0.35),
    bottom = textGrob(F1D3$labels$x, hjust = 0.8, vjust = -0.5))
  )

#############################################################################################################################

#Fig 1E: heatmap (leiden clusters)

#Select clusters of interest
# bb_genebubbles(
#   cds_main,
#   genes = c("CD14","CD4","CD8A", "CD3E", "CD79A", "MS4A1", "CD19"
#   ),
#   cell_grouping = c("leiden", "leiden_assignment_1"),
#   return_value = "data"
# ) |> 
#   ggplot(mapping = aes(x = leiden, 
#                        y = gene_short_name, 
#                        color = expression,
#                        size = proportion)) +
#   geom_point() +
#   scale_size_area() +
#   scale_color_viridis_c() +
#   facet_wrap(~leiden_assignment_1, scales = "free_x", ) +
#   theme_minimal_grid(font_size = 8) +
#   theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
#   theme(axis.text.y = element_text(face = "italic")) +
#   labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")

F1E1 <- bb_var_umap(ln_cds, "leiden", overwrite_labels = T
  ) + labs(title = "RT LN")
F1E1

#p3 <- bb_var_umap(ln_cds, "leiden", overwrite_labels = T)+labs(title = "RT LN")

#Heatmap: 
###filter for clusters of interst: c('3','11','6','2','5','9','1')

#Top markers should be used to look at top 50genes in clusters of interest.
LN_B_leiden_Top50_tm <- monocle3::top_markers(ln_cds, group_cells_by = "leiden", genes_to_test_per_group = 50, cores = 12)
#T_tables_out <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Fig1 human RT data/leiden clustering/LN_B_leiden_Top50_tm.csv"
#write_csv(LN_B_leiden_Top50_tm, file = file.path(T_tables_out, "LN_B_leiden_Top50_tm.csv"))
#LN_B_leiden_Top50_tm <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Fig1 human RT data/leiden clustering/LN_B_leiden_Top50_tm.csv")
#bb_var_umap(ln_cds,"leiden", overwrite_labels = T) /p3

#All leiden_assignment_1 B cell clusters
markers <- LN_B_leiden_Top50_tm |> filter(cell_group %in% c("3","11","6","5","1","2","9")) |> pull(gene_short_name)


f1_mat <- bb_aggregate(obj = filter_cds(ln_cds,
                                        cells = bb_cellmeta(ln_cds) |>
                                          filter(leiden %in% c('3','11','6','2','5','9','1')),
                                        genes = bb_rowmeta(ln_cds) |>
                                          filter(gene_short_name %in% markers)),
                       cell_group_df = bb_cellmeta(ln_cds) |>
                         select(cell_id, leiden)) |>
  t() |>
  scale() |>
  t()

rownames(f1_mat) <- tibble(feature_id = rownames(f1_mat)) |>
  left_join(bb_rowmeta(ln_cds) |>
              select(feature_id, gene_short_name)) |>
  pull(gene_short_name)
#f1_mat
f1_colfun = circlize::colorRamp2(breaks = c(min(f1_mat),
                                            0,
                                            max(f1_mat)),
                                 colors = heatmap_3_colors)

#Mouse markers discussed in prev manuscript version: c("Ccr7", "Cdk4", "Cxcr5", "Birc5", "Il4","Npm1","Jun","Junb", "Fos","Fosb","Atf3","Atf4","Myc","Cd69","Il10", "Top2a", "Hmgb1", "Hmgb2", "Cd83","Ube2a", "Tubb5", "Tuba1b")
#Human markers discussed "PRMT5", "MKI67", "TOP2A","PCNA","CALM2","CALM3","HMGB1","HMGB2", "UBE2C", "UBE2S", "TUBA1A","TUBA1B", "TUBA1C", "BCL21A", "IL4l1", "CD83", "TCL1A", "CCL3", "CCL4"
#######human markers in top 50/clust = MKI67, TOP2A, HMGB2, UBE2C, UBE2S, TUBB4B(not TUBA1A or TUBA1C), TCL1A, CD83, TCL1A 
###added to list due to presence in mouse data CDK4, BIRC5,NPM1,Junb (up)
###other interesting genes CDK1, RRM2, CCNA2, CCNB1, CCNB2, CDC20,	CDCA3, CDCA5, CDCA8, MS4A1, (KIF proteins)

highlights <- c("MKI67", "TOP2A","HMGB2", "UBE2C", "UBE2S","TUBA1B", "TUBA1C", "TUBB", "TCL1A", "BIRC5", "NPM1", "CDK1", "JUNB")#, "AURKB", "PLK1", "CDK4","TUBB4B","CD83")

#Additional lymphoma associated genes
# lymphoma_genes<- readxl::read_excel("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/disgenet.org_LymphomaGenes_C0024299_disease_gda_summary.xlsx")[3]
# filt<- LN_B_leiden_Top50_tm |> filter(LN_B_leiden_Top50_tm$gene_short_name %in% lymphoma_genes$Gene)
# tm_gois <- filter(filt, cell_group %in% c('11','3'))[["gene_short_name"]]
# highlights_plus <- unique(c(highlights, tm_gois))

#BCR associated genes
# highlights <- #BCR gene list
#   readxl::read_excel(
#     "~/network/T/Labs/EHL/Rosa/Ethan/10X/Deepa_CLL_study/BCR genelist.xlsx",
#     skip = 1,
#   ) |> as.vector()

  
fig1_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(f1_mat) %in% highlights),
  labels = rownames(f1_mat)[rownames(f1_mat) %in% highlights],
  labels_gp = gpar(fontsize = 5),
  padding = unit(2.8, "mm")
  # link_width = unit(5, "mm"),
  # link_height = link_width,
  #link_gp = unit(10, "mm"),
  # extend = unit(0, "mm)
))


F1E2 <- grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(f1_mat,
                        col = f1_colfun,
                        name = "Expression", 
                        show_row_names = F, 
                        right_annotation = fig1_anno,
                        row_dend_width = unit(4, "mm"),
                        column_dend_height = unit(2, "mm"),
                        heatmap_legend_param = list(legend_direction = "vertical",
                                                    #legend_width = unit(1, "mm"),
                                                    title_position = "topleft-rot", 
                                                    title_gp = gpar(fontsize = 6)
                                               ))))
#F1E2 <- grid.grabExpr(draw(F1E2, heatmap_legend_side = "bottom"))


F1E <- F1E1 / F1E2 +
  plot_layout(heights = c(1, 2.5))
F1E
as_ggplot(F1E2)
