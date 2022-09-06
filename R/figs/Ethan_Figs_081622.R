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
colData(cds_main)$partition_assignment_1 <-
  recode(
    colData(cds_main)$partition,
    "1" = "B",
    "2" = "B",
    "3" = "T",
    "4" = "T",
    "5" = "Mono",
    "6" = "B",
    "7" = "B"
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
                   colData(cds_main)$clonotype_id %in% "clonotype1"]


#supplemental figs - to show clustering pattern is due to clustering of both pts
S1A <- bb_var_umap(cds_main, "patient", facet_by = "value", legend_pos = "none", foreground_alpha = 0.2)
#V(D)J

S1B <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(patient == "Pt 1")), "clonotype_id", facet_by = "disease_tissue")

# possible alt S1A
# bb_var_umap(cds_main, "patient", value_to_highlight = "pt_2712", legend_pos = "none", plot_title = "pt_2712") +
# bb_var_umap(cds_main, "patient", value_to_highlight = "pt_1245", legend_pos = "none", plot_title = "pt_1245")

S1C <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "Pt 1")), "partition")

#possible alt S1B -----
#S1B <- bb_var_umap(cds_main, var = "partition", facet_by = "pt_recode")

S1D <- bb_genebubbles(
  cds_main,
  genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19"
 ),
  cell_grouping = c("partition", "partition_assignment_1"),
  return_value = "data"
) |> 
  ggplot(mapping = aes(x = partition, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~partition_assignment_1, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
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

S1E <- as_ggplot(grid.arrange(patchworkGrob(S1EP1/S1EP2), left = S1EP1$labels$y, bottom = textGrob(S1EP1$labels$x, hjust = 0.9)))

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
    "partition_assignment_1",
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
  labs(fill = "Cluster")+ 
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
F1A_optional_dotplot <-bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "Pt 1" &
             colData(cds_main)$clonotype_id %in% "clonotype1" &
             colData(cds_main)$partition_assignment_1 %in% "B"],
  markers = c("PRMT5", "MYC", "MKI67"),
  group_cells_by = "disease_tissue",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)

# bb_genebubbles(
#   obj = filter_cds(
#     cds_main,
#     cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")
#     ),
#   genes = c("MYC", "PRMT5", "MKI67"),
#   cell_grouping = "partition_assignment_1") + facet_wrap(~disease_tissue)


#Fig1D
F1D0 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                    "partition_assignment_1", facet_by = "f1d_lbl") +
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
    bottom = textGrob(F1D3$labels$x, hjust = 0.8, vjust = -1.2))
  )


# Fig1E 
# COMBAK finish heatmap

# TODO determine number of genes to look into and 
LN_Top50_tm <-monocle3::top_markers(ln_cds, group_cells_by = "partition", genes_to_test_per_group = 50, cores = 10)
#write_csv(LN_Top100_tm, file = file.path(tables_out, "fig1_heatmap_LNpart_Top100_tm.csv"))
#LN_Top500_tm <-monocle3::top_markers(ln_cds, group_cells_by = "partition", genes_to_test_per_group = 500, cores = 10)
#write_csv(LN_Top500_tm, file = file.path(WalkerTables, "LN_Top500_tm.csv"))
#LN_Top2000_tm <-monocle3::top_markers(ln_cds, group_cells_by = "partition", genes_to_test_per_group = 2000, cores = 12)
#write_csv(LN_Top2000_tm, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Fig1 human RT data", "LN_Top2000_tm.csv"))

# TODO determine markers discussed in human RT and those discussed in mouse data for heatmap annotation + any additions
#FIXME bb_gene_umap(cds_main[,colData(cds_main)$disease_tissue == "RT LN"],gene_or_genes = "PRMT5")

#bb_var_umap(cds_main, "partition", facet_by = "patient")
#bb_var_umap(cds_main, "partition_assignment_1", facet_by = "patient")
markers <- LN_Top50_tm |> filter(cell_group %in% c("1","2","6")) |> pull(gene_short_name)

f1_mat <- bb_aggregate(obj = filter_cds(ln_cds, 
                                          cells = bb_cellmeta(ln_cds) |> 
                                            filter(partition %in% c("1","2","6")),
                                          genes = bb_rowmeta(ln_cds) |> 
                                            filter(gene_short_name %in% markers)), 
                         cell_group_df = bb_cellmeta(ln_cds) |> 
                           select(cell_id, partition)) |> 
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

highlights <- as.list(LN_Top50_tm |> filter(cell_group %in% c("2")) |>
            head(sort(LN_leiden_Top100_tm$marker_score, decreasing =
                        TRUE), n = 20) |> pull(gene_short_name))
#TODO fill in remaining genes from paper
#paper_highlights <- list(c("Myc","MKI67","PRMT5",))

#highlights <- c("Rel","PRMT5","Myc","CD83", "TUBB5", "TUBA1B","WNT10A", "CCR7", "CDK1", "BIRC5","JUNB","ATF4","CD69","NFKBIA","NFKBID","REL","UBC","CCNA2", "BLNK")

#NOTE: Mouse markers discussed in prev manuscript version: c("Ccr7", "Cdk4", "Cxcr5", "Birc5", "Il4","Npm1","Jun","Junb", "Fos","Fosb","Atf3","Atf4","Myc","Cd69","Il10", "Top2a", "Hmgb1", "Hmgb2", "Cd83","Ube2a", "Tubb5", "Tuba1b")

fig1_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(f1_mat) %in% highlights),
  labels = rownames(f1_mat)[rownames(f1_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 0
))

F1E<-ComplexHeatmap::Heatmap(f1_mat, 
                        col = f1_colfun, 
                        name = "Expression", show_row_names = F, right_annotation = fig1_anno)



# QUESTION use leiden clustering on heatmap and all previous analysises
#Leiden heatmap clusters
#TEMP LN_leiden_Top100_tm <- monocle3::top_markers(ln_cds, group_cells_by = "leiden", genes_to_test_per_group = 100, cores = 10)

#bb_var_umap(ln_cds,"partition_assignment_1", overwrite_labels = T)
#bb_var_umap(ln_cds,"leiden", overwrite_labels = T)

#QUESTION should markers only come from clusters 3&11?
#TEMP markers <- LN_leiden_Top100_tm |> filter(cell_group %in% c("3","11","6","5","1","9")) |> pull(gene_short_name)


#TEMP f1_mat <- bb_aggregate(obj = filter_cds(ln_cds, 
#                                         cells = bb_cellmeta(ln_cds) |> 
#                                           filter(leiden %in% c("3","11","6","5","1","9")),
#                                         genes = bb_rowmeta(ln_cds) |> 
#                                           filter(gene_short_name %in% markers)), 
#                        cell_group_df = bb_cellmeta(ln_cds) |> 
#                          select(cell_id, leiden)) |> 
#   t() |> 
#   scale() |> 
#   t()
# 
# rownames(f1_mat) <- tibble(feature_id = rownames(f1_mat)) |> 
#   left_join(bb_rowmeta(ln_cds) |> 
#               select(feature_id, gene_short_name)) |> 
#   pull(gene_short_name)
#f1_mat
#TEMP f1_colfun = circlize::colorRamp2(breaks = c(min(f1_mat),
#                                             0,
#                                             max(f1_mat)),
#                                  colors = heatmap_3_colors)

# QUESTION determine genes to highlight &
#   how to come up with that list (nadeu et al - lymphma / leukemia gens &CLL / RT gense, paper mentioned genes - both mus&homo)

# TEMP highlights <-
#   as.list(LN_leiden_Top100_tm |> filter(cell_group %in% c("3")) |>
#             head(sort(LN_leiden_Top100_tm$marker_score, decreasing =
#                         TRUE), n = 20) |> pull(gene_short_name)) 
                              #|>unique() ## highlights[[!duplicated(highlights$gene_short_name),]]


#highlights <- c("Rel","PRMT5","Myc","CD83", "TUBB5", "TUBA1B","WNT10A", "CCR7", "CDK1", "BIRC5","JUNB","ATF4","CD69","NFKBIA","NFKBID","REL","UBC","CCNA2", "BLNK")
#Mouse markers discussed in prev manuscript version: c("Ccr7", "Cdk4", "Cxcr5", "Birc5", "Il4","Npm1","Jun","Junb", "Fos","Fosb","Atf3","Atf4","Myc","Cd69","Il10", "Top2a", "Hmgb1", "Hmgb2", "Cd83","Ube2a", "Tubb5", "Tuba1b")

#TEMP fig1_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
#   at = which(rownames(f1_mat) %in% highlights),
#   labels = rownames(f1_mat)[rownames(f1_mat) %in% highlights],
#   labels_gp = gpar(fontsize = 8),
#   padding = 0
# ))

# TEMP ComplexHeatmap::Heatmap(f1_mat, 
#                         col = f1_colfun, 
#                         name = "Expression", show_row_names = F, right_annotation = fig1_anno)
