#Create Disease_tissue column in cds
colData(cds_main)$disease_tissue <- paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
# Reordering group factor levels - For following graphing order
cds_main$disease_tissue <- factor(cds_main$disease_tissue,      
                                  levels = c("CLL PBMC", "RT PBMC", "RT LN"))
#partition assignment
bb_genebubbles(
  obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
  genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8"), cell_grouping = "partition") + labs(x = "Partition Cluster", y = NULL) 
colData(cds_main)$partition_assignment_1 <- recode(colData(cds_main)$partition, "1" = "B", "2" = "B", "3" = "T", "4" = "T", "5" = "Mono", "6" = "B", "7" = "B")

#supp Fig 1A
suppfig_1A <- bb_genebubbles(
  cds_main,
  genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8"
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
#suppfig_1A
#ggsave("suppfig_1A.pdf", path = figures_out, width = 4.95, height = 2.45)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)), "partition_assignment_1") + facet_grid(~patient)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)), "partition") +
  facet_wrap(~patient)

# LN_B_clst2 <-cds_main[,colData(cds_main)$partition == "2"]
# colData(LN_B_clst2)

#Top Markers Partition Cluster 2
# Fig1_tm <- monocle3::top_markers(cds_main,
#                      group_cells_by = "partition", genes_to_test_per_group = 300, cores = 10)
# write_csv(Fig1_tm, file = file.path(WalkerTables, "Fig1_300tm.csv"))


#Fig1A
bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "pt_2712")), "partition_assignment_1") +
  facet_wrap(~disease_tissue)

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient = "pt_2712")), gene_or_genes = "PRMT5"
) + 
  facet_wrap(~disease_tissue)

#Compare to Lit data - Nadeu et al
nadeu_11b <- readxl::read_excel("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx", sheet = "Supplementary Table 11b", skip = 5, col_names = c("feature_id", "gene_short_name", "mean", "l2fc", "se", "p", "padj", "direction"))
view(nadeu_11b)

cds_main <- nadeu_11b |> 
  filter(direction == "Up") |> 
  filter(padj < 0.05) |> 
  mutate(feature_id = str_remove(feature_id, "\\..*")) |> 
  mutate(nadeu_RT_gene = TRUE) |> 
  bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)

#Density Mapping. 
##Can sub "log_local_n" for "density"
bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
  "log_local_n") +
  facet_grid(row = vars(patient),col = vars(disease_tissue))

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)), "partition_assignment_1") +
  facet_grid(row = vars(patient), col = (vars(disease_tissue)))

<<<<<<< Updated upstream
#Supp Fig - Nadeu et al RT UP aggregate gene expression mapping 
=======
#Nadeu et al RT UP aggregate gene expression mapping 
>>>>>>> Stashed changes
bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
) + 
  facet_grid(row = vars(patient), col = (vars(disease_tissue)))

<<<<<<< Updated upstream
####Supp Fig - Plot Nadeu et al CLL genes/RT downreg genes
=======
####Plot Nadeu et al CLL genes/RT downreg genes
>>>>>>> Stashed changes
#rm(cds_main)
#colData(cds_main)$disease_tissue <- paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
#colData(cds_main)$partition_assignment_1 <- recode(colData(cds_main)$partition, "1" = "B", "2" = "B", "3" = "T", "4" = "T", "5" = "Mono", "6" = "B", "7" = "B")

#cds_main <- nadeu_11b |> 
#  filter(direction == "Down") |> 
#  filter(padj < 0.05) |> 
#  mutate(feature_id = str_remove(feature_id, "\\..*")) |> 
#  mutate(nadeu_RT_gene = TRUE) |> 
#  bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)
# 
# bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
# ) + 
#   facet_grid(row = vars(patient), col = (vars(disease_tissue)))

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "partition_assignment_1"
<<<<<<< Updated upstream
) + facet_wrap(~disease_tissue)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_1245")), "partition_assignment_1"
) + facet_wrap(~disease_tissue)

 bb_gene_umap(
   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), "PRMT5"
 ) + facet_grid(row = vars(patient), col = (vars(disease_tissue)))

#Gene dot plot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
=======
  ) + facet_wrap(~disease_tissue)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_1245")), "partition_assignment_1"
 ) + facet_wrap(~disease_tissue)

# bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), "PRMT5"
# ) + facet_grid(row = vars(patient), col = (vars(disease_tissue)))

#Gene dot plot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_1245" &
>>>>>>> Stashed changes
             colData(cds_main)$clonotype_id %in% "clonotype1" &
             colData(cds_main)$partition_assignment_1 %in% "B"],
  markers = c("PRMT5", "MYC", "MKI67"),
  group_cells_by = "disease_tissue",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)

# bb_genebubbles(
#   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")),
#   genes = c("MYC", "PRMT5", "MKI67"), cell_grouping = "partition_assignment_1") + facet_grid(row = vars(patient), col = (vars(disease_tissue)))

<<<<<<< Updated upstream
=======


>>>>>>> Stashed changes
view(unique(cds_main$sample))

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B") |> filter(clonotype_id == "clonotype1") |> filter(disease_tissue == "RT LN")
<<<<<<< Updated upstream
  ), gene_or_genes = c("PRMT5", "MYC", "MKI67")
) + facet_grid(row = vars(patient), col = vars(disease_tissue))

#####Fig 1A dotplot -> scatter plot
#Still working on this...

#Fig1D
bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(disease_tissue == "RT LN")), "partition_assignment_1") +
  facet_wrap(~patient)

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(clonotype_id == "clonotype1") |> filter(partition_assignment_1 == "B")), gene_or_genes = c("PRMT5", "MYC", "MKI67")
) + facet_grid(row = vars(patient), col = vars(disease_tissue))

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(clonotype_id == "clonotype1") |> filter(partition_assignment_1 == "B")), gene_or_genes = "PRMT5"
) + facet_grid(row = vars(patient), col = vars(disease_tissue))

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(clonotype_id == "clonotype1") |> filter(partition_assignment_1 == "B")), gene_or_genes = "MYC"
) + facet_grid(row = vars(patient), col = vars(disease_tissue))

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(clonotype_id == "clonotype1") |> filter(partition_assignment_1 == "B")), gene_or_genes = "MYC"
) + facet_grid(row = vars(patient), col = vars(disease_tissue))), gene_or_genes = c("PRMT5", "MYC", "MKI67")
  ) + facet_grid(row = vars(patient), col = vars(disease_tissue))

# bb_gene_dotplot(
#   cds_main[, colData(cds_main)$patient == "pt_2712" &
#              colData(cds_main)$partition_assignment_1 == "B"],
#   markers = c("PRMT5"),
#   group_cells_by = "disease_tissue",
#   group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
#   colorscale_name = "Expression",
#   sizescale_name = "Proportion\nExpressing",
# ) + labs(x = NULL, y = NULL)

#### violin plot of PRMT5 expression
# bb_gene_violinplot(cds_main[, colData(cds_main)$patient == "pt_2712" &
#                               colData(cds_main)$clonotype_id %in% "clonotype1"], 
#                    variable = "sample", 
#                    genes_to_plot = "PRMT5")
#####Fig 1A dotplot -> scatter plot
library(Seurat)
FeatureScatter(object = cds_main, feature1 = 'CD9', feature2 = 'CD3E')
FeatureScatter(object = cds_main, feature1 = 'PRMT5', feature2 = 'MYC')


#clonotype1 <- cds_main[, colData(cds_main)$clonotype_id %in% "clonotype1"]
# bb_gene_umap(
  # filter_cds(clonotype1, cells = bb_cellmeta(clonotype1) |> filter(patient == "pt_1245")), gene_or_genes = c("PRMT5", "MYC", "MKI67")) + facet_wrap(~disease_tissue)

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(clonotype_id == "clonotype1") |> filter(partition_assignment_1 == "B")), gene_or_genes = c("PRMT5", "MYC", "MKI67")
  ) + facet_grid(row = vars(patient), col = vars(disease_tissue))

FeaturePlot(cds_main, features = "MS4A1")
