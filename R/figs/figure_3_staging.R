# make some clusters assignments

colData(mouse_cds_list[[1]])$kmeans_10_clusters = factor(colData(mouse_cds_list[[1]])$kmeans_10_clusters, levels = as.character(1:10))


colData(mouse_cds_list[[1]])$k_10_assignment = recode(colData(mouse_cds_list[[1]])$kmeans_10_clusters, 
                                                    "1" = "Cd8+ T",
                                                    "2" = "B",
                                                    "3" = "B",
                                                    "4" = "B",
                                                    "5" = "Cd4+ T",
                                                    "6" = "B",
                                                    "7" = "Neutrophils",
                                                    "8" = "B",
                                                    "9" = "Monocytes",
                                                    "10" = "Low Quality",
                                                    )

# make logical values for CD19+CD5+ cells
mat <- monocle3::exprs(mouse_cds_list[[1]])

cd19_tbl <- colnames(mat[ ,mat["ENSMUSG00000030724", ] > 0]) |> as_tibble() |> mutate(Cd19_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[1]] <- bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = cd19_tbl)

cd5_tbl <- colnames(mat[ ,mat["ENSMUSG00000024669", ] > 0]) |> as_tibble() |> mutate(Cd5_pos = TRUE) |> rename(cell_id = value)
mouse_cds_list[[1]] <- bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = cd5_tbl)


colData(mouse_cds_list[[1]])$cd19_cd5_pos <- colData(mouse_cds_list[[1]])$Cd19_pos & colData(mouse_cds_list[[1]])$Cd5_pos  

# get top markers for kmeans clusters
fig3_kmeans_10_tm <- monocle3::top_markers(mouse_cds_list[[1]], group_cells_by = "kmeans_10_clusters")

dir.create("data")
save(mouse_cds_list, file = "data/mouse_cds_list.rda", compress = "bzip2")
save(fig3_kmeans_10_tm, file = "data/fig3_kmeans_10_tm.rda", compress = "bzip2")

mouse_cds_list[[1]] |> bb_cellmeta()
mouse_cds_list[[1]] |> bb_cellmeta() |> group_by(sample_id, genotype, tissue) |> summarise()
mouse_cds_list[[2]] |> bb_cellmeta() |> group_by(sample_id, genotype, tissue) |> summarise()
mouse_cds_list[[3]] |> bb_cellmeta() |> group_by(sample_id, genotype, tissue) |> summarise()
mouse_cds_list[[4]] |> bb_cellmeta() |> group_by(sample_id, genotype, tissue) |> summarise()
bb_var_umap(mouse_cds_list[[2]], "sample_id", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
bb_var_umap(mouse_cds_list[[1]], "genotype")
test <- bb_align(mouse_cds_list[[1]], "sample")

bb_var_umap(test, "sample_id")
bb_var_umap(test, "genotype", value_to_highlight = "PRMT5")
bb_var_umap(test, "genotype", value_to_highlight = "TCL1")
bb_var_umap(test, "genotype", facet_by = "value")
bb_gene_umap(test, "Cd19")

# figure 3A
bb_var_umap(mouse_cds_list[[1]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")

bb_var_umap(mouse_cds_list[[1]], "k_10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F)

# supplemental figure
bb_genebubbles(mouse_cds_list[[1]], genes = c("Cd19", "Cd3d", "Cd4", "Cd14", "Itgam", "Cd8a", "Cd177", "Pdcd1", "Foxp3"), cell_grouping = "kmeans_10_clusters")


bb_var_umap(mouse_cds_list[[1]], "kmeans_10_clusters", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T)


# figure 3B
bb_var_umap(mouse_cds_list[[1]], "cd19_cd5_pos", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", value_to_highlight = TRUE)

# figure 3D

markers <- kmeans_10_tm |> 
  filter(cell_group %in% c("2", "3", "4", "6", "8")) |> 
  pull(gene_short_name)

fig3_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[1]], 
                              cells = bb_cellmeta(mouse_cds_list[[1]]) |> filter(kmeans_10_clusters %in% c("2", "3", "4", "6", "8")),
                              genes = bb_rowmeta(mouse_cds_list[[1]]) |> filter(gene_short_name %in% markers)), cell_group_df = bb_cellmeta(mouse_cds_list[[1]]) |> select(cell_id, kmeans_10_clusters)) |> 
  t() |> 
  scale() |> 
  t()
rownames(fig3_mat) <- tibble(feature_id = rownames(fig3_mat)) |> left_join(bb_rowmeta(mouse_cds_list[[1]]) |> select(feature_id, gene_short_name)) |> pull(gene_short_name)

fig3_colfun = circlize::colorRamp2(breaks = c(min(fig3_mat),
                                              0,
                                              max(fig3_mat)),
                                   colors = heatmap_3_colors)

highlights <- c("Wnt10a")
fig3_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig3_mat) %in% highlights),
  labels = rownames(fig3_mat)[rownames(fig3_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2
))

ComplexHeatmap::Heatmap(fig3_mat, 
                        col = fig3_colfun, 
                        name = "Expression", show_row_names = F, right_annotation = fig3_anno)
