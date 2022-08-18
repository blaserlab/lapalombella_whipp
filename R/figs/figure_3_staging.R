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
