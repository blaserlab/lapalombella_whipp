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


# figure 3A
bb_cellmeta(mouse_cds_list[[1]])

bb_var_umap(mouse_cds_list[[1]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")

bb_var_umap(mouse_cds_list[[1]], "k_10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F)

# supplemental figure
bb_genebubbles(mouse_cds_list[[1]], genes = c("Cd19", "Cd3d", "Cd4", "Cd14", "Itgam", "Cd8a", "Cd177", "Pdcd1", "Foxp3"), cell_grouping = "kmeans_10_clusters")

kmeans_10_tm <- monocle3::top_markers(mouse_cds_list[[1]], group_cells_by = "kmeans_10_clusters")

bb_var_umap(mouse_cds_list[[1]], "kmeans_10_clusters", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T)


# figure 3B

mat <- monocle3::exprs(mouse_cds_list[[1]])
bb_rowmeta(mouse_cds_list[[1]]) |> filter(gene_short_name == "Cd5")
bb_cellmeta(mouse_cds_list[[1]]) |> glimpse()
cd19_tbl <- colnames(mat[ ,mat["ENSMUSG00000030724", ] > 0]) |> as_tibble() |> mutate(Cd19_pos = TRUE) |> rename(cell_id = value)

mouse_cds_list[[1]] <- bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = cd19_tbl)

cd5_tbl <- colnames(mat[ ,mat["ENSMUSG00000024669", ] > 0]) |> as_tibble() |> mutate(Cd5_pos = TRUE) |> rename(cell_id = value)

mouse_cds_list[[1]] <- bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = cd5_tbl)

c(TRUE, TRUE) & c(TRUE, TRUE)
TRUE & FALSE

c(1, 1, 1) == c(1, 2, 3)

colData(mouse_cds_list[[1]])$cd19_cd5_pos <- colData(mouse_cds_list[[1]])$Cd19_pos & colData(mouse_cds_list[[1]])$Cd5_pos  

bb_var_umap(mouse_cds_list[[1]], "cd19_cd5_pos", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")




