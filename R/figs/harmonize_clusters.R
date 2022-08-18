mouse_cds_list[[1]] |> bb_cellmeta()

fig3_5_intersect_sample_id <- intersect(mouse_cds_list[[1]] |> bb_cellmeta() |> pull(sample_id) |> unique(),
          mouse_cds_list[[2]] |> bb_cellmeta() |> pull(sample_id) |> unique())

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


mouse_cds_list[[1]] |> 
  bb_cellmeta() |> 
  filter(sample_id %in% fig3_5_intersect_sample_id) |> 
  select(cell_id, clust = kmeans_10_clusters) |> 
  mutate(clust = paste0("3.", clust))

mouse_cds_list[[2]] |> 
  bb_cellmeta() |> 
  filter(sample_id %in% fig3_5_intersect_sample_id) |> 
  select(cell_id, clust = kmeans_10_clusters) |> 
  mutate(clust = paste0("5.", clust))




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
