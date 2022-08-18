fig3_5_intersect_sample_id <- intersect(mouse_cds_list[[1]] |> bb_cellmeta() |> pull(sample_id) |> unique(),
          mouse_cds_list[[2]] |> bb_cellmeta() |> pull(sample_id) |> unique())

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


fig3_fig5 <- bind_rows(
  mouse_cds_list[[1]] |>
    bb_cellmeta() |>
    filter(sample_id %in% fig3_5_intersect_sample_id) |>
    select(cell_id, clust = kmeans_10_clusters) |>
    mutate(clust = paste0("3.", clust)),
  
  mouse_cds_list[[2]] |>
    bb_cellmeta() |>
    filter(sample_id %in% fig3_5_intersect_sample_id) |>
    select(cell_id, clust = kmeans_10_clusters) |>
    mutate(clust = paste0("5.", clust))
) |> 
  group_by(clust) |> 
  nest()

harmony_cluster_matrix <- map_dfr(.x = fig3_fig5$clust,
    .f = \(x, dat = fig3_fig5) {
      query_barcodes <- dat |> 
        filter(clust == x) |> 
        unnest(data) |> 
        pull(cell_id)
      map2_dfr(.x = dat$data,
              .y = dat$clust,
          .f = \(x,
                 y,
                 q = query_barcodes) {
            x <- x |> pull(cell_id)
            jac <- jaccard(x, q)
            tibble(jac = jac, subject = y)
            
          }) |> mutate(query = x) 
      
    }) |>
  pivot_wider(names_from = "subject", values_from = "jac") |> 
  bb_tbl_to_matrix()

harmony_colfun <- circlize::colorRamp2(breaks = c(0, 1), colors = c("grey80", "red"))

ComplexHeatmap::Heatmap(harmonized_cluster_matrix, 
                        col = harmony_colfun, 
                        name = "Jaccard")
