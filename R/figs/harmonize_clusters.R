#Fig3 (PRMT5 vs TCL1) & Fig5 (PRMT5xTCL1 vs TCL1) Spleen Cluster Intersection
#Only pulls TCL1 mice from both cds -> modified to include all mice
# fig3_5_intersect_sample_id <- intersect(mouse_cds_list[[1]] |> bb_cellmeta() |> pull(mouse) |> unique(), 
#           mouse_cds_list[[2]] |> bb_cellmeta() |> pull(mouse) |> unique())
fig3_5_intersect_sample_id <- unique(c(mouse_cds_list[[1]] |> bb_cellmeta() |> pull(mouse) |> unique(),
                                        mouse_cds_list[[2]] |> bb_cellmeta() |> pull(mouse) |> unique()))

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


fig3_fig5 <- bind_rows(
  mouse_cds_list[[1]] |>
    bb_cellmeta() |>
    filter(mouse %in% fig3_5_intersect_sample_id) |>
    select(cell_id, clust = kmeans_10_clusters) |>
    mutate(clust = paste0("3.", clust)),
  
  mouse_cds_list[[2]] |>
    bb_cellmeta() |>
    filter(mouse %in% fig3_5_intersect_sample_id) |>
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

ComplexHeatmap::Heatmap(harmony_cluster_matrix, 
                        col = harmony_colfun, 
                        name = "Jaccard")

##Supp5 (PRMT5xTCL1 vs TCL1: Lymphnode) & Fig5 (PRMT5xTCL1 vs TCL1: Spleen) Cluster Intersection
#No overlapping sample IDs -> modified to include all mice
# figS5_5_intersect_sample_id <- intersect(mouse_cds_list[[2]] |> bb_cellmeta() |> pull(sample_id) |> unique(), 
#                                         mouse_cds_list[[4]] |> bb_cellmeta() |> pull(sample_id) |> unique())

figS5_5_intersect_sample_id <- unique(c(mouse_cds_list[[2]] |> bb_cellmeta() |> pull(mouse) |> unique(),
                                          mouse_cds_list[[4]] |> bb_cellmeta() |> pull(mouse) |> unique()))
 
########Jaccard Clusters By Genotype
# colData(mouse_cds_list[[2]])$genotype <- recode(colData(mouse_cds_list[[2]])$genotype,
#                                                 "PRMT5" = "PRMT5/TCL1",
#                                                 "TCL1" = "TCL1",
#                                                 "P/T" = "PRMT5/TCL1")
# colData(mouse_cds_list[[4]])$genotype <- recode(colData(mouse_cds_list[[4]])$genotype,
#                                                 "PRMT5" = "PRMT5/TCL1",
#                                                 "TCL1" = "TCL1",
#                                                 "P/T" = "PRMT5/TCL1")
# figS5_5_intersect_sample_id <- unique(c(filter_cds(
#   mouse_cds_list[[2]],
#   cells = bb_cellmeta(mouse_cds_list[[2]]) |>
#     filter(genotype == "PRMT5/TCL1")) |> bb_cellmeta() |> pull(mouse) |> unique(),
#   filter_cds(
#     mouse_cds_list[[4]],
#     cells = bb_cellmeta(mouse_cds_list[[4]]) |>
#       filter(genotype == "PRMT5/TCL1")) |> bb_cellmeta() |> pull(mouse) |> unique()))
#unique(mouse_cds_list[[4]]$genotype); unique(mouse_cds_list[[2]]$genotype)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


figS5_fig5 <- bind_rows(
  mouse_cds_list[[4]] |>
    bb_cellmeta() |>
    filter(mouse %in% figS5_5_intersect_sample_id) |>
    select(cell_id, clust = kmeans_10_clusters) |>
    mutate(clust = paste0("S5.", clust)),
  
  mouse_cds_list[[2]] |>
    bb_cellmeta() |>
    filter(mouse %in% figS5_5_intersect_sample_id) |>
    select(cell_id, clust = kmeans_10_clusters) |>
    mutate(clust = paste0("5.", clust))
) |> 
  group_by(clust) |> 
  nest()

harmony_cluster_matrix2 <- map_dfr(.x = figS5_fig5$clust,
                                  .f = \(x, dat = figS5_fig5) {
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

harmony_colfun2 <- circlize::colorRamp2(breaks = c(0, 1), colors = c("grey80", "red"))

ComplexHeatmap::Heatmap(harmony_cluster_matrix2, 
                        col = harmony_colfun2, 
                        name = "Jaccard")
