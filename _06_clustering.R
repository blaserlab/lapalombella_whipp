source("00_packages_functions.R")

# cluster cells
cds_aligned <- cluster_cells(cds_aligned, cluster_method = "leiden")
cds_louvain <- cluster_cells(cds_aligned, cluster_method = "louvain")
  
colData(cds_aligned)$leiden <- clusters(cds_aligned)
colData(cds_aligned)$partition <- partitions(cds_aligned)
colData(cds_aligned)$louvain <- clusters(cds_louvain)

rm(cds_louvain)

custom_variable_plot(cds_aligned, var = "partition")
custom_variable_plot(cds_aligned, var = "leiden")
custom_variable_plot(cds_aligned, var = "louvain")

marker_test_res_partition <-
  top_markers(
    cds_aligned,
    group_cells_by = "partition",
    reference_cells = 1000,
    cores = 39
  )

marker_test_res_leiden <-
  top_markers(
    cds_aligned,
    group_cells_by = "leiden",
    reference_cells = 1000,
    cores = 39
  )

marker_test_res_louvain <-
  top_markers(
    cds_aligned,
    group_cells_by = "louvain",
    reference_cells = 1000,
    cores = 39
  )


