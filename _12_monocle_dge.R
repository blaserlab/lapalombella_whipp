source("00_packages_functions.R")

#monocle_dge_bcells_timepoint <-
#  top_markers(
#    cds = cds_final[, colData(cds_final)$seurat.celltype.l1 == "B"],
#    group_cells_by = "timepoint",
#    genes_to_test_per_group = 50,
#    reference_cells = 1000,
#    cores = 39,
#  )
#monocle_dge_bcells_timepoint %>% arrange(cell_group) %>% rename(timepoint = cell_group) %>% write_csv("data_out/timepoint_top_markers.csv")

monocle_dge_12_bysample <-
  top_markers(
    cds = cds_final[, colData(cds_final)$partition %in% c("1","2")],
    group_cells_by = "partition",
    genes_to_test_per_group = 500,
    reference_cells = 1000,
    cores = 39,
  )
view(monocle_dge_12_bysample)
colData(cds_final)
view(cds_main_top_markers)
monocle_dge_12_bysample %>% arrange(cell_group) %>% write_csv("data_out/top_markers.csv")

monocle_dge_bcell_bypartition <-
  top_markers(
    cds = cds_final[, colData(cds_final)$seurat.celltype.l1 == "B"],
    group_cells_by = "partition",
    genes_to_test_per_group = 500,
    reference_cells = 1000,
    cores = 39,
  )

view(monocle_dge_bcell_bypartition)

monocle_dge_leiden_bysample <-
  top_markers(
    cds = cds_final[, colData(cds_final)$leiden %in% c("1","9", "5", "3", "11")],
    group_cells_by = "leiden",
    genes_to_test_per_group = 500,
    reference_cells = 1000,
    cores = 39,
  )
view(monocle_dge_leiden_bysample)
write.table(monocle_dge_leiden_bysample, "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/monocle_dge_leiden_bysample.txt", sep="\t", row.names = F)
