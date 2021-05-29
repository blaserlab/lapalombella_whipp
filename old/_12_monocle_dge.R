source("00_packages_functions.R")

monocle_dge_bcells_timepoint <-
  top_markers(
    cds = cds_final[, colData(cds_final)$predicted.celltype.l1 == "B"],
    group_cells_by = "timepoint",
    genes_to_test_per_group = 50,
    reference_cells = 1000,
    cores = 39,
  )
monocle_dge_bcells_timepoint %>% arrange(cell_group) %>% rename(timepoint = cell_group) %>% write_csv("data_out/timepoint_top_markers.csv")
