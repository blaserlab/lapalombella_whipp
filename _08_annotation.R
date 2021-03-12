source("00_packages_functions.R")

colData(cds_final)$refUMAP_1 <- seurat_data$refUMAP_1
colData(cds_final)$refUMAP_2 <- seurat_data$refUMAP_2
colData(cds_final)$predicted.celltype.l1 <- seurat_data$predicted.celltype.l1
colData(cds_final)$predicted.celltype.l2 <- seurat_data$predicted.celltype.l2
colData(cds_final)$sanity_check <- seurat_data$barcode
sum(rownames(colData(cds_final))!=colData(cds_final)$sanity_check)
colData(cds_final)$sanity_check <- NULL

colData(cds_aligned)$refUMAP_1 <- seurat_data$refUMAP_1
colData(cds_aligned)$refUMAP_2 <- seurat_data$refUMAP_2
colData(cds_aligned)$predicted.celltype.l1 <- seurat_data$predicted.celltype.l1
colData(cds_aligned)$predicted.celltype.l2 <- seurat_data$predicted.celltype.l2
colData(cds_aligned)$sanity_check <- seurat_data$barcode
sum(rownames(colData(cds_aligned))!=colData(cds_aligned)$sanity_check)
colData(cds_aligned)$sanity_check <- NULL

marker_test_res_partition_anno <- marker_test_res_partition %>%
  left_join(
    colData(cds_aligned) %>%
      as_tibble() %>%
      group_by(partition, predicted.celltype.l1) %>%
      summarise(n = n()) %>%
      top_n(1) %>%
      select(cell_group = partition, predicted.celltype.l1)
  ) %>%
  left_join(
    colData(cds_aligned) %>%
      as_tibble() %>%
      group_by(partition, predicted.celltype.l2) %>%
      summarise(n = n()) %>%
      top_n(1) %>%
      select(cell_group = partition, predicted.celltype.l2)
  ) %>%
  write_csv("data_out/marker_test_res_partition_anno.csv")

marker_test_res_leiden_anno <- marker_test_res_leiden %>%
  left_join(
    colData(cds_aligned) %>%
      as_tibble() %>%
      group_by(leiden, predicted.celltype.l1) %>%
      summarise(n = n()) %>%
      top_n(1) %>%
      select(cell_group = leiden, predicted.celltype.l1)
  ) %>%
  left_join(
    colData(cds_aligned) %>%
      as_tibble() %>%
      group_by(leiden, predicted.celltype.l2) %>%
      summarise(n = n()) %>%
      top_n(1) %>%
      select(cell_group = leiden, predicted.celltype.l2)
  ) %>%
  write_csv("data_out/marker_test_res_leiden_anno.csv")

