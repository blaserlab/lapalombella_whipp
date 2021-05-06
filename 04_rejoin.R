source("00_packages_functions.R")

cds_list_rejoined <- pmap(
  .l = list(
    cds = cds_list,
    qc_data = map(ind_qc_res,1),
    doubletfinder_data = doubletfinder_list
  ),
  .f = function(cds, qc_data, doubletfinder_data) {
    # qc_data <- qc_data, 1)
    cds_tbl <- as_tibble(colData(cds))
    cds_tbl <- left_join(cds_tbl, qc_data)
    cds_tbl <- left_join(cds_tbl, doubletfinder_data)
    cds_df <- as.data.frame(cds_tbl)
    row.names(cds_df) <- cds_df$barcode
    cds <- new_cell_data_set(expression_data = cds@assays@data$counts,
                             cell_metadata = cds_df, gene_metadata = rowData(cds))
    return(cds)}
) %>% set_names(names(cds_list)) 

# combine the cds's-------------------------------

cds <- combine_cds(cds_list = cds_list_rejoined, keep_all_genes = TRUE)
