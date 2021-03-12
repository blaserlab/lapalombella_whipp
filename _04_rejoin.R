source("00_packages_functions.R")

cds_list_rejoined <- pmap(
  .l = list(
    cds = cds_list,
    qc_data = purrr::map(ind_qc_res, 1),
    doubletfinder_data = doubletfinder_list
  ),
  .f = join_metadata
)


# combine the cds's-------------------------------

cds <- combine_cds(cds_list = cds_list_rejoined, keep_all_genes = TRUE)

