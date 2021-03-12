source("00_packages_functions.R")

anticipated_doublet_rate <- summarized_sequencing_metrics %>%
  mutate(anticipated_doublet_rate = cds_dim_cells/100000) %>%
   pull(anticipated_doublet_rate)

doubletfinder_list <-
  pmap(
    .l = list(
      cds = cds_list,
      doublet_prediction = anticipated_doublet_rate,
      qc_table = map(ind_qc_res,1)
    ),
    .f = find_homotypic_doublets
  )

names(doubletfinder_list) <- names(cds_list)

