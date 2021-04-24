doubletfinder_list <-
  pmap(
    .l = list(
      cds = cds_list,
      doublet_prediction = read_csv("data_out/summarized_gex_sequencing_metrics.csv") %>%
        mutate(anticipated_doublet_rate = `Estimated Number of Cells`/100000) %>%
        pull(anticipated_doublet_rate),
      qc_table = map(ind_qc_res,1)
    ),
    .f = bb_doubletfinder
  )

names(doubletfinder_list) <- names(cds_list)

