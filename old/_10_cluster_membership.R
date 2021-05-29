source("00_packages_functions.R")

cluster_timepoint_membership <- colData(cds_aligned) %>%
  as_tibble() %>%
  group_by(predicted.celltype.l1,timepoint) %>%
  summarise(n = n()) %>% 
  pivot_wider(names_from = predicted.celltype.l1, values_from = n) %>%
  tbl_to_matrix()



