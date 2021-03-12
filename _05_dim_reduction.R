source("00_packages_functions.R")

# remove cells flagged with qc.any == TRUE
cds_clean <- cds[,colData(cds)$qc.any == FALSE]

# preprocess the cds
cds_clean <- preprocess_cds(cds = cds_clean, num_dim = 100)

# reduce dimensions
cds_clean <- reduce_dimension(cds_clean, cores = 39)

#previz to look at doublets
# custom_variable_plot(cds_clean, var = "doubletfinder_low_conf", foreground_alpha = 0.2 ,cell_size = 1) 
# custom_variable_plot(cds_clean, var = "doubletfinder_high_conf", foreground_alpha = 0.2 ,cell_size = 1) 
# remove the high confidence homotypic doublets
cds_final <- cds_clean[,colData(cds_clean)$doubletfinder_high_conf=="Singlet"]

#remove columns we don't need
colData(cds_final)$doubletfinder_low_conf <- NULL

# recalculate PCAs with cds_final
cds_final <- preprocess_cds(cds = cds_final, num_dim = 100)

# reduce dimensions with cds_final
cds_final <- reduce_dimension(cds_final, cores = 39)

# align by patient id
cds_aligned <- align_cds(cds = cds_final, alignment_group = "pt")

# reduce dimensions with cds_aligned
cds_aligned <- reduce_dimension(cds_aligned, cores = 39)




