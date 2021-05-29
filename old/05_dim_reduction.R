# remove cells flagged with qc.any == TRUE
cds_clean <- cds[,colData(cds)$qc.any == FALSE]

# preprocess the cds
cds_clean <- preprocess_cds(cds = cds_clean, num_dim = 100)

# reduce dimensions
cds_clean <- reduce_dimension(cds_clean, cores = 39)

#previz to look at doublets
bb_var_umap(cds_clean, var = "doubletfinder_low_conf", foreground_alpha = 0.2 ,cell_size = 1)
bb_var_umap(cds_clean, var = "doubletfinder_high_conf", foreground_alpha = 0.2 ,cell_size = 1)
# remove the high confidence homotypic doublets
cds_full <- cds_clean[,colData(cds_clean)$doubletfinder_high_conf=="Singlet"]

#remove columns we don't need
colData(cds_full)$doubletfinder_low_conf <- NULL
colData(cds_full)$doubletfinder_high_conf <- NULL
colData(cds_full)$qc.any <- NULL

# recalculate PCAs with cds_full
cds_full <- preprocess_cds(cds = cds_full, num_dim = 100)

# reduce dimensions with cds_full
cds_full <- reduce_dimension(cds_full, cores = 39)

# make a new column for patient ID
colData(cds_full)$patient <- recode(colData(cds_full)$sample, 
                                    "L33_19972712RTLN" = "2712",
                                    "L34_19972712RTPBMC" = "2712",
                                    "L35_19972712CLLPBMC" = "2712",
                                    "L36_19971245RTLN" = "1245",
                                    "L37_19971245RTPBMC" = "1245",
                                    "L38_19971245CLLPBMC" = "1245",
                                    "L39_19973183RTBM" = "3183",
                                    "L40_19973183RTPBMC" = "3183",
                                    )


# align by patient id
cds_aligned <- align_cds(cds = cds_full, alignment_group = "patient")

# reduce dimensions with cds_aligned
cds_aligned <- reduce_dimension(cds_aligned, cores = 39)
