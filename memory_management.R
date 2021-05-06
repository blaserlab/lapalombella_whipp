# cached_items_1 <- c(
#   "cds",
#   "cds_clean",
#   "cds_list",
#   "cds_list_rejoined"
# )
# 
# save.pigz(list = cached_items_1, file = "rdata/cached_items_1.RData", n.cores = 8)
# rm(list = cached_items_1)

# save cds full
saveRDS.pigz(object = cds_full, file = "rdata/cds_full.rds", n.cores = 8)
# rm(cds_full)

# save cds aligned
saveRDS.pigz(object = cds_aligned, file = "rdata/cds_aligned.rds", n.cores = 8)
# rm(cds_aligned)

save.image.pigz("rdata/active_workspace.RData",n.cores = 8)
