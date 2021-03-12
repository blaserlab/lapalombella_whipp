source("00_packages_functions.R")

# make a new composite column

colData(cds_final)$pt_timepoint <- paste0(colData(cds_final)$pt, "_timepoint_",colData(cds_final)$timepoint)

pseudobulk_res <- map(
  .x =
    list(cds_final[, colData(cds_final)$predicted.celltype.l1 == "B" &
                     colData(cds_final)$timepoint %in% c("1", "2")],
         cds_final[, colData(cds_final)$predicted.celltype.l1 == "B" &
                     colData(cds_final)$timepoint %in% c("1", "3")],
         cds_final[, colData(cds_final)$predicted.celltype.l1 == "B" &
                     colData(cds_final)$timepoint %in% c("2", "3")]),
  .f = function(x, rv = "pt_timepoint", cv = "timepoint") {
    return(pseudobulk_dge(
      cds_deseq = x,
      replicate_variable = rv,
      class_variable = cv
    ))
  }
)


# timepoint 2 vs 1, positive l2fc indicates up at timepoint 2
pseudobulk_res[[1]][[2]] %>% filter(padj<0.05)

# timepoint 3 vs 1, positive l2fc indicates up at timepoint 3
pseudobulk_res[[2]][[2]] %>% filter(padj<0.05)

# timepoint 3 vs 2, positive l2fc indicates up at timepoint 3
pseudobulk_res[[3]][[2]] %>% filter(padj<0.05)
