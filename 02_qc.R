ind_qc_res <- pmap(.l = list(cds = cds_list,
                             cds_name = names(cds_list),
                             genome = rep("human", times = length(cds_list))),
                   .f = function(cds, cds_name, genome) {
                     return_list <- bb_qc(cds = cds, cds_name = cds_name, genome = genome)
                     return(return_list)
                   })


