source("00_packages_functions.R")

ind_qc_res <- pmap(.l = list(cds = cds_list, 
                             genome = rep("human", times = length(cds_list))),
                   .f = function(cds, genome) {
                     return_list <- qc_func(cds = cds, genome = genome)
                     return(return_list)
                   })
