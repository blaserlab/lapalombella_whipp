source("00_packages_functions.R")

# identify gene modules ####--------------------------------------------------------------
# warning!! Slow!! Destructive!!
# pr_graph_test_res <-
#   graph_test(
#     cds = cds_final,
#     neighbor_graph = "knn",
#     cores = 16,
#     verbose = TRUE
#   )
# 
# pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
# 
# gene_module_df <-
#   find_gene_modules(cds_final[pr_deg_ids,], cores = 39)
# 
# gene_module_df_anno <-
#   left_join(gene_module_df, as_tibble(rowData(cds_final))) %>%
#   write_csv("data_out/gene_module_df_anno.csv")
# 
# # print out the individual gene modules
# dir.create("data_out/individual_modules")
# 
# lapply(
#   X = gene_module_df_anno %>% pull(module) %>% unique(),
#   FUN = function(x) {
#     gene_module_df_anno %>% filter(module == x) %>% write_csv(paste0("data_out/individual_modules/module_", x, ".csv"))
#   }
# )

# make the gene module figures ####--------------------------------------------------------
# make a composite variable of cluster_timepont

colData(cds_final)$cluster_timepoint <- paste0(colData(cds_final)$predicted.celltype.l1, "_", colData(cds_final)$timepoint)

cell_group_df <-
  tibble::tibble(cell = row.names(colData(cds_final)),
                 cell_group = colData(cds_final)$cluster_timepoint)
agg_mat <-
  aggregate_gene_expression(cds = cds_final, 
                            gene_group_df = gene_module_df, 
                            cell_group_df = cell_group_df)

row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))

# GO analysis of gene module genes
modules_of_interest <- c("1","5","7","9")
module_titles <- paste0("Module ",modules_of_interest)
module_pvals <- c(0.01,0.01,0.01,0.01)

module_summaries <- pmap(
  .l = list(module_of_interest = modules_of_interest,
            module_title = module_titles,
            module_pval = module_pvals),
  .f = possibly(function(module_of_interest, module_title, module_pval, moddata = gene_module_df_anno) {
    module_genes_named <- moddata %>%
      select(module, gene_short_name) %>%
      mutate(selected = ifelse(module == module_of_interest, 1, 0)) %>%
      pull(selected)
    names(module_genes_named) <-moddata %>%
      select(module, gene_short_name) %>%
      mutate(selected = ifelse(module == module_of_interest, 1, 0)) %>%
      pull(gene_short_name)
    module_genes_named
    
    sampleGOdata <- new(
      "topGOdata",
      description = "Simple session",
      ontology = "BP",
      allGenes = module_genes_named,
      geneSel = selector,
      nodeSize = 10,
      annot = annFUN.org,
      mapping = "org.Hs.eg.db",
      ID = "symbol"
    )
    
    resultFisher <-
      runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    
    res_table <- GenTable(
      sampleGOdata,
      classicFisher = resultFisher,
      orderBy = "classicFisher",
      ranksOf = "classicFisher",
      topNodes = 100
    ) %>%
      as_tibble(rownames = "Rank")
    
    resultFisher_tbl <-
      tibble(goterm = names(resultFisher@score),
             pval = resultFisher@score)
    
    mod_mat <-
      simplifyEnrichment::GO_similarity(
        resultFisher_tbl %>% filter(pval < module_pval) %>% pull(goterm),
        ont = "BP",
        db = "org.Hs.eg.db",
      )
    
    plot <-
      simplifyEnrichment::simplifyGO(mod_mat, column_title = module_title)
    
    return(list(sampleGOdata, resultFisher, res_table, plot))
  }, 
  otherwise = NULL,
  quiet = F
)) %>% set_names(nm = module_titles)

module_rrvgo <- map(
  .x = module_summaries,
  .f = function(x) {
    simMatrix <-
      calculateSimMatrix(x = x[[3]]$GO.ID,
                         ont = "BP",
                         orgdb = "org.Hs.eg.db")
    scores <- setNames(-log10(ifelse(is.na(
      as.numeric(x[[3]]$classicFisher)),
      1e-30,
      as.numeric(x[[3]]$classicFisher)
    )),
    x[[3]]$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = 0.9,
                                    orgdb = "org.Hs.eg.db")
    returnlist <- list(simMatrix, scores, reducedTerms)
    names(returnlist) <- c("simMatrix", "scores", "reducedTerms")
    return(returnlist)
  }
) %>% set_names(module_titles)

View(module_rrvgo$`Module 7`$reducedTerms)
