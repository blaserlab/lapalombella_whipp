
bb_var_umap(cds_aligned[,colData(cds_aligned)$partition %in% c("1","2","5","7","9")], var = "disease") + facet_grid(rows = vars(patient), cols = vars(tissue)) + theme(panel.background = element_rect(color = "grey80"))
bb_var_umap(cds_aligned, var = "patient")
bb_var_umap(cds_aligned, var = "partition",overwrite_labels = T)
bb_var_umap(cds_aligned, var = "leiden",overwrite_labels = T)
bb_gene_umap(cds = cds_aligned, gene_or_genes = "CD79A")

# transfer cluster columns to cds_full
colData(cds_full)$aligned_partition <- colData(cds_aligned)$partition
colData(cds_full)$aligned_leiden <- colData(cds_aligned)$leiden
colData(cds_full)$aligned_louvain <- colData(cds_aligned)$louvain

# pseudobulk for blood vs ln in Richter's only
pseudobulk_tissue_res <- map(
  .x =
    list(cds_full[, colData(cds_full)$aligned_leiden %in% c("1","2") &
                     colData(cds_full)$patient %in% c("1245","2712") &
                     colData(cds_full)$disease == "RT"]),
  .f = function(x, rv = "patient", cv = "tissue") {
    return(bb_pseudobulk(
      cds_deseq = x,
      replicate_variable = rv,
      class_variable = cv
    ))
  }
)


# leiden 3 and 13 vs 4 and 11
colData(cds_full)$leiden_subcluster <- recode(colData(cds_full)$aligned_leiden, "3" = "3_13", "13" = "3_13", "4" = "4_11", "11" = "4_11")

pseudobulk_leiden_sublcuster_res <- map(
  .x =
    list(cds_full[, colData(cds_full)$leiden_subcluster %in% c("3_13","4_11") &
                     colData(cds_full)$patient %in% c("1245","2712") &
                     colData(cds_full)$disease == "RT"]),
  .f = function(x, rv = "patient", cv = "leiden_subcluster") {
    return(bb_pseudobulk(
      cds_deseq = x,
      replicate_variable = rv,
      class_variable = cv
    ))
  }
)

# CLL vs richters

pseudobulk_cll_rt_res <- map(
  .x =
    list(cds_aligned[, colData(cds_aligned)$leiden %in% c("1","8") &
                     colData(cds_aligned)$patient %in% c("1245","2712") &
                     colData(cds_aligned)$tissue == "PBMC"]),
  .f = function(x, rv = "patient", cv = "disease") {
    return(bb_pseudobulk(
      cds_deseq = x,
      replicate_variable = rv,
      class_variable = cv
    ))
  }
)
pseudobulk_cll_rt_res[[1]][[2]] %>% View()









