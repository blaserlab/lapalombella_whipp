bb_genebubbles(
   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
   genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8A","CD4"), cell_grouping = "leiden") + labs(x = "Leiden Cluster", y = NULL) 

#leiden_clusters<- bb_var_umap(cds_main, "leiden", overwrite_labels = T)
#ggsave("leiden_clusters.pdf", path = WalkerTables)

#bb_var_umap(cds_main, "leiden", overwrite_labels = T, facet_by= "disease_tissue")
#ln_cds<-cds_main[,colData(cds_main)$disease_tissue == "RT LN"]


colData(cds_main)$leiden_assignment_1 <-
  recode(
    colData(cds_main)$leiden,
    "1" = "CLL_B",
    "2" = "Int_B",
    "3" = "RT_B",
    "4" = "T",
    "5" = "Int_B",
    "6" = "Int_B",
    "7" = "T",
    "8" = "T",
    "9" = "CLL_B",
   "10" = "T",
   "11" = "RT_B",
   "12" = "T",
   "13" = "Mono",
   "14" = "T",
    "15" = "Int_B",
    "16" = "T",
    "17" = "Int_B",
    "18" = "B"
  )
exp_design <- 
  bb_cellmeta(cds_main) |>
  group_by(disease_tissue, leiden_assignment_1) |>
  summarise()
#exp_design

pseudobulk_res <-
  bb_pseudobulk_mf(cds = cds_main,
                   pseudosample_table = exp_design, 
                   design_formula = "~ leiden_assignment_1",
                   result_recipe = c("leiden_assignment_1", "RT_B", "CLL_B"))

pseudobulk_res$Header

# Differential expression results.  Positive L2FC indicates up in B vs T
#upregulated
leiden_3n11_vs_1n9_up <- pseudobulk_res$Result %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj)

RTup<- leiden_3n11_vs_1n9_up[1:100,2]
RTup<- as.vector(RTup$gene_short_name)
#write_csv(leiden_3n11_vs_1n9_upregulated, file = file.path(WalkerTables, "leiden_3n11_vs_1n9_upregulated.csv"))

#downregulated
leiden_3n11_vs_1n9_downregulated<- pseudobulk_res$Result %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj)
#write_csv(leiden_3n11_vs_1n9_downregulated, file = file.path(WalkerTables, "leiden_3n11_vs_1n9_downregulated.csv"))

#GO
bb_goenrichment(RTup, as_tibble(rowData(cds_main)))
blaseRtools:::bb_goenrichment

colData(rowData(cds_main))
