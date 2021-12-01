colData(cds_aligned) %>% as_tibble()
bb_var_umap(cds_aligned, "density", facet_by = c("tissue", "disease"), rows = (vars(tissue)), cols = vars(disease)) + 
  theme(panel.background = element_rect(color = "grey80"))
bb_var_umap(cds_aligned, "seurat_l1")
bb_gene_umap(cds_aligned, "CD3E")

bb_var_umap(cds_aligned, "density", facet_by = c("sample")) + 
  theme(panel.background = element_rect(color = "grey80"))
