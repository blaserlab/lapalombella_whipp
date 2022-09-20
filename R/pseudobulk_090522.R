bb_genebubbles(
   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
   genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8A","CD4"), cell_grouping = "leiden") + labs(x = "Leiden Cluster", y = NULL) 

#leiden_clusters<- bb_var_umap(cds_main, "leiden", overwrite_labels = T)
#ggsave("leiden_clusters.pdf", path = WalkerTables)

#bb_var_umap(cds_main, "leiden", overwrite_labels = T, facet_by= "disease_tissue")
#ln_cds<-cds_main[,colData(cds_main)$disease_tissue == "RT LN"]


colData(cds_main)$leiden_assignment_2 <-
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
bb_var_umap(cds_main, "leiden_assignment_1", overwrite_labels = T, facet_by= "disease_tissue")
bb_var_umap(cds_main, "leiden_assignment_2", overwrite_labels = T, facet_by= "disease_tissue")
bb_var_umap(cds_main, "seurat_celltype_l1", overwrite_labels = T, facet_by= "disease_tissue")
bb_var_umap(cds_main, "seurat_celltype_l2", overwrite_labels = T, facet_by= "disease_tissue")



exp_design <- 
  bb_cellmeta(cds_main) |>
  group_by(disease_tissue, leiden_assignment_2) |>
  summarise()
#exp_design

pseudobulk_res <-
  bb_pseudobulk_mf(cds = cds_main,
                   pseudosample_table = exp_design, 
                   design_formula = "~ leiden_assignment_2",
                   result_recipe = c("leiden_assignment_2", "RT_B", "CLL_B"))

pseudobulk_res$Header

# Differential expression results.  Positive L2FC indicates up in B vs T
#upregulated
genes_to_highlight <- c("MKI67", "PRMT5","FOXM1", "BIRC5")

volcano_data_RTvCLL <- pseudobulk_res$Result %>%
  mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))

volcano_plot_RTvCLL <-
  ggplot(
    volcano_data_RTvCLL,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21,
             size = 0.5,
             alpha = 0.4) +
  geom_text_repel(color = "black",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  size = 3,
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log<sub>2</sub> fold change") +
  ylab("-log<sub>10</sub> adjusted p-value") +
  theme(axis.title.x =  element_markdown()) +
  theme(axis.title.y = element_markdown()) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in CLL\nUp in RT \U21D2",title = "Pseudobulk:RT clusters 3 & 11 vs CLL clusters 1 & 9")+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-1.1*max(abs(range(volcano_data_RT %>% dplyr::filter(!is.na(padj)) %>% pull(log2FoldChange)))), 1.1*max(abs(range(volcano_data_RT %>% filter(!is.na(padj)) %>% pull(log2FoldChange))))))
volcano_plot_RTvCLL

pseudobulk_RT <- leiden_3n11_vs_1n9_up |> filter(padj < 0.01)

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
