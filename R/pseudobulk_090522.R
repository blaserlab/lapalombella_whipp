# 
# S1.1 <- bb_genebubbles(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(leiden_assignment_1 == "B")),
#   genes = c("PRMT5", "MYC", "MS4A1"
#   ),
#   cell_grouping = c("leiden", "leiden_assignment_1"),
#   return_value = "data"
# ) |> 
#   ggplot(mapping = aes(x = leiden, 
#                        y = gene_short_name, 
#                        color = expression,
#                        size = proportion)) +
#   geom_point() +
#   scale_size_area() +
#   scale_color_viridis_c() +
#   facet_wrap(~leiden_assignment_1, scales = "free_x", ) +
#   theme_minimal_grid(font_size = 8) +
#   theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
#   theme(axis.text.y = element_text(face = "italic")) +
#   labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
# S1.1
#leiden_clusters<- bb_var_umap(cds_main, "leiden", overwrite_labels = T)
#ggsave("leiden_clusters.pdf", path = WalkerTables)

#ln_cds<-cds_main[,colData(cds_main)$disease_tissue == "RT LN"]

bb_var_umap(cds_main, "leiden", overwrite_labels = T, facet_by= "disease_tissue")

colData(cds_main)$leiden_assignment_2 <-
  recode(
    colData(cds_main)$leiden,
    "1" = "CLL_B",
    "2" = "Int_PBMC_B",
    "3" = "RT_B",
    "4" = "T",
    "5" = "Int_LN_B",
    "6" = "Int_LN_B",
    "7" = "T",
    "8" = "T",
    "9" = "CLL_B",
   "10" = "T",
   "11" = "RT_B",
   "12" = "T",
   "13" = "Mono",
   "14" = "T",
    "15" = "B",
    "16" = "T",
    "17" = "B",
    "18" = "B"
  )
bb_var_umap(cds_main, "leiden", overwrite_labels = T, facet_by= "disease_tissue")
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

# Differential expression results.  Positive L2FC indicates up in B vs T upregulated
genes_to_highlight <- unique(c("FOSB", "JUN","PRMT5","FOXM1", F1_highlights))
genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(pseudobulk_res$Result, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]


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
                  box.padding = 0.23, #0.5
                  point.padding = 0.1, #0.25
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
  coord_cartesian(xlim = c(-1.0*max(abs(range(volcano_data_RTvCLL %>% dplyr::filter(!is.na(padj)) %>% pull(log2FoldChange)))), 1.0*max(abs(range(volcano_data_RTvCLL %>% filter(!is.na(padj)) %>% pull(log2FoldChange))))))
volcano_plot_RTvCLL

#Gene Ontology
pseudo3n11_goenrichment <-
  bb_goenrichment(
    query = dplyr::filter(pseudobulk_res$Result, padj < 0.1 &
                            log2FoldChange >= 0.58) |> pull(gene_short_name),
    reference = bb_rowmeta(cds_main),
    go_db = "org.Hs.eg.db"
  )
pseudo3n11_gosummary_0.9 <- bb_gosummary(x = pseudo3n11_goenrichment, 
                                 reduce_threshold = 0.9,
                                 go_db = "org.Hs.eg.db")
pseudobulk_Clust3n11_GO_PCA <-
  bb_goscatter(simMatrix = pseudo3n11_gosummary_0.9$simMatrix,
               reducedTerms = pseudo3n11_gosummary_0.9$reducedTerms)
#pseudoGO barplot
pseudo3n11_goenrichment$res_table$classicFisher[pseudo3n11_goenrichment$res_table$classicFisher=="< 1e-30"]<-"1.0e-30"
pseudo3n11_goenrichment$res_table$Rank <- as.numeric(as.character(pseudo3n11_goenrichment$res_table$Rank))
pseudob_top25 <- filter(pseudo3n11_goenrichment$res_table, as.numeric(pseudo3n11_goenrichment$res_table$Rank) <= 25) |> mutate(neg_log10_pval = -log10(as.numeric(classicFisher))) |> rename(GO_Term = Term, Genes_Mapped = Significant)
pseudob_top25

pseudob_3n11vs1n9_GObp<- ggplot(data=pseudob_top25, aes(reorder(x= GO_Term, y= neg_log10_pval, neg_log10_pval), y= neg_log10_pval, fill = Genes_Mapped)) +
  geom_bar(stat="identity") + 
  coord_flip() + scale_fill_viridis_c() + labs(x = "GO Terms", y = "-log(pval)")#theme(axis.title.x = element_text("GO Terms"), axis.title.y = element_text("-log(pval)"))
pseudob_3n11vs1n9_GObp
ggsave("pseudob_3n11vs1n9_GObp.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1")

#Clust 3 GO:
C3_goenrichment <-
  bb_goenrichment(
    query = dplyr::filter(LN_B_leiden_Top50_tm, cell_group %in% '3')[["gene_short_name"]],
    reference = bb_rowmeta(cds_main),
    go_db = "org.Hs.eg.db"
  )
C3_gosummary_0.9 <- bb_gosummary(x = C3_goenrichment, 
                                   reduce_threshold = 0.9,
                                   go_db = "org.Hs.eg.db")
S1_Clust3_GO_PCA <- bb_goscatter(simMatrix = C3_gosummary_0.9$simMatrix,
             reducedTerms = C3_gosummary_0.9$reducedTerms)
#barplot
C3_goenrichment$res_table$classicFisher[C3_goenrichment$res_table$classicFisher=="< 1e-30"]<-"1.0e-30"
C3_goenrichment$res_table$Rank <- as.numeric(as.character(C3_goenrichment$res_table$Rank))
C3_top25 <- filter(C3_goenrichment$res_table, as.numeric(C3_goenrichment$res_table$Rank) <= 25) |> mutate(neg_log10_pval = -log10(as.numeric(classicFisher))) |> rename(GO_Term = Term, Genes_Mapped = Significant)
C3_top25

C3_GObp<- ggplot(data=C3_top25, aes(reorder(x= GO_Term, y= neg_log10_pval, neg_log10_pval), y= neg_log10_pval, fill = Genes_Mapped)) +
  geom_bar(stat="identity") + 
  coord_flip() + scale_fill_viridis_c() + labs(x = "GO Terms", y = "-log(pval)")#theme(axis.title.x = element_text("GO Terms"), axis.title.y = element_text("-log(pval)"))
C3_GObp
ggsave("C3_GObp.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1")
#ggsave("S1_Clust3_GO_PCA.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1")

#Clust 11 GO:
C11_goenrichment <- bb_goenrichment(query = dplyr::filter(LN_B_leiden_Top50_tm, cell_group %in% '11')[["gene_short_name"]], 
                                   reference = bb_rowmeta(cds_main),
                                   go_db = "org.Hs.eg.db")
C11_gosummary_0.9 <- bb_gosummary(x = C11_goenrichment, 
                                  reduce_threshold = 0.9,
                                  go_db = "org.Hs.eg.db")
S1_Clust11_GO_PCA <- bb_goscatter(simMatrix = C11_gosummary_0.9$simMatrix,
                                 reducedTerms = C11_gosummary_0.9$reducedTerms)
#ggsave("S1_Clust11_GO_PCA.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1")

C11_goenrichment$res_table$classicFisher[C11_goenrichment$res_table$classicFisher=="< 1e-30"]<-"1.0e-30"
C11_goenrichment$res_table$Rank <- as.numeric(as.character(C11_goenrichment$res_table$Rank))
C11_top25 <- filter(C11_goenrichment$res_table, as.numeric(C11_goenrichment$res_table$Rank) <= 25) |> mutate(neg_log10_pval = -log10(as.numeric(classicFisher))) |> rename(GO_Term = Term, Genes_Mapped = Significant)
C11_top25

C11_GObp<- ggplot(data=C11_top25, aes(reorder(x= GO_Term, y= neg_log10_pval, neg_log10_pval), y= neg_log10_pval, fill = Genes_Mapped)) +
  geom_bar(stat="identity") + 
  coord_flip() + scale_fill_viridis_c() + labs(x = "GO Terms", y = "-log(pval)")#theme(axis.title.x = element_text("GO Terms"), axis.title.y = element_text("-log(pval)"))
C11_GObp
ggsave("C11_GObp.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1")


