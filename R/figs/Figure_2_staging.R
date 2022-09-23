T_tables_out <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/Bulk RNAseq"
#Supp 2 - CLEAR RNAseq -Additional TCL1 mice added
CLEARseq_data <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/Bulk RNAseq/CLEAR_v2_sheet_prmt5_vs_tscl1.2.3.4.csv")[,c(1,3,7)]
#Genes highlighted in F3 or F5 heatmaps
# F3_F5_highlights <- c("CCR7", "BIRC5", "JUNB", "ATF4", "CD69", "TOP2A", "HMGB1", "HMGB2", "CD83", "TUBB5", "TUBA1B", "BLNK", "HSPA5", 
#                                "CCDC34", "UBE2C", "IL2RG", "CKS1B", "SMC2", "STMN1", "TYMS", "UBC", "GAPDH", "APOE", "XRCC1",
#                                "CD79A","ZFP36","PPP1R15A", "LDHA", "MKI67","CDK1", "GPX4", "BTG1", "ITM2B", "XBP1", "REL", "AURKB",
#                                "EIF5A","C1QBP", "CD79B", "TK1","NFKBIA","HSP90AA1", "CENPM", "NR4A1","PFDN5","EIF4A2", "PIM1", "NPM1", 
#                                "MARCKSL1", "ACTB", "NAPSA","CD37", "CD19", "CORO1A", "KLF2", "HNRNPK","RAC2", "MS4A1") |>str_to_title()
#F3_F5_highlights <- str_to_title(F3_F5_highlights)
#Lymphome Genes
# lymphoma_genes<- readxl::read_excel("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/disgenet.org_LymphomaGenes_C0024299_disease_gda_summary.xlsx")[3]
# lymphoma_genes <- str_to_title(lymphoma_genes$Gene)
# filt<- CLEARseq_data |> filter(CLEARseq_data$Gene %in% lymphoma_genes) |> filter(padj <0.05 & abs(log2FoldChange)>=1)
#write_csv(PRMT5vTCL1_CLEAR_filtered, file = file.path(T_tables_out, "PRMT5vTCL1_CLEAR_filtered.csv"))

#IPA_top_path<-c("AURKA","BAK1","BCL2L1","BCL2L11","BRAF","BRCA1","CDC25B","CDK9","CDKN1A","CDKN2D","CHEK1","E2F1","E2F8","EP300","FOS","FYN","HDAC6","HDAC7","HIF1A","ITGB2","JUN","MAPK1","MAPK14","NBN","NFKBIA","NRAS","PIK3CB","PRKAR1A","RAC2","RASGRP1","RHOH","TCF3","TGFBR1","ZBTB17")
#IPA_B <- c("CD274","CXCR4","CD19","CD79A","CD79B","FOS","FYN","JUN","JUND","LYN","MAPK1","NFKBIA","NRAS","PIK3CB","PLCG2","PTEN","SYK", "BRAF", "AURKA", "BCL2L1")      
F2_highlights_lg2FCcutoff1 <-c("CD274", "CD79A","CD79B","FOS","FOSB","FYN","JUN","JUND","LYN","NFKBIA","NFKBIZ","PIK3CB","SYK", "BRAF", "BCL2L1","BRCA1","XPO1","BRCA2","CXCL16","ZAP70","FOSB","Zfp36","Ubc","RAC2","PLK2","Il1r2","Atf3","Arg2","MMP24","Hspa1l","Elane","Stat4")      
F2_highlights_lg2FCcutoff2 <-c("CD79B","FOS","FOSB","FYN","JUN","JUND","CXCL16","ZAP70","FOSB","Zfp36","Ubc","PLK2","Il1r2","Atf3","Arg2","MMP24","Hspa1l","Elane","Stat4")      

F2_highlights<- F2_highlights_lg2FCcutoff1 |> str_to_title()


#####VOLCANO PLOT
CLEAR_volcano_data <- CLEARseq_data |>
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 2) |>
  mutate(text_label = ifelse(Gene %in% F2_highlights, Gene, ""))

volcano_CLEAR <-
  ggplot(
    CLEAR_volcano_data,
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
                  box.padding = 0.1, #0.5
                  point.padding = 0.0, #0.25
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  size = 3,
                  segment.size = 0.15,
                  force = 1,
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
  labs(title = "PRMT5 vs TCL1")+
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-0.75*max(abs(range(CLEAR_volcano_data |> dplyr::filter(!is.na(padj)) |> pull(log2FoldChange)))), 1*max(abs(range(CLEAR_volcano_data |> filter(!is.na(padj)) |> pull(log2FoldChange))))))
volcano_CLEAR
volcano_CLEAR_log2FCcutoff2_padjcutoff0.05 <-volcano_CLEAR
ggsave("volcano_CLEAR_log2FCcutoff2_padjcutoff0.05.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/Bulk RNAseq", width = 7.99, height = 5.82)



