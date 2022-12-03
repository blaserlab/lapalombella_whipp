#TODO modifications to CDS. Incorporate into data package -----------

#Create Disease_tissue column in cds
colData(cds_main)$disease_tissue <-
  paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
# Reordering group factor levels - For following graphing order
cds_main$disease_tissue <- factor(cds_main$disease_tissue,
                                  levels = c("CLL PBMC", "RT PBMC", "RT LN"))

colData(cds_main)$leiden_assignment_1 <-
  recode(
    colData(cds_main)$leiden,
    "1" = "B",
    "2" = "B",
    "3" = "B",
    "4" = "T",
    "5" = "B",
    "6" = "B",
    "7" = "T",
    "8" = "T",
    "9" = "B",
    "10" = "T",
    "11" = "B",
    "12" = "T",
    "13" = "Mono",
    "14" = "T",
    "15" = "B",
    "16" = "T",
    "17" = "B",
    "18" = "B"
  )

colData(cds_main)$f1d_lbl <- paste0(colData(cds_main)$disease_tissue, "-", recode(colData(cds_main)$patient, "pt_2712" = "Pt 1", "pt_1245" = "Pt 2"))
colData(cds_main)$patient <- recode(colData(cds_main)$patient, "pt_2712" = "Pt 1", "pt_1245" = "Pt 2")

#Compare to Lit data - Nadeu et al
nadeu_11b <-
  readxl::read_excel(
    "~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx",
    sheet = "Supplementary Table 11b",
    skip = 5,
    col_names = c(
      "feature_id",
      "gene_short_name",
      "mean",
      "l2fc",
      "se",
      "p",
      "padj",
      "direction"
    )
  )

cds_main <- nadeu_11b |> 
  filter(direction == "Up") |> 
  filter(padj < 0.05) |> 
  mutate(feature_id = str_remove(feature_id, "\\..*")) |> 
  mutate(nadeu_RT_gene = TRUE) |> 
  bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)

cds_main <- nadeu_11b |>
 filter(direction == "Down") |>
 filter(padj < 0.05) |>
 mutate(feature_id = str_remove(feature_id, "\\..*")) |>
 mutate(nadeu_CLL_gene = TRUE) |>
 bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)

#leiden_assignment_2 assigned w/S1E & F1E1
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
###############################################################################

#Figure 1A
F1AP1 <-
  bb_var_umap(
    filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
    "leiden_assignment_1", overwrite_labels = T,
    facet_by = "disease_tissue"
  ) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1AP2 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
                     "log_local_n",
                     facet_by = "disease_tissue") +
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "Log<sub>10</sub><br>Cells") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.4,"cm"))+
  theme(legend.title = ggtext::element_markdown())
F1AP3 <- bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
                      gene_or_genes = "PRMT5") +
  facet_wrap( ~ disease_tissue) +
  theme(strip.text = element_blank()) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*PRMT5*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.4,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
  
F1A <- F1AP1/F1AP2/F1AP3

#Figure 1D
F1D0 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                    "leiden_assignment_1", facet_by = "f1d_lbl", overwrite_labels = T) +
    theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(panel.background = element_rect(color = "black"))+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D1<- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                   "density", facet_by = "f1d_lbl") + 
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*Density*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D2<- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("PRMT5") 
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*PRMT5*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D3<- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("MYC")
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*MYC*") +
  theme(legend.title=element_text(size=9))+
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D4 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("MKI67")
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*MKI67*") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.key.size = unit(0.3,"cm"))+
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

F1D <-
  as_ggplot(grid.arrange(
    patchworkGrob(F1D0 / F1D1 / F1D2 / F1D3 / F1D4),
    left = textGrob(F1D3$labels$y, rot=90, vjust = 1.5, hjust=0.35),
    bottom = textGrob(F1D3$labels$x, hjust = 0.8, vjust = -0.5))
  )

#Figure 1E
F1E1 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
  gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
) +
  facet_wrap( ~ disease_tissue) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*RT Gene*<br>Expression") +
  theme(legend.title = ggtext::element_markdown())

F1E2 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "Pt 1")),
  gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_CLL_gene)
) +
  facet_wrap( ~ disease_tissue) +
  theme(strip.text = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*CLL Gene*<br>Expression") +
  theme(legend.title = ggtext::element_markdown())

F1E <-
  as_ggplot(grid.arrange(
    patchworkGrob(F1E1 / F1E2),
    left = F1E1$labels$y,
    bottom = textGrob(F1E1$labels$x, hjust = 0.85, vjust = -1)
  ))

#Figure 1F
F1F <- bb_var_umap(filter_cds(cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(disease_tissue == "RT LN")), "leiden", overwrite_labels = T
  ) + labs(title = "Combined RT LN")

#Figure 1G
#Top markers
LN_B_leiden_Top50_tm <- monocle3::top_markers(
  filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |>
      filter(disease_tissue == "RT LN") |>
      filter(leiden_assignment_1 == "B")),
    group_cells_by = "leiden",
    genes_to_test_per_group = 50,
    cores = 12)

markers <- LN_B_leiden_Top50_tm |> filter(cell_group %in% c("3","11")) |> pull(gene_short_name)

#subsetting
ln_cds<-cds_main[,colData(cds_main)$disease_tissue == "RT LN" &
                   colData(cds_main)$leiden_assignment_1 == "B"]

f1_mat <- bb_aggregate(obj = filter_cds(ln_cds,
                                        cells = bb_cellmeta(ln_cds) |>
                                          filter(leiden %in% c('3','11','6','2','5','9','1')),
                                        genes = bb_rowmeta(ln_cds) |>
                                          filter(gene_short_name %in% markers)),
                       cell_group_df = bb_cellmeta(ln_cds) |>
                         select(cell_id, leiden)) |>
  t() |>
  scale() |>
  t()

rownames(f1_mat) <- tibble(feature_id = rownames(f1_mat)) |>
  left_join(bb_rowmeta(ln_cds) |>
              select(feature_id, gene_short_name)) |>
  pull(gene_short_name)

f1_colfun = circlize::colorRamp2(breaks = c(min(f1_mat),
                                            0,
                                            max(f1_mat)),
                                 colors = heatmap_3_colors)

F1_highlights <- c("CCNA2", "HMGB2", "NPM1", "TUBB", "CDK1", "MKI67", "TUBA1B", "TUBA1C", "TCL1A", "TOP2A", "BIRC5", "JUNB", "UBE2S", "UBE2C", "STMN1", "CKS1B", "MS4A1", "GAPDH", "PFDN5", "BTG1", "UBC", "CD19", "AURKB", "CCR7", "KLF2", "ZFP36", "CD37", "XBP1", "RAC2", "GAS5", "BUB1", "HMMR", "CDKN2A", "ANP32B", "IGHG3", "PLK1", "TPX2", "AURKA")

fig1_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(f1_mat) %in% F1_highlights),
  labels = rownames(f1_mat)[rownames(f1_mat) %in% F1_highlights],
  labels_gp = gpar(fontsize = 6),
  padding = unit(4.5, "mm")
  ))

F1G <- grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(f1_mat,
                        col = f1_colfun,
                        name = "Expression", 
                        show_row_names = F, 
                        right_annotation = fig1_anno,
                        row_dend_width = unit(4, "mm"),
                        column_dend_height = unit(2, "mm"),
                        heatmap_legend_param = list(legend_direction = "vertical",
                                                    #legend_width = unit(1, "mm"),
                                                    title_position = "topleft", 
                                                    title_gp = gpar(fontsize = 6)
                                               ))))

F1FG <- F1F / F1G +
  plot_layout(heights = c(1, 2.5))
F1FG

#Figure 1H: Pseudobulk Volcano
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

# Differential expression results.
genes_to_highlight <- unique(c("FOSB", "PRMT5","FOXM1", F1_highlights))
genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(pseudobulk_res$Result, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]

volcano_data_RTvCLL <- pseudobulk_res$Result %>%
  mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))

library(ggtext)
volcano_pseudob_RTvCLL <-
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
                  box.padding = 0.22, #0.5
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
volcano_pseudob_RTvCLL
