library(gridExtra)
library(patchwork)

#Create Disease_tissue column in cds
colData(cds_main)$disease_tissue <- paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
# Reordering group factor levels - For following graphing order
cds_main$disease_tissue <- factor(cds_main$disease_tissue,      
                                  levels = c("CLL PBMC", "RT PBMC", "RT LN"))
#partition assignment
# bb_genebubbles(
#   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
#   genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8"), cell_grouping = "partition") + labs(x = "Partition Cluster", y = NULL) 
colData(cds_main)$partition_assignment_1 <- recode(colData(cds_main)$partition, "1" = "B", "2" = "B", "3" = "T", "4" = "T", "5" = "Mono", "6" = "B", "7" = "B")

#Compare to Lit data - Nadeu et al
nadeu_11b <- readxl::read_excel("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx", sheet = "Supplementary Table 11b", skip = 5, col_names = c("feature_id", "gene_short_name", "mean", "l2fc", "se", "p", "padj", "direction"))

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

colData(cds_main)$pt_recode <- recode(colData(cds_main)$patient, "pt_2712" = "Pt 1", "pt_1245" = "Pt 2")

colData(cds_main)$f1d_lbl <- paste0(colData(cds_main)$disease_tissue, "-", colData(cds_main)$pt_recode)


#supplemental figs
S1A <- bb_var_umap(cds_main, "pt_recode", facet_by = "value", legend_pos = "none", foreground_alpha = 0.2)

# possible alt S1A
# bb_var_umap(cds_main, "patient", value_to_highlight = "pt_2712", legend_pos = "none", plot_title = "pt_2712") +
# bb_var_umap(cds_main, "patient", value_to_highlight = "pt_1245", legend_pos = "none", plot_title = "pt_1245")

S1B <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "pt_2712")), "partition")

#possible alt S1B -----
#S1B <- bb_var_umap(cds_main, var = "partition", facet_by = "pt_recode")

S1C <- bb_genebubbles(
  cds_main,
  genes = c("CD14", "CD3E", "CD79A", "MS4A1", "CD19","CD8"
 ),
  cell_grouping = c("partition", "partition_assignment_1"),
  return_value = "data"
) |> 
  ggplot(mapping = aes(x = partition, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~partition_assignment_1, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")

#Supp Fig - Nadeu et al RT UP aggregate gene expression mapping 
#Nadeu et al RT UP aggregate gene expression mapping 
# bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
# ) + 
#   facet_grid(row = vars(patient), col = (vars(disease_tissue)))

S1DP1 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
) + 
  facet_wrap(~disease_tissue)+
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
S1DP2 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_CLL_gene)
) +
  facet_wrap(~disease_tissue)+
  theme(strip.text = element_blank()) +
   theme(
  #   axis.line.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   axis.text.y = element_blank(),
     axis.title.x = element_blank(),
     axis.title.y = element_blank()
   ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*CLL Gene*<br>Expression") +
  theme(legend.title = ggtext::element_markdown())
#S1D <- S1DP1/S1DP2

S1D <- as_ggplot(grid.arrange(patchworkGrob(S1DP1/S1DP2), left = S1DP1$labels$y, bottom = textGrob(S1DP1$labels$x, hjust = 0.9)))

#ggsave("suppfig_1A.pdf", path = figures_out, width = 4.95, height = 2.45)

# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main)|>filter(patient == "pt_2712")), "partition", facet_by = "disease_tissue") 

# LN_B_clst2 <-cds_main[,colData(cds_main)$partition == "2"]
# colData(LN_B_clst2)

#Top Markers Partition Cluster 2
# Fig1_tm <- monocle3::top_markers(cds_main,
#                      group_cells_by = "partition", genes_to_test_per_group = 300, cores = 10)
# write_csv(Fig1_tm, file = file.path(WalkerTables, "Fig1_300tm.csv"))


# main Fig 1

F1AP1 <-
  bb_var_umap(
    filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
    "partition_assignment_1",
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
  labs(fill = "Cluster")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1AP2 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
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
  theme(legend.title = ggtext::element_markdown())
F1AP3 <- bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
                      gene_or_genes = "PRMT5") +
  facet_wrap( ~ disease_tissue) +
  theme(strip.text = element_blank()) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*PRMT5*") +
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
  
F1A <- F1AP1/F1AP2/F1AP3





# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "partition_assignment_1"
# ) + facet_wrap(~disease_tissue)
# 
# bb_var_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_1245")), "partition_assignment_1"
# ) + facet_wrap(~disease_tissue)
# 
#  bb_gene_umap(
#    filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")), "PRMT5"
#  ) + facet_grid(row = vars(patient), col = (vars(disease_tissue)))

#Gene dot plot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
             colData(cds_main)$clonotype_id %in% "clonotype1" &
             colData(cds_main)$partition_assignment_1 %in% "B"],
  markers = c("PRMT5", "MYC", "MKI67"),
  group_cells_by = "disease_tissue",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)

# bb_genebubbles(
#   obj = filter_cds(
#     cds_main,
#     cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")
#     ),
#   genes = c("MYC", "PRMT5", "MKI67"),
#   cell_grouping = "partition_assignment_1") + facet_wrap(~disease_tissue)


#Fig1D
F1D0 <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
                    "partition_assignment_1", facet_by = "f1d_lbl") +
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
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
F1D4 <- bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
  gene_or_genes = c("MKI67")
) + facet_wrap(~f1d_lbl)+
  theme(strip.text = element_blank()) +
  theme(
    # axis.line.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(panel.background = element_rect(color = "black")) +
  labs(color = "*MKI67*") +
  theme(legend.title = ggtext::element_markdown())+ 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# F1D5<- bb_gene_umap(
#   filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")), 
#   gene_or_genes = c("BIRC5")
# ) + facet_wrap(~f1d_lbl)+
#   theme(strip.text = element_blank()) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )+
#   theme(panel.background = element_rect(color = "black")) +
#   labs(color = "*BIRC5*") +
#   theme(legend.title = ggtext::element_markdown())

F1D <-
  as_ggplot(grid.arrange(
    patchworkGrob(F1D0 / F1D1 / F1D2 / F1D3 / F1D4),
    left = textGrob(F1D3$labels$y, rot=90, vjust = 1.5, hjust=0.35),
    bottom = textGrob(F1D3$labels$x, hjust = 0.8, vjust = -1.2))
  )
