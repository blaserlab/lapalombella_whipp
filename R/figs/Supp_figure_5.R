source("R/dependencies.R")
source("R/configs.R")

#Supplemental Figure 5A
S5A1<-bb_var_umap(mouse_cds_list[[4]], "k10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")+
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank())

S5A2 <- bb_var_umap(mouse_cds_list[[4]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2") + 
  theme(strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.key.size = unit(3, 'mm'))

S5A <- (S5A1/S5A2)
S5A <- grid.arrange(patchworkGrob(S5A1/S5A2), left = textGrob("UMAP 2", rot = 90, vjust = 1.5), bottom = textGrob("UMAP 1", vjust = -1))

#Supplemental Figure 5B
S5B <- 
  bb_genebubbles(
    mouse_cds_list[[4]],
    genes = c(
      "Ms4a1",
      "Cd79a",
      "Cd19",
      "Cd3d",
      "Cd4",
      "Cd14",
      "Itgam",
      "Cd8a",
      "Cd177",
      "Pdcd1",
      "Foxp3"
    ),
    cell_grouping = c("kmeans_10_harmonized", "k10_assignment"),
    return_value = "data"
  ) |> 
  ggplot(mapping = aes(x = kmeans_10_harmonized, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~k10_assignment, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
S5B

#Supplemental Figure 5C
S5C <- bb_var_umap(
  mouse_cds_list[[4]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  labs(y = "UMAP 2", x = "UMAP 1") +
  theme(legend.justification = "center")
S5C

#Supplemental Figure 5D
S5D_plotlist <- map(.x = c("Myc","Ly6a","Mki67","Ccr7","Cxcr5","Il10","Ctla4","Egr1","Cd274","Cd93"),
                    .f = \(x, dat = mouse_cds_list[[4]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) +
                        facet_wrap( ~ genotype) + scale_x_continuous(breaks = c(-10,0, 10)) + 
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white")) +
                        theme(axis.line = element_blank()) + 
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(legend.position = "none") + theme(strip.text = element_blank())
                      p 
                    })

S5D1 <-
  (S5D_plotlist[[1]]/plot_spacer()/ bb_gene_violinplot(        #+theme(axis.text.y = element_text())
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression")) + plot_layout(heights = c(1, -0.1 ,1))

S5D2 <- S5D_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                      cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                        filter(k10_assignment == "B")), variable = "genotype",
                                                           genes_to_plot = "Ly6a", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D3 <- S5D_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                      cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                        filter(k10_assignment == "B")), variable = "genotype",
                                                           genes_to_plot = "Mki67", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D4 <-S5D_plotlist[[4]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                       filter(k10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Ccr7", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D5 <-(S5D_plotlist[[5]]+theme(legend.position = "right",
                                legend.key.size = unit(4, "mm"),
                                legend.margin = margin(c(0, -7, 0, -6.5)),
                                legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                                                                             cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                                                                               filter(k10_assignment == "B")), variable = "genotype",
                                                                                                                  genes_to_plot = "Cxcr5", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

S5D6 <-
  (S5D_plotlist[[6]]+theme(axis.text.y = element_text())) / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[4]],
      cells = bb_cellmeta(mouse_cds_list[[4]]) |>
        filter(k10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1 , 1))
S5D7 <-S5D_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                       filter(k10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Ctla4", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D8 <-S5D_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                       filter(k10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Egr1", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D9 <-S5D_plotlist[[9]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                       filter(k10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Cd274", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S5D10 <-(S5D_plotlist[[10]])/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[4]], 
                                                                         cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                                                           filter(k10_assignment == "B")), variable = "genotype",
                                                              genes_to_plot = "Cd93", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

S5D_1 <- S5D1|S5D2|S5D3|S5D4|S5D5
S5D_1
S5D_2 <-S5D6|S5D7|S5D8|S5D9|S5D10
S5D_2

#Pseudobulk - size effect (log2FC on Supp Figure 5D violin plots)
s5_k10_B_cds<- filter_cds(mouse_cds_list[[4]],
                          cells = bb_cellmeta(mouse_cds_list[[4]]) |>
                            filter(k10_assignment == "B"))

exp_design <- 
  bb_cellmeta(s5_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()

rowData(s5_k10_B_cds)$id <- rownames(rowData(s5_k10_B_cds))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = s5_k10_B_cds,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PRMT5/TCL1", "TCL1"))

pseudobulk_res$Result |> filter(gene_short_name == "Myc") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Ly6a") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Mki67") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Ccr7")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cxcr5")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Il10")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Ctla4")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Egr1") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cd274")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cd93")|> select(id, gene_short_name, log2FoldChange)

#Supplemental Figure 5E: Heatmap
#Top Markers
S5_LN_Bclust_Top50markers<- monocle3::top_markers(
  filter_cds(
    mouse_cds_list[[4]],
    cells = bb_cellmeta(mouse_cds_list[[4]]) |> filter(k10_assignment == "B")
  ),
  group_cells_by = "kmeans_10_harmonized",
  genes_to_test_per_group = 50,
  cores = 12
)

#####S5C -Heatmap -RT Clust2&7 and TCL1 Clust4
S5E1<-bb_var_umap(mouse_cds_list[[4]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")+labs(x="UMAP 1", y="UMAP 2")
S5E1  
markers <- S5_LN_Bclust_Top50markers |> 
  filter(cell_group %in% c("LN.2", "LN.7","LN.4")) |> 
  pull(gene_short_name)

S5_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[4]], 
                                        cells = bb_cellmeta(mouse_cds_list[[4]]) |> 
                                          filter(kmeans_10_harmonized %in% c("LN.2","LN.4","LN.5","LN.7")),
                                        genes = bb_rowmeta(mouse_cds_list[[4]]) |> 
                                          filter(gene_short_name %in% markers)), 
                       cell_group_df = bb_cellmeta(mouse_cds_list[[4]]) |> 
                         select(cell_id, kmeans_10_harmonized)) |> 
  t() |> 
  scale() |> 
  t()

rownames(S5_mat) <- tibble(feature_id = rownames(S5_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[4]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

S5_colfun = circlize::colorRamp2(breaks = c(min(S5_mat),
                                            0,
                                            max(S5_mat)),
                                 colors = heatmap_3_colors)

highlights <-
  c(
    "Rac2",
    "Coro1a",
    "Ms4a1",
    "Cd37",
    "Actb",
    "Cd83",
    "Cd79b",
    "Napsa",
    "Cks1b",
    "Cenpm",
    "Mki67",
    "Ube2c",
    "Cdk1",
    "Birc5",
    "Hmgb1",
    "Tubb5",
    "Smc2",
    "Aurkb",
    "Stmn1",
    "Top2a",
    "Hmgb2",
    "Tk1",
    "Tyms",
    "Tuba1b"
  )

S5_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(S5_mat) %in% highlights),
  labels = rownames(S5_mat)[rownames(S5_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2.5
))

S5E2 <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(S5_mat, 
                                                   col = S5_colfun, 
                                                   name = "Expression", 
                                                   show_row_names = F, 
                                                   right_annotation = S5_anno,
                                                   row_dend_width = unit(4, "mm"),
                                                   column_dend_height = unit(4, "mm"),
                                                   heatmap_legend_param = list(legend_direction = "vertical",
                                                                               #legend_width = unit(1, "mm"),
                                                                               title_position = "lefttop-rot", 
                                                                               title_gp = gpar(fontsize = 6)
                                                   ))))

S5E <- S5E1 / S5E2 + plot_layout(heights = c(1, 2))
S5E
