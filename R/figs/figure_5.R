source("R/dependencies.R")
source("R/configs.R")
#Figure 5E
F5E1 <-
  bb_var_umap(
    mouse_cds_list[[2]],
    "k_10_assignment",
    alt_dim_x = "aggr_UMAP_1",
    alt_dim_y = "aggr_UMAP_2",
    overwrite_labels = T,
    facet_by = "genotype"
  )+ labs(x = "UMAP 1", y = "UMAP 2") 

F5E2 <- bb_var_umap(
    mouse_cds_list[[2]],
    "density",
    facet_by = "genotype",
    alt_dim_x = "aggr_UMAP_1",
    alt_dim_y = "aggr_UMAP_2"
  )+ labs(x = "UMAP 1", y = "UMAP 2")
F5E <- F5E1/F5E2
F5E

#Figure 5F
F5F<- 
  bb_var_umap(mouse_cds_list[[2]], var = "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) +
  facet_wrap(~genotype)+labs(x = "UMAP 1", y = "UMAP 2")
F5F

#Figure 5G: Heatmap
#top markers
F5_k10_Top50_tm <-
  monocle3::top_markers(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(kmeans_10_harmonized %in% c("5.1", "5.3", "5.5", "5.6", "5.9"))
    ),
    group_cells_by = "kmeans_10_harmonized",
    genes_to_test_per_group = 50,
    cores = 10
  )

markers <- F5_k10_Top50_tm |> 
  filter(cell_group %in% c("5.1","5.6")) |> 
  pull(gene_short_name)

fig5_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[2]], 
                                          cells = bb_cellmeta(mouse_cds_list[[2]]) |> 
                                            filter(kmeans_10_harmonized %in% c("5.1", "5.3", "5.5", "5.6", "5.8","5.9")),
                                          genes = bb_rowmeta(mouse_cds_list[[2]]) |> 
                                            filter(gene_short_name %in% markers)), 
                         cell_group_df = bb_cellmeta(mouse_cds_list[[2]]) |> 
                           select(cell_id, kmeans_10_harmonized)) |> 
  t() |> 
  scale() |> 
  t()

rownames(fig5_mat) <- tibble(feature_id = rownames(fig5_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[2]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

fig5_colfun = circlize::colorRamp2(breaks = c(min(fig5_mat),
                                              0,
                                              max(fig5_mat)),
                                   colors = heatmap_3_colors)

highlights <-
  c(
    "Tubb5",
    "Tuba1b",
    "Npm1",
    "Marcksl1",
    "Actb",
    "Gapdh",
    "Napsa",
    "Cd37",
    "Ldha",
    "Cd19",
    "Coro1a",
    "Klf2",
    "Eif5a",
    "Hnrnpk",
    "Rac2",
    "Ms4a1",
    "Rpl3"
  )

fig5_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig5_mat) %in% highlights),
  labels = rownames(fig5_mat)[rownames(fig5_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2.5
))

F5G <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(fig5_mat, 
                                                  col = fig5_colfun, 
                                                  name = "Expression", 
                                                  show_row_names = F, 
                                                  right_annotation = fig5_anno,
                                                  row_dend_width = unit(4, "mm"),
                                                  column_dend_height = unit(4, "mm"),
                                                  heatmap_legend_param = list(legend_direction = "vertical",
                                                                              #legend_width = unit(1, "mm"),
                                                                              title_position = "lefttop-rot", 
                                                                              title_gp = gpar(fontsize = 6)
                                                  ))))

F5FG <- F5F / F5G + plot_layout(heights = c(1, 1.75))
F5FG

#Figure 5H
F5H_plotlist <- map(.x = c("Myc", "Mki67", "Egr1", "Cxcr5", "Ccr7", "Il10", "Ctla4", "Cd274"),
                    .f = \(x, dat = mouse_cds_list[[2]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25 #adjusted cell size - default is 0.5
                      ) +
                        scale_color_distiller(palette = "Oranges", #?scale_color_distiller(...)
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) + #fix scale limit 
                        facet_wrap( ~ genotype) + scale_y_continuous(breaks = c(-10,0, 10)) + 
                        scale_x_continuous(breaks = c(-5,5, 15))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white")) +
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        #labs(x = NULL, y = x) +
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(strip.text = element_blank()) + theme(legend.position = "none") #+
                        #theme(legend.title = element_blank()) 
                      #if (x != "Myc") p <- p + theme(strip.text = element_blank())
                      p 
                    })

F5H1<-
  F5H_plotlist[[1]]/ plot_spacer() / bb_gene_violinplot(
    filter_cds(mouse_cds_list[[2]],
               cells = bb_cellmeta(mouse_cds_list[[2]]) |>
                 filter(k_10_assignment == "B")),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1 , 1))

F5H2 <- F5H_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Mki67",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1, -0.1 ,1))
  
F5H3 <-F5H_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Egr1",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1, -0.1 ,1))
  
F5H4 <- (F5H_plotlist[[4]]+theme(legend.position = "right",
                               legend.key.size = unit(4, "mm"),
                               legend.margin = margin(c(0, -7, 0, -6.5)),
                               legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cxcr5",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

F5H5 <-
  F5H_plotlist[[5]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Ccr7",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1, 1))
F5H6 <- F5H_plotlist[[6]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
F5H7 <- F5H_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Ctla4",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1, 1))
F5H8 <- F5H_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[2]],
      cells = bb_cellmeta(mouse_cds_list[[2]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cd274",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))

F5H_1 <- 
  F5H1 |
  F5H2 |
  F5H3 |
  F5H4

F5H_2 <- 
  F5H5 |
  F5H6 |
  F5H7 |
  F5H8
F5H_1
F5H_2

#Pseudobulk - size effect (log2FC on Figure 5H violin plots)
f5_k10_B_cds<- filter_cds(mouse_cds_list[[2]],
                          cells = bb_cellmeta(mouse_cds_list[[2]]) |>
                            filter(k_10_assignment == "B"))
unique(colData(f5_k10_B_cds)$kmeans_10_harmonized)

exp_design <- 
  bb_cellmeta(f5_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()

rowData(f5_k10_B_cds)$id <- rownames(rowData(f5_k10_B_cds))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = f5_k10_B_cds,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PRMT5/TCL1", "TCL1"))

pseudobulk_res$Result |> filter(gene_short_name == "Myc") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Mki67") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Egr1") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cxcr5")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Ccr7")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Il10")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Ctla4")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cd274")|> select(id, gene_short_name, log2FoldChange)
