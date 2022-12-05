#Supplemental Figure 2E
S2E <- bb_genebubbles(
  mouse_cds_list[[1]],
  genes = c(
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
  cell_grouping = c("kmeans_10_harmonized", "k_10_assignment"),
  return_value = "data"
) |> 
  ggplot(mapping = aes(x = kmeans_10_harmonized, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~k_10_assignment, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
S2E

#Supplemental Figure 2F
S2F_plotlist <- map(.x = c("Ccr7", "Il4", "Cd69", "Cd93", "Cxcr5", "Myc", "Il10", "Mki67", "Npm1"),
                    .f = \(x, dat = mouse_cds_list[[1]]) {
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
                        #theme(axis.text.y = element_blank()) +
                        labs(x = NULL, y = NULL, title = x) +
                        theme(axis.title.y = element_text(face = "italic")) +
                        theme(legend.position = "none") + theme(strip.text = element_blank())
                      p 
                    })
S2F1 <-
  (S2F_plotlist[[1]]+theme(axis.text.y = element_text()))/plot_spacer()/ bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Ccr7",
    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1, -0.1 ,1))

S2F2 <-
  S2F_plotlist[[2]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il4",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
S2F3 <-
  S2F_plotlist[[3]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cd69",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
S2F4 <- (
  S2F_plotlist[[4]] + theme(
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.margin = margin(c(0,-7, 0,-6.5)),
    legend.title = element_blank()
  )
) / plot_spacer() / bb_gene_violinplot(
  filter_cds(
    mouse_cds_list[[1]],
    cells = bb_cellmeta(mouse_cds_list[[1]]) |>
      filter(k_10_assignment == "B")
  ),
  variable = "genotype",
  genes_to_plot = "Cd93",
  pseudocount = 0,
  jitter_fill = "transparent",
  violin_alpha = 0.55,
  jitter_alpha = .1,
  include_jitter = TRUE
) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))

S2F5 <-
  (S2F_plotlist[[5]] + theme(axis.text.y = element_text())) / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Cxcr5",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
    t = 0,
    r = -1.5,
    b = 0,
    l = -7
  ))) + labs(y = "B cell Expression") + plot_layout(heights = c(1,-0.1 , 1))

S2F6 <-
  S2F_plotlist[[6]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Myc",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
S2F7 <-
  S2F_plotlist[[7]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Il10",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
S2F8 <-
  S2F_plotlist[[8]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Mki67",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))
S2F9 <-
  S2F_plotlist[[9]] / plot_spacer() / bb_gene_violinplot(
    filter_cds(
      mouse_cds_list[[1]],
      cells = bb_cellmeta(mouse_cds_list[[1]]) |>
        filter(k_10_assignment == "B")
    ),
    variable = "genotype",
    genes_to_plot = "Npm1",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = .1,
    include_jitter = TRUE
  ) + theme(strip.text = element_blank()) + theme(axis.title.y = element_blank()) + plot_layout(heights = c(1,-0.1 , 1))

S2F_1 <- S2F1 | S2F2 | S2F3 | S2F4
S2F_1
S2F_2 <- S2F5 | S2F6 | S2F7 | S2F8 | S2F9
S2F_2

#Pseudobulk - size effect (log2FC on Supplemantal Figure 2F:violin plots)
s2_k10_B_cds<- filter_cds(mouse_cds_list[[1]],
                          cells = bb_cellmeta(mouse_cds_list[[1]]) |>
                            filter(k_10_assignment == "B"))

exp_design <- 
  bb_cellmeta(s2_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()

rowData(s2_k10_B_cds)$id <- rownames(rowData(s2_k10_B_cds))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = s2_k10_B_cds,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PRMT5", "TCL1"))

pseudobulk_res$Result |> filter(gene_short_name == "Ccr7") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Il4") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cd69") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cd93") |> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Cxcr5")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Myc")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Il10")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Mki67")|> select(id, gene_short_name, log2FoldChange)
pseudobulk_res$Result |> filter(gene_short_name == "Npm1")|> select(id, gene_short_name, log2FoldChange)
