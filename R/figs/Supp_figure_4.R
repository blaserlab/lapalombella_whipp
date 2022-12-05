#Supplemental Figure 4E
S4E <- 
  bb_genebubbles(
    mouse_cds_list[[2]],
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
S4E

#Supplemental Figure 4F
S4F <- bb_var_umap(
  mouse_cds_list[[2]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +labs(x = "UMAP 1", y= "UMAP 2") +
  theme(legend.justification = "center")
S4F

#Supplemental Figure 4G
S4G_plotlist <- map(.x = c("Cd93", "Il4"),
                    .f = \(x, dat = mouse_cds_list[[2]]) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.25 #adjusted cell size - default is 0.5
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,1)) + #set fixed scale (limits)
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
                        theme(strip.text = element_blank()) + theme(legend.position = "none")
                      p 
                    })

S4G1 <- S4G_plotlist[[1]]/plot_spacer()/bb_gene_violinplot(
  filter_cds(
    mouse_cds_list[[2]],
    cells = bb_cellmeta(mouse_cds_list[[2]]) |>
      filter(k_10_assignment == "B")),
  variable = "genotype",
  genes_to_plot = "Cd93",
  pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE
) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
  t = 0,
  r = -1.5,
  b = 0,
  l = -7
))) + labs(y = "B cell Expression") +plot_layout(heights = c(1,-0.1,1)) #+ stat_compare_means(method = "t.test")
S4G2 <- (S4G_plotlist[[2]]+theme(legend.position = "right",
                                 legend.key.size = unit(4, "mm"),
                                 legend.margin = margin(c(0, -7, 0, -6.5)),
                                 legend.title = element_blank(), legend.text = ))/plot_spacer()/bb_gene_violinplot(
                                   filter_cds(
                                     mouse_cds_list[[2]],
                                     cells = bb_cellmeta(mouse_cds_list[[2]]) |>
                                       filter(k_10_assignment == "B")
                                   ),
                                   variable = "genotype",
                                   genes_to_plot = "Il4",
                                   pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) +plot_layout(heights = c(1,-0.1,1))

S4G<- S4G1|S4G2
S4G
