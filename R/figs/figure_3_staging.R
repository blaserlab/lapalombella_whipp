colData(mouse_cds_list[[1]])$kmeans_10_harmonized <- recode(colData(mouse_cds_list[[1]])$kmeans_10_clusters, 
                                                       "1" = "3.1",
                                                       "2" = "3.2",
                                                       "3" = "3.3", 
                                                       "4" = "3.4",
                                                       "5" = "3.5",
                                                       "6" = "3.6", 
                                                       "7" = "3.7",
                                                       "8" = "3.8",
                                                       "9" = "3.9", 
                                                       "10" = "3.10")

colData(mouse_cds_list[[1]])$kmeans_10_harmonized <- factor(colData(mouse_cds_list[[1]])$kmeans_10_harmonized, 
                                                            levels = paste0("3.", 1:10))

colData(mouse_cds_list[[1]])$cd19_cd5_label

colData(mouse_cds_list[[1]])$cd19_cd5_pos
mouse_cds_list[[1]] <- bb_cellmeta(mouse_cds_list[[1]]) |> 
  filter(cd19_cd5_pos) |> 
  select(cell_id) |> 
  mutate(cd19_cd5_label = "CD19+/CD5+ cells") |> 
  bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = _)


# figure 3A
bb_var_umap(mouse_cds_list[[1]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")

bb_var_umap(mouse_cds_list[[1]], "k_10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T)

# supplemental figure
bb_genebubbles(
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
  
          


bb_var_umap(mouse_cds_list[[1]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T)


# figure 3B
bb_var_umap(
  mouse_cds_list[[1]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  labs(x = "UMAP 1", y= "UMAP 2") +
  theme(legend.justification = "center")
# figure 3C
fig3c_plotlist <- map(.x = c("Ighm", "Pax5"),
                      .f = \(x, dat = mouse_cds_list[[1]]) {
                        p <- bb_gene_umap(
                          dat,
                          gene_or_genes = x,
                          alt_dim_x = "aggr_UMAP_1",
                          alt_dim_y = "aggr_UMAP_2"
                        ) +
                          scale_color_distiller(palette = "Oranges",
                                                direction = 1,
                                                na.value = "grey80") +
                          facet_wrap( ~ genotype) +
                          theme(panel.background = element_rect(color = "black")) +
                          theme(axis.line = element_blank()) +
                          theme(axis.ticks = element_blank()) +
                          theme(axis.text = element_blank()) +
                          labs(x = NULL, y = x) +
                          theme(axis.title.y = element_text(face = "italic")) +
                          theme(legend.position = "none")
                        if (x != "Ighm") p <- p + theme(strip.text = element_blank())
                        p 
                      })

fig3c_plotlist[[1]]/fig3c_plotlist[[2]]  
  


# figure 3D



markers <- fig3_kmeans_10_tm |> 
  filter(cell_group %in% c("2", "3", "4", "6", "8")) |> 
  pull(gene_short_name)

# write_csv(fig3_kmeans_10_tm, file = file.path(tables_out, "fig3_kmeans_10_tm.csv"))

fig3_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[1]], 
                              cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                filter(kmeans_10_clusters %in% c("2", "3", "4", "6", "8")),
                              genes = bb_rowmeta(mouse_cds_list[[1]]) |> 
                                filter(gene_short_name %in% markers)), 
                         cell_group_df = bb_cellmeta(mouse_cds_list[[1]]) |> 
                           select(cell_id, kmeans_10_harmonized)) |> 
  t() |> 
  scale() |> 
  t()

rownames(fig3_mat) <- tibble(feature_id = rownames(fig3_mat)) |> 
  left_join(bb_rowmeta(mouse_cds_list[[1]]) |> 
              select(feature_id, gene_short_name)) |> 
  pull(gene_short_name)

fig3_colfun = circlize::colorRamp2(breaks = c(min(fig3_mat),
                                              0,
                                              max(fig3_mat)),
                                   colors = heatmap_3_colors)

highlights <- c("Wnt10a")

fig3_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig3_mat) %in% highlights),
  labels = rownames(fig3_mat)[rownames(fig3_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 2
))

ComplexHeatmap::Heatmap(fig3_mat, 
                        col = fig3_colfun, 
                        name = "Expression", show_row_names = F, right_annotation = fig3_anno)
