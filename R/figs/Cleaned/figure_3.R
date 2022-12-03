source("R/dependencies.R")
source("R/configs.R")
#TODO incorporate cds modifications into datapkg

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

colData(mouse_cds_list[[1]])$k_10_assignment <- recode(colData(mouse_cds_list[[1]])$k_10_assignment, "Low Quality" = "Other")

mouse_cds_list[[1]] <- bb_cellmeta(mouse_cds_list[[1]]) |> 
  filter(cd19_cd5_pos) |> 
  select(cell_id) |> 
  mutate(cd19_cd5_label = "CD19+/CD5+ cells") |> 
  bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = _)
##########################################

#Figure 3A
F3A1 <- bb_var_umap(mouse_cds_list[[1]], "k_10_assignment", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) + 
theme(
     axis.line.x = element_blank(),
     axis.ticks.x = element_blank(),
     axis.title.x = element_blank(),
     axis.text.x = element_blank(),
     axis.title.y = element_blank())

F3A2 <- bb_var_umap(mouse_cds_list[[1]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2") + 
  theme(strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.key.size = unit(3, 'mm'))

F3A <- (F3A1/F3A2) 
F3A

#Figure 3B
F3B <- bb_var_umap(
  mouse_cds_list[[1]],
  "cd19_cd5_label",
  facet_by = "genotype",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  value_to_highlight = "CD19+/CD5+ cells",
  palette = "#ed718d", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) + labs(y = "UMAP 2", x = "UMAP 1") +
  theme(legend.justification = "center")
F3B

#Figure 3C
fig3c_plotlist <-
  map(
    .x = c(
      "Ighm",
      "Pax5",
      "Ighd",
      "Ighe",
      "Ebf1",
      "Cd93",
      "Cd69",
      "Spn",
      "Myc",
      "Mki67"
    ),
    .f = \(x, dat = mouse_cds_list[[1]]) {
      p <- bb_gene_umap(
        dat,
        gene_or_genes = x,
        alt_dim_x = "aggr_UMAP_1",
        alt_dim_y = "aggr_UMAP_2",
        cell_size = 0.25 #adjusted cell size - default is 0.5
      ) +
        scale_color_distiller(
          palette = "Oranges",
          direction = 1,
          na.value = "grey80",
          limits = c(0, 3)
        ) + #set fixed scale (limits)
        facet_wrap(~ genotype) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        theme(panel.background = element_rect(color = "black", fill = "white")) +
        theme(axis.line = element_blank()) +
        scale_x_continuous(breaks = c(-10, 0, 10)) + scale_y_continuous(breaks = c(-10, 0, 10)) +
        labs(x = NULL, y = x) +
        theme(axis.title.y = element_text(face = "italic")) +
        theme(legend.position = "none")
      if (x != "Ighm")
        p <- p + theme(strip.text = element_blank())
      p
    }
  )

F3C <-
  fig3c_plotlist[[1]] / plot_spacer() / fig3c_plotlist[[2]] / plot_spacer() /
  (
    fig3c_plotlist[[3]] + theme(
      legend.position = "right",
      legend.key.size = unit(4, "mm"),
      legend.margin = margin(c(0, -10, 0, -7))
    ) + theme(legend.title = element_blank())
  ) / plot_spacer() / fig3c_plotlist[[4]] / plot_spacer() /
  fig3c_plotlist[[5]] / plot_spacer() / fig3c_plotlist[[6]] / plot_spacer() /
  fig3c_plotlist[[7]] / plot_spacer() / fig3c_plotlist[[8]] / plot_spacer() /
  fig3c_plotlist[[9]] / plot_spacer() / fig3c_plotlist[[10]] + plot_layout(heights = c(
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1,-0.5,
    1
  ))
F3C
#ggsave("F3C.pdf", width = 4.5, height = 8.25)

#Figure 3D
F3D <-
  bb_var_umap(
    mouse_cds_list[[1]],
    "kmeans_10_harmonized",
    alt_dim_x = "aggr_UMAP_1",
    alt_dim_y = "aggr_UMAP_2",
    overwrite_labels = T,
    facet_by = "genotype"
  ) + labs(y = "UMAP 2", x = "UMAP 1")

#Figure 3E: Heatmap
##subset cds
k10_Bclust <- filter_cds(mouse_cds_list[[1]],
                         cells = bb_cellmeta(mouse_cds_list[[1]]) |>
                           filter(kmeans_10_clusters %in% c("2", "3", "4", "6", "8")))
##Top markers
F3_kmeans10_tm_Top50 <-
  monocle3::top_markers(
    k10_Bclust,
    group_cells_by = "kmeans_10_harmonized",
    genes_to_test_per_group = 50,
    cores = 10
  )

markers <- F3_kmeans10_tm_Top50 |> 
  filter(cell_group %in% c('3.3', '3.4', '3.6')) |> 
  pull(gene_short_name)

fig3_mat <- bb_aggregate(obj = filter_cds(mouse_cds_list[[1]], 
                              cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                filter(kmeans_10_harmonized %in% c("3.3", "3.4", "3.6","3.2","3.8")),
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

highlights <-
  c("Hmgb1",
    "Cd69",
    "Hmgb2",
    "Junb",
    "Top2a",
    "Ccr7",
    "Birc5",
    "Cd83",
    "Atf4",
    "Tuba1b",
    "Tubb5",
    "Hspa5",
    "Ccdc34",
    "Ube2c",
    "Il2rg",
    "Cks1b",
    "Smc2",
    "Stmn1",
    "Tyms",
    "Ubc",
    "Gapdh",
    "Apoe",
    "Xrcc1",
    "Cd79a",
    "Zfp36",
    "Ppp1r15a",
    "Ldha",
    "Mki67",
    "Cdk1",
    "Gpx4",
    "Btg1",
    "Itm2b",
    "Xbp1",
    "Rel",
    "Aurkb",
    "Eif5a",
    "C1qbp",
    "Cd79b",
    "Tk1",
    "Nfkbia",
    "Hsp90aa1",
    "Cenpm",
    "Nr4a1",
    "Pfdn5",
    "Eif4a2",
    "Pim1",
    "Cdca3",
    "Pclaf",
    "Ccnb2",
    "Rrm2"
  )

fig3_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig3_mat) %in% highlights),
  labels = rownames(fig3_mat)[rownames(fig3_mat) %in% highlights],
  labels_gp = gpar(fontsize = 5),
  padding = 1.2
))

F3E<- grid.grabExpr(draw(
  ComplexHeatmap::Heatmap(fig3_mat,
                          col = fig3_colfun,
                          name = "Expression", 
                          show_row_names = F, 
                          right_annotation = fig3_anno,
                          row_dend_width = unit(3, "mm"),
                          column_dend_height = unit(3, "mm"),
                          heatmap_legend_param = list(legend_direction = "vertical",
                                                      title_position = "lefttop-rot", 
                                                      title_gp = gpar(fontsize = 6)
                          ))))

F3DE <- F3D / F3E + plot_layout(heights = c(1, 5))
F3DE
