#mouse_cds_list
  #File1 - Spleen - PRMT5 vs TCL1
  #File2 - Spleen - PRMT5xTCL1 vs TCL1
  #File3 - LN - PRMT5 vs TCL1
  #File4 - LN - PRMT5xTCL1 vs TCL1

# unique(mouse_cds_list[[2]]$mouse)
# #Spleens
# ####PRMT5xTCL1:M0955-RT, M0942, M1040
# ####TCL1:M0244, M0228, M0229, M0322
# ####PRMT5: M0980
# #Mistake, recode M0980 as Eu-PRMT5/TCL1
colData(mouse_cds_list[[2]])$genotype <- recode(colData(mouse_cds_list[[2]])$genotype,
                                                "PRMT5" = "PRMT5/TCL1",
                                                "TCL1" = "TCL1",
                                                "P/T" = "PRMT5/TCL1")
unique(mouse_cds_list[[2]]$genotype)

# unique(mouse_cds_list[[4]]$mouse)
# #LN
# ###Mice:"M0244"    "M0955-RT" "M0228"   "M0322"    "M0942"    "M1040"    "M0980"
# #same mice besides M0229

#Fig5 Figs
colData(mouse_cds_list[[2]])
# cds_sub <- mouse_cds_list[[2]]
# cds_sub <- cds_sub[,colData(cds_sub)$genotype == "Eμ-PRMT5/TCL1" |
#                      colData(cds_sub)$genotype == "Eμ-TCL1"]
# unique(cds_sub$genotype)
# unique(colData(cds_sub)$mouse)


#filter for PRMT5xTCL1 & TCL1 mice only


#recode/harmoize clusters
colData(mouse_cds_list[[2]])$kmeans_10_harmonized <- recode(colData(mouse_cds_list[[2]])$kmeans_10_clusters, 
                                                            "1" = "5.1",
                                                            "2" = "5.2",
                                                            "3" = "5.3", 
                                                            "4" = "5.4",
                                                            "5" = "5.5",
                                                            "6" = "5.6", 
                                                            "7" = "5.7",
                                                            "8" = "5.8",
                                                            "9" = "5.9", 
                                                            "10" = "5.10")

colData(mouse_cds_list[[2]])$kmeans_10_harmonized <- factor(colData(mouse_cds_list[[2]])$kmeans_10_harmonized, 
                                                            levels = paste0("5.", 1:10))

fig5_k10_harmonized<- bb_var_umap(mouse_cds_list[[2]], var = "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) +
  facet_wrap(~genotype)
#ggsave("fig5_k10_harmonized.pdf", path = figures_out, width = 7.8, height = 4.7)
fig5_density_supp <- bb_var_umap(mouse_cds_list[[2]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
#ggsave("fig5_density_supp.pdf", path = figures_out, width = 7.8, height = 4.7)

#celltype calling
cds_sub <- mouse_cds_list[[2]]

# cds_sub$kmeans_10_harmonized <- factor(cds_sub$kmeans_10_harmonized,
#                                                    levels = c("5.1", "5.2", "5.3", "5.4", "5.5", "5.6","5.7","5.8","5.9","5.10"))
# levels(cds_sub$kmeans_10_harmonized) <- c("5.1", "5.2", "5.3", "5.4", "5.5", "5.6","5.7","5.8","5.9","5.10")

fig5_marker_genebubbles_supp <- bb_genebubbles(
  obj = filter_cds(cds_sub, cells = bb_cellmeta(cds_sub)),
  genes = c("Cd19", "Ms4a1", "Cd79a", "Cd3d","Cd4","Cd14","Itgam","Cd8a","Cd177","Pdcd1","Foxp3"), cell_grouping = "kmeans_10_harmonized") 
ggsave("fig5_marker_genebubbles_supp.pdf", path = figures_out, width = 4.6, height = 2.9)


#cluster assignment
cds_sub <- mouse_cds_list[[2]]
colData(cds_sub)$k_10_assignment <- recode(colData(cds_sub)$kmeans_10_harmonized, "5.1" = "B", "5.2" = "Cd8+ T", "5.3" = "B", "5.4" = "Cd4+ T", "5.5" = "B", "5.6" = "B", "5.7" = "Neutrophils", "5.8" = "B", "5.9" = "B", "5.10" = "Monocytes")

fig5d <- bb_var_umap(cds_sub, "k_10_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype")
ggsave("fig5d.pdf", path = figures_out, width = 7.8, height = 4.7, labeler = )

colData(cds_sub)

#Gene dotplot
#fig5_marker_genebubbles_supp2 <- 
bb_genebubbles(
  cds_sub,
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
#ggsave("fig5_marker_genebubbles_supp2.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs", height = 6.1, width = 5.8)

#UMAPs
#"fig5_mouse_facet_supp" <- 
  bb_var_umap(mouse_cds_list[[2]], "kmeans_10_harmonized", facet_by = "mouse", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
#ggsave("fig5_mouse_facet_supp.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs", height = 5.8, width = 5.3)


#Heatmap/topmarkers
  #top markers
fig5_kmeans_10_tm_Top100 <-monocle3::top_markers(mouse_cds_list[[2]], group_cells_by = "kmeans_10_harmonized", genes_to_test_per_group = 100, cores = 10)
fig5_kmeans_10_tm_Top100
#########################################################################################
Spln_Bcell_Top100 <- monocle3::top_markers(filter_cds(cds_sub, bb_cellmeta(cds_sub), cds_sub$k_10_assignment == "B"),
                                           group_cells_by = "kmeans_10_harmonized", genes_to_test_per_group = 100, cores = 10)


  write_csv(fig5_kmeans_10_tm_Top100, file = file.path(tables_out, "fig5_kmeans_100_tm.csv"))

markers <- fig5_kmeans_10_tm_Top100 |> 
  filter(cell_group %in% c("5.1", "5.3", "5.5", "5.6", "5.8","5.9")) |> 
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

#highlights <- c("Cd83", "Tubb5", "Tuba1b","Wnt10a", "Ccr7", "Cdk1", "Birc5","Junb","Atf4","Cd69","Nfkbia","Nfkbid","Rel","Ubc","Ccna2", "Blnk")
highlights <- c("Cd73")

fig5_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig5_mat) %in% highlights),
  labels = rownames(fig5_mat)[rownames(fig5_mat) %in% highlights],
  labels_gp = gpar(fontsize = 8),
  padding = 0
))

ComplexHeatmap::Heatmap(fig5_mat, 
                        col = fig5_colfun, 
                        name = "Expression", show_row_names = F, right_annotation = fig5_anno)

fig5_kmeans_10_tm_Top100 |> filter(cell_group == "5.6") |> arrange(desc(marker_score))
fig5_kmeans_10_tm_Top100 |> filter(cell_group == "5.6") |> arrange(desc(mean_expression))



         