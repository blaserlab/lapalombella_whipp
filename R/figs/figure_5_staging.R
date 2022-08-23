#mouse_cds_list
  #File1 - Spleen - PRMT5 vs TCL1
  #File2 - Spleen - PRMT5xTCL1 vs TCL1
  #File3 - LN - PRMT5 vs TCL1
  #File4 - LN - PRMT5xTCL1 vs TCL1

unique(mouse_cds_list[[2]]$mouse)
#Spleens
####PRMT5xTCL1:M0955-RT, M0942, M1040
####TCL1:M0244, M0228, M0229, M0322
####PRMT5: M0980
#Mistake, recode M0980 as Eu-PRMT5/TCL1
colData(mouse_cds_list[[2]])$genotype <- recode(colData(mouse_cds_list[[2]])$genotype,
                                                "PRMT5" = "Eμ-PRMT5/TCL1",
                                                "TCL1" = "Eμ-TCL1",
                                                "P/T" = "Eμ-PRMT5/TCL1")

unique(mouse_cds_list[[4]]$mouse)
#LN
###Mice:"M0244"    "M0955-RT" "M0228"   "M0322"    "M0942"    "M1040"    "M0980"
#same mice besides M0229

#Fig5 Figs
colData(mouse_cds_list[[2]])
unique(mouse_cds_list[[2]]$genotype)
cds_sub <- mouse_cds_list[[2]]

#filter for PRMT5xTCL1 & TCL1 mice only
cds_sub <- cds_sub[,colData(cds_sub)$genotype == "Eμ-PRMT5/TCL1" |
                     colData(cds_sub)$genotype == "Eμ-TCL1"]
unique(cds_sub$genotype)
unique(colData(cds_sub)$mouse)

#recode clusters
colData(mouse_cds_list[[2]])$genotype <- recode(colData(mouse_cds_list[[2]])$genotype,
                                                "PRMT5" = "Eμ-PRMT5/TCL1",
                                                "TCL1" = "Eμ-TCL1",
                                                "P/T" = "Eμ-PRMT5/TCL1")
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

bb_var_umap(cds_sub, var = "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T) +
  facet_wrap(~genotype)

#cluster assignment
bb_genebubbles(
  obj = filter_cds(cds_sub, cells = bb_cellmeta(cds_sub)),
  genes = c("Cd19", "Ms4a1", "Cd79a", "Cd3d","Cd4","Cd14","Itgam","Cd8a","Cd177","Pdcd1","Foxp3"), cell_grouping = "k_10_assigment") 

#check if assigned 
colData(cds_sub)$k_10_assignment <- recode(colData(cds_sub)$kmeans_10_harmonized, "5.1" = "B", "5.2" = "Cd8+ T", "5.3" = "B", "5.4" = "Cd4+ T", "5.5" = "B", "5.6" = "B", "5.7" = "Neutrophils", "5.8" = "Low Quality", "5.9" = "B", "5.10" = "Monocytes")
colData(cds_sub)

bb_genebubbles(
  cds_sub,
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

bb_var_umap(cds_sub, var = "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", facet_by = "genotype", overwrite_labels = T)
bb_var_umap(cds_sub, var = "density", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", facet_by = "genotype")

bb_var_umap(mouse_cds_list[[2]], "kmeans_10_harmonized", facet_by = "mouse", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
