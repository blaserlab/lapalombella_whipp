WalkerAccess <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs"
WalkerTables <- "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables"

#LN
colData(mouse_cds_list[[4]])$genotype <- recode(colData(mouse_cds_list[[4]])$genotype,
                                                "PRMT5" = "PRMT5/TCL1",
                                                "TCL1" = "TCL1",
                                                "P/T" = "PRMT5/TCL1")
unique(mouse_cds_list[[4]]$genotype)

S5A2 <- bb_var_umap(mouse_cds_list[[4]], "density", facet_by = "genotype", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
cds_ln <- mouse_cds_list[[4]]

LN_fig_marker_genebubbles_supp <- bb_genebubbles(
  obj = filter_cds(cds_ln, cells = bb_cellmeta(cds_ln)),
  genes = c("Cd19", "Ms4a1", "Cd79a", "Cd3d","Cd4","Cd14","Itgam","Cd8a","Cd177","Pdcd1","Foxp3"), cell_grouping = "graphclust")

colData(mouse_cds_list[[4]])$kmeans_10_harmonized <- recode(colData(mouse_cds_list[[4]])$kmeans_10_clusters, 
                                                            "1" = "LN.1",
                                                            "2" = "LN.2",
                                                            "3" = "LN.3", 
                                                            "4" = "LN.4",
                                                            "5" = "LN.5",
                                                            "6" = "LN.6", 
                                                            "7" = "LN.7",
                                                            "8" = "LN.8",
                                                            "9" = "LN.9", 
                                                            "10" = "LN.10")

colData(mouse_cds_list[[2]])$kmeans_10_harmonized <- factor(colData(mouse_cds_list[[2]])$kmeans_10_harmonized, 
                                                            levels = paste0("5.", 1:10))
#ggsave("fig5_graphbased_genebubbles.pdf", path = WalkerAccess, width = 8, height = 5.5)

# Number of ea genotype
unique(colData(mouse_cds_list[[4]])$genotype)
length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "TCL1"])) #n=3
length(unique(colData(mouse_cds_list[[4]])$mouse[colData(mouse_cds_list[[4]])$genotype == "PRMT5/TCL1"])) #n=4

#TODO switch to kmeans
colData(cds_ln)$gbased_assignment <- recode(colData(cds_ln)$graphclust, "1" = "Cd8+ T", "2" = "B", "3" = "Cd8+ T",
                                             "4" = "Cd8+ T", "5" = "Cd8+ T", "6" = "B","7" = "B", "8" = "Cd4+ T",
                                             "9" = "B", "10" = "Cd8+ T", "11" = "B", "12" = "B",
                                             "13" = "B", "14" = "B", "15" = "B", "16" = "B", "17" = "B",
                                             "18" = "Cd4+ T", "19" = "Monocytes", "20" = "Neutrophils", "21" = "T")
colData(cds_ln)
bb_var_umap(cds_ln, "gbased_assignment", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F, facet_by = "genotype")

bb_var_umap(cds_ln, "graphclust", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F, facet_by = "genotype")
#ggsave(".pdf", path = WalkerAccess, width = 10, height = 3.9)

bb_var_umap(cds_ln, "kmeans_10_clusters", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = F, facet_by = "genotype")
colData(cds_ln)
bb_gene_umap(cds_ln, gene_or_genes = "Mki67", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)
bb_gene_umap(cds_ln, gene_or_genes = "Myc", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+facet_wrap(~genotype)

#Top Markers of all B cells by genotype (top 100)
# LnB<-cds_ln[,colData(cds_ln)$gbased_assignment == "B"]
# 
# monocle3::top_markers(LnB, group_cells_by = "genotype", genes_to_test_per_group = 100, cores = 10)
