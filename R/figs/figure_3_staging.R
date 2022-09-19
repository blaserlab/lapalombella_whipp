#unique(colData(mouse_cds_list[[1]])$tissue)
#unique(colData(mouse_cds_list[[1]])$genotype)
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

#unique(colData(mouse_cds_list[[1]])$k_10_assignment)
colData(mouse_cds_list[[1]])$k_10_assignment <- recode(colData(mouse_cds_list[[1]])$k_10_assignment, "Low Quality" = "Other")

colData(mouse_cds_list[[1]])$cd19_cd5_label

colData(mouse_cds_list[[1]])$cd19_cd5_pos
mouse_cds_list[[1]] <- bb_cellmeta(mouse_cds_list[[1]]) |> 
  filter(cd19_cd5_pos) |> 
  select(cell_id) |> 
  mutate(cd19_cd5_label = "CD19+/CD5+ cells") |> 
  bb_tbl_to_coldata(mouse_cds_list[[1]], min_tbl = _)

k10_Bclust <- filter_cds(mouse_cds_list[[1]], 
                         cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                           filter(kmeans_10_clusters %in% c("2", "3", "4", "6", "8")))

# figure 3A
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
F3A <- grid.arrange(patchworkGrob(F3A1/F3A2), left = textGrob("UMAP 2", rot = 90, vjust = 1.5), bottom = textGrob("UMAP 1", vjust = -1))

# low_q<- mouse_cds_list[[1]]
# low_q <- low_q[,colData(low_q)$k_10_assignment == "Low Quality"]
# monocle3::top_markers(low_q, group_cells_by = "genotype", genes_to_test_per_group = 100, cores = 10)

#bb_var_umap(mouse_cds_list[[1]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T)


# figure 3B
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
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  labs(y = "UMAP 2", x = "UMAP 1") +
  theme(legend.justification = "center")
# figure 3C
fig3c_plotlist <- map(.x = c("Ighm","Pax5", "Ighd", "Ighe", "Ebf1", "Cd93", "Cd69", "Spn","Myc", "Mki67"),
                      .f = \(x, dat = mouse_cds_list[[1]]) {
                        p <- bb_gene_umap(
                          dat,
                          gene_or_genes = x,
                          alt_dim_x = "aggr_UMAP_1",
                          alt_dim_y = "aggr_UMAP_2", cell_size = 0.25 #adjusted cell size - default is 0.5
                        ) +
                          scale_color_distiller(palette = "Oranges",
                                                direction = 1,
                                                na.value = "grey80", limits = c(0,3)) + #set fixed scale (limits)
                          facet_wrap( ~ genotype) + 
                          theme_minimal() + 
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
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

F3C <-
fig3c_plotlist[[1]] / fig3c_plotlist[[2]]/ (fig3c_plotlist[[3]]+ theme(
  legend.position = "right",
  legend.key.size = unit(3, "mm"),
  legend.margin = margin(c(0, -10, 0, -7)) #c(top,right,bottom,left)
) + theme(legend.title = element_blank())) / fig3c_plotlist[[4]] /
  fig3c_plotlist[[5]] / fig3c_plotlist[[6]] / fig3c_plotlist[[7]] / fig3c_plotlist[[8]] /
  fig3c_plotlist[[9]] / fig3c_plotlist[[10]] 

#ggsave("fig3c.pdf", path = figures_out, width = 6, height = 15)

# figure 3D
F3D1 <- bb_var_umap(mouse_cds_list[[1]], "kmeans_10_harmonized", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2", overwrite_labels = T, facet_by = "genotype") +labs(y= "UMAP 2", x = "UMAP 1")
#F3A1/F3D1
F3_kmeans10_tm_Top50 <-monocle3::top_markers(k10_Bclust, group_cells_by = "kmeans_10_harmonized", genes_to_test_per_group = 50, cores = 10)
#write_csv(F3_kmeans10_tm_Top50, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap Tables/Fig3_PRMT5_vs_TCL1", "F3_kmeans10_tm_Top50.csv"))
#F3_kmeans10_tm_Top50 <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap Tables/Fig3_PRMT5_vs_TCL1/F3_kmeans10_tm_Top50.csv")

markers <- F3_kmeans10_tm_Top50 |> 
  filter(cell_group %in% c('3.3', '3.4', '3.6')) |> 
  pull(gene_short_name)
# write_csv(fig3_kmeans_10_tm, file = file.path(tables_out, "fig3_kmeans_10_tm.csv"))

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
#fig3_mat
fig3_colfun = circlize::colorRamp2(breaks = c(min(fig3_mat),
                                              0,
                                              max(fig3_mat)),
                                   colors = heatmap_3_colors)

#highlights <- c("Cd83", "Tubb5", "Tuba1b","Wnt10a", "Ccr7", "Cdk1", "Birc5","Junb","Atf4","Cd69","Nfkbia","Nfkbid","Rel","Ubc","Ccna2", "Blnk")
mus_paper_highlights <- c("Ccr7", "Cdk4", "Cxcr5", "Birc5", "Il4","Npm1","Jun","Junb", "Fos","Fosb","Atf3","Atf4","Myc","Cd69","Il10", "Top2a", "Hmgb1", "Hmgb2", "Cd83","Ube2a", "Tubb5", "Tuba1b", "S100a8", "S100a9")
 highlights3a <- as.vector(F3_kmeans10_tm_Top50 |> filter(F3_kmeans10_tm_Top50$gene_short_name %in% mus_paper_highlights)|> pull(gene_short_name))

 #Additional lymphoma associated genes: DisGeNet.org - Lymphoma; CUI: C0024299
lymphoma_genes<- readxl::read_excel("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/disgenet.org_LymphomaGenes_C0024299_disease_gda_summary.xlsx")
lymphoma_genes$Gene <- str_to_title(lymphoma_genes$Gene)
filt<- F3_kmeans10_tm_Top50 |> filter(F3_kmeans10_tm_Top50$gene_short_name %in% lymphoma_genes$Gene)
lymphoma_gois <- filter(filt, cell_group %in% c('3.3','3.4','3.6'))[["gene_short_name"]]
#Human top markers
RTLN_tm <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap Tables/Fig1 human RT data/leiden clustering/LN_B_leiden_Top50_tm.csv")
RTLN_tm$gene_short_name <- str_to_title(RTLN_tm$gene_short_name)
RTLN_tm <- filter(RTLN_tm, cell_group == c('3','11'))
human_overlap<- F3_kmeans10_tm_Top50 |> filter(F3_kmeans10_tm_Top50$gene_short_name %in% RTLN_tm$gene_short_name)
human_overlap <- filter(human_overlap, cell_group %in% c('3.3','3.4','3.6'))[["gene_short_name"]]
highlights <- unique(c(highlights3a,lymphoma_gois, human_overlap))


fig3_anno <- ComplexHeatmap::rowAnnotation(link =  anno_mark(
  at = which(rownames(fig3_mat) %in% highlights),
  labels = rownames(fig3_mat)[rownames(fig3_mat) %in% highlights],
  labels_gp = gpar(fontsize = 5),
  padding = 0.8
))

F3D2<- grid.grabExpr(draw(
  ComplexHeatmap::Heatmap(fig3_mat,
                          col = fig3_colfun,
                          name = "Expression", 
                          show_row_names = F, 
                          right_annotation = fig3_anno,
                          row_dend_width = unit(3, "mm"),
                          column_dend_height = unit(3, "mm"),
                          heatmap_legend_param = list(legend_direction = "vertical",
                                                      #legend_width = unit(1, "mm"),
                                                      title_position = "lefttop-rot", 
                                                      title_gp = gpar(fontsize = 6)
                          ))))
F3D <- F3D1 / F3D2 + plot_layout(heights = c(1, 3))
#F3D

# supplemental figure 2
S2_gene_dotplot <- bb_genebubbles(
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
#ggsave("S2_gene_dotplot.pdf", path = T_Figs, width = 5, height = 4.3)

S2E_plotlist <- map(.x = c("Ccr7", "Il4", "Cd69", "Cd93", "Cxcr5", "Myc", "Il10", "Mki67", "Npm1"),
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
                          labs(x = x, y = NULL) +
                          theme(axis.title.y = element_text(face = "italic")) +
                          theme(legend.position = "none") + theme(strip.text = element_blank())
                        p 
                      })

S2E1 <- S2E_plotlist[[1]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                        cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                          filter(k_10_assignment == "B")), variable = "genotype",
                                             genes_to_plot = "Ccr7", pseudocount = 0)+ theme(strip.text = element_blank())

S2E2 <- S2E_plotlist[[2]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                  cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                    filter(k_10_assignment == "B")), variable = "genotype",
                                       genes_to_plot = "Il4", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
S2E3 <- S2E_plotlist[[3]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                  cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                    filter(k_10_assignment == "B")), variable = "genotype",
                                       genes_to_plot = "Cd69", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
S2E4 <-S2E_plotlist[[4]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                  cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                    filter(k_10_assignment == "B")), variable = "genotype",
                                       genes_to_plot = "Cd93", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
# S2E_1 <- (S3E1 | S3E2 | S3E3 | S3E4)
# ggsave("S2E_1.pdf", path = T_Figs, width = 7.5, height = 9.75)


S2E5 <- S2E_plotlist[[5]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                  filter(k_10_assignment == "B")), variable = "genotype",
                                     genes_to_plot = "Cxcr5", pseudocount = 0)+ theme(strip.text = element_blank())
S2E6 <-S2E_plotlist[[6]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                  filter(k_10_assignment == "B")), variable = "genotype",
                                     genes_to_plot = "Myc", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
S2E7 <-S2E_plotlist[[7]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                  filter(k_10_assignment == "B")), variable = "genotype",
                                     genes_to_plot = "Il10", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
S2E8 <-S2E_plotlist[[8]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                  filter(k_10_assignment == "B")), variable = "genotype",
                                     genes_to_plot = "Mki67", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())
S2E9 <-S2E_plotlist[[9]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                  filter(k_10_assignment == "B")), variable = "genotype",
                                     genes_to_plot = "Npm1", pseudocount = 0)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())

# S2E_2 <- (S2E5 | S2E6 | S2E7 | S2E8 | S2E9) 
# ggsave("S2E_2.pdf", path = T_Figs, width = 7.5, height = 3)

##################################################################################################################################################
#TODO Supp 2F - Previously - Upregulated genes in 3.3, 3.4, 3.6
### Supp2F gene modules

#bb_gene_modules(obj = mouse_cds_list[[1]], n_cores = )
#F3_cds <- mouse_cds_list[[1]]

pr_graph_test_res <-
  graph_test(
    cds = mouse_cds_list[[1]],
    neighbor_graph = "knn",
    cores = 8,
    verbose = TRUE
  )
#bb_gene_umap(mouse_cds_list[[1]], gene_or_genes = "Birc5", alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2") + facet_wrap(~genotype)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

#pr_deg_ids_fix = intersect(pr_deg_ids, rownames(mouse_cds_list[[1]]@preprocess_aux$gene_loadings))
preprocess_mat <- F3_cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_v %*% diag(F3_cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_sdev)

gene_module_df <-
  monocle3::find_gene_modules(mouse_cds_list[[1]][pr_deg_ids,], cores = 8) |>
  dplyr::rename(feature_id = id)

rowData(cds_main)

##TODO Supp2F top markers
F3_kmeans10_tm_Top50 <- readr::read_csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Tables/Heatmap Tables/Fig3_PRMT5_vs_TCL1/F3_kmeans10_tm_Top50.csv")
#Cluster 3.3
query_C3.3 <- dplyr::filter(F3_kmeans10_tm_Top50, cell_group %in% '3.3')[["gene_short_name"]]
query_C3.3

C3.3_goenrichment <- bb_goenrichment(query = query_C3.3, 
                                     reference = bb_rowmeta(mouse_cds_list$file_1),
                                     go_db = "org.Mm.eg.db")
C3.3_goenrichment
C3.3_gosummary_0.8 <- bb_gosummary(x = C3.3_goenrichment, 
                                   reduce_threshold = 0.8,
                                   go_db = "org.Mm.eg.db")
C3.3_gosummary_0.9 <- bb_gosummary(x = C3.3_goenrichment, 
                                   reduce_threshold = 0.9,
                                   go_db = "org.Mm.eg.db")
C3.3_gosummary_0.98 <- bb_gosummary(x = C3.3_goenrichment, 
                                   reduce_threshold = .98,
                                   go_db = "org.Mm.eg.db")

bb_goscatter(simMatrix = C3.3_gosummary_0.8$simMatrix,
             reducedTerms = C3.3_gosummary_0.8$reducedTerms)
bb_goscatter(simMatrix = C3.3_gosummary_0.9$simMatrix,
             reducedTerms = C3.3_gosummary_0.9$reducedTerms)
bb_goscatter(simMatrix = C3.3_gosummary_0.98$simMatrix,
             reducedTerms = C3.3_gosummary_0.98$reducedTerms)
#Supp 2F Cluster 3.4
query_C3.4 <- dplyr::filter(F3_kmeans10_tm_Top50, cell_group %in% '3.4')[["gene_short_name"]]
C3.4_goenrichment <- bb_goenrichment(query = query_C3.4, 
                                     reference = bb_rowmeta(mouse_cds_list$file_1),
                                     go_db = "org.Mm.eg.db")
C3.4_goenrichment
C3.4_gosummary_0.8 <- bb_gosummary(x = C3.4_goenrichment, 
                                   reduce_threshold = 0.8,
                                   go_db = "org.Mm.eg.db")
C3.4_gosummary_0.9 <- bb_gosummary(x = C3.4_goenrichment, 
                                   reduce_threshold = 0.9,
                                   go_db = "org.Mm.eg.db")
C3.4_gosummary_0.95 <- bb_gosummary(x = C3.4_goenrichment, 
                                   reduce_threshold = 0.95,
                                   go_db = "org.Mm.eg.db")

bb_goscatter(simMatrix = C3.4_gosummary_0.8$simMatrix,
             reducedTerms = C3.4_gosummary_0.8$reducedTerms)
bb_goscatter(simMatrix = C3.4_gosummary_0.9$simMatrix,
             reducedTerms = C3.4_gosummary_0.9$reducedTerms)
bb_goscatter(simMatrix = C3.4_gosummary_0.95$simMatrix,
             reducedTerms = C3.4_gosummary_0.95$reducedTerms)
#Supp 2F Cluster 3.6
query_C3.6 <- dplyr::filter(F3_kmeans10_tm_Top50, cell_group %in% '3.6')[["gene_short_name"]]
C3.6_goenrichment <- bb_goenrichment(query = query_C3.6, 
                                     reference = bb_rowmeta(mouse_cds_list$file_1),
                                     go_db = "org.Mm.eg.db")
C3.6_goenrichment
C3.6_gosummary_0.8 <- bb_gosummary(x = C3.6_goenrichment, 
                                   reduce_threshold = 0.8,
                                   go_db = "org.Mm.eg.db")
C3.6_gosummary_0.9 <- bb_gosummary(x = C3.6_goenrichment, 
                                   reduce_threshold = 0.9,
                                   go_db = "org.Mm.eg.db")

bb_goscatter(simMatrix = C3.6_gosummary_0.8$simMatrix,
             reducedTerms = C3.6_gosummary_0.8$reducedTerms)
bb_goscatter(simMatrix = C3.6_gosummary_0.9$simMatrix,
             reducedTerms = C3.6_gosummary_0.9$reducedTerms)

#TODO Supp 2G - Previously - Depleted genes in clusters 3.4 & 3.6  


# ####################################################################################################################################################
# S3E1a <- S3E_plotlist[[1]]|S3E_plotlist[[2]]|S3E_plotlist[[3]]|S3E_plotlist[[4]]
# S3E1a2 <- bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
#                               cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
#                                 filter(k_10_assignment == "B")), variable = "genotype",
#                    genes_to_plot = "Ccr7", pseudocount = 0) | bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
#                               cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
#                                 filter(k_10_assignment == "B")), variable = "genotype",
#                    genes_to_plot = "Il4", pseudocount = 0) +theme(axis.title.y = element_blank()) | bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
#                                         cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
#                                           filter(k_10_assignment == "B")), variable = "genotype",
#                              genes_to_plot = "Cd69", pseudocount = 0) +theme(axis.title.y = element_blank()) | bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
#                                         cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
#                                           filter(k_10_assignment == "B")), variable = "genotype",
#                              genes_to_plot = "Cd93", pseudocount = 0) +theme(axis.title.y = element_blank())
# # S3E1aOne <- S3E1a1 | S3E1a2 | S3E1a3 | S3E1a4
# # S3E1a / S3E1aOne + plot_layout(widths = c(5,1))
# 
# S3E1b <- S3E_plotlist[[5]]|S3E_plotlist[[6]]|S3E_plotlist[[7]]|S3E_plotlist[[8]]|S3E_plotlist[[9]]
# b2 <- bb_gene_violinplot(
#   filter_cds(
#     mouse_cds_list[[1]],
#     cells = bb_cellmeta(mouse_cds_list[[1]]) |>
#       filter(k_10_assignment == "B")
#   ),
#   variable = "genotype",
#   genes_to_plot = "Cxcr5",
#   pseudocount = 0
# )+ theme(axis.title.y = element_blank()) | bb_gene_violinplot(
#   filter_cds(
#     mouse_cds_list[[1]],
#     cells = bb_cellmeta(mouse_cds_list[[1]]) |>
#       filter(k_10_assignment == "B")
#   ),
#   variable = "genotype",
#   genes_to_plot = "Myc",
#   pseudocount = 0
# ) + theme(axis.title.y = element_blank()) |
#   bb_gene_violinplot(
#     filter_cds(
#       mouse_cds_list[[1]],
#       cells = bb_cellmeta(mouse_cds_list[[1]]) |>
#         filter(k_10_assignment == "B")
#     ),
#     variable = "genotype",
#     genes_to_plot = "Il10",
#     pseudocount = 0
#   ) + theme(axis.title.y = element_blank()) |
#   bb_gene_violinplot(
#     filter_cds(
#       mouse_cds_list[[1]],
#       cells = bb_cellmeta(mouse_cds_list[[1]]) |>
#         filter(k_10_assignment == "B")
#     ),
#     variable = "genotype",
#     genes_to_plot = "Mki67",
#     pseudocount = 0
#   ) + theme(axis.title.y = element_blank()) |
#   bb_gene_violinplot(
#     filter_cds(
#       mouse_cds_list[[1]],
#       cells = bb_cellmeta(mouse_cds_list[[1]]) |>
#         filter(k_10_assignment == "B")
#     ),
#     variable = "genotype",
#     genes_to_plot = "Npm1",
#     pseudocount = 0
#   ) + theme(axis.title.y = element_blank())
# 
# S3E1b / S3E1b2
