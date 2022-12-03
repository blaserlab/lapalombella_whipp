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

#TODO add limits and reduce cell size
S2E_plotlist <- map(.x = c("Ccr7", "Il4", "Cd69", "Cd93", "Cxcr5", "Myc", "Il10", "Mki67", "Npm1"),
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
# S2E_plotlist[[1]]/S2E_plotlist[[2]]/S2E_plotlist[[3]]/S2E_plotlist[[4]]/S2E_plotlist[[5]]
# S2E_plotlist[[6]]/S2E_plotlist[[7]]/S2E_plotlist[[8]]/S2E_plotlist[[9]]

###############################
#Pseudobulk - for size effect (log2FC &padj on violin plots)
#B cell CDS
#unique(colData(mouse_cds_list[[1]])$genotype)

s2_k10_B_cds<- filter_cds(mouse_cds_list[[1]],
                          cells = bb_cellmeta(mouse_cds_list[[1]]) |>
                            filter(k_10_assignment == "B"))

exp_design <- 
  bb_cellmeta(s2_k10_B_cds) |>
  group_by(sample, genotype) |>
  summarise()
exp_design

#exp_design
s2_k10_B_cds_copy <- s2_k10_B_cds 

rowData(s2_k10_B_cds_copy)$id <- rownames(rowData(s2_k10_B_cds_copy))

pseudobulk_res <-
  bb_pseudobulk_mf(cds = s2_k10_B_cds_copy,
                   pseudosample_table = exp_design, 
                   design_formula = "~genotype",
                   result_recipe = c("genotype", "PRMT5", "TCL1"))

#less conservative approach (pseudobulk is a very conservative approach)
#bb_monocle_regression(cds = f5_k10_B_cds, gene_or_genes = "Cd93", form = "~genotype")

pseudobulk_res$Header 

pseudobulk_res$Result |> filter(gene_short_name == "Cd93")
pseudobulk_res$Result |> filter(gene_short_name == "Il4")
Fig3_PRMT5vsTCL1_spleenB_pseudobulk<- pseudobulk_res$Result
write.csv(Fig3_PRMT5vsTCL1_spleenB_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Fig3_PRMT5vsTCL1_spleenB_pseudobulk.csv")

S2E1 <-
  (S2E_plotlist[[1]]+theme(axis.text.y = element_text()))/plot_spacer()/ bb_gene_violinplot(
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

# S2E1 <- S2E_plotlist[[1]]/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
#                                                         cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
#                                                           filter(k_10_assignment == "B")), variable = "genotype",
#                                              genes_to_plot = "Ccr7", pseudocount = 0)+ theme(strip.text = element_blank())

S2E2 <- S2E_plotlist[[2]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                      cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                        filter(k_10_assignment == "B")), variable = "genotype",
                                                           genes_to_plot = "Il4", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S2E3 <- S2E_plotlist[[3]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                      cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                        filter(k_10_assignment == "B")), variable = "genotype",
                                                           genes_to_plot = "Cd69", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S2E4 <-(S2E_plotlist[[4]]+theme(legend.position = "right",
                                legend.key.size = unit(4, "mm"),
                                legend.margin = margin(c(0, -7, 0, -6.5)),
                                legend.title = element_blank()))/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                                                             cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                                                               filter(k_10_assignment == "B")), variable = "genotype",
                                                                                                  genes_to_plot = "Cd93", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

# S2F_1 <- 
#   (S2E1 |plot_spacer()|
#      S2E2 |plot_spacer()|
#      S2E3 |plot_spacer()|
#      S2E4) + plot_layout(widths = c(1,-0.1,1,-0.1,1,-0.1,1))
S2F_1 <- S2E1|S2E2|S2E3|S2E4
ggsave("S2F_1.pdf", path = T_Figs, width = 11.6, height = 3.3)

S2E5 <- (S2E_plotlist[[5]]+theme(axis.text.y = element_text()))/plot_spacer()/ bb_gene_violinplot(
  filter_cds(
    mouse_cds_list[[1]],
    cells = bb_cellmeta(mouse_cds_list[[1]]) |>
      filter(k_10_assignment == "B")
  ),
  variable = "genotype",
  genes_to_plot = "Cxcr5",
  pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE
) + theme(strip.text = element_blank()) + theme(axis.title.y = element_text(margin = margin(
  t = 0,
  r = -1.5,
  b = 0,
  l = -7
))) + labs(y = "B cell Expression") + plot_layout(heights = c(1, -0.1 ,1))

S2E6 <-S2E_plotlist[[6]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                       filter(k_10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Myc", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S2E7 <-S2E_plotlist[[7]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                       filter(k_10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Il10", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S2E8 <-S2E_plotlist[[8]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                       filter(k_10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Mki67", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))
S2E9 <-S2E_plotlist[[9]]/plot_spacer()/bb_gene_violinplot(filter_cds(mouse_cds_list[[1]], 
                                                                     cells = bb_cellmeta(mouse_cds_list[[1]]) |> 
                                                                       filter(k_10_assignment == "B")), variable = "genotype",
                                                          genes_to_plot = "Npm1", pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = .1, include_jitter = TRUE)+ theme(strip.text = element_blank()) +theme(axis.title.y = element_blank())+ plot_layout(heights = c(1, -0.1 ,1))

S2F_2 <- (S2E5 | S2E6 | S2E7 | S2E8 | S2E9) 
ggsave("S2F_2.pdf", path = T_Figs, width = 11.6, height = 3.3)

# S2F_2 <- 
#   (S2E5 |plot_spacer()|
#      S2E6 |plot_spacer()|
#      S2E7 |plot_spacer()|
#      S2E8 |plot_spacer()|
#      S2E9) + plot_layout(widths = c(1,-0.1,1,-0.1,1,-0.1,1,-0.1,1))
# ggsave("S2F_2.pdf", path = T_Figs, width = 11.6, height = 4.4)


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
