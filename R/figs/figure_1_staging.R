# modifications to CDS.  Should incorporate into data package -----------
# aggr_umap_tbl <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/cellranger_aggr/lapalombella_whipp_aggr_20211130/outs/count/analysis/umap/2_components/projection.csv", col_names = c("cell_id", "aggr_UMAP_1", "aggr_UMAP_2"), skip = 1) |> 
#   mutate(barcode_truncated = str_remove(cell_id, "[[:digit:]]")) |> 
#   mutate(sample_num = str_extract(cell_id, "[[:digit:]]")) |> 
#   mutate(sample_name = recode(sample_num, 
#                               "1" = "L33_19972712RTLN", 
#                               "2" = "L34_19972712RTPBMC", 
#                               "3" = "L35_19972712CLLPBMC", 
#                               "4" = "L36_19971245RTLN", 
#                               "5" = "L37_19971245RTPBMC", 
#                               "6" = "L38_19971245CLLPBMC", 
#                               )) |> 
#   mutate(cell_id = paste0(barcode_truncated, "1_", sample_name)) |> count(sample_name)
#   select(cell_id, aggr_UMAP_1, aggr_UMAP_2)
# 
# aggr_cluster_tbl <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/cellranger_aggr/lapalombella_whipp_aggr_20211130/outs/count/analysis/clustering/kmeans_9_clusters/clusters.csv", col_names = c("cell_id", "cluster"), skip = 1) |> 
#   mutate(barcode_truncated = str_remove(cell_id, "[[:digit:]]")) |> 
#   mutate(sample_num = str_extract(cell_id, "[[:digit:]]")) |> 
#   mutate(sample_name = recode(sample_num, 
#                               "1" = "L33_19972712RTLN", 
#                               "2" = "L34_19972712RTPBMC", 
#                               "3" = "L35_19972712CLLPBMC", 
#                               "4" = "L36_19971245RTLN", 
#                               "5" = "L37_19971245RTPBMC", 
#                               "6" = "L38_19971245CLLPBMC", 
#                               )) |> 
#   mutate(cell_id = paste0(barcode_truncated, "1_", sample_name)) |>
#   select(cell_id, aggr_cluster = cluster) |> 
#   mutate(aggr_cluster = as.character(aggr_cluster))
# 
# 
# cds_main <- bb_tbl_to_coldata(obj = cds_main, min_tbl = aggr_umap_tbl)
# cds_main <- bb_tbl_to_coldata(obj = cds_main, min_tbl = aggr_cluster_tbl)
# 
colData(cds_main)$disease_tissue <- paste0(colData(cds_main)$disease, " ", colData(cds_main)$tissue)
colData(cds_main)$partition_assignment_1 <- recode(colData(cds_main)$partition, "1" = "B", "2" = "B", "3" = "T", "4" = "T", "5" = "Mono", "6" = "B")

nadeu_11b <- readxl::read_excel("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx", sheet = "Supplementary Table 11b", skip = 5, col_names = c("feature_id", "gene_short_name", "mean", "l2fc", "se", "p", "padj", "direction"))

cds_main <- nadeu_11b |> 
  filter(direction == "Up") |> 
  filter(padj < 0.05) |> 
  mutate(feature_id = str_remove(feature_id, "\\..*")) |> 
  mutate(nadeu_RT_gene = TRUE) |> 
  bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)
  

  
# figs ---------------------------------------

bb_var_umap(cds_main, "density", facet_by = "patient", alt_dim_y = "aggr_UMAP_1", alt_dim_x = "aggr_UMAP_2")
bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
  "log_local_n",
  facet_by = "disease_tissue"
)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
  "partition"
)
bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "PRMT5"
) + facet_wrap(~disease_tissue)

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "CD79A"
)
bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "CD3E"
)


bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "CD14"
)

bb_genebubbles(
  obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")),
  genes = c("CD14", "CD3E", "CD79A"), cell_grouping = "partition_assignment_1"
)

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), "partition_assignment_1"
  
)

bb_gene_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(patient == "pt_2712")), gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
) + 
  facet_wrap(~disease_tissue)



bb_cellmeta(cds_main) |> glimpse()
bb_cellmeta(cds_main) |> count(patient, disease, tissue)
plot_cells(cds_main, group_cells_by = "partition")

bb_cellmeta(cds_main) |> count(sample)

