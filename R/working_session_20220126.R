# create a new metadata column with better labels for a figure
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

# violin plot of PRMT5 expression
bb_gene_violinplot(cds_main[, colData(cds_main)$patient == "pt_2712" &
                              colData(cds_main)$clonotype_id %in% "clonotype1"], 
                   variable = "sample", 
                   genes_to_plot = "PRMT5")

# umap plot of PRMT5 expression
bb_gene_umap(cds_main[, colData(cds_main)$patient == "pt_2712" &
                        colData(cds_main)$clonotype_id %in% "clonotype1"], 
             gene_or_genes = "TP53")

# gene dotplot for pt_2712 in malignant cells
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
             colData(cds_main)$clonotype_id %in% "clonotype1"],
  markers = c("PRMT5", "MYC"),
  group_cells_by = "nice_label",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)

# density plot of size factors by sample
bb_cellmeta(cds_main) %>%
  ggplot(mapping = aes(x = Size_Factor, color = sample)) +
  geom_density()

bb_gene_modules(cds_main)

colData(cds_main)
rowData(cds_main)

#bb_gene_umap(cds_main[, colData(cds_main)$patient == "pt_2712" &
#                        colData(cds_main)$clonotype_id %in% "clonotype1"], 
#             gene_or_genes = "PRMT5")

#rowData(cds_main)$module

 rowdata <- fData(cds_main)
view(rowdata)

rowData(cds_main)$supermodule == 1
blaseRtools::bb_gene_modules

cds_subset <- cds_main[, colData(cds_main)$patient == "pt_2712" &
                         colData(cds_main)$clonotype_id %in% "clonotype1"]
rowData(cds_subset)

bb_gene_pseudotime(cds_subset)

order_cells(cds_subset)
learn_graph(cds_subset)

cluster_cells(cds_subset, reduction_method = "UMAP")
bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset, reduction_method = "UMAP"))))

cds_subset <- choose_cells(cds_subset)

blaseRtools::bb_gene_pseudotime