source("00_packages_functions.R")

# load the reference dataset

seurat_reference <- LoadH5Seurat("~/workspace_pipelines/OSUBlaserlab_analyses/brad/scrnaseq_analysis/preprocessing/references/pbmc_multimodal.h5seurat")


exprs_renamed_long <- broom::tidy(exprs(cds_final)) %>%
  as_tibble() %>%
  left_join(.,rowData(cds_final) %>% as_tibble(), by = c("row" = "id")) %>%
  select(gene_short_name, column, value) %>%
  select(row = gene_short_name, column, value)
exprs_renamed_long_rna <- exprs_renamed_long 
exprs_renamed_df <- as.data.frame(exprs_renamed_long_rna)
exprs_renamed_rows <- exprs_renamed_df %>% select(row) %>% unique() %>% data.frame()
exprs_renamed_cols <- exprs_renamed_df %>% select(column) %>% unique() %>% data.frame()
exprs_renamed_df$rowIdx <- match(exprs_renamed_df$row, exprs_renamed_rows$row)
exprs_renamed_df$colIdx <- match(exprs_renamed_df$column, exprs_renamed_cols$column)



exprs_spars_rna <- sparseMatrix(
  i = exprs_renamed_df$rowIdx,
  j = exprs_renamed_df$colIdx,
  x = as.integer(exprs_renamed_df$value),
  dims = c(nrow(exprs_renamed_rows),nrow(exprs_renamed_cols)),
  dimnames = list(exprs_renamed_rows$row,exprs_renamed_cols$column)
) 


seurat_cds <- CreateSeuratObject(counts = exprs_spars_rna,assay = "RNA")

seurat_cds <- SCTransform(seurat_cds)


anchors <- FindTransferAnchors(
  reference = seurat_reference,
  query = seurat_cds,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50, 
  features = intersect(rownames(x = seurat_reference), VariableFeatures(object = seurat_cds)),
  reference.assay = "SCT",
  query.assay = "SCT",
  verbose = T
)


seurat_cds <- TransferData(
  anchorset = anchors,
  reference = seurat_reference,
  query = seurat_cds,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  )
)

seurat_cds <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = seurat_reference,
  query = seurat_cds,
  new.reduction.name = "ref.spca"
)

seurat_cds <- ProjectUMAP(
  query = seurat_cds,
  query.reduction = "ref.spca",
  reference = seurat_reference,
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

DimPlot(object = seurat_reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

seurat_metadata <- as_tibble(seurat_cds@meta.data, rownames = "barcode")
seurat_dims <- as_tibble(seurat_cds@reductions[["ref.umap"]][[1:dim(seurat_cds@reductions[["ref.umap"]])[1]]], rownames = "barcode")
seurat_data <- left_join(seurat_dims,seurat_metadata)


rm(
  seurat_reference,
  seurat_cds,
  anchors,
  exprs_renamed_cols,
  exprs_renamed_df,
  exprs_renamed_long,
  exprs_renamed_long_rna,
  exprs_renamed_rows,
  exprs_spars_rna,
  seurat_dims,
  seurat_metadata
)

