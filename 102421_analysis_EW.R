library("tidyverse")
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
#devtools::install_github('barkasn/fastSave')
library(fastSave)
#install.packages("lazyData")
library(lazyData)

library(blaseRtools)

#Install most recent tar.gz analysis
#install.packages("~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/datapkg/lapalombella.whipp.datapkg_0.0.0.9001.tar.gz", repos = NULL, type = "source")

#only run requireData function once at the start of each session
requireData("lapalombella.whipp.datapkg")
library(lapalombella.whipp.datapkg)

#bcr_data
colData(cds_main)
cds_main
CDS_new <- cds_main

pt2712 <- CDS_new[,patient] 

bb_gene_umap(cds = cds_main, gene_or_genes = "CD79A") + facet_wrap(vars(sample))
#bb_var_umap(cds_main, var = "specimen_clonotype", legend_pos = "none") + facet_wrap(vars(sample))
bb_var_umap(cds_main, var = "clonotype_id") + facet_wrap(vars(sample))

bb_var_umap(cds_main, var = "cdr3s_aa", legend_pos = "none") + facet_wrap(vars(sample))
bb_var_umap(cds_main[,colData(cds_main)$cdr3s_aa %in% c("IGH:CATDHDDSSGYYHSWGYFDLW;IGK:CQQRSNWPLTF","IGH:CARHLWFGEYRAMGYFDYW;IGL:CQSADSSGTYYVF")], var = "disease") + facet_wrap(vars(sample))
bb_var_umap(cds_main[,colData(cds_main)$cdr3s_aa %in% c("IGH:CATDHDDSSGYYHSWGYFDLW;IGK:CQQRSNWPLTF")], var = "disease") + facet_wrap(rows = vars(patient), cols = vars(tissue)) + theme(panel.background = element_rect(color = "grey80"))
#No clonal pop identified between CLL & RT in sample L38
bb_var_umap(cds_main[,colData(cds_main)$cdr3s_aa %in% c("IGH:CATDHDDSSGYYHSWGYFDLW;IGK:CQQRSNWPLTF")], var = "disease") + facet_grid(rows = vars(patient), cols = vars(tissue)) + theme(panel.background = element_rect(color = "grey80"))

#UMAP
bb_var_umap(cds_main, var = "partition", facet_by = "sample")
bb_var_umap(cds_main, var = "leiden", facet_by = "sample", overwrite_labels = TRUE)
bb_gene_umap(cds = cds_main, gene_or_genes = c("CD8", "CD3E")) + facet_wrap(vars(sample))
bb_gene_umap(cds = cds_main, gene_or_genes = c("PEBP1")) + facet_wrap(vars(sample))
bb_var_umap(cds_main[,colData(cds_main)$partition %in% c("1","2")], var = "disease") + facet_grid(rows = vars(patient), cols = vars(tissue)) + theme(panel.background = element_rect(color = "grey80"))

#Seurat cell type calling within clusters
bb_var_umap(cds_main, var = "seurat_celltype_l1") + facet_wrap(vars(sample))
bb_var_umap(cds_main, var = "seurat_celltype_l2") + facet_wrap(vars(sample))

view(cds_main_top_markers)
view(joined_bcr_data)
bb_gene_umap(cds = cds_main, gene_or_genes = "PDPN") + facet_wrap(vars(sample))
bb_gene_umap(cds = cds_main, gene_or_genes = "RBL2") + facet_wrap(vars(sample))

colData(cds_main)
####################
bb_gene_violinplot(cds = cds_main, variable = "seurat_celltype_l1", genes_to_plot = c("PRMT5"))
bb_gene_violinplot(cds = cds_main, variable = "sample", genes_to_plot = c("PRMT5" %in% patient = "pt_2712")) + facet_wrap(vars(patient="pt_2712", clonotype_id="clonotype1"))
###Subset CDS on 1997-2712#####
CDS_new <- cds_main[,colData(cds_main)$patient == "pt_2712" & colData(cds_main)$clonotype_id %in% "clonotype1"]
bb_gene_violinplot(cds = (CDS_new[,colData(CDS_new)$clonotype_id %in% "clonotype1"]), 
                   variable = "sample", genes_to_plot = c("PRMT5"), ytitle = "PRMT5   Expression", pseudocount = 1) ###+ facet_wrap(vars(tissue, disease))
#CDS_1_2 <- CDS_new[,colData(CDS_new)$partition %in% c('1','2')]
#bb_gene_violinplot(cds = (CDS_1_2[,colData(CDS_1_2)$clonotype_id %in% "clonotype1"]), 
#                   variable = "sample", genes_to_plot = c("PRMT5"), ytitle = "PRMT5   Expression", pseudocount = 1) ###+ facet_wrap(vars(tissue, disease))
colData(CDS_new)
########2712 dot plot
#&Rename/order samples
colData(cds_main)$nice_label <- recode(colData(cds_main)$sample, "L34_19972712RTPBMC" = "RT PBMC", "L33_19972712RTLN" = "RT LN", "L35_19972712CLLPBMC" = "CLL PBMC")
bb_gene_dotplot(cds_main[, colData(cds_main)$patient == "pt_2712" & colData(cds_main)$clonotype_id %in% "clonotype1"],
                markers = c("PRMT5", "MYC"), group_cells_by = "nice_label", group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"), 
                sizescale_name = NULL, max.size = 17.5)

colData(CDS_1_2)
bb_gene_umap(cds = CDS_new, gene_or_genes = c("PRMT5")) + facet_wrap(vars(sample))
#Evaluate a function
blaseRtools::bb_gene_dotplot

#PRMT5_pos <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
#PRMT5_pos <- 
#  row.names(subset(fData(CDS_new), gene_short_name == "PRMT5", mean_expression > 0))
#CDS_new2 <- CDS_new %>% subset(disp_table, mean_expression >= 0.1)

#choose_cells(CDS_new)
#pData(CDS_new)
#fData(CDS_new) %>% row.names(fData(CDS_new), gene_short_name)

#cds_subset <- CDS_new[colData(CDS_New)$gene_short_name %in% "PRMT5"]

#####Subset CDS on 1997-1245
CDS_1245 <- cds_main[,colData(cds_main)$patient %in% "pt_1245"]
bb_gene_violinplot(cds = (CDS_1245[,colData(CDS_1245)$clonotype_id %in% "clonotype1"]), 
                   variable = "sample",genes_to_plot = "PRMT5"), scale_y_log10()) ####+ facet_wrap(vars(tissue, disease))

bb_gene_umap(cds = cds_main, gene_or_genes = c("PRMT5")) + facet_wrap(vars(sample))
bb_gene_umap(cds = cds_main, gene_or_genes = c("TIGIT","CD8A")) + facet_wrap(vars(sample))

bb_gene_dotplot(cds = CDS_new[,colData(CDS_new)$clonotype_id %in% "clonotype1"], markers = c("PRMT5", "CD79A"), group_cells_by = "partition")

bb_gene_dotplot(cds = CDS_new, markers = c("PRMT5", "CD79A"), group_cells_by = "partition")
