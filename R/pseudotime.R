renv::restore()
blaseRtemplates::git_easy_branch(branch = "ethan_working")               #

library(blaseRtemplates)
library(blaseRtools)
devtools::install_github("blaserlab/blaseRdata")
bb_renv_datapkg("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/datapkg")

#library(lazyData)
requireData("lapalombella.whipp.datapkg")
library("lapalombella.whipp.datapkg")

usethis::edit_r_profile()

colData(cds_main)

# create a new metadata column with better labels for a figure
#for use in group_cells_by = "nice_label"
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

#generate cds subset by pt & B cells (via clonotype_id)
cds_subset2712 <- cds_main[, colData(cds_main)$patient == "pt_2712" &
                         colData(cds_main)$clonotype_id %in% "clonotype1"]

#cds_subset1245 <- cds_main[, colData(cds_main)$patient == "pt_1245" &
#                         colData(cds_main)$clonotype_id %in% "clonotype1"]

bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP"))))
#blaseRtools::bb_gene_pseudotime
cds_subset2712<- order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP")))
colData(cds_subset2712)
cds_subset2712 <- learn_graph(cds_subset2712)

#colData(cds_subset)
plot_cells(cds_subset2712,
           color_cells_by = "clonotype_id",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds_subset2712 <- order_cells(cds_subset2712)

#manually selected root nodes
plot_cells(cds_subset2712,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

gtest <- graph_test(cds_subset2712, neighbor_graph="principal_graph", cores=4)
#monocle3::graph_test

pr_deg_ids <- row.names(subset(gtest, q_value < 0.05))
cds_subset <- choose_cells(cds_subset)



#state_1_cells <- subset(pData(cds), pData(cds)$State == 1)
#state_1_cell_names <- rownames(state_1_cells)
#state_1_exprs <- subset(exprs(cds), colnames(exprs(cds)) %in% state_1_cell_names)
#state_1_genes <- rownames(state_1_exprs)