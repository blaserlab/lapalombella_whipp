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


bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset, reduction_method = "UMAP"))))

cds_subset <- choose_cells(cds_subset)

blaseRtools::bb_gene_pseudotime