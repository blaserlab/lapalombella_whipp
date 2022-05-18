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

cds_subset1245 <- cds_main[, colData(cds_main)$patient == "pt_1245" &
                         colData(cds_main)$clonotype_id %in% "clonotype1"]

#bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP"))))
#blaseRtools::bb_gene_pseudotime

cds_subset2712<- order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP")))
colData(cds_subset2712)

cds_subset2712 <- cluster_cells(cds_subset2712)
cds_subset2712 <- learn_graph(cds_subset2712)
colData(cds_subset2712)
plot_cells(cds_subset2712,
           color_cells_by = "sample",
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

monocle3::graph_test
#use trace and substitute Matrix::rBind with rbind
trace('calculateLW', edit = T, where = asNamespace("monocle3"))

gtest <- monocle3::graph_test(cds_subset2712, neighbor_graph="principal_graph", cores=4)

pr_deg_ids <- row.names(subset(gtest, q_value < 0.05))
view(pr_deg_ids)
pr_deg_ids_q_0.01 <- row.names(subset(gtest, q_value < 0.01))
pr_deg_ids_q_0.01
view(gtest)
write.csv(gtest)
pr_deg_ids_q_0.0001 <- row.names(subset(gtest, q_value < 0.0001))
gtestOK <- filter(gtest, status == "OK")
#q0.gtest <- filter(gtest, q_value < 4.407045e-307)
gene_module_df <- find_gene_modules(cds_subset2712[pr_deg_ids_q_0.0001,], resolution=c(0,10^seq(-6,-1)))

q0.gtest4 <- filter(gtest, q_value < 4.407045e-307 & module == '4')
view(gtest)
write.csv(q0.gtest4)
view(pr_deg_ids_q_0)
write.table(pr_deg_ids_q_0)
#packageVersion("monocle3")
#cds_subset <- choose_cells(cds_subset)

plot_cells(cds_subset2712, genes=c("LINC01011","FOXO3B","LGMN", "CDK1","AURKB","AURKA","NME1","PLK1","CCL3", "CCL4","CD38", "CD3E"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
bb_gene_dotplot(cds_subset2712, gene_or_genes)
bb_var_umap(cds_subset1245, var = "sample")
bb_gene_umap(cds_subset1245, gene_or_genes = c("CCL3", "CCL4"))
bb_gene_umap(cds_main, gene_or_genes = c("CD57", "CD3"))

#nice labeling
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

#dotplot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
             colData(cds_main)$clonotype_id %in% "clonotype1"],
  markers = c("CCL3", "CCL4", "AURKB","AURKA","NME1","CDK1"),
  group_cells_by = "nice_label",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)

