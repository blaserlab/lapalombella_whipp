# uncomment and run this to install packages indicated in lock file
# renv::restore()

# blaseRtools and additional dependencies you may have to install 
# renv::install("/usr/lib/R/site-library/blaseRtools")
# renv::install("/usr/lib/R/site-library/DESeq2")
# renv::install("/usr/lib/R/site-library/genefilter")
# renv::install("/usr/lib/R/site-library/annotate")
# renv::install("/usr/lib/R/site-library/AnnotationDbi")
# renv::install("/usr/lib/R/site-library/KEGGREST")
# renv::install("/usr/lib/R/site-library/Biostrings")
# renv::install("/usr/lib/R/site-library/geneplotter")
# renv::install("/usr/lib/R/site-library/DoubletFinder")
# renv::install("/usr/lib/R/site-library/Seurat")
# renv::install("/usr/lib/R/site-library/SeuratDisk")
# renv::install("/usr/lib/R/site-library/rrvgo")
# renv::install("/usr/lib/R/site-library/GO.db")
# renv::install("/usr/lib/R/site-library/GOSemSim")
# renv::install("/usr/lib/R/site-library/scater")
# renv::install("/usr/lib/R/site-library/topGO")
# renv::install("/usr/lib/R/site-library/fastSave")
# renv::install("/usr/lib/R/site-library/fgsea")
# renv::install("/usr/lib/R/site-library/fastmatch")
# renv::install("/usr/lib/R/site-library/tinytex")
# renv::install("/usr/lib/R/site-library/reshape2")
# renv::install("/usr/lib/R/site-library/pheatmap")
# renv::install("/usr/lib/R/site-library/Rcpp")
# renv::install("/usr/lib/R/site-library/pander")


# load core packages for the analysis
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("pander"))


# run this to update the data package in renv
bb_renv_datapkg("~/network/X/Labs/Blaser/collaborators/lapalombella_whipp_network/datapkg")


# load the data set into a hidden environment
requireData("lapalombella.whipp.datapkg")
