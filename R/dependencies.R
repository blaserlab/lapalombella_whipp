# renv --------------------------------------------------------------------

# set up the renv from scratch

# renv::init(bioconductor = TRUE)

# restore the renv from the lockfile

# renv::restore()



# package installation ----------------------------------------------------

# # Try this first...it's faster:
# blaseRtemplates::easy_install("<package name>", how = "link_from_cache")
# blaseRtemplates::easy_install("ggtext", how = "link_from_cache")

# # If you need a new package or an update, try this:
#blaseRtemplates::easy_install("blaserlab/blaseRtools", how = "new_or_update")
# blaseRtemplates::easy_install("rr1859/R.4Cker", "new_or_update")

# # If you are installing from a "tarball", use this:
# blaseRtemplates::easy_install("/path/to/tarball.tar.gz")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages


# load core packages for the analysis -------------------------------------
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("blaseRtemplates"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("blaseRtools"))
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
suppressPackageStartupMessages(library("SingleCellExperiment"))
library("topGO")

# run this to update the data package in renv
blaseRtemplates::bb_renv_datapkg("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/datapkg")


# load the data set into a hidden environment
requireData("lapalombella.whipp.datapkg")
