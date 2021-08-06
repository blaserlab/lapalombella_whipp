library("blaseRtools")# devtools::install_github("git@github.com:blaserlab/blaseRtools.git")
library("monocle3")
library("tidyverse")
library("cowplot")
library("circlize")
library("ComplexHeatmap")
library("fastSave")
library("readxl")

theme_set(theme_cowplot(font_size = 10))

experimental_group_palette <- c(
  "group" = "color"# "#3C5488",#blue; "#DC0000",#red
)

jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#313695","white","#A50026")

add_cds_factor_columns<-function(cds, columns_to_add){
  for (i in 1:length(columns_to_add)) {
    colData(cds)$new<-unname(columns_to_add[i])
    names(colData(cds))[names(colData(cds)) == "new"] <- names(columns_to_add[i])
  }
  return(cds)
}
