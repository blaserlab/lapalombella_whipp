setwd(
  "/workspace/workspace_pipelines/lapalombella_whipp"
)

#libraries
library("ggrepel")
library("htmltools")
library("ggpubr")
library("pheatmap")
library("ggplotify")
library("RColorBrewer")
library("assertthat")
library("data.table")
library("monocle3")
library("tidyverse")
library("cowplot")
library("ggsci")
library("scales")
library("patchwork")
library("scater")
library("reticulate")
library("DoubletFinder")#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library("Seurat")
library("Matrix")
library("rhdf5")
library("glue")
library("broom")
library("Matrix.utils")
#library("DESeq2")
library("Seurat")
library("SeuratDisk")
library("patchwork")
library("SeuratData")
library("circlize")
library("ComplexHeatmap")
library("rstatix")
library("WRS2")
library("Cairo")
library("tidyselect")
library("simplifyEnrichment")
library("topGO")
library("Rgraphviz")
library("GOsummaries")
library("gprofiler2")
library("magick")
library("prodlim")
library("ViSEAGO")
library("rrvgo")
library("biclust")
library("fgsea")
library("readxl")
library("viridis")
library("ggplotify")
library("GenomicFeatures")
library("DescTools")
library("fastSave")
library("cicero")
library("Gviz")
library("GenomicRanges")
library("RMariaDB")
library("biomaRt")
library("BSgenome")
library("BSgenome.Drerio.UCSC.danRer11")
library("lattice")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# fix file path header
fix_file_path <- function(x) {
  x <- str_replace_all(x,"\\\\","/")
  x <- str_replace(x,"X:/","~/network/X/")
  # x <- str_replace(x, " ", "\\\\ ")
}
 


#unmask important functions####
filter <- dplyr::filter
mutate <- dplyr::mutate
group_by <- dplyr::group_by
select <- dplyr::select
#custom operators####
`%notin%` <- Negate(`%in%`)

#helper functions####

rank_in_group <- function(data) {
  
  data %>%
    dplyr::mutate(constcol = 1) %>%
    dplyr::mutate(rank = cumsum(.data$constcol)) %>%
    dplyr::select(-.data$constcol)
  
}

tbl_to_matrix <- function(data) {
  data <- data %>%
    as.data.frame()
  rownames(data) <- data[,1]
  return(as.matrix(data[,-1]))
}

remove_dupes <- function(data,column) {
  data <- data %>%
    group_by(!!sym(column)) %>%
    mutate(duplicate_flag = n() > 1) %>%
    filter(!duplicate_flag) %>% # removes all rows where zf_id is not unique  
    select(-duplicate_flag)
  return(data)
}

se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))

data_summary <- function(x) {
  m <- median(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary1 <- function(x) {
  m <- median(x)
  ymin <- m - (IQR(x) / 2)
  ymax <- m + (IQR(x) / 2)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary2 <- function(x) {
  m <- mean(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

data_summary3 <- function(x) {
  m <- median(x)
  ymin <- m - mad(x)
  ymax <- m + mad(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

# selector function for topGO
selector <- function(theScore) {
  return (theScore == 1)
}
    

# standard 2color jitter geom
standard_2color_jitter <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_jitter(
      shape = 21,
      width = 0.2,
      size = 1,
      stroke = 0.25,
      alpha = 0.4
    )
}

# blank plot to fill in plot_grid
empty_data <- data.frame(
  horiz = c(1,2),
  vert = c(3,4)
)
blank_plot<-ggplot(data = empty_data, aes(x = horiz, y = vert))+
  geom_point(color = "transparent")+
  theme(axis.ticks = element_line(color = "transparent"),
        axis.title = element_text(color = "transparent"),
        axis.line = element_line(color = "transparent"),
        axis.text = element_text(color = "transparent"))

# print stat report to text file

print_full_stats <- function(data, class_var, class1, class2, value_var, out = NULL) {
	if (!is.null(out)) {
	  sink(out) 
	}
	print(data %>%
	      group_by(!!sym(class_var)) %>%
	      summarise(n = n(),
			mean = mean(!!sym(value_var)),
			sd = sd(!!sym(value_var)),
			se = se(!!sym(value_var)),
			median = median(!!sym(value_var)),
			IQR = IQR(!!sym(value_var)),
			shapiro_p = shapiro.test(!!sym(value_var))[[2]]))
	class1_data <- data %>%
		filter(!!sym(class_var) == class1) %>%
		pull(!!sym(value_var))
	class2_data <- data %>%
		filter(!!sym(class_var) == class2) %>%
		pull(!!sym(value_var))
	cat("_________________________________")
	cat("\n\n")
	cat(paste0("Class 1 is ",class1," and Class 2 is ",class2,"\n"))
	cat("_________________________________")
	cat("\n")
	print(shapiro.test(class1_data))
	cat("_________________________________")
	cat("\n")
	print(shapiro.test(class2_data))
	cat("_________________________________")
	cat("\n")
	print(t.test(class1_data, class2_data, alternative = "two.sided", var.equal = T))
	cat("_________________________________")
	cat("\n")
	print(t.test(class1_data, class2_data, alternative = "two.sided", var.equal = F))
	cat("_________________________________")
	cat("\n")
	print(wilcox.test(class1_data, class2_data))
	if (!is.null(out)){
	  sink()
	}
}


# custom scrnaseq functions

plot_cells_alt <-
  function (cds,
            gene_or_genes,
            cell_size = 1,
            alpha = 1 ,
            ncol = NULL,
            plot_title = NULL,
            color_legend_title = NULL
            ) {
    
    data <- plot_cells(cds = cds, genes = gene_or_genes)[["data"]]
    # data$feature_label <-
    #   factor(data$feature_label, levels = gene_or_genes)
    background_data <- data %>% filter(is.na(value))
    foreground_data <- data %>% filter(!is.na(value))
    p <- ggplot() +
      geom_point(
        data = background_data,
        aes(x = data_dim_1, y = data_dim_2),
        color = "grey80",
        shape = 1,
        size = cell_size,
        stroke = 0.25
      ) +
      geom_point(
        data = foreground_data,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          color = log10(value)
        ),
        shape = 16,
        size = cell_size,
        alpha = alpha
      ) +
      scale_color_viridis_c(end = 0.8) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        color = color_legend_title,
        title = plot_title
      ) +
      facet_wrap(facets = vars(feature_label), ncol = ncol) +
      theme(strip.background = element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }


cumulative_max_expr<-function(extracted_df,gene_list){
  return(extracted_df %>% 
           tbl_df() %>% 
           select(-feature_id,-gene_short_name,-id) %>% 
           pivot_wider(names_from = feature_label, values_from = value) %>% 
           replace(., is.na(.),0) %>%
           mutate(max_val = do.call(pmax, c(select(., one_of(gene_list))))) %>%
           mutate(max_val = na_if(max_val,0)))
}

add_cds_factor_columns<-function(cds, columns_to_add){
  for (i in 1:length(columns_to_add)) {
    colData(cds)$new<-unname(columns_to_add[i])
    names(colData(cds))[names(colData(cds)) == "new"] <- names(columns_to_add[i])
  }
  return(cds)
}

custom_variable_plot<-function(cds,
                               var, 
                               value_to_highlight = NULL, 
                               foreground_alpha = 1, 
                               legend_pos = "right", 
                               cell_size = 0.5, 
                               legend_title = NULL,
                               plot_title = NULL,
                               palette = NULL,
                               ref_dim_x = NULL,
                               ref_dim_y = NULL,
                               overwrite_labels = FALSE,
                               group_label_size = 3,
                               alt_label_col = NULL) {
  data<-plot_cells(cds)[["data"]]
  data_long<-data %>% pivot_longer(cols = (!!sym(var)), names_to = "var")
  dim_x <- ifelse(is.null(ref_dim_x),"data_dim_1",ref_dim_x)
  dim_y <- ifelse(is.null(ref_dim_y),"data_dim_2",ref_dim_y)
  
  # generate text data frame for variable labels if you are going to use them
  if (is.null(alt_label_col)) {
    text_df <- data_long %>% group_by(value)
  } else {
    text_df <- data %>% pivot_longer(cols = !!sym(alt_label_col), names_to = "var") %>% group_by(value)
  }
  median_coord_df <-
    text_df %>%
    summarise(
      fraction_of_group = n(),
      text_x = median(!!sym(dim_x)),
      text_y = median(!!sym(dim_y))
    )
  text_df <- left_join(text_df, median_coord_df) %>%
    mutate(label = value)
  text_df <-
    text_df %>% group_by(label,text_x,text_y) %>% summarise()

  # make the main plot
  plot<-ggplot()
  if(!is.null(value_to_highlight)){
    data_background<-data_long %>% filter(value %notin% value_to_highlight)
    data_long<-data_long %>% filter(value %in% value_to_highlight)
    plot<-plot+
      geom_point(data = data_background,
                 aes(x = !!sym(dim_x),
                     y = !!sym(dim_y)),
                 stroke = 0.25,
                 shape = 1,
                 size = cell_size,
                 color = "grey80")
  }
  plot<-plot+
    geom_point(data = data_long,
               aes(x = !!sym(dim_x),
                   y = !!sym(dim_y),
                   fill = value,
                   color = value
                   ),
               stroke = 0.25,
               shape = 21,
               alpha = foreground_alpha,
               size = cell_size)
  if (class(data_long$value)=="numeric") {
    plot<-plot+
      scale_fill_viridis_c(guide = "colorbar",na.value = "transparent")+
      scale_color_viridis_c(guide = F, na.value = "grey80")
  } else if (length(palette) == 1 && palette == "viridis") {
    plot<-plot+
      scale_fill_viridis_d(begin = 0.1,end = 0.9)+
      scale_color_viridis_d(begin = 0.1,end = 0.9, guide = F)
  } else if (length(palette) == 1 && palette == "rcolorbrewer") {
    plot<-plot+scale_color_brewer(palette = "Paired",guide = F)+
      scale_fill_brewer(palette = "Paired")
  } else if (!is.null(palette)) {
    plot<-plot+scale_color_manual(values = palette, guide = F)+
      scale_fill_manual(values = palette)
  } else {
    plot<-plot+scale_color_discrete(guide = F)+
      scale_fill_discrete()
  }
  if(class(data_long$value)!= "numeric") {
    plot<-plot+guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1, color = "transparent")))
  }
  plot<-plot+labs(x = ifelse(is.null(ref_dim_x),"UMAP 1",ref_dim_x),
                  y = ifelse(is.null(ref_dim_y),"UMAP 2", ref_dim_y),
                  title = plot_title,
                  fill = legend_title) +
    theme(plot.title = element_text(hjust = 0.5))
  plot<-plot +
    theme(legend.position = legend_pos)#+coord_fixed()

  #option to overwrite labels
  if (overwrite_labels == T) {
    plot <- plot +
      theme(legend.position = "none") +
      geom_text_repel(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,
        min.segment.length = 1)
  }

  return(plot)
 # return(text_df) 
}



custom_cp_plot <- function(cds,
                           var = NULL,
                           alpha = 1,
                           cp ,
                           overwrite_labels = T,
                           legend_pos = "none",
                           cell_size = 0.5,
                           legend_title = NULL,
                           plot_title = NULL,
                           group_label_size = 3,
                           value_to_highlight = NULL,
                           alt_color_var = NULL,
                           facet_by = NULL,
                           palette = NULL) {
    #extract the data from the input cds
  data <- plot_cells(cds)[["data"]]
  #convert to long format
  data_long <- pivot_longer(data = data,
                            cols = cp,
                            names_to = "cp")
  # generate text data frame for cluster/partition labels
  text_df <- data_long %>% group_by(value)
  median_coord_df <-
    data_long %>% group_by(value) %>% summarise(
      fraction_of_group = n(),
      text_x = median(data_dim_1),
      text_y = median(data_dim_2)
    )
  text_df <- left_join(text_df, median_coord_df) %>%
    mutate(label = value)
  text_df <-
    text_df %>% group_by(label) %>% summarise(text_x = dplyr::first(text_x),
                                              text_y = dplyr::first(text_y))
  # if highlighting a categorical variable, generate background data and keep data_long as foreground
  if (!is.null(value_to_highlight)) {
    background_data_long <- data_long %>% filter((!!sym(var)) != value_to_highlight)
    data_long <- data_long %>% filter((!!sym(var)) %in% value_to_highlight)
  }
  #initialize the plot
  plot <- ggplot()
  # lay down the background plot if using
  if(!is.null(value_to_highlight)){
    plot<-plot+
      geom_point(data = background_data_long,
                 aes(x = data_dim_1, y = data_dim_2),
                 stroke = 0.25,
                 shape = 1,
                 size = cell_size,
                 color = "grey80")
  }

  # make the main colored plot from data_long
  if (!is.null(alt_color_var)){
    plot <- plot +
      geom_point(
        data = data_long,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          fill = (!!sym(alt_color_var)),
          color = (!!sym(alt_color_var))),
        stroke = 0.25,
        shape = 21,
        alpha = alpha,
        size = cell_size
        )
  } else {
    plot <- plot +
      geom_point(
        data = data_long,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          fill = value,
          color = value),
        stroke = 0.25,
        shape = 21,
        alpha = alpha,
        size = cell_size)
  }

  plot<-plot+
    scale_color_discrete(guide = F) +
    scale_fill_discrete()+
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = plot_title,
      fill = legend_title) +
    theme(plot.title = element_text(hjust = 0.5))
  # overwrite labels if you want to
  if (overwrite_labels == T) {
    plot <- plot +
      theme(legend.position = legend_pos) +
      ggrepel::geom_text_repel(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,min.segment.length = 1) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))
  } else {
    plot <- plot +
      theme(legend.position = legend_pos) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))}
  #option to facet
  if (!is.null(facet_by)) {
    plot <- plot +
      facet_wrap(facets = vars(!!sym(facet_by)),) +
      theme(strip.background = element_blank())
  }
  if (!is.null(palette) && palette[1] =="rcolorbrewer") {
    colourCount = length(unique(data_long$value))
    getPalette = colorRampPalette(brewer.pal(12, "Paired"))
    plot<-plot+
      scale_color_manual(values = getPalette(colourCount), guide = F)+
      scale_fill_manual(values = getPalette(colourCount))
  }

  if (!is.null(palette) && length(palette)>1) {
    plot <- plot +
      scale_color_manual(values = palette, guide = F) +
      scale_fill_manual(values = palette)
  }

    return(plot)
}




custom_violin_plot <-
  function(cds,
           variable,
           genes_to_plot,
           pseudocount = 1,
           include_jitter = FALSE,
           ytitle = "Log10(Expr+1)",
           plot_title = NULL,
           rows = 1,
           show_x_label = TRUE,
           legend_pos = "none",
           comparison_list = NULL,
           palette = NULL, 
           violin_alpha = 1,
           jitter_alpha = 1,
           facet_scales = "fixed",
           order_genes = TRUE
           #sig_lab_y = 1,
           #yplotmax
  ) {
    my_comparisons <-
      comparison_list#(list(c(comparator1,comparator2),c(comparator1,comparator3)...))
    data_to_plot <-
      plot_genes_violin(cds_subset = cds[rowData(cds)$gene_short_name %in% genes_to_plot,], group_cells_by = variable)[["data"]]
    if(order_genes) {
      data_to_plot <- 
        data_to_plot %>% mutate(gene_short_name = factor(gene_short_name, levels = genes_to_plot))
    }
    p1 <-
      ggplot(data = data_to_plot, aes(
        x = !!as.name(variable),
        y = log10(expression +
                    pseudocount)
      )) #expression already normalized when data extracted by violin plot function
    if (include_jitter == TRUE) {
      p1 <-
        p1 + geom_jitter(
          shape = 21,
          size = 0.5,
          color = "black",
          alpha = jitter_alpha,
          width = 0.2
        )
    }
    p1 <- p1 +
      geom_violin(
        scale = "width",
        color = "black",
        trim = T,
        size = 0.5,
        aes(fill = !!as.name(variable)),
        draw_quantiles = 0.5
      )
    # p1 <- p1 +
    #   ylim(0,yplotmax)
    # if (!is.null(comparison_list)) {
    #   p1<-p1+stat_compare_means(
    #     comparisons = my_comparisons,
    #     method = "wilcox.test",
    #     size = 2,
    #     label = "p.signif",
    #     hide.ns = F,
    #     label.y = sig_lab_y
    #   )
    # } 
    p1 <- p1+
      theme(legend.position = legend_pos) +
      theme(legend.direction = "horizontal") +
      theme(legend.justification = "center") +
      labs(
        x = "",
        y = ytitle,
        title = plot_title,
        fill = NULL
      ) 
    if(is.null(palette)) {
      p1 <- p1 + 
        scale_fill_viridis_d(
          alpha = 0.6,
          begin = 0.1,
          end = 0.9
          )
      } else {
      p1 <- p1 + 
        scale_fill_manual(
          values = alpha(palette, violin_alpha)
          )
    }
    p1 <- p1 +
      theme(plot.title = element_text(hjust = 0.5)) +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~ gene_short_name, nrow = rows, scales = facet_scales) +
      theme(strip.background = element_rect(fill = "transparent"))
    if (show_x_label == F) {
      p1 <- p1 + theme(axis.text.x = element_blank())
    }
    return(p1)
  }

plot_genes_in_pseudotime_alt<-function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
                                        ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                        trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                        vertical_jitter = NULL, horizontal_jitter = NULL) 
{
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  tryCatch({
    pseudotime(cds_subset)
  }, error = function(x) {
    stop(paste("No pseudotime calculated. Must call order_cells first."))
  })
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                      "plotted."))
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)), 
                          msg = paste("pseudotime must be a column in", "colData. Please run order_cells", 
                                      "before running", "plot_genes_in_pseudotime."))
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  assertthat::assert_that(!is.null(size_factors(cds_subset)))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                      "plotted."))
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", 
                     by.y = "row.names")
  cds_exprs$adjusted_expression <- cds_exprs$expression
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
  model_expectation <- model_predictions(model_tbl, new_data = colData(cds_subset))
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                             function(x) {
                               data.frame(expectation = model_expectation[x$f_id, 
                                                                          x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
  if (!is.null(color_cells_by)) {
    q <- q + geom_point(aes_string(color = color_cells_by), 
                        size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
    if (class(colData(cds_subset)[, color_cells_by]) == "numeric") {
      q <- q + viridis::scale_color_viridis(option = "C")
    }
  }
  else {
    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
  }
  q <- q + geom_line(aes(x = pseudotime, y = expectation), 
                     data = cds_exprs)
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  q <- q + ylab("Expression")
  q <- q + xlab("pseudotime")
  q <- q + monocle_theme_opts()
  q
}

# qc function ####
qc_func <- function(cds, genome = c("human", "mouse", "zfish")) {
  if (genome == "human") {
    mito_pattern <-"^MT-"
  }
  if (genome == "mouse") {
    mito_pattern <- "^mt-"
  }
  if (genome == "zfish") {
    mito_pattern <- "^mt-"
  }
  cds <-
    scater::addPerCellQC(cds,subsets=list(Mito=grep(mito_pattern, rowData(cds)$gene_short_name)))
  cds_tbl <- as_tibble(colData(cds)) %>%
    mutate(
      qc.detected = isOutlier(
        detected,
        log = TRUE,
        type = "lower",
        nmads = 2
      ),
      qc.mito = isOutlier(
        subsets_Mito_percent,
        type = "higher",
        nmads = 2
      ),
      qc.any = qc.detected | qc.mito,
      pre_post = "pre"
    )
  qc_detected_thresh <-
    attr(
      isOutlier(
        cds_tbl$detected,
        log = TRUE,
        type = "lower",
        nmads = 2
      ),
      "thresholds"
    )
  qc_mito_thresh <-
    attr(isOutlier(cds_tbl$subsets_Mito_percent, type = "higher", nmads = 2),
         "thresholds")
  cds_tbl_filtered_any <-
    cds_tbl %>% filter(qc.any == FALSE) %>% mutate(pre_post = "post")
  cds_tbl_filtered_detected <-
    cds_tbl %>% filter(qc.detected == FALSE) %>% mutate(pre_post = "post")
  cds_tbl_filtered_mito <-
    cds_tbl %>% filter(qc.mito == FALSE) %>% mutate(pre_post = "post")
  cds_to_plot_any <- bind_rows(cds_tbl, cds_tbl_filtered_any)
  cds_to_plot_detected <-
    bind_rows(cds_tbl, cds_tbl_filtered_detected)
  cds_to_plot_mito <- bind_rows(cds_tbl, cds_tbl_filtered_mito)
  plot_any <-
    ggplot(cds_to_plot_any, aes(x = fct_rev(pre_post), y = log10(detected))) +
    geom_violin() + 
    labs(y = "log10(detected features)", title = "qc.any", x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  plot_mito <-
    ggplot(cds_to_plot_mito, aes(x = fct_rev(pre_post), y = subsets_Mito_percent)) +
    geom_violin() +
    geom_hline(yintercept = qc_mito_thresh[[2]]) +
    scale_y_continuous(breaks = c(qc_mito_thresh[[2]], seq(
      0, max(cds_to_plot_mito$subsets_Mito_percent), by = 20
    ))) +
    labs(y = "Percent Mitochondrial", title = "qc.mito", x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  plot_detected <-
    ggplot(cds_to_plot_detected, aes(
      x = fct_rev(pre_post),
      y = log10(detected)
    )) +
    geom_violin() +
    geom_hline(yintercept = log10(qc_detected_thresh[[1]])) +
    scale_y_continuous(breaks = c(log10(qc_detected_thresh[[1]]), seq(0, max(
      log10(cds_to_plot_detected$detected)
    ), by = 1))) +
    labs(y = "Log10(detected features)", title = "qc.detected", x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  cds_return <- cds_tbl %>% select(barcode, qc.any)
  return_list <-
    list(cds_return,
         plot_any,
         plot_mito,
         plot_detected,
         qc_detected_thresh,
         qc_mito_thresh)
  return(return_list)
}


# run doubletfinder #### 
find_homotypic_doublets <-
  function(cds,
           doublet_prediction,
           qc_table
           ) {
    # system(paste0("gunzip -k ", directory, "/*"))
    keepers <- qc_table %>% filter(qc.any == FALSE) %>% pull(barcode)
    #seu_data <- Read10X(data.dir = directory)
    seu_object <- CreateSeuratObject(exprs(cds))
    seu_object@meta.data$barcode <- rownames(seu_object@meta.data)
    seu_object <- subset(seu_object, subset = barcode %in% keepers)
    seu_object <- NormalizeData(seu_object)
    seu_object <-
      FindVariableFeatures(seu_object,
                           selection.method = "vst",
                           nfeatures = 2000)
    all.genes  <- rownames(seu_object)
    seu_object <- ScaleData(seu_object, features = all.genes)
    seu_object <-
      RunPCA(seu_object, features = VariableFeatures(object = seu_object))
    seu_object <- FindNeighbors(seu_object)
    seu_object <- FindClusters(seu_object)
    seu_object <- RunUMAP(seu_object, dims = 1:10)
    
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list <-
      paramSweep_v3(seu_object, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seu_object@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(doublet_prediction * length(annotations))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seu_object_lowconf <-
      doubletFinder_v3(
        seu_object,
        PCs = 1:10,
        pN = 0.25,
        pK = 0.09,
        nExp = nExp_poi,
        reuse.pANN = FALSE,
        sct = FALSE
      )
    seu_object_highconf <-
      doubletFinder_v3(
        seu_object,
        PCs = 1:10,
        pN = 0.25,
        pK = 0.09,
        nExp = nExp_poi.adj,
        reuse.pANN = FALSE,
        sct = FALSE
      )
    
    ## extract the barcode and singlet/doublet calls
    solm <-
      as_tibble(seu_object_lowconf@meta.data) %>% select(barcode, doubletfinder_low_conf = contains("classifications"))
    sohm <-
      as_tibble(seu_object_highconf@meta.data) %>% select(barcode, doubletfinder_high_conf = contains("classifications"))
    
    ## join together and export as one
    
    return(full_join(solm, sohm))
  }


custom_gene_dotplot <- function (cds, markers, group_cells_by = "cluster", reduction_method = "UMAP", 
                                 norm_method = c("log", "size_only"), lower_threshold = 0, 
                                 max.size = 10, ordering_type = c("cluster_row_col", "maximal_on_diag", 
                                                                  "none"), axis_order = c("group_marker", "marker_group"), 
                                 flip_percentage_mean = FALSE, pseudocount = 1, scale_max = 3, 
                                 scale_min = -3, colorscale_name = NULL, sizescale_name = NULL) 
{
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  if (!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", 
                                                  "partition") | group_cells_by %in% names(colData(cds)), 
                            msg = paste("group_cells_by must be a column in", 
                                        "the colData table."))
  }
  norm_method = match.arg(norm_method)
  gene_ids = as.data.frame(fData(cds)) %>% tibble::rownames_to_column() %>% 
    dplyr::filter(rowname %in% markers | gene_short_name %in% 
                    markers) %>% dplyr::pull(rowname)
  if (length(gene_ids) < 1) 
    stop(paste("Please make sure markers are included in the gene_short_name\",\n               \"column of the fData!"))
  if (flip_percentage_mean == FALSE) {
    major_axis <- 1
    minor_axis <- 2
  }
  else if (flip_percentage_mean == TRUE) {
    major_axis <- 2
    minor_axis <- 1
  }
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression")
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  if (group_cells_by == "cluster") {
    cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    cell_group <- tryCatch({
      partitions(cds, reduction_method = reduction_method)
    }, error = function(e) {
      NULL
    })
  }
  else {
    cell_group <- colData(cds)[, group_cells_by]
  }
  if (length(unique(cell_group)) < 2) {
    stop(paste("Only one type in group_cells_by. To use plot_genes_by_group,", 
               "please specify a group with more than one type. "))
  }
  names(cell_group) = colnames(cds)
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>% 
    dplyr::summarize(mean = mean(log(Expression + pseudocount)), 
                     percentage = sum(Expression > lower_threshold)/length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, 
                        ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, 
                        ExpVal$mean)
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, "gene_short_name"]
  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 + 
                                                                                     major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id
  if (ordering_type == "cluster_row_col") {
    row_dist <- stats::as.dist((1 - stats::cor(t(res)))/2)
    row_dist[is.na(row_dist)] <- 1
    col_dist <- stats::as.dist((1 - stats::cor(res))/2)
    col_dist[is.na(col_dist)] <- 1
    ph <- pheatmap::pheatmap(res, useRaster = T, cluster_cols = TRUE, 
                             cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
                             clustering_distance_cols = col_dist, clustering_distance_rows = row_dist, 
                             clustering_method = "ward.D2", silent = TRUE, filename = NA)
    ExpVal$Gene <- factor(ExpVal$Gene, levels = colnames(res)[ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group, levels = row.names(res)[ph$tree_row$order])
  }
  else if (ordering_type == "maximal_on_diag") {
    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for (i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
    if (major_axis == 1) {
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), 
                                            max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene, levels = dimnames(res)[[2]][max_ind_vec])
    }
    else {
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)), 
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group, levels = dimnames(res)[[1]][max_ind_vec])
    }
  }
  else if (ordering_type == "none") {
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }
  if (flip_percentage_mean) {
    g <- ggplot(ExpVal, aes(y = Gene, x = Group)) + geom_point(aes(colour = percentage, 
                                                                   size = mean)) + viridis::scale_color_viridis(name = ifelse(is.null(colorscale_name),"proportion",colorscale_name)) + 
      scale_size(name = ifelse(is.null(sizescale_name),"log(mean + 0.1)",sizescale_name), range = c(0, 
                                                     max.size))
  }
  else {
    g <- ggplot(ExpVal, aes(y = Gene, x = Group)) + geom_point(aes(colour = mean, 
                                                                   size = percentage)) + viridis::scale_color_viridis(name = ifelse(is.null(colorscale_name),"log(mean + 0.1)",colorscale_name)) + 
      scale_size(name = ifelse(is.null(sizescale_name),"proportion",sizescale_name), range = c(0, max.size))
  }
  if (group_cells_by == "cluster") {
    g <- g + xlab("Cluster")
  }
  else if (group_cells_by == "partition") {
    g <- g + xlab("Partition")
  }
  else {
    g <- g + xlab(group_cells_by)
  }
  g <- g + ylab("Gene") + theme(axis.text.x = element_text(angle = 30, 
                                                                                  hjust = 1))
  if (axis_order == "marker_group") {
    g <- g + coord_flip()
  }
  g
}

# a function to join up the qc and doubletfinder data------------------------------
join_metadata <- function(cds, qc_data, doubletfinder_data) {
  cds_tbl <- as_tibble(colData(cds))
  cds_tbl <- left_join(cds_tbl, qc_data)
  cds_tbl <- left_join(cds_tbl, doubletfinder_data)
  cds_df <- as.data.frame(cds_tbl)
  row.names(cds_df) <- cds_df$barcode
  cds <- new_cell_data_set(
    expression_data = cds@assays@data$counts,
    cell_metadata = cds_df,
    gene_metadata = rowData(cds)
  )
  
  
  return(cds)
}

# scatter or bubble plot for plotting go term summaries-------------------------------------
rrvgo_scatter <-
  function (simMatrix,
            reducedTerms,
            size = "score",
            addLabel = TRUE,
            labelSize = 4) {
    x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE,
                  k = 2)
    df <-
      cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points),
                                                        reducedTerms$go), c("term", "parent", "parentTerm", "size", "score")])
    p <-
      ggplot(df, aes(x = V1, y = V2, color = parentTerm)) +
      geom_point(aes(size = !!sym(size)), alpha = 0.5) +
      scale_color_discrete(guide = FALSE) +
      scale_size_continuous(guide = FALSE, range = c(0, 10)) +
      geom_text_repel(
        aes(label = parentTerm),
        segment.size = 0.25,
        data = subset(df, parent == rownames(df)),
        box.padding = grid::unit(0.5, "lines"),
        size = labelSize,
        color = "black",
        max.overlaps = 100,
        force = 2,
        seed = 1234,
        segment.curvature = -0.1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        min.segment.length = 0
      ) +
      labs(x = "PCoA 1", y = "PCoA 2")
   
   return(p)
  }

# a function to summarize go terms

summarize_go <- function(x, reduce_threshold) {
    simMatrix <-
      calculateSimMatrix(x = x[[3]]$GO.ID,
                         ont = "BP",
                         orgdb = "org.Dr.eg.db")
    scores <- setNames(-log10(ifelse(is.na(
      as.numeric(x[[3]]$classicFisher)),
      1e-30,
      as.numeric(x[[3]]$classicFisher)
    )),
    x[[3]]$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = reduce_threshold,
                                    orgdb = "org.Dr.eg.db")
    returnlist <- list(simMatrix, scores, reducedTerms)
    names(returnlist) <- c("simMatrix", "scores", "reducedTerms")
    return(returnlist)
}
