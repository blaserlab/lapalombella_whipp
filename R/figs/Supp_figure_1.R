#Supplemental Figure 1A
S1A <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(patient == "Pt 1")), "clonotype_id", facet_by = "disease_tissue")+ theme(legend.position = "bottom")
S1A

#Supplemental Figure 1B
S1B <- bb_genebubbles(
  cds_main,
  genes = c("CD14","CD4","CD8A", "CD3E", "CD79A", "MS4A1", "CD19"
  ),
  cell_grouping = c("leiden", "leiden_assignment_1"),
  return_value = "data"
) |> 
  ggplot(mapping = aes(x = leiden, 
                       y = gene_short_name, 
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~leiden_assignment_1, scales = "free_x", ) +
  theme_minimal_grid(font_size = 8) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
S1B

#Supplemental Figure 1E
S1E <-
  (
    bb_var_umap(
      cds_main,
      "patient",
      value_to_highlight = "Pt 1",
      legend_pos = "none",
      plot_title = "All Pt 1 cells",
      palette = "#EF8A62"
    ) +
      bb_var_umap(
        cds_main,
        "patient",
        value_to_highlight = "Pt 2",
        legend_pos = "none",
        plot_title = "All Pt 2 cells",
        palette = "#67A9CF"
      )
  ) /
  (
    bb_var_umap(
      filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(tissue == "LN")),
      "patient",
      value_to_highlight = "Pt 1",
      legend_pos = "none",
      plot_title = "RT LN-Pt 1",
      palette = "#EF8A62"
    ) +
      bb_var_umap(
        filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(tissue == "LN")),
        "patient",
        value_to_highlight = "Pt 2",
        legend_pos = "none",
        plot_title = "RT LN-Pt 2",
        palette = "#67A9CF"
      ) 
  )
S1E

#Supplemental Figure 1F

S1F1 <-bb_var_umap(
  filter_cds(cds_main, cells= bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
  "prmt5_label",
  #facet_by = "disease_tissue",
  value_to_highlight = "PRMT5+ cells",
  palette = "#481F70FF", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  #labs(x = "UMAP 1", y= "UMAP 2") + 
  labs(y=NULL,x=NULL) +theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(color = "black", fill = "white"))

S1F2<- bb_var_umap(
  filter_cds(cds_main, cells= bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
  "prmt5_myc_label",
  #facet_by = "disease_tissue",
  value_to_highlight = "PRMT5+/MYC+ cells",
  palette = "#365D8DFF", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  labs(x = "UMAP 1", y= "UMAP 2") + labs(y=NULL,x=NULL) +theme(legend.position = "right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(color = "black", fill = "white"))

S1F3<-bb_var_umap(
  filter_cds(cds_main, cells= bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
  "prmt5_mki67_label",
  #facet_by = "disease_tissue",
  value_to_highlight = "PRMT5+/MKI67+ cells",
  palette = "#21908CFF", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  labs(x = "UMAP 1", y= "UMAP 2") + labs(y=NULL,x=NULL) +theme(legend.position = "right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(color = "black", fill = "white"))

S1F4<-bb_var_umap(
  filter_cds(cds_main, cells= bb_cellmeta(cds_main) |> filter(disease_tissue == "RT LN")),
  "prmt5_myc_mki67_label",
  #facet_by = "disease_tissue",
  value_to_highlight = "PRMT5+/MYC+/MKI67+ cells",
  palette = "#47C16EFF", 
  legend_pos = "bottom", 
  foreground_alpha = 0.6
) +
  #theme(axis.text = element_blank()) +
  #theme(axis.ticks = element_blank()) +
  labs(x = "UMAP 1", y= "UMAP 2") + labs(y=NULL,x=NULL) +theme(legend.position = "right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(color = "black", fill = "white"))
S1F1/S1F2/S1F3/S1F4

#Supplemental Figure 1G
S1G<- bb_var_umap(cds_main, "leiden", overwrite_labels = T, facet_by= "disease_tissue")
S1G
