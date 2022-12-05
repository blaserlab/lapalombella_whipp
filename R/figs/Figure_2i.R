#TODO add RNAseq to datapkg?

#Supplemental Figure 2i - CLEAR RNAseq: PRMT5 vs. TCL1 mice
CLEARseq_data <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/CLEAR Bulk RNAseq/CLEAR_v2_sheet_prmt5_vs_tscl1.2.3.4.csv")[,c(1,3,7)]
####################################################################

S2I_highlights <-
  c(
    "CD79B",
    "FOS",
    "FOSB",
    "FYN",
    "JUN",
    "JUND",
    "CXCL16",
    "ZAP70",
    "FOSB",
    "Zfp36",
    "Ubc",
    "PLK2",
    "Il1r2",
    "Atf3",
    "Arg2",
    "MMP24",
    "Hspa1l",
    "Elane",
    "Stat4"
  ) |> str_to_title()

#####VOLCANO PLOT
CLEAR_volcano_data <- CLEARseq_data |>
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 2) |>
  mutate(text_label = ifelse(Gene %in% S2I_highlights, Gene, ""))

volcano_CLEAR <-
  ggplot(
    CLEAR_volcano_data,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21,
             size = 0.5,
             alpha = 0.4) +
  geom_text_repel(color = "black",
                  fontface = "italic",
                  box.padding = 0.1, #0.5
                  point.padding = 0.0, #0.25
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  size = 3,
                  segment.size = 0.15,
                  force = 1,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log<sub>2</sub> fold change") +
  ylab("-log<sub>10</sub> adjusted p-value") +
  theme(axis.title.x =  element_markdown()) +
  theme(axis.title.y = element_markdown()) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(title = "PRMT5 vs TCL1")+
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-0.75*max(abs(range(CLEAR_volcano_data |> dplyr::filter(!is.na(padj)) |> pull(log2FoldChange)))), 1*max(abs(range(CLEAR_volcano_data |> filter(!is.na(padj)) |> pull(log2FoldChange))))))
volcano_CLEAR
