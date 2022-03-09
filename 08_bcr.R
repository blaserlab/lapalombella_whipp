#Pull Analysis Configs
analysis_configs <-
  read_excel(
    "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/analysis_configs.xlsx",
    sheet = "ingest"
  ) %>%
  mutate(directory = bb_fix_file_path(`10X_dir`)) %>%
  dplyr::select(-`10X_dir`) %>%
  left_join(
    read_excel(
      "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/analysis_configs.xlsx",
      sheet = "metadata"
    )
  )
# load the bcr data
  map_dfr(
   .x = analysis_configs %>% pull(specimen) %>% unique(),
  .f = function(x) {
    read_csv(paste0(analysis_configs %>% filter(lib_type == "bcr", specimen == x) %>% pull(directory),
                   "/outs/vdj_b/filtered_contig_annotations.csv"),
             col_types = "clcldcccccllccddcc") #%>%
       #dplyr::select(barcode, raw_clonotype_id) %>%
       #filter(!is.na(raw_clonotype_id)) %>%
       #mutate(sample_clonotype = paste0(str_replace(x,"_.*",""),"_",raw_clonotype_id)) #%>%
       #mutate(barcode_sample = paste0(barcode, "_", x)) %>%
       #select(barcode_sample, sample_clonotype) %>%
       #distinct()
    
  }
) %>% View()
  # right_join(colData(cds_aligned) %>% as_tibble(rownames = "barcode_sample")) %>% View()
  # pull(sample_clonotype)

#Clonotype Filtering

bcr <- left_join(
  colData(cds_main) %>% as_tibble(rownames = "barcode_sample"),
  map_dfr(
    .x = analysis_configs %>% pull(specimen) %>% unique(),
    .f = function(x) {
      read_csv(paste0(analysis_configs %>% filter(lib_type == "bcr", specimen == x) %>% pull(directory),
                      "/outs/vdj_b/filtered_contig_annotations.csv"),
               col_types = "clcldcccccllccddcc") %>%
      dplyr::select(barcode, raw_clonotype_id) %>%
      filter(!is.na(raw_clonotype_id)) %>%
      mutate(sample_clonotype = paste0(str_replace(x,"_.*",""),"_",raw_clonotype_id)) %>%
      mutate(barcode_sample = paste0(barcode, "_", x)) %>%
      dplyr::select(-barcode) %>%
      distinct()  
    }
  ))

bcr$sample_clonotype
colData(cds_main)$sample_clonotype <- bcr$sample_clonotype
colData(cds_aligned)$santity_check <- bcr$barcode_sample

#sanity check to determine how many times the barcode_sample we added is not equal to rownames in original cds_aligned (returns 0)
sum(colData(cds_aligned)$sanity_check != rownames(colData(cds_aligned)))
#remove sanity_check
colData(cds_aligned)$sanity_check <- NULL

colData(cds_main)
bb_gene_umap(cds = cds_main, gene_or_genes = "CD19") + facet_wrap(vars(sample))
bb_var_umap(cds_main, var = "specimen_clonotype", legend_pos = "none") + facet_wrap(vars(sample))
bb_var_umap(cds_main, var = "partition", facet_by = "sample")


#colData(cds_aligned)$clonotype <- str_replace(colData(cds_aligned)$sample_clonotype,".*_","")
