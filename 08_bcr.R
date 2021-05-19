# load the bcr data
  map_dfr(
  .x = analysis_configs %>% pull(specimen) %>% unique(),
  .f = function(x) {
    read_csv(paste0(analysis_configs %>% filter(lib_type == "bcr", specimen == x) %>% pull(directory),
                   "/outs/vdj_b/filtered_contig_annotations.csv"),
             col_types = "clcldcccccllccddcc") # %>%
      # select(barcode, raw_clonotype_id) %>%
      # filter(!is.na(raw_clonotype_id)) %>%
      # mutate(sample_clonotype = paste0(str_replace(x,"_.*",""),"_",raw_clonotype_id)) %>%
      # mutate(barcode_sample = paste0(barcode, "_", x)) %>%
      # select(barcode_sample, sample_clonotype) %>%
      # distinct()
    
  }
) %>% View()
  # right_join(colData(cds_aligned) %>% as_tibble(rownames = "barcode_sample")) %>% View()
  # pull(sample_clonotype)

colData(cds_aligned)$clonotype <- str_replace(colData(cds_aligned)$sample_clonotype,".*_","")
