analysis_configs <-
  read_excel(
    "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/analysis_configs.xlsx",
    sheet = "ingest"
  ) %>%
  mutate(directory = bb_fix_file_path(`10X_dir`)) %>%
  select(-`10X_dir`) %>%
  left_join(
    read_excel(
      "~/network/X/Labs/Blaser/single_cell/lapalombella_whipp/analysis_configs.xlsx",
      sheet = "metadata"
    )
  )

cds_list <-
  map(
    .x = analysis_configs %>% filter(lib_type == "gex") %>% pull(directory),
    .f = function(x, 
                  data = analysis_configs %>% filter(lib_type == "gex")) {
      #directory <- data %>% filter(directory == x) %>% pull(directory)
      cds <- bb_load_multi_counts(x)
      cds <-
        add_cds_factor_columns(
          cds = cds,
          columns_to_add = c(
            "disease" = data %>% filter(directory == x) %>% pull(disease),
            "tissue" = data %>% filter(directory == x) %>% pull(tissue),
            "sample_date" = data %>% filter(directory == x) %>% pull(sample_date),
            "myc" = data %>% filter(directory == x) %>% pull(myc)
            
          )
        )
      return(cds)
    }
  ) %>% set_names(analysis_configs %>% filter(lib_type == "gex") %>% pull(specimen))

dir.create("data_out")

map2_dfr(.x = paste0(analysis_configs %>% 
                      filter(lib_type == "gex") %>%
                      select(specimen, directory) %>%
                      filter(complete.cases(.)) %>% 
                      pull(directory), "/outs/count/metrics_summary.csv"),
         .y = analysis_configs %>% 
                       filter(lib_type == "gex") %>%
                       select(specimen, directory) %>%
                       filter(complete.cases(.)) %>% 
                       pull(specimen),
        .f = function(x, y) {
          read_csv(x) %>%
            mutate(specimen = y) %>%
            select(
              specimen,
              `Estimated Number of Cells`,
              `Mean Reads per Cell`,
              `Median Genes per Cell`,
              `Fraction Reads in Cells`
            ) 
          
        }
        ) %>% write_csv("data_out/summarized_gex_sequencing_metrics.csv")


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

colData(cds_aligned)



bb_var_umap(cds = cds_aligned,var = "density", facet_by = "sample")

bb_gene_umap(cds_aligned, gene_or_genes = rowData(cds_aligned) %>% 
               as_tibble() %>% 
               filter(gene_short_name %in% c("CD79A", "CD5")) %>% 
               select(id) %>% 
               mutate(gene_grouping = "Ethan Genes") %>%
               as.data.frame())
