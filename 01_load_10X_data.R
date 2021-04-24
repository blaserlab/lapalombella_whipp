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

