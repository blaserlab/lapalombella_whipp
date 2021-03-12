source("00_packages_functions.R")

analysis_configs <-
  read_excel(
    "~/network/X/Labs/Blaser/single_cell/rosa_mar2021/analysis_configs.xlsx",
    sheet = "ingest"
  ) %>%
  mutate(directory = fix_file_path(`10X_dir`)) %>%
  select(-`10X_dir`) %>%
  left_join(
    read_excel(
      "~/network/X/Labs/Blaser/single_cell/rosa_mar2021/analysis_configs.xlsx",
      sheet = "metadata"
    )
  )


# cds_list <-
#   map(
#     .x = analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names),
#     .f = function(x, data = analysis_configs %>% filter(lib_type == "gex")) {
#       directory <- data %>% filter(pipestance_names == x) %>% pull(directory)
#       cds <- load_cellranger_data(directory)
#       cds <-
#         add_cds_factor_columns(
#           cds = cds,
#           columns_to_add = c(
#             "pt" = data %>% filter(pipestance_names == x) %>% pull(patient),
#           )
#         )
#       return(cds)
#     }
#   ) %>% set_names(analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names))

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

