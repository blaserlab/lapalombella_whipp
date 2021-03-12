source("00_packages_functions.R")

analysis_configs <- read_excel("",sheet = "ingest") %>% 
  mutate(directory = fix_file_path(`10X_dir`)) %>%
  select(-`10X_dir`) %>%
  left_join(read_excel("",sheet = "metadata"))
  

cds_list <-
  map(
    .x = analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names),
    .f = function(x, data = analysis_configs %>% filter(lib_type == "gex")) {
      directory <- data %>% filter(pipestance_names == x) %>% pull(directory)
      cds <- load_cellranger_data(directory)
      cds <-
        add_cds_factor_columns(
          cds = cds,
          columns_to_add = c(
            "pt" = data %>% filter(pipestance_names == x) %>% pull(patient),
          )
        )
      return(cds)
    }
  ) %>% set_names(analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names))

dir.create("data_out")
summarized_sequencing_metrics <-
  tibble(pipestance_path = analysis_configs %>% filter(lib_type == "gex") %>% pull(directory)) %>%
  mutate(metrics_summary_path = paste0(pipestance_path, "/outs/metrics_summary.csv")) %>%
  mutate(cds_dim_cells = unname(sapply(X = cds_list, FUN = dim)[2, ])) %>%
  mutate(cds_name = names(sapply(X = cds_list, FUN = dim)[2, ])) %>%
  left_join(.,
            bind_rows(lapply(
              X = .$metrics_summary_path, FUN = read_csv
            )),
            by = c("cds_dim_cells" = "Estimated Number of Cells")) %>% # this is your sanity check.  Joining on cell number derived from the CDS object and the metrics summary
  select(
    cds_name,
    cds_dim_cells,
    `Mean Reads per Cell`,
    `Median Genes per Cell`,
    `Fraction Reads in Cells`
  ) %>%
  write_csv("data_out/summarized_sequencing_metrics.csv")

