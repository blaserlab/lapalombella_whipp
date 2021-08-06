# load the bcr data
troubleshoot <- left_join(
  colData(cds_aligned) %>% as_tibble(rownames = "barcode_sample"),
  map_dfr(
  .x = analysis_configs %>% pull(specimen) %>% unique(),
  .f = function(x) {
    read_csv(
      paste0(
        analysis_configs %>% filter(lib_type == "bcr", specimen == x) %>% pull(directory),
        "/outs/vdj_b/filtered_contig_annotations.csv"
      ),
      col_types = "clcldcccccllccddcc"
    ) %>%
    select(barcode, raw_clonotype_id) %>%
      filter(!is.na(raw_clonotype_id)) %>%
      mutate(sample_clonotype = paste0(str_replace(x, "_.*", ""), "_", raw_clonotype_id)) %>% 
      mutate(barcode_sample = paste0(barcode, "_", x)) %>% 
      select(-barcode) %>%
      distinct()
      
  }
))
troubleshoot$sample_clonotype
colData(cds_aligned)$sample_clonotype <- troubleshoot$sample_clonotype
colData(cds_aligned)$sanity_check <- troubleshoot$barcode_sample

# sanity check good
sum(colData(cds_aligned)$sanity_check != rownames(colData(cds_aligned)))
colData(cds_aligned)$sanity_check <- NULL

bb_var_umap(cds_aligned, var = "sample_clonotype", legend_pos = "none") + facet_wrap(vars(sample_short))

bb_gene_umap(cds_aligned[,colData(cds_aligned)$sample_short == "L37"], gene_or_genes = "CD79A")

colData(cds_aligned[,colData(cds_aligned)$sample_short == "L37"])

troubleshoot %>% filter(str_detect(sample, "L37"))












bb_var_umap(cds = cds_aligned,var = "density", facet_by = "sample")

bb_gene_umap(cds_aligned, gene_or_genes = rowData(cds_aligned) %>% 
               as_tibble() %>% 
               filter(gene_short_name %in% c("CD79A", "CD5")) %>% 
               select(id) %>% 
               mutate(gene_grouping = "Ethan Genes") %>%
               as.data.frame())
