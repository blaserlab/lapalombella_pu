bb_var_umap(cds = cds_main, var = "specimen", value_to_highlight = "P9")
bb_var_umap(cds_main, "leukemia_phenotype")
bb_var_umap(cds_main[, colData(cds_main)$leukemia_phenotype == "AML"], "partition")
bb_var_umap(cds_main[, colData(cds_main)$leukemia_phenotype == "AML"], "leiden")
bb_var_umap(cds_main[, colData(cds_main)$leukemia_phenotype %in% c("AML", "No leukemia")], "partition") +
  facet_grid(rows = vars(tissue), cols = vars(leukemia_phenotype)) +
  theme(panel.background = element_rect(color = "grey80"))

bb_gene_umap(cds_main[, colData(cds_main)$leukemia_phenotype %in% c("AML", "No leukemia")], c("Cd14", "Ly6a", "Kit", "Cd3e", "Cd19"))



bb_var_umap(cds = cds_main[, colData(cds_main)$leukemia_phenotype == "AML"], var = "specimen") +
  facet_grid(cols = vars(primary_or_engraftment), rows = vars(tissue)) +
  theme(panel.background = element_rect(color = "grey80"))


colData(cds_main) %>%
  as_tibble() %>%
  # filter(specimen %in% c("P2", "P14_S1303_SPN_AML")) %>%
  ggplot(mapping = aes(x = log(Size_Factor), color = specimen)) +
  geom_density()


sequencing_metrics

write_csv(sequencing_metrics, str_glue("{tables_out}/sequencing_metrics.csv"))
