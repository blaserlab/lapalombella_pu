bb_var_umap(cds = cds_main, var = "specimen")
bb_var_umap(cds = cds_main, var = "leukemia_phenotype") + facet_grid(cols = vars(primary_or_engraftment), rows = vars(tissue))
