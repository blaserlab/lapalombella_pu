colData(cds_main)
cds_main
bb_var_umap(cds_main, "partition", overwrite_labels = T)
bb_var_umap(cds_main, "leiden")
bb_var_umap(cds_main, "louvain")

colData(cds_main) %>%
  as_tibble() %>%
  ggplot(mapping = aes(x = log2(Size_Factor), color = specimen)) +
  geom_density()

bb_var_umap(cds_main, "Size_Factor")
bb_gene_umap(cds_main, "Cd3e")
bb_gene_umap(cds_main, "Icos")
bb_gene_umap(cds_main, "Kit")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Kit")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Hoxa9")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Mpo")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Spi1")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Cebpa")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Myb")
bb_gene_umap(cds_main[,colData(cds_main)$partition == "9"], "Runx1")
bb_gene_umap(cds = cds_main[,colData(cds_main)$partition == "9"], gene_or_genes = "Ly6a")


prealignment_top_markers %>% as_tibble()


partition_top_markers <- top_markers(cds_main, group_cells_by = "partition", genes_to_test_per_group = 50, reference_cells = 1000, cores = 39)

prealignment_top_markers <- read_csv("/workspace/workspace_pipelines/lapalombella.pu.datapkg/data/prealignment_top_markers.csv")
prealignment_top_markers %>%
  filter(cluster_method == "partition") %>%
  View()

# partition 9 looks like myeloid
bb_var_umap(cds_main[,colData(cds_main)$partition == 9], "leiden", overwrite_labels = T)

prealignment_top_markers %>%
  filter(cell_group %in% c("leiden 16", "leiden 36", "leiden 17", "leiden 22", "leiden 41", "leiden 27", "leiden 33") ) %>%
  arrange(cell_group) %>%
  View()
