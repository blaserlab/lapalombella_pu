
# good
leiden_auc_mat <-
  t(regulonAUCmat) |>
  as_tibble(rownames = "cell_id") |>
  pivot_longer(-cell_id) |>
  left_join(bb_cellmeta(cds_p568) |>
              select(cell_id, leiden)) |>
  group_by(leiden, name) |>
  summarise(meanAUC = mean(value)) |>
  pivot_wider(names_from = name, values_from = meanAUC) |>
  bb_tbl_to_matrix()
auc_hm <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(t(scale(leiden_auc_mat)),
                                                     name = "Scaled AUC",
                                                     column_title = "Leiden Cluster",
                                                     row_title = "Regulon",
                                                     column_title_side = "bottom",
                                                     row_names_gp = gpar(fontsize = 8),
                                                     column_dend_height = unit(3, "mm"))), wrap = TRUE)

# good
rss_6_rankplot <- SCENIC::plotRSS_oneSet(rss, setName = "6") + labs(title = "Leiden Cluster 6") + theme(plot.title = element_text(hjust = 0.5)) # cluster ID
rss_9_rankplot <- SCENIC::plotRSS_oneSet(rss, setName = "9") + labs(title = "Leiden Cluster 9") + theme(plot.title = element_text(hjust = 0.5)) # cluster ID

# good
rss_mat <- rss_data |> select(-Z) |> pivot_wider(names_from = Topic, values_from = RSS) |> bb_tbl_to_matrix()
rss_hm <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      t(scale(rss_mat)),
      name = "Scaled RSS",
      column_title = "Leiden Cluster",
      row_title = "Regulon",
      column_title_side = "bottom",
      row_names_gp = gpar(fontsize = 8),
      column_dend_height = unit(3, "mm")
    )
  ))




umap_aml_highlight <- bb_var_umap(cds_p568, "leukemia_phenotype", value_to_highlight = "AML")
umap_pvr <- bb_gene_umap(cds_p568, "Pvr")
umap_leiden <- bb_var_umap(cds_p568, "leiden", overwrite_labels = TRUE, plot_title = "Leiden Clusters")
bb_cellmeta(cds_p568) |>
  count(leiden, leukemia_phenotype)

regulon_tbl <- regulons |> enframe(name = "regulon", value = "gene_short_name") |> unnest(gene_short_name)
regulon_module_summary <- regulons |>
  enframe(name = "regulon", value = "gene_short_name") |>
  unnest(gene_short_name) |>
  left_join(bb_rowmeta(cds_p568)) |>
  count(regulon, module_labeled)

regulons_of_interest <- c("Irf8", "Irf8_extended", "Irf5_extended", "Klf4_extended", "Mafb_extended", "Tcf7l2", "Tcf7l2_extended", "Ctcf_extended", "Foxp1_extended", "Nr3c1_extended")

regulon_module_tables <- map(.x = regulons_of_interest, .f = \(x, dat = regulon_module_summary) {
  filter(dat, regulon == x) |>
    select(-regulon)
}) |> set_names(regulons_of_interest)


umap_regulons <- map(.x = regulons_of_interest,
    .f = \(x,
           dat = cds_p568,
           regulon_genes = regulon_tbl) {
      regulon_genes <- regulon_genes |> filter(regulon == x) |> select(gene_short_name, regulon)
      dat <- filter_cds(dat, genes = bb_rowmeta(dat) |> filter(gene_short_name %in% regulon_genes$gene_short_name))
      bb_gene_umap(dat, gene_or_genes = regulon_genes)
    }) |> set_names(regulons_of_interest)

umap_regulons$Irf8
umap_regulons$Irf8_extended
umap_regulons$Irf5_extended
umap_regulons$Klf4_extended
umap_regulons$Mafb_extended
umap_regulons$Tcf7l2
umap_regulons$Tcf7l2_extended
umap_regulons$Ctcf_extended
umap_regulons$Foxp1_extended
umap_regulons$Nr3c1_extended



umap_regulons_wt_aml_ball <- map(.x = regulons_of_interest,
    .f = \(x,
           dat = cds_WT_AML_bALL,
           regulon_genes = regulon_tbl) {
      regulon_genes <- regulon_genes |> filter(regulon == x) |> select(gene_short_name, regulon)
      dat <- filter_cds(dat, genes = bb_rowmeta(dat) |> filter(gene_short_name %in% regulon_genes$gene_short_name))
      bb_gene_umap(dat, gene_or_genes = regulon_genes)
    }) |> set_names(regulons_of_interest)


umap_pvr_wt_aml_ball <- bb_gene_umap(cds_WT_AML_bALL, "Pvr")

umap_regulons_wt_aml_ball$Irf8
umap_regulons_wt_aml_ball$Irf8_extended
umap_regulons_wt_aml_ball$Irf5_extended
umap_regulons_wt_aml_ball$Klf4_extended
umap_regulons_wt_aml_ball$Mafb_extended
umap_regulons_wt_aml_ball$Tcf7l2
umap_regulons_wt_aml_ball$Tcf7l2_extended
umap_regulons_wt_aml_ball$Ctcf_extended
umap_regulons_wt_aml_ball$Foxp1_extended
umap_regulons_wt_aml_ball$Nr3c1_extended

umap_leiden_consensus_wt_aml_ball <- bb_var_umap(cds_WT_AML_bALL, "leiden_best_consensus_all", overwrite_labels = FALSE, palette = brewer.pal(n = 11, name = "Set3"))
umap_aml_highlight_wt_aml_ball <- bb_var_umap(cds_WT_AML_bALL, "leukemia_phenotype", value_to_highlight = "AML")
umap_all_highlight_wt_aml_ball <- bb_var_umap(cds_WT_AML_bALL, "leukemia_phenotype", value_to_highlight = "pre-B ALL")

write_csv(regulon_tbl, fs::path(tables_out, "regulon_tbl.csv"))
