# the new cds, cds_WT_AML_bALL has all of the no_leukemia and leukemia samples in it.
# check out ?cds_WT_AML_bALL to see how it was made
# here are some example plots you may find useful
bb_var_umap(cds_WT_AML_bALL, "leiden_best_consensus_all", overwrite_labels = TRUE)

bb_var_umap(filter_cds(cds_WT_AML_bALL,
                       cells = bb_cellmeta(cds_WT_AML_bALL) |>
                         filter(leukemia_phenotype == "No leukemia")),
            "density",
            facet_by = "leukemia_phenotype") +
bb_var_umap(filter_cds(cds_WT_AML_bALL,
                       cells = bb_cellmeta(cds_WT_AML_bALL) |>
                         filter(leukemia_phenotype == "AML")),
            "density",
            facet_by = "leukemia_phenotype") +
bb_var_umap(filter_cds(cds_WT_AML_bALL,
                       cells = bb_cellmeta(cds_WT_AML_bALL) |>
                         filter(leukemia_phenotype == "pre-B ALL")),
            "density",
            facet_by = "leukemia_phenotype")



bb_plot_genes_in_pseudotime(
  cds = filter_cds(
    cds_WT_AML_bALL,
    cells = bb_cellmeta(cds_WT_AML_bALL) |> filter(
      leiden_best_consensus_all %in% c("tusi_MPP",
                                       "stumpf_Monoblasts",
                                       "stumpf_Monocytes")
    )
  ),
  gene_or_genes = c("Cd34", "Cd68"),
  pseudotime_dim = "pseudotime_leiden_best_consensus_all_tusi_MPP",
  color_cells_by = "leiden_best_consensus_all",
  legend_title = "Cluster",
  vertical_jitter = 0.1
) +
  cowplot::panel_border()
