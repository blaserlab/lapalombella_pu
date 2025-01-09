conflicts_prefer(base::as.data.frame)
bb_var_umap(cds_main_human, "genotype", facet_by = "value")
bb_var_umap(cds_main_human, "density", facet_by = "genotype")
bb_var_umap(cds_main_human, "density", facet_by = "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2") +
  lims(x = c(-20, 20), y = c(-20, 20))
bb_var_umap(cds_main_human, "celltype.l2_ref", overwrite_labels = TRUE)
bb_var_umap(cds_main_human, "celltype.l1_ref", overwrite_labels = TRUE)

bb_var_umap(filter_cds(cds_main_human, cells = bb_cellmeta(cds_main_human) |> filter(celltype.l1_ref %in% c("CD8 T", "CD4 T", "other T"))), "density", facet_by = "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")
bb_cellmeta(cds_main_human) |> glimpse()
bb_var_umap(cds_main_human, "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2") +
  lims(x = c(-20, 20), y = c(-20, 20))
bb_var_umap(cds_main_human, "genotype")
bb_rowmeta(cds_main_human, experiment_type = "Antibody Capture") |> View()

mdsc_abs <- c("CD11b",
  "CD14",
  "CD33",
  "HLA-DR")

mdsc_gene_bub_dat <- bb_genebubbles(cds_main_human, genes = mdsc_abs, cell_grouping = "leiden", experiment_type = "Antibody Capture", return_value = "data")

ggplot(mdsc_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid()

# gmdsc:  11b+ 33 dim/- DRlow/neg 14 neg
bb_var_umap(cds_main_human, "leiden", value_to_highlight = c("11", "25","37", "49"))

# mmdsc:  llb+ 14+ DRlow/neg
bb_var_umap(cds_main_human, "leiden", value_to_highlight = c("14", "22", "35", "36", "40", "59"))


bb_cite_umap(cds_main_human, "CD33")
bb_cite_umap(cds_main_human, "CD14")
bb_cite_umap(cds_main_human, "HLA-DR")


bb_var_umap(cds_main_human, "genotype", facet_by = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")+
  lims(x = c(-20, 20), y = c(-20, 20))

bb_var_umap(cds_main_human, "genotype", facet_by = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")+
  lims(x = c(-20, 20), y = c(-20, 20))

bb_var_umap(cds_main_human, "genotype", overwrite_labels = TRUE, alt_label_col = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2", foreground_alpha = 0.1)+
  lims(x = c(-18, 18), y = c(-18, 18)) +
  theme(legend.position = "right")

p1 <- bb_var_umap(cds_main_human, "genotype", value_to_highlight = "WT", overwrite_labels = TRUE, alt_label_col = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2", foreground_alpha = 0.1)+
  lims(x = c(-18, 18), y = c(-18, 18)) +
  theme(legend.position = "right")

p2 <- bb_var_umap(cds_main_human, "genotype", value_to_highlight = "tp53", overwrite_labels = TRUE, alt_label_col = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2", foreground_alpha = 0.1)+
  lims(x = c(-18, 18), y = c(-18, 18)) +
  theme(legend.position = "right")

p3 <- bb_var_umap(cds_main_human, "genotype", value_to_highlight = "tet2", overwrite_labels = TRUE, alt_label_col = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2", foreground_alpha = 0.1)+
  lims(x = c(-18, 18), y = c(-18, 18)) +
  theme(legend.position = "right")

p4 <- bb_var_umap(cds_main_human, "genotype", value_to_highlight = "comutant", overwrite_labels = TRUE, alt_label_col = "celltype.l1_ref", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2", foreground_alpha = 0.1)+
  lims(x = c(-18, 18), y = c(-18, 18)) +
  theme(legend.position = "right")

p1+p2+p3+p4



bb_var_umap(cds_combined, "data_set")
bb_var_umap(cds_)
bb_gene_umap(cds_combined, "CD3E")
bb_gene_umap(cds_combined, "CD14")
bb_gene_umap(cds_combined, "CD19")
bb_var_umap(cds_combined, "genotype", facet_by = c("value", "data_set"))
bb_var_umap(cds_combined, "density", facet_by = c("genotype", "data_set"))
bb_cellmeta(cds_combined) |> glimpse()
bb_var_umap(cds_combined, "partition", facet_by = "data_set", value_to_highlight = c("3", "6"))
bb_var_umap(cds_combined, "partition", facet_by = "data_set", overwrite_labels = TRUE)

bb_var_umap(cds_main_human_unaligned, "genotype")
bb_var_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(celltype.l1_ref == "CD8 T")), "celltype.l2_ref", overwrite_labels = FALSE)

bb_cellmeta(cds_main_human_unaligned) |> glimpse()
bb_cellmeta(cds_main_human_unaligned) |> count(celltype.l1_ref)

cds_combined <- bb_rowmeta(cds_combined) |>
  mutate(exhaustion = ifelse(gene_short_name %in% some_vector_of_genes), TRUE, FALSE) |>
  select(feature_id, exhaustion) |>
  bb_tbl_to_rowdata(obj = cds_combined, min_tbl = _)

tempdir()

# how to fork a project from github
# blaseRtemplates::initialize_github("blaserlab/lapalombella.pu.datapkg", dest = "/workspace/ethan_workspace")
