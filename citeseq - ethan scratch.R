bb_var_umap(cds_main_human_unaligned, "genotype", facet_by = "value")
bb_var_umap(cds_main_human_unaligned, "density", facet_by = "genotype")
bb_var_umap(cds_main_human, "density", facet_by = "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2") +
  lims(x = c(-20, 20), y = c(-20, 20))
bb_var_umap(cds_main_human, "celltype.l2_ref", overwrite_labels = TRUE)
bb_var_umap(cds_main_human, "celltype.l1_ref", overwrite_labels = TRUE)

bb_var_umap(filter_cds(cds_main_human, cells = bb_cellmeta(cds_main_human) |> filter(celltype.l1_ref %in% c("CD8 T", "CD4 T", "other T"))), "density", facet_by = "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")
bb_cellmeta(cds_main_human) |> glimpse()
bb_var_umap(cds_main_human, "genotype", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2") +
  lims(x = c(-20, 20), y = c(-20, 20))
bb_var_umap(cds_main_human_unaligned, "genotype") +
  lims(x = c(-20, 20), y = c(-20, 20))
bb_var_umap(cds_main_human, "genotype")
bb_rowmeta(cds_main_human_unaligned, experiment_type = "Antibody Capture") |> View()

mdsc_abs <- c("CD11b",
              "CD14",
              "CD33",
              "HLA-DR")

mdsc_gene_bub_dat <-
  bb_genebubbles(
    cds_main_human,
    genes = mdsc_abs,
    cell_grouping = "leiden",
    experiment_type = "Antibody Capture",
    return_value = "data"
  )

ggplot(mdsc_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid()

# gmdsc:  11b+ 33 dim/- DRlow/neg 14 neg
bb_var_umap(
  cds_main_human,
  "leiden",
  alt_dim_x = "prealignment_dim1",
  alt_dim_y = "prealignment_dim2",
  value_to_highlight = c("11", "25", "37", "49")
)

# mmdsc:  llb+ 14+ DRlow/neg
bb_var_umap(
  cds_main_human,
  "leiden",
  value_to_highlight = c("14", "22", "35", "36", "40", "59")
)


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

############################################################################################################################
#Ethan scratch
bb_cellmeta(cds_main_human_unaligned) |> glimpse()
unique(colData(cds_main_human_unaligned)$pid)
unique(colData(cds_main_human_unaligned)$leiden)
head(colData(cds_main_human_unaligned))

view(cds_main_human_unaligned_top_markers)

bb_var_umap(cds_main_human_unaligned, "genotype")#, facet_by = "value")
bb_var_umap(cds_main_human_unaligned, "density", facet_by = "genotype")

#clustering viewing
bb_var_umap(cds_main_human_unaligned, "partition", overwrite_labels = T, foreground_alpha = 0.1)

bb_var_umap(cds_main_human, "leiden", foreground_alpha = 0.1)

#by patient sample
bb_var_umap(cds_main_human_unaligned, "pid", foreground_alpha = 0.1) +
  lims(x = c(-18, 18), y = c(-18, 18))
#by geno
bb_var_umap(cds_main_human_unaligned, "genotype", foreground_alpha = 0.1) +
  lims(x = c(-20, 20), y = c(-20, 20))

bb_rowmeta(cds_main_human_unaligned, experiment_type = "Antibody Capture") |> View()

mdsc_abs <- c("CD11b",
              "CD14",
              "CD33",
              "HLA-DR")

mdsc_gene_bub_dat <-
  bb_genebubbles(
    filter_cds(cds_main_human, cells = bb_cellmeta(cds_main_human) |>
                 filter(!(celltype.l1_ref %in% c("CD8 T", "CD4 T", "other T", "B", "NK")))),
    genes = mdsc_abs,
    cell_grouping = "leiden",
    experiment_type = "Antibody Capture",
    return_value = "data"
  )
ggplot(mdsc_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c() +
  theme_minimal_grid()

#gene level mdsc markers
mdsc_markers <- c("ITGAM",
              "CD14",
              "CD33",
              "HLA-DRB1",
              "HLA-DRB3",
              "HLA-DRB4",
              "HLA-DRB5",
              "FUT4")

bb_genebubbles(
    cds_main_human_unaligned,
    genes = mdsc_markers,
    cell_grouping = "leiden",
    experiment_type = "Gene Expression",
    return_value = "data"
  ) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "ITGAM" ~ "CD11b (ITGAM)",
    gene_short_name == "FUT4" ~ "CD15 (FUT4)",
    TRUE ~ gene_short_name
  )) |>
ggplot(mapping = aes(x = leiden,
                               y = gene_short_name,
                               size = proportion,
                               fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c() +
  theme_minimal_grid() +
  labs(title = "MDSC Markers - transcipts")

#lymphoid & dendritic lineages
lineage_abs <- c("CD45", "CD3", "CD19", "CD56", "CD11c", "CD123", "HLA-DR")
#CD56 is a NK marker
#CD123 & 11c are dendritic cell markers

lin_gene_bub_dat <-
  bb_genebubbles(
    cds_main_human_unaligned,
    genes = lineage_abs,
    cell_grouping = "leiden",
    experiment_type = "Antibody Capture",
    return_value = "data"
  )

ggplot(lin_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c() +
  theme_minimal_grid()

bb_cite_umap(cds_main_human_unaligned, "CD33")
bb_cite_umap(cds_main_human_unaligned, "CD14")
bb_gene_umap(cds_main_human_unaligned, gene_or_genes = "CD14")

bb_cite_umap(cds_main_human_unaligned, "HLA-DR")
bb_cite_umap(cds_main_human_unaligned, "CD11b")
bb_gene_umap(cds_main_human_unaligned, gene_or_genes = "ITGAM") + facet_wrap(~genotype)

bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref")
bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref", facet_by = "genotype")
bb_var_umap(cds_main_human_unaligned, "celltype.l2_ref", facet_by = "genotype")

###################################################################################################
#CD155
unique(colData(cds_main_human_unaligned)$genotype)
bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "WT")), "CD155 (PVR)", plot_title = "WT - CD155 protein")
bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "tet2")), "CD155 (PVR)", plot_title = "Tet2 - CD155 protein")
bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "tp53")), "CD155 (PVR)", plot_title = "TP53 - CD155 protein")
bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "comutant")), "CD155 (PVR)", plot_title = "Comutant - CD155 protein")
