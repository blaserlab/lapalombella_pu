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


# the informative plots here
bb_var_umap(cds_combined, "partition", facet_by = "data_set", value_to_highlight = c("3", "6"))
bb_var_umap(cds_combined, "leiden.1", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "partition", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "genotype", facet_by = "data_set")

# how to count human cells per sample
bb_cellmeta(cds_combined) |>
  count(pid, genotype)

# how to actually sample the numan cells
bb_cellmeta(cds_combined) |>
  filter(data_set == "human") |>
  filter(!sample %in% c("the ones to exclude")) |>
  slice_sample(n = 1873, by = pid)

# mouse similar but use specimen variable
bb_cellmeta(cds_combined) |> glimpse()
bb_cellmeta(cds_combined) |>
  count(specimen)

# want to 1. Define regions of enrichment by human genotype.  Do this by lumping leiden.1 clusters together, using recode if you want.  2. Normalize human cell number per sample by downsampling to validate that you lumped the clusters together fairly.  Then show number of human cells of a particular genotype are really enriched as you say they are.  3.  Normalize the mouse cell numbers per specimen.  Show that the number of cells coming from clusters 3 and 6 is higher in the human comutant regions.

bb_var_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(celltype.l1_ref == "CD8 T")), "celltype.l2_ref", overwrite_labels = FALSE)

bb_cellmeta(cds_main_human_unaligned) |> glimpse()
bb_cellmeta(cds_main_human_unaligned) |> count(celltype.l1_ref)

cds_combined <- bb_rowmeta(cds_combined) |>
  mutate(exhaustion = ifelse(gene_short_name %in% some_vector_of_genes), TRUE, FALSE) |>
  select(feature_id, exhaustion) |>
  bb_tbl_to_rowdata(obj = cds_combined, min_tbl = _)


# how to fork a project from github


bb_cellmeta(cds_main_human_unaligned) |> glimpse()
bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref", overwrite_labels = TRUE)

human_aml_only_cds <- filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(celltype.l1_ref %in% c("Mono", "DC", "other")))
bb_var_umap(human_aml_only_cds, "genotype")
bb_var_umap(human_aml_only_cds, "density")
human_aml_only_cds_tm <- monocle3::top_markers(human_aml_only_cds, group_cells_by = "genotype", genes_to_test_per_group = 100, cores = 20)

View(human_aml_only_cds_tm)

human_aml_only_cds_tm |>
  count(gene_short_name, cell_group) |>
  arrange(desc(n))


orthos <- read_tsv("https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt") |>
  mutate(`DB Class Key` = as.integer(`DB Class Key`))

duplicates <- orthos |>
  filter(`Common Organism Name` %in% c("mouse, laboratory", "human")) |>
  # filter(!`DB Class Key` %in% c(48325866, 48332781)) |>
  select(`DB Class Key`, `Common Organism Name`, `Symbol`) |>
  dplyr::summarise(n = dplyr::n(), .by = c(`DB Class Key`, `Common Organism Name`)) |>
  filter(n >1L) |>
  pull(`DB Class Key`)


human_ensembl <- bb_rowmeta(cds_main_human_unaligned) |> select(human_ensembl = id, gene_short_name)
mouse_ensembl <- bb_rowmeta(cds_main) |> select(mouse_ensembl = id, gene_short_name)


human_mouse_orthos <- orthos |>
  filter(`Common Organism Name` %in% c("mouse, laboratory", "human")) |>
  filter(!`DB Class Key` %in% c(48325866, 48332781)) |>
  filter(!`DB Class Key` %in% duplicates) |>
  select(`DB Class Key`,`Common Organism Name`, `Symbol`) |>
  pivot_wider(names_from = `Common Organism Name`, values_from = Symbol) |>
  select(-`DB Class Key`) |>
  rename(mouse = `mouse, laboratory`) |>
  left_join(human_ensembl, by = c("human" = "gene_short_name")) |>
  left_join(mouse_ensembl, by = c("mouse" = "gene_short_name"))
human_mouse_orthos

left_join(human_aml_only_cds_tm, human_mouse_orthos, by = c("gene_short_name" = "human")) |>
  as_tibble() |>
  select(feature_id = mouse_ensembl, human_genotype = cell_group) |>
  filter(!is.na(feature_id)) |>
  group_by(feature_id) |>
  mutate(n = n()) |>
  filter(n == 1) |>
  ungroup()
  count(feature_id) |> arrange(desc(n))


cds_main <- left_join(human_aml_only_cds_tm, human_mouse_orthos, by = c("gene_short_name" = "human")) |>
  as_tibble() |>
  select(feature_id = mouse_ensembl, human_genotype = cell_group) |>
  filter(!is.na(feature_id)) |>
  group_by(feature_id) |>
  mutate(n = n()) |>
  filter(n == 1) |>
  ungroup() |>
  select(-n) |>
  bb_tbl_to_rowdata(cds_main, min_tbl = _)

bb_cellmeta(cds_main)
bb_rowmeta(cds_main)
bb_var_umap(cds_main, "partition", value_to_highlight = c("3", "6"))
bb_gene_umap(cds_main, gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, human_genotype) |> filter(!is.na(human_genotype)))

bb_cellmeta(cds_combined) |> glimpse()
bb_var_umap(cds_combined, "partition", facet_by = "data_set", value_to_highlight = c("3", "6"))




