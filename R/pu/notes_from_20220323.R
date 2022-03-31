# get the top markers from the whole dataset
prealignment_top_markers |>
  group_by(cell_group) |>
  filter(cluster_method == "partition") |>
  summarise()

# plot all the cells from whole dataset
bb_var_umap(cds = cds_main, "partition")

# subset to include only wt and aml cells with old dimensions
olddims_wt_aml <- blaseRtools::filter_cds(cds = cds_main,
                        cells = bb_cellmeta(cds_main) |>
                          filter(leukemia_phenotype %in% c("AML", "No leukemia")))

# plot the wt aml cells with old dimensions
bb_var_umap(olddims_wt_aml, "partition")

colData(olddims_wt_aml)$partition_assignment <-
  recode(colData(olddims_wt_aml)$partition,
         "1" = "pre-Neu2/3",
         "2" = "Neu",
         "3" = "immNeu",
         "4" = "pre-Neu1",
         "5" = "HSC/Prog",
         "6" = "Blast-like",
         "7" = "Mono",
         "8" = "cMoP/DC",
         "9" = "B",
         "10" = "T/NK 1",
         "11" = "PC",
         "12" = "T/NK 2",
         "13" = "Th1",
         "14" = "Th2",
         "15" = "Th3",
         "16" = "plasma cells"
  )

bb_var_umap(olddims_wt_aml, "partition_assignment")

write_csv(olddims_wt_aml, file = str_glue("{tables_out}/olddims_wt_aml.csv"))

# use bb_cell_anno() using a new reference cds or seurat object to get new annotations if you want.

bb_var_umap(cds_wt_aml_marrow, "leiden_assignment", overwrite_labels = T, facet_by = "genotype")

# use the monocle3 instructions for pseudotime.https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
#  monocle3::learn_graph() and monocle3::order_cells() are the key functions.


bb_var_umap(cds_wt_aml_marrow, "partition")

