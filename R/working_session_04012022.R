# look at each specimen to see what is overlapping and not
bb_var_umap(cds_p568, "specimen")
bb_var_umap(cds_p568, "partition")
cds_p568_tm |> filter(cell_group == "partition 2") |> View()
cds_p568_tm |> filter(cell_group == "partition 4") |> View()

# calculate the differential representation of specimens into partitions
bb_cluster_representation(
  cds_p568[, colData(cds_p568)$specimen %in% c("P6", "P8")],
  cluster = "partition",
  class_var = "specimen",
  experimental_class = "P6",
  control_class = "P8",
  return_value = "table"
)

bb_cluster_representation(
  cds_p568[, colData(cds_p568)$specimen %in% c("P5", "P8")],
  cluster = "partition",
  class_var = "specimen",
  experimental_class = "P5",
  control_class = "P8",
  return_value = "table"
)

# drill down into smaller popluations
bb_var_umap(cds_p568, "leiden")
bb_var_umap(cds_p568, "leiden")
bb_var_umap(cds_p568, "louvain")
cds_subset <- filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |>
    filter(specimen %in% c("P6", "P8")) |>
    filter(partition == "2")
)

bb_var_umap(cds_subset, "louvain", overwrite_labels = T)
bb_var_umap(cds_subset, "density", facet_by = "specimen")

# inspect the top markers
cds_p568_tm |> filter(cell_group == "louvain 37") |> View()

# visualze selected top markers
bb_gene_dotplot(cds_subset,
                markers = c("Il7r", "Myb"),
                group_cells_by = "louvain")

plot_cells(cds_subset,
           color_cells_by = "louvain",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


# calculate the differential representation of specimens into partitions (P5 and 8)
bb_cluster_representation(
  cds_p568[, colData(cds_p568)$specimen %in% c("P5", "P8")],
  cluster = "partition",
  class_var = "specimen",
  experimental_class = "P5",
  control_class = "P8",
  return_value = "table"
)

cds_subset_1 <- filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |>
    filter(specimen %in% c("P5", "P8")) |>
    filter(partition == "3")
)

bb_var_umap(cds_subset_1, "louvain", overwrite_labels = T)
bb_var_umap(cds_subset_1, "density", facet_by = "specimen")

# inspect the top markers
cds_p568_tm |> filter(cell_group == "louvain 19") |> View()
cds_p568_tm |> filter(cell_group == "louvain 26") |> View()


# visualze selected top markers
bb_gene_dotplot(cds_subset,
                markers = c("Il7r", "Myb"),
                group_cells_by = "louvain")

plot_cells(cds_subset,
           color_cells_by = "louvain",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# visualze selected top markers
bb_gene_dotplot(cds_subset_1,
                markers = c("Klf4", "Itga4", "Fcgr4", "Flt3"),
                group_cells_by = "louvain")

plot_cells(cds_subset,
           color_cells_by = "louvain",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,cell_size = 1,
           trajectory_graph_segment_size = 0.5,  graph_label_size = 2)


cds_subset <- order_cells(cds_subset)
plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,cell_size = 1)

cds_subset_1 <- order_cells(cds_subset_1)
plot_cells(cds_subset_1,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,cell_size = 1)

genes_p568_Ball <- c("Il7r", "Myb")
lineage_cds_p568_Ball <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes_p568_Ball,
                   colData(cds_subset)$specimen %in% c("P6", "P8")]

plot_genes_in_pseudotime(lineage_cds_p568_Ball,
                         color_cells_by="louvain",
                         min_expr=0.5)

genes_p568_aml <- c("Flt3", "Klf4", "Gata2", "Spi1")
lineage_cds_p568_aml <- cds_subset_1[rowData(cds_subset_1)$gene_short_name %in% genes_p568_aml,
                          colData(cds_subset_1)$specimen %in% c("P5", "P8")]

plot_genes_in_pseudotime(lineage_cds_p568_aml,
                         color_cells_by="louvain",
                         min_expr=0.5)

cds_p568ann<-bb_cds_anno(query_cds=cds_p568, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )

colData(cds_p568ann)
bb_var_umap(cds_p568ann, "predicted.leiden_assignment")
bb_var_umap(cds_p568ann, "leukemia_phenotype")

