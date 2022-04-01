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
