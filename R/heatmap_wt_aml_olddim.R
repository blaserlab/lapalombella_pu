# I added the assignment for top20
top20 <- tm_wt_aml_marrow |>
  filter(cluster_method == "partition") |>
  group_by(cell_group) |>
  slice_min(order_by = marker_test_q_value, n = 20) |>
  pull(gene_short_name)
bb_cellmeta(olddims_wt_aml)
colData(olddims_wt_aml)

top20

# put in a new row metadata column
# you didnt define your top20 genes here
rowData(olddims_wt_aml)$top20 <- ifelse(rowData(olddims_wt_aml)$gene_short_name %in% top20,
                                  "yes",
                                  "no")

# filter the cds and pipe into aggregate gene expression
agg_mat_wt_aml_olddim <-
  olddims_wt_aml |>
  filter_cds(genes = bb_rowmeta(olddims_wt_aml) |>
               filter(top20 == "yes")) |>
  aggregate_gene_expression(cell_group_df = bb_cellmeta(olddims_wt_aml) |>
                              select(cell_id, partition_assignment))


max(agg_mat_wt_aml_olddim)

min(agg_mat_wt_aml_olddim)


# convert from sparse to regular matrix
agg_mat_wt_aml_olddim <- as.matrix(agg_mat_wt_aml_olddim)
max(agg_mat_wt_aml_olddim)

agg_mat_wt_aml_olddim


# fix the rownames
rownames(agg_mat_wt_aml_olddim) <-
  left_join(tibble(feature_id = rownames(agg_mat_wt_aml_olddim)),
            bb_rowmeta(olddims_wt_aml)) |>
  pull(gene_short_name)
max(agg_mat_wt_aml_olddim)
# transpose and then put all of the genes (columns) on the same scale
agg_mat_wt_aml_olddim <- scale(t(agg_mat_wt_aml_olddim))



# make a list of genes you want to point out
# put whatever genes you want from top20 here
heatmap_highlights <- c(
  "S100a8",
  "S100a9",
  "Mpo",
  "Ctsg",
  "Elane",
  "Ybx1",
  "Gzmb",
  "Tox",
  "Cd8b1",
  "Vpreb1",
  "Cd3d"
)
# make the heatmap color scale
# see configs.R for what these colors are
col_fun_heatmap <-
  colorRamp2(breaks = c(min(agg_mat_wt_aml_olddim),
                        0,
                        max(agg_mat_wt_aml_olddim)),
             colors = heatmap_3_colors)


# make the annotation object
heatmap_anno_df <-
  map(
    .x = heatmap_highlights,
    .f = function(x) {
      index <- which(colnames(agg_mat_wt_aml_olddim) == x)
      return(index)
    }
  ) %>% set_names(heatmap_highlights) %>%
  bind_cols() %>%
  pivot_longer(everything()) %>%
  as.data.frame()

heatmap_gene_anno <- HeatmapAnnotation(
  foo = anno_mark(
    at = heatmap_anno_df$value,
    labels = heatmap_anno_df$name,
    labels_gp = gpar(fontsize = 8),
    padding = 1.5,
    labels_rot = 45
  ),
  which = "column"
)


# make the heatmap finally
partition_heatmap <- grid.grabExpr(draw(
  Heatmap(
    matrix = agg_mat_wt_aml_olddim, # you had the wrong matrix here
    col = col_fun_heatmap,
    name = "Expression",
    heatmap_legend_param = list(
      title_gp = gpar(fontface = "plain", fontsize = 9),
      grid_width = unit(0.14, "in"),
      labels_gp = gpar(fontsize = 8)
    ),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    column_dend_side = "bottom",
    show_row_names = T,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = F,
    top_annotation = heatmap_gene_anno,
    row_dend_gp = gpar(lwd = 0.5),
    column_dend_gp = gpar(lwd = 0.5),
    row_title = "Partition",
    column_title = "Top 20 Genes"
  )
), wrap = T)

# plot the heatmap
plot_grid(partition_heatmap)

# save the heatmap
save_plot(
  plot_grid(partition_heatmap),
  filename = file.path(figs_out, "wt_aml_heatmap.png"),
  base_width = 5.5,
  base_height = 3.5
)
