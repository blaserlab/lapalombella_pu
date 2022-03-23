# get the top20 of all partitions, ranked by q

top20 <- prealignment_top_markers |>
  filter(cluster_method == "partition") |>
  group_by(cell_group) |>
  slice_min(order_by = marker_test_q_value, n = 20) |>
  pull(gene_short_name)

# put in a new row metadata column

rowData(cds_main)$top20 <- ifelse(rowData(cds_main)$gene_short_name %in% top20,
                                  "yes",
                                  "no")

# filter the cds and pipe into aggregate gene expression
agg_mat_partition <- cds_main |>
  filter_cds(genes = bb_rowmeta(cds_main) |>
               filter(top20 == "yes")) |>
  aggregate_gene_expression(cell_group_df = bb_cellmeta(cds_main) |>
                              select(cell_id, partition))

# convert from sparse to regular matrix
agg_mat_partition <- as.matrix(agg_mat_partition)

# fix the rownames
rownames(agg_mat_partition) <-
  left_join(tibble(feature_id = rownames(agg_mat_partition)),
            bb_rowmeta(cds_main)) |>
  pull(gene_short_name)


# transpose and then put all of the genes (columns) on the same scale
agg_mat_partition <- scale(t(agg_mat_partition))

# make a list of genes you want to point out
# put whatever genes you want from top20 here
heatmap_highlights <- c(
  "Mpeg1",
  "Ikzf2",
  "Cd3g",
  "Cd34",
  "Gata2"
)


# make the heatmap color scale
# see configs.R for what these colors are
col_fun_heatmap <-
  colorRamp2(breaks = c(min(agg_mat_partition),
                        0,
                        max(agg_mat_partition)),
             colors = heatmap_3_colors)

# make the annotation object
heatmap_anno_df <-
  map(
    .x = heatmap_highlights,
    .f = function(x) {
      index <- which(colnames(agg_mat_partition) == x)
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
    matrix = agg_mat_partition,
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

