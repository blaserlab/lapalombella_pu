
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
ComplexHeatmap::Heatmap(t(scale(leiden_auc_mat)))

# good
SCENIC::plotRSS_oneSet(rss, setName = "6") # cluster ID

# good
rss_mat <- rss_data |> select(-Z) |> pivot_wider(names_from = Topic, values_from = RSS) |> bb_tbl_to_matrix()
ComplexHeatmap::Heatmap(t(scale(rss_mat)))
