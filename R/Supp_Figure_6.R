#Supplemental Figure 6
#Supp Figs D-H

# MDSC proteins of interest
prots <- c("CD33", "CD11b", "CD14", "HLA-DR")

mdsc_prot_dat <-
  bb_genebubbles(cds_main_human_unaligned,
                 genes = prots,
                 cell_grouping = c("leiden"),
                 experiment_type = "Antibody Capture",
                 return_value = "data"
  )

set.seed(123)

# Loop through each protein
for (prot in prots) {
  # Filter for the current protein and run k-means
  prot_data <- mdsc_prot_dat |>
    filter(gene_short_name == prot) |>
    select(expression, proportion)

  # Run k-means
  kmeans_result <- kmeans(prot_data, centers = 3)

  # Add cluster assignments back to the dataset
  clustered_data <- mdsc_prot_dat |>
    filter(gene_short_name == prot) |>
    mutate(cluster = as.factor(kmeans_result$cluster))

  # Plot the protein expression/proportion clustering results
  plot <- ggplot(clustered_data, aes(x = expression, y = proportion, label = leiden)) +
    geom_point(aes(color = cluster), size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    labs(
      title = paste("Expression vs. Proportion for", prot, "with K-means Clusters"),
      x = "Expression",
      y = "Proportion",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(legend.position = "right")

  # geom_text_repel()# Print the plot
  print(plot)

  # Calculate x-axis range for each cluster
  cluster_ranges <- clustered_data |>
    group_by(cluster) |>
    summarise(
      min_expression = min(expression, na.rm = TRUE),
      max_expression = max(expression, na.rm = TRUE),
      leiden_values = paste(unique(leiden), collapse = ", ") # Combine unique Leiden values
    )

  # Order clusters by max_expression and assign new cluster names
  cluster_ranges <- cluster_ranges |>
    arrange(max_expression) |>
    mutate(
      cluster = case_when(
        row_number() == 1 ~ "negative",
        row_number() == 2 ~ "low",
        row_number() == 3 ~ "high"
      )
    )

  # create df name based on the prot
  assign(
    paste0(prot, "_cluster_ranges"),
    cluster_ranges
  )
  # Print the ranges
  print(paste("X-axis ranges for", prot))
  print(cluster_ranges)
}

head(CD33_cluster_ranges)

#CD33 low and high expression clusters
cd33_clusts <- mdsc_prot_dat |>
  filter(gene_short_name == "CD33", expression > CD33_cluster_ranges |>
           filter(cluster == "low") |>
           pull(min_expression)
  ) |> pull(leiden) |> as.character()
#CD11b low and high expression clusters
cd11b_clusts <- mdsc_prot_dat |>
  filter(gene_short_name == "CD11b", expression > CD11b_cluster_ranges |>
           filter(cluster == "low") |>
           pull(min_expression)
  ) |> pull(leiden) |> as.character()
#HLA-DR negative and low expression clusters
hladr_clusts <- mdsc_prot_dat |>
  filter(gene_short_name == "HLA-DR", expression < `HLA-DR_cluster_ranges` |>
           filter(cluster == "high") |>
           pull(min_expression)
  ) |> pull(leiden) |> as.character()

# Intersection between vectors - potential mdsc leiden cluster list
mdsc_clusts <- Reduce(intersect, list(cd33_clusts, cd11b_clusts, hladr_clusts))

# is the mdsc population present in comutant patients--------------------

cds_main_human_unaligned <- bb_cellmeta(cds_main_human_unaligned) |>
  select(cell_id, leiden) |>
  mutate(mdsc_clust = ifelse(leiden %in% mdsc_clusts, "MDSC", "other")) |>
  select(-leiden) |>
  bb_tbl_to_coldata(obj = cds_main_human_unaligned, min_tbl = _)

Supp_Fig_6DE <-
  bb_var_umap(
    cds_main_human_unaligned,
    "mdsc_clust",
    cell_size = 0.1,
    foreground_alpha = 0.05,
    palette = experimental_group_palette_1,
    rasterize = F
  ) + labs(title = "MDSC Clusters") +
  bb_var_umap(
    cds_main_human_unaligned,
    "genotype",
    cell_size = 0.1,
    foreground_alpha = 0.05,
    overwrite_labels = F,
    palette = experimental_group_palette_1,
    rasterize = T
  ) + labs(title = "Genotype")


save_plot(
  filename = fs::path(figs_out, "Supp_Fig_6DE.pdf"),
  Supp_Fig_6DE,
  base_width = 7,
  base_height = 2.5
)

#mdsc cluster protein bubbles -----------------------------
mdsc_protein_dat <- bb_genebubbles(
  cds_main_human_unaligned,
  genes = c("HLA-DR", "CD33", "CD11b", "CD14"),
  cell_grouping = "mdsc_clust",
  experiment_type = "Antibody Capture",
  scale_expr = FALSE,
  expression_threshold = 0.6,
  return_value = "data"
)

Supp_Fig_6F <- ggplot(mdsc_protein_dat,
                      aes(x = mdsc_clust,
                          y = gene_short_name,
                          size = proportion,
                          fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() +
  labs(x = NULL, y = NULL, fill = "Binding", size = "Proportion") +
  theme(legend.box = "horizontal")

save_plot(
  # filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_6F.pdf"),
  plot = Supp_Fig_6F,
  base_width = 4.0,
  base_height = 2.0
)

#mdsc cluster transcript bubbles--------------
mdsc_markers <- c("ITGAM",
                  "CD14",
                  "CD33",
                  "HLA-DRA",
                  "FUT4"
                  # "CD34",
                  # "ARG1
)
mdsc_transcript_dat <- bb_genebubbles(
  cds_main_human_unaligned,
  genes = mdsc_markers,
  cell_grouping = "mdsc_clust",
  experiment_type = "Gene Expression",
  scale_expr = FALSE,
  #expression_threshold = 0.6,
  return_value = "data"
) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "ITGAM" ~ "CD11b (ITGAM)",
    gene_short_name == "FUT4" ~ "CD15 (FUT4)",
    gene_short_name == "CEACAM8" ~ "CD66b (CEACAM8)",
    TRUE ~ gene_short_name
  ))

Supp_Fig_6G <- ggplot(mdsc_transcript_dat,
                      aes(x = mdsc_clust,
                          y = gene_short_name,
                          size = proportion,
                          fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "D") +
  theme_minimal_grid() +
  labs(x = NULL, y = NULL, fill = "Expression", size = "Proportion") +
  theme(legend.box = "horizontal")
Supp_Fig_6G

save_plot(
  # filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_6G.pdf"),
  plot = Supp_Fig_6G,
  base_width = 4.75,
  base_height = 2.0
)

cellgeno_sums <- bb_cellmeta(cds_main_human_unaligned) |>
  count(genotype, name = "geno_sum")
cellcount2 <- bb_cellmeta(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                                       filter((leiden %in% mdsc_clusts)))) |>
  group_by(genotype, pid) |>
  summarise(n = n()) |>
  left_join(cellgeno_sums) |>
  mutate(total = sum(cellgeno_sums$geno_sum)) |>
  mutate(ratio = geno_sum/total) |>
  mutate(normalized_cell_frac = n/ratio)

#factor levels
cellcount2$genotype
cellcount2$genotype <- factor(cellcount2$genotype,
                              levels = c("WT", "tet2","tp53", "comutant"))
#stacked mdsc bar plot
Supp_Fig_6H <- ggplot(cellcount2,
                      aes(x = 1, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.9) +
  labs( x = NULL,
        y = "MDSC Cell Fraction (Normalized)",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = experimental_group_palette_1)

Supp_Fig_6H

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_6H.pdf"),
  plot = Supp_Fig_6H,
  base_width = 3.5,
  base_height = 3.0
)

# fisher exact test for mdsc--------------------
values <- bb_cellmeta(cds_main_human_unaligned) |>
  mutate(genotype_binary = ifelse(genotype == "comutant", "comutant", "other")) |>
  count(mdsc_clust, genotype_binary) |>
  pull(n)

mat <- matrix(values, nrow = 2)
colnames(mat) <- c("mdsc_yes", "mdsc_no")
rownames(mat) <- c("comutant", "other")
mat

fisher.test(mat[c(2:1), c(2:1)])
