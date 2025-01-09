source("~/whipp_workspace/lapalombella_pu/R/dependencies.R")
source("~/whipp_workspace/lapalombella_pu/R/configs.R")

###################################################################################################################
#Ethan scratch
# bb_cellmeta(cds_main_human_unaligned) |> glimpse()
# unique(colData(cds_main_human_unaligned)$pid)
# unique(colData(cds_main_human_unaligned)$leiden)
# head(colData(cds_main_human_unaligned))
#
# view(cds_main_human_unaligned_top_markers)
#
# bb_var_umap(cds_main_human_unaligned, "genotype", foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "density", facet_by = "genotype")
#
# #clustering viewing
# bb_var_umap(cds_main_human, "partition", overwrite_labels = T, foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "partition", overwrite_labels = T, foreground_alpha = 0.1)
#
# bb_var_umap(cds_main_human, "leiden", foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "leiden", foreground_alpha = 0.1)
#
# bb_var_umap(cds_main_human_unaligned, "louvain", foreground_alpha = 0.1)
#
# #by patient sample
# bb_var_umap(cds_main_human_unaligned, "pid", foreground_alpha = 0.1) +
#   lims(x = c(-18, 18), y = c(-18, 18))
# #by geno
# bb_var_umap(cds_main_human_unaligned, "genotype", foreground_alpha = 0.1) +
#   lims(x = c(-20, 20), y = c(-20, 20))
#
# #reference aligned cell annotations
# bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref", overwrite_labels = F, foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "celltype.l2_ref", overwrite_labels = F, , foreground_alpha = 0.1)
# bb_var_umap(cds_main_human, "celltype.l3_ref", overwrite_labels = T, foreground_alpha = 0.1)
# ###################################################################################################################
#
# #Protein expression eval
# bb_rowmeta(cds_main_human_unaligned, experiment_type = "Antibody Capture") |> View()
# #mdsc markers (M-MDSC: CD11b+CD33+CD14+HLA-DR-lo/neg, PMN-MDSC: CD11b+CD33+CD15+CD66b+/-HLA-DR-lo/neg)
# mdsc_abs <- c("CD11b","CD33", "HLA-DR", "CD14")
#
# unique(colData(cds_main_human_unaligned)$celltype.l1_ref)
#
# mdsc_gene_bub_dat <-
#   bb_genebubbles(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
#                               filter(!(celltype.l1_ref %in% c("B", "other T", "CD4 T", "CD8 T", "NK")))),
#                  genes = mdsc_abs,
#                  cell_grouping = c("celltype.l1_ref", "leiden"),
#                  experiment_type = "Antibody Capture",
#                  return_value = "data"
#   )
#
# #Plot MDSC protein markers across leiden clusters
# ggplot(mdsc_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
#   geom_point(pch = 21) +
#   scale_size_area() +
#   scale_fill_viridis_c(option = "A") +
#   theme_minimal_grid() + labs(title = "MDSC Markers - protein")
#
# head(mdsc_gene_bub_dat)
#
# cd33_clusts <- mdsc_gene_bub_dat |>
#    filter(gene_short_name == "CD33", expression > 0.5) |> pull(leiden) |> as.character()
# cd11b_clusts <- mdsc_gene_bub_dat |>
#   filter(gene_short_name == "CD11b", expression > 0.5) |> pull(leiden) |> as.character()
# hladr_neg_clusts <- mdsc_gene_bub_dat |>
#   filter(gene_short_name == "HLA-DR", expression < 0.5) |> pull(leiden) |> as.character()
#
# # Intersection between vectors - potential mdsc leiden cluster list
# conflicts_prefer(base::intersect)
# mdsc_clusts <- Reduce(intersect, list(cd33_clusts, cd11b_clusts, hladr_neg_clusts))
#
# bb_var_umap(cds_main_human_unaligned, "leiden", value_to_highlight = mdsc_clusts, foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "pid", foreground_alpha = 0.1, overwrite_labels = T) +
#   lims(x = c(-18, 18), y = c(-18, 18))
#
# #filtered for potential mdsc clusters
# ggplot(mdsc_gene_bub_dat |>
#          filter(leiden %in% mdsc_clusts),
#        aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
#   geom_point(pch = 21) +
#   scale_size_area() +
#   scale_fill_viridis_c(option = "A") +
#   theme_minimal_grid() + labs(title = "MDSC Markers - protein")
#
# #CD33 gene expression in Mono, DC, other, and NA cell types
# bb_gene_violinplot(
#   filter_cds(
#     cds_main_human_unaligned,
#     cells = bb_cellmeta(cds_main_human_unaligned) |>
#       filter(!(
#         celltype.l1_ref %in% c("B", "other T", "CD4 T", "CD8 T", "NK")
#       ))
#   ),
#   variable = "leiden",
#   experiment_type = "Gene Expression",
#   genes_to_plot = "CD33",
#   pseudocount = 0,
#   jitter_fill = "transparent",
#   violin_alpha = 0.55,
#   jitter_alpha = 0.1,
#   include_jitter = TRUE
# ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1)) +
#   labs(title = "CD33", x = "Leiden clusters", y = "Gene Expression")
#
# # #lymphoid & dendritic lineages
# # lineage_abs <- c("CD45", "CD3", "CD19", "CD56", "CD11c", "CD123", "HLA-DR")
# # #CD56 is a NK marker
# # #CD123 & 11c are dendritic cell markers
# #
# # lin_gene_bub_dat <-
# #   bb_genebubbles(
# #     cds_main_human_unaligned,
# #     genes = lineage_abs,
# #     cell_grouping = "leiden",
# #     experiment_type = "Antibody Capture",
# #     return_value = "data"
# #   )
# #
# # ggplot(lin_gene_bub_dat, aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
# #   geom_point(pch = 21) +
# #   scale_size_area() +
# #   scale_fill_viridis_c() +
# #   theme_minimal_grid()
#
# bb_cite_umap(cds_main_human_unaligned, "CD33")
# bb_cite_umap(cds_main_human_unaligned, "CD14")
# bb_gene_umap(cds_main_human_unaligned, gene_or_genes = "CD14")
#
# bb_cite_umap(cds_main_human_unaligned, "HLA-DR")
# bb_cite_umap(cds_main_human_unaligned, "CD11b")
# bb_gene_umap(cds_main_human_unaligned, gene_or_genes = c("ITGAM", "CD33")) + facet_wrap(~genotype)
#
# bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref")
# bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref", facet_by = "genotype")
# bb_var_umap(cds_main_human_unaligned, "celltype.l2_ref", facet_by = "genotype")

###################################################################################################
#CD155
# unique(colData(cds_main_human_unaligned)$genotype)
# bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "WT")), "CD155 (PVR)", plot_title = "WT - CD155 protein")
# bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "tet2")), "CD155 (PVR)", plot_title = "Tet2 - CD155 protein")
# bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "tp53")), "CD155 (PVR)", plot_title = "TP53 - CD155 protein")
# bb_cite_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(genotype %in% "comutant")), "CD155 (PVR)", plot_title = "Comutant - CD155 protein")

###################################################################################################
#12.31.24 MDSC identification
bb_cellmeta(cds_main_human_unaligned) |> glimpse()

# MDSC proteins of interest
prots <- c("CD33", "CD11b", "CD14", "HLA-DR")

#cell unaligned data eval
mdsc_prot_dat <-
  bb_genebubbles(#filter_cds(cds_main_human, cells = bb_cellmeta(cds_main_human) |>
                  #            filter(!(celltype.l1_ref %in% c("B", "other T", "CD4 T", "CD8 T", "NK")))),
    cds_main_human_unaligned,
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
    select(expression, proportion) # Data for clustering

  # Run k-means
  kmeans_result <- kmeans(prot_data, centers = 3) # Adjust centers as needed

  # Add cluster assignments back to the dataset
  clustered_data <- mdsc_prot_dat |>
    filter(gene_short_name == prot) |>
    mutate(cluster = as.factor(kmeans_result$cluster))

  # Plot the protein expression/proportion clustering results
  plot <- ggplot(clustered_data, aes(x = expression, y = proportion, label = leiden)) +
    geom_point(aes(color = cluster), size = 3) + # K-means cluster coloration
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) + # Add leiden labels
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
conflicts_prefer(base::intersect)
mdsc_clusts <- Reduce(intersect, list(cd33_clusts, cd11b_clusts, hladr_clusts))


mdsc_leiden_clusts <-
  bb_var_umap(
    cds_main_human_unaligned,
    "leiden",
    value_to_highlight = mdsc_clusts,
    foreground_alpha = 0.1
  ) + labs(title = "MDSC leiden clusters") +
  bb_var_umap(
    cds_main_human_unaligned,
    "genotype",
    foreground_alpha = 0.1,
    overwrite_labels = F
  ) + labs(title = "Genotype")

mdsc_leiden_clusts

#mdsc cluster protein bubbles
ggplot(mdsc_prot_dat |>
         filter(leiden %in% mdsc_clusts),
       aes(x = leiden, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() + labs(title = "MDSC Markers - protein")

#mdsc cluster transcript bubbles
mdsc_markers <- c("ITGAM",
                  "CD14",
                  "CD33",
                  "HLA-DRA",
                  "FUT4",
                  "CD34",
                  "ARG1",
                  "CEACAM8")

bb_genebubbles(
  filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
               filter((leiden %in% mdsc_clusts))),
  genes = mdsc_markers,
  cell_grouping = "leiden",
  experiment_type = "Gene Expression",
  return_value = "data"
) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "ITGAM" ~ "CD11b (ITGAM)",
    gene_short_name == "FUT4" ~ "CD15 (FUT4)",
    gene_short_name == "CEACAM8" ~ "CD66b (CEACAM8)",
    TRUE ~ gene_short_name
  )) |>
  ggplot(mapping = aes(x = leiden,
                       y = gene_short_name,
                       size = proportion,
                       fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() +
  labs(title = "MDSC Markers - transcipts")

# Stacked bar chart: mdsc population genotype composition
#normalized to cell quantity contributed from each leukemia phenotype
colData(cds_main_human_unaligned)

cellgeno_sums <- bb_cellmeta(cds_main_human_unaligned) |>
  count(genotype, name = "geno_sum")

cellcount <- bb_cellmeta(cds_main_human_unaligned) |>
  group_by(leiden, genotype) |>
  summarise(n = n()) |> filter(leiden %in% c(1, 2, 21, 26, 32, 39, 43, 60, 67, 70, 93, 95)) |>
  left_join(cellgeno_sums) |>
  mutate(total = sum(cellgeno_sums$geno_sum)) |>
  mutate(ratio = geno_sum/total) |>
  mutate(normalized_cell_frac = n/ratio/4)

#write.csv(cellcount, file = "~/network/T/Labs/EHL/Rosa/Ethan/10x/Tet2_P53/cellcount.csv")

#factor levels
cellcount$leiden <- factor(cellcount$leiden,
                              levels = c(1, 2, 21, 26, 32, 39, 43, 60, 67, 70, 93, 95)) #for all

cellcount$genotype <- factor(cellcount$genotype,
                                       levels = c("WT", "tet2","tp53", "comutant"))
#plot
mdsc_bp <-
  ggplot(cellcount,
         aes(x = leiden, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.9)

mdsc_bp

#mdsc proportions by genotype

cellcount2 <- bb_cellmeta(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                                      filter((leiden %in% mdsc_clusts)))) |>
  group_by(genotype) |>
  summarise(n = n()) |>
  left_join(cellgeno_sums) |>
  mutate(total = sum(cellgeno_sums$geno_sum)) |>
  mutate(ratio = geno_sum/total) |>
  mutate(normalized_cell_frac = n/ratio/4)

#factor levels
cellcount2$genotype <- factor(cellcount$genotype,
                             levels = c("WT", "tet2","tp53", "comutant"))
#plot
mdsc_bp <-
  ggplot(cellcount2,
         aes(x = 1, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.9) +
  scale_y_continuous(expand = c(0, 0)) +  # Removes extra spacing at the bottom
  labs(
    x = NULL,
    y = "Normalized Cell Fraction",
    title = "MDSC Cluster Genotype Fractions",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Removes x-axis labels
    axis.ticks.x = element_blank()  # Removes x-axis ticks
  )

mdsc_bp
head(cellcount2)


################################################################################
#01.02.24
# bb_cellmeta(cds_main_human_unaligned) |> glimpse()
#
# # make logical values for CD33+CD11b+ cells
# mat <- monocle3::exprs(cds_main_human_unaligned)
#
# str(mat)
#
# matrix_dir <- as.matrix(mat@matrix)
#
# gene_expression <- mat["ENSG00000105383", , drop = FALSE]
#
# gene_expression_values <- matrix_dir["ENSG00000105383", , drop = FALSE]
#
# cols_above_zero <- colnames(gene_expression_values)[gene_expression_values > 0]
#
#
# str(gene_expression)
#
# cols_above_zero <- colnames(gene_expression)[gene_expression > 0]
#
# # Check if it's a sparse matrix or a different structure
# if (inherits(gene_expression, "CsparseMatrix") || inherits(gene_expression, "RsparseMatrix")) {
#   # If it's sparse, convert to a dense vector using as.vector() specifically for sparse matrices
#   gene_expression <- as.vector(gene_expression)
# }
#
# # Find the columns where gene expression is greater than 0
# selected_columns <- colnames(mat)[gene_expression > 0]
# # Find the columns where gene expression is greater than 0
# selected_columns <- colnames(mat)[gene_expression > 0]
#
# cd33_tbl <- colnames(mat[ ,mat["ENSG00000105383", ] > 0]) |> as_tibble() |> mutate(CD33_pos = TRUE) |> rename(cell_id = value)
# mouse_cds_list[[4]] <- bb_tbl_to_coldata(mouse_cds_list[[4]], min_tbl = cd19_tbl)
