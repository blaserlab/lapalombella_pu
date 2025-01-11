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
cellcount2$genotype
cellcount2$genotype <- factor(cellcount2$genotype,
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
#generatung new logical col data for aggregate gene expression
# nadeu_11b <- readxl::read_excel("~/network/X/Labs/Blaser/share/collaborators/lapalombella_whipp_network/queries/41591_2022_1927_MOESM3_ESM.xlsx", sheet = "Supplementary Table 11b", skip = 5, col_names = c("feature_id", "gene_short_name", "mean", "l2fc", "se", "p", "padj", "direction"))
#
# #murine DKO partition cluster 3/6 signature
# cds_m <- nadeu_11b |>
#   filter(direction == "Up") |>
#   filter(padj < 0.05) |>
#   mutate(feature_id = str_remove(feature_id, "\\..*")) |>
#   mutate(nadeu_RT_gene = TRUE) |>
#   bb_tbl_to_rowdata(obj = cds_main, min_tbl = _)
#
# bb_gene_umap(
#   filter_cds(
#     cds_main,
#     cells = bb_cellmeta(cds_main) |> filter(partition_assignment_1 == "B")
#   ),
#   gene_or_genes = bb_rowmeta(cds_main) |> select(feature_id, nadeu_RT_gene)
# ) +
#   facet_grid(row = vars(patient), col = (vars(disease_tissue)))

################################################################################
#T cell exhaustion
bb_var_umap(cds_main_human_unaligned, "genotype", foreground_alpha = 0.1)
bb_var_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((
        celltype.l1_ref %in% c("CD8 T", "CD4 T", "other T")
      ))
  ),
  "genotype",
  overwrite_labels = F,
  foreground_alpha = 0.1
) + labs(title = "T cells")

bb_var_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((
        celltype.l1_ref %in% c("CD8 T")
      ))
  ),
  "genotype",
  foreground_alpha = 0.1
) + labs(title = "CD8 T cells") +theme_minimal()

bb_cite_umap(cds_main_human_unaligned, "CD8") + theme_minimal() +labs(title = "CD8 Antibody Capture")

#CD8 T cell proportions - normalized
# Stacked bar chart: CD8 T population genotype composition
#normalized to cell quantity contributed from each leukemia genotype
cellgeno_sums <- bb_cellmeta(cds_main_human_unaligned) |>
  count(genotype, name = "geno_sum") #geno_sum - total cell count for each genotype

total_geno_sum <- sum(cellgeno_sums$geno_sum)  # total cells across all genotypes

cellcount <- bb_cellmeta(cds_main_human_unaligned) |>
  group_by(celltype.l1_ref, genotype, pid) |>
  summarise(n = n()) |> filter(celltype.l1_ref %in% c("CD8 T")) |>
  left_join(cellgeno_sums) |>
  mutate(ratio = geno_sum/total_geno_sum) |> #ratio- proportion of cells contributed from ea geno
  mutate(normalized_cell_frac = (n/ratio)/4) #scaled count of cells for each celltype.l1_ref & genotype

#factor levels
cellcount$genotype <- factor(cellcount$genotype,
                             levels = c("comutant","tet2","tp53", "WT"))

#CD8 T cell proportions stacked bar plot
t_prop_bp <-
  ggplot(cellcount,
         aes(x = celltype.l1_ref, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill",
           stat = "identity",
           width = 0.9) + labs(title = "CD8 T cell proportions", subtitle = "cell type & genotype scaled cell counts") +
  theme_minimal() +
  theme(#axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

t_prop_bp

#CD8 T cell proportions - normalized
# Stacked bar chart: CD8 T population sample composition
#normalized to cell quantity contributed from each sample

# Calculate the total cell count for each genotype
geno_sums <- bb_cellmeta(cds_main_human_unaligned) |>
  group_by(genotype) |>
  summarise(total_geno_sum = n(), .groups = "drop")  # Total cell count for each genotype

# Calculate the total cell count for each patient (pid)
pid_sums <- bb_cellmeta(cds_main_human_unaligned) |>
  group_by(pid) |>
  summarise(pid_sum = n(), .groups = "drop")  # Total cell count for each patient

# Calculate normalized cell fractions
count <- bb_cellmeta(cds_main_human_unaligned) |>
  group_by(celltype.l1_ref, genotype, pid) |>
  summarise(n = n(), .groups = "drop") |>
  filter(celltype.l1_ref %in% c("CD8 T")) |>
  left_join(geno_sums, by = "genotype") |>
  left_join(pid_sums, by = "pid") |>
  mutate(ratio = pid_sum / total_geno_sum) |> # Ratio: proportion of cells contributed by the patient to the genotype
  mutate(normalized_cell_frac = n / ratio)  # Scaled count of CD8 T cells by genotype and pid

# Scatter plot of normalized_cell_frac for each pid, grouped by genotype
ggplot(count, aes(x = genotype, y = normalized_cell_frac, color = pid)) +
  geom_point(position = position_dodge(width = 0.4)) +
  # geom_bar(position = "fill",
  #          stat = "identity",
  #          width = 0.9) +
  labs(
    title = "CD8 T cell Count Normalized by genotype and pid",
    x = "Genotype",
    y = "Normalized Cell #",
    color = "Patient ID"
  ) +
  theme_minimal() +  # Apply minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

#Exhaustion marker eval
#CD8 T cell exhaustion marker gene bubble - protein
prot_exhaust <- c(
  "CD279 (PD-1)",
  "CD152 (CTLA-4)",
  "TIGIT (VSTM3)",
  "CD244 (2B4)",
  "CD272 (BTLA)",
  "CD223 (LAG-3)"
)

mdsc_prot_dat <-
  bb_genebubbles(cds_main_human_unaligned,
                 genes = prots,
                 cell_grouping = c("leiden"),
                 experiment_type = "Antibody Capture",
                 return_value = "data"
  )


cd8_prot_bub <- bb_genebubbles(
  filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
               filter((celltype.l1_ref %in% "CD8 T"))),
  genes = prot_exhaust,
  cell_grouping = "genotype",
  experiment_type = "Antibody Capture",
  return_value = "data"
) |>
ggplot(mapping = aes(x = genotype, y = gene_short_name, size = proportion, fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() + labs(title = "CD8 T cell Exhaustion Marker Expression - protein")
cd8_prot_bub
#CD8 T exhaustion marker genes - upregulated inhibitory receptors & TF TOX
bb_gene_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref %in% c("CD8 T")))
  ),
  gene_or_genes = c(
    "PDCD1",
    "CTLA4",
    "LAG3",
    "TIGIT",
    "HAVCR2",
    "TOX")) +
  facet_wrap(~ genotype) +
  labs(title = "CD8 T cell Exhaustion Marker Aggregate Transcript Expression",
       subtitle = "Marker Genes: PD-1 (PDCD1), CTLA4, TIGIT, TIM3 (HAVCR2),TOX") +
  theme_minimal()

#CD8 T cell exhaustion marker gene bubble - transcripts
exhaust_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2","TOX")

bb_genebubbles(
  filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
               filter((celltype.l1_ref %in% "CD8 T"))),
  genes = exhaust_genes,
  cell_grouping = "genotype",
  experiment_type = "Gene Expression",
  return_value = "data"
) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "HAVCR2" ~ "TIM3 (HAVCR2)",
    gene_short_name == "PDCD1" ~ "PD-1 (PDCD1)",
    TRUE ~ gene_short_name
  )) |>
  ggplot(mapping = aes(x = genotype,
                       y = gene_short_name,
                       size = proportion,
                       fill = expression)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() +
  labs(title = "CD8 T cell Exhaustion Marker Expression - transcipts")

# bb_var_umap(cds_main_human_unaligned, "genotype", foreground_alpha = 0.1)
# bb_var_umap(cds_main_human_unaligned, "celltype.l1_ref")
# bb_var_umap(cds_main_human_unaligned, "celltype.l2_ref")
# bb_var_umap(cds_main_human_unaligned, "celltype.l3_ref")

#protein level exhaustion markers
bb_rowmeta(cds_main_human_unaligned, experiment_type = "Antibody Capture") |>
  View()

#T cell exhaustion - inhibitory receptors
bb_cite_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref == "CD8 T" & genotype == "comutant"))
  ),
  antibody = c(
    "CD279 (PD-1)",
    "CD152 (CTLA-4)",
    "TIGIT (VSTM3)",
    "CD244 (2B4)",
    "CD272 (BTLA)",
    "CD223 (LAG-3)"),
  cell_size = 0.5, plot_title = "comutant") + theme_minimal() +
  lims(x = c(-2.5, 18), y = c(-10, 15))

bb_cite_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref == "CD8 T" & genotype == "tet2"))
  ),
  antibody = c(
    "CD279 (PD-1)",
    "CD152 (CTLA-4)",
    "TIGIT (VSTM3)",
    "CD244 (2B4)",
    "CD272 (BTLA)",
    "CD223 (LAG-3)"),
  cell_size = 0.5) +
  lims(x = c(-2.5, 18), y = c(-10, 15))

bb_cite_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref == "CD8 T" & genotype == "tp53"))
  ),
  antibody = c(
    "CD279 (PD-1)",
    "CD152 (CTLA-4)",
    "TIGIT (VSTM3)",
    "CD244 (2B4)",
    "CD272 (BTLA)",
    "CD223 (LAG-3)"),
  cell_size = 0.5) +
  lims(x = c(-2.5, 18), y = c(-10, 15))

bb_cite_umap(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref == "CD8 T" & genotype == "WT"))
  ),
  antibody = c(
    "CD279 (PD-1)",
    "CD152 (CTLA-4)",
    "TIGIT (VSTM3)",
    "CD244 (2B4)",
    "CD272 (BTLA)",
    "CD223 (LAG-3)"),
  cell_size = 0.5) +
  lims(x = c(-2.5, 18), y = c(-10, 15))

blaseRtools::bb_cite_umap
#Treg markers
bb_gene_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                          filter((celltype.l1_ref == "CD8 T"))), gene_or_genes = c("FOXP3", "IL2RA"))
#genes down regulated in exhausted CD8 T
bb_gene_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                          filter((celltype.l1_ref == "CD8 T" & genotype == "comutant"))), gene_or_genes = c("GZMB", "PRF1", "IL2", "TNF", "IFNG", "HNF1A")) +
  lims(x = c(-2.5, 18), y = c(-10, 15)) + labs(title = "comutant")+theme_minimal()

bb_gene_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                          filter((celltype.l1_ref == "CD8 T" & genotype == "tet2"))), gene_or_genes = c("GZMB", "PRF1", "IL2", "TNF", "IFNG", "HNF1A")) +
  lims(x = c(-2.5, 18), y = c(-10, 15)) + labs(title = "tet2")+theme_minimal()

bb_gene_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                          filter((celltype.l1_ref == "CD8 T" & genotype == "tp53"))), gene_or_genes = c("GZMB", "PRF1", "IL2", "TNF", "IFNG", "HNF1A")) +
  lims(x = c(-2.5, 18), y = c(-10, 15)) + labs(title = "tp53")+theme_minimal()

bb_gene_umap(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                          filter((celltype.l1_ref == "CD8 T" & genotype == "WT"))), gene_or_genes = c("GZMB", "PRF1", "IL2", "TNF", "IFNG", "HNF1A")) +
  lims(x = c(-2.5, 18), y = c(-10, 15)) + labs(title = "WT")+theme_minimal()

####################
#Exhaustion Scoring#
####################
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "TOX")

#extract expression mat
expr_mat <- exprs(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                               filter((celltype.l1_ref == "CD8 T"))))

#calc mean expression across the genes for each cell
exhaustion_scores <- colMeans(expr_mat, na.rm = TRUE)

# Add the exhaustion score to colData
cd8_cds <- filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |>
                        filter((celltype.l1_ref == "CD8 T")))
colData(cd8_cds)$Exhaustion_Score <- exhaustion_scores

colData(cd8_cds)

bb_var_umap(cd8_cds, "Exhaustion_Score", facet_by = "genotype")+
  labs(title = "Exhaustion Score across CD8 T Cells")
####################
bb_gene_violinplot(
    cd8_cds,
    variable = "leiden",
    experiment_type = "Gene Expression",
    genes_to_plot = "CD33",
    pseudocount = 0,
    jitter_fill = "transparent",
    violin_alpha = 0.55,
    jitter_alpha = 0.1,
    include_jitter = TRUE
  ) + theme(axis.title.y = element_blank()) + theme(strip.text = element_blank()) + plot_layout(heights = c(1,-0.1 , 1)) +
    labs(title = "CD33", x = "Leiden clusters", y = "Gene Expression")

###############################################################################
#Do the authors see the same patterns in human patients carrying TET2+TP53 mutations?
#Similar to clusters 3 and 6 observed in mice

bb_var_umap(cds_combined, "partition", facet_by = "data_set", value_to_highlight = c("3", "6"))
bb_var_umap(cds_combined, "leiden.1", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "partition", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "genotype", facet_by = "data_set")


#Pull leiden.1 cluster values from mouse cells from partitions 3 & 6
leiden.1_ms3.6_equiv <-
  bb_cellmeta(filter_cds(
    cds_combined,
    cells = bb_cellmeta(cds_combined) |> filter(data_set %in% c("mouse") &
                                                  partition == c("3", "6"))
  )) |> count(leiden.1) |>
  mutate(percentage = n / sum(n) * 100)

leiden.1_ms3.6_equiv <- leiden.1_ms3.6_equiv |> filter(percentage > 0.5) |> pull(leiden.1)

#Plot highlighting these clusters
bb_var_umap(cds_combined, "leiden.1", facet_by = "data_set", value_to_highlight = leiden.1_ms3.6_equiv) + theme_minimal()
bb_var_umap(cds_combined, "genotype", facet_by = "data_set")+theme_minimal()
bb_gene_umap(cds_combined, gene_or_genes = c("CD34", "CD33")) + facet_wrap(~ data_set) +theme_minimal() + labs(title = "CD33/CD34 aggregate gene expression")


# count human cells per sample
bb_cellmeta(cds_combined) |>
  count(pid, genotype) |> arrange(n)

# Down sample the human cells
set.seed(123)

filter_cds(cds_combined, bb_cellmeta(cds_combined) |>
  filter(data_set == "human") |>
  filter(!pid %in% c("U18-6524", "U22-0332", "U18-3620", "U17-2250")) |>
  slice_sample(n = 1873, by = pid))

#Before
bb_var_umap(filter_cds(cds_combined, bb_cellmeta(cds_combined) |>
                         filter(data_set == "human") |>
                         filter(!pid %in% c("U18-6524", "U22-0332", "U18-3620", "U17-2250"))),
            "leiden.1",
            facet_by = "data_set",
            value_to_highlight = leiden.1_ms3.6_equiv) + theme_minimal()
#After down sampling
bb_var_umap(filter_cds(cds_combined, bb_cellmeta(cds_combined) |>
                         filter(data_set == "human") |>
                         filter(!pid %in% c("U18-6524", "U22-0332", "U18-3620", "U17-2250")) |>
                         slice_sample(n = 1873, by = pid)),
            "leiden.1",
            facet_by = "data_set",
            value_to_highlight = leiden.1_ms3.6_equiv) + theme_minimal()

#Show number of human cells of a particular genotype are really enriched as you say they are.
bb_cellmeta(filter_cds(cds_combined, bb_cellmeta(cds_combined) |>
                         filter(data_set == "human") |>
                         filter(!pid %in% c("U18-6524", "U22-0332", "U18-3620", "U17-2250")) |>
                         slice_sample(n = 1873, by = pid))) |>
  filter(leiden.1 %in% leiden.1_ms3.6_equiv)

# mouse similar but use specimen variable
bb_cellmeta(cds_combined) |> glimpse()
bb_cellmeta(cds_combined) |> filter(data_set == "mouse") |>
  count(specimen, genotype) |> arrange(n)


# want to 1. Define regions of enrichment by human genotype.
#Do this by lumping leiden.1 clusters together, using recode if you want.
#2. Normalize human cell number per sample by downsampling to validate that you lumped the clusters
#together fairly.
#Then show number of human cells of a particular genotype are really enriched as you say they are.
#3.  Normalize the mouse cell numbers per specimen.
#Show that the number of cells coming from clusters 3 and 6 is higher in the human comutant regions.

#TODO check to make sure this was performed correctly

#Down sample the human cells
set.seed(123)
hm_aligned_count <- bb_cellmeta(filter_cds(cds_combined, bb_cellmeta(cds_combined) |>
                                       filter(data_set == "human") |>
                                       filter(!pid %in% c("U18-6524", "U22-0332", "U18-3620", "U17-2250")) |>
                                       slice_sample(n = 1873, by = pid))) |>
  filter(leiden.1 %in% leiden.1_ms3.6_equiv) |>
  group_by(genotype, pid) |>
  summarise(n = n()) |>
  ggplot(mapping = aes(x = genotype, y = n, fill = genotype)) +
  geom_violin(trim = TRUE, alpha = 0.6) +  # Violin plot layer
  geom_jitter(aes(color = pid), width = 0.2, alpha = 0.8) +  # Add jitter for individual points
  labs(
    title = "dKO AML Mouse Clusters 3 & 6 Aligned Human Cells",
    subtitle = "Down sampled for comparison",
    x = "Genotype",
    y = "Mouse Aligned Normalized Human Cell #",
    fill = "Genotype",
    color = "Patient ID"
  ) +
  theme_minimal() +  # Apply minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "right"  # Position the legend
  )

hm_aligned_count
