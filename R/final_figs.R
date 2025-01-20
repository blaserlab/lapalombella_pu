###################################################################################################
#12.31.24 MDSC identification
#bb_cellmeta(cds_main_human_unaligned) |> glimpse()

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
mdsc_clusts <- Reduce(intersect, list(cd33_clusts, cd11b_clusts, hladr_clusts))

# is the mdsc population present in comutant patients--------------------

cds_main_human_unaligned <- bb_cellmeta(cds_main_human_unaligned) |>
  select(cell_id, leiden) |>
  mutate(mdsc_clust = ifelse(leiden %in% mdsc_clusts, "MDSC", "other")) |>
  select(-leiden) |>
  bb_tbl_to_coldata(obj = cds_main_human_unaligned, min_tbl = _)

mdsc_leiden_clusts <-
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


save_plot(filename = fs::path(figs_out, "mdsc_clust_umap.pdf"), mdsc_leiden_clusts, base_width = 7, base_height = 2.5)

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

mdsc_proteinbubble_plot <- ggplot(mdsc_protein_dat,
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
  filename = fs::path(figs_out, "mdsc_proteinbubble_plot.pdf"),
  plot = mdsc_proteinbubble_plot,
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

mdsc_genebubble_plot <- ggplot(mdsc_transcript_dat,
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
mdsc_genebubble_plot

save_plot(
  # filename = "temp.pdf",
  filename = fs::path(figs_out, "mdsc_genebubble_plot.pdf"),
  plot = mdsc_genebubble_plot,
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
mdsc_bp_stack <- ggplot(cellcount2,
                        aes(x = 1, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.9) +  # Change position to "stack"
  labs( x = NULL,
    y = "MDSC Cell Fraction (Normalized)",
      ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Removes x-axis labels
    axis.ticks.x = element_blank()  # Removes x-axis ticks
  ) +
  scale_fill_manual(values = experimental_group_palette_1)

mdsc_bp_stack

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "mdsc_geno_proportions_barplot.pdf"),
  plot = mdsc_bp_stack,
  base_width = 3.5,
  base_height = 3.0
)

# fisher exact test for mdsc--------------------
values <- bb_cellmeta(cds_main_human_unaligned) |>
  mutate(genotype_binary = ifelse(genotype == "comutant", "comutant", "other")) |>
  count(mdsc_clust, genotype_binary) |>
  pull(n)

mat <- matrix(values, nrow = 2)
colnames(mat) <- c("mdsc_no", "mdsc_yes")
rownames(mat) <- c("comutant", "other")
mat

fisher.test(mat[c(2:1), c(2:1)])

################################################################################
#T cell exhaustion

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
    #title = "CD8 T cell Count Normalized by genotype and pid",
    x = "Genotype",
    y = "CD8 T Cell # (Normalized)",
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



cd8_prot_dat <- bb_genebubbles(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref %in% "CD8 T"))
  ),
  genes = prot_exhaust,
  cell_grouping = "genotype",
  experiment_type = "Antibody Capture",
  return_value = "data",
  scale_expr = TRUE,
  expression_threshold = 0.8
)

cd8_exhaustion_prot_bubble <-
  ggplot(
    cd8_prot_dat,
    mapping = aes(
      x = genotype,
      y = gene_short_name,
      size = proportion,
      fill = expression
    )
  ) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid() +
  labs(x = NULL, y = NULL, size = "Proportion", fill = "Binding")
cd8_exhaustion_prot_bubble

save_plot(
  # filename = "temp.pdf",
  filename = fs::path(figs_out, "cd8_exhaustion_prot_bubble.pdf"),
  plot = cd8_exhaustion_prot_bubble,
  base_width = 5,
  base_height = 3
)

#CD8 T cell exhaustion marker gene bubble - transcripts
exhaust_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2","TOX")

cd8_exh_gene_dat <- bb_genebubbles(
  filter_cds(
    cds_main_human_unaligned,
    cells = bb_cellmeta(cds_main_human_unaligned) |>
      filter((celltype.l1_ref %in% "CD8 T"))
  ),
  genes = exhaust_genes,
  cell_grouping = "genotype",
  experiment_type = "Gene Expression",
  return_value = "data"
) |>
  dplyr::mutate(
  gene_short_name = dplyr::case_when(
    gene_short_name == "HAVCR2" ~ "TIM3 (HAVCR2)",
    gene_short_name == "PDCD1" ~ "PD-1 (PDCD1)",
    TRUE ~ gene_short_name
  )
)

cd8_exh_genebubble <- ggplot(cd8_exh_gene_dat, mapping = aes(
  x = genotype,
  y = gene_short_name,
  size = proportion,
  fill = expression
)) +
  geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "D") +
  theme_minimal_grid() +
  labs(x = NULL, y = NULL, size = "Proportion", fill = "Expression")
cd8_exh_genebubble

save_plot(
  # filename = "temp.pdf",
  filename = fs::path(figs_out, "cd8_exh_genebubble.pdf"),
  plot = cd8_exh_genebubble,
  base_width = 5,
  base_height = 3
)



####################
#Exhaustion Scoring#
####################
# exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "TOX")
#
# cds_main_human_unaligned <- bb_rowmeta(cds_main_human_unaligned) |>
#   select(feature_id, gene_short_name) |>
#   mutate(exhaustion = ifelse(gene_short_name %in% exhaust_genes, "exhaustion_yes", "other")) |>
#   bb_tbl_to_rowdata(cds_main_human_unaligned, min_tbl = _)
#
# agg_exh_gene_score <- bb_aggregate(filter_cds(cds_main_human_unaligned, cells = bb_cellmeta(cds_main_human_unaligned) |> filter(celltype.l2_ref %in% c("CD8 TEM", "CD8 TCM", "CD8 Naive"))), gene_group_df = bb_rowmeta(cds_main_human_unaligned) |> select(feature_id, exhaustion)) |> t()
#
# left_join(
#   bb_cellmeta(cds_main_human_unaligned) |> filter(celltype.l2_ref %in% c("CD8 TEM", "CD8 TCM", "CD8 Naive")),
#   as_tibble(agg_exh_gene_score) |>
#     mutate(cell_id = rownames(agg_exh_gene_score))
# ) |>
#   group_by(pid, genotype, celltype.l2_ref) |>
#   summarise(mean_exhaustion_score = mean(exhaustion_yes)) |>
#   ggplot(aes(x = genotype, y = mean_exhaustion_score, fill = genotype)) +
#   geom_jitter(aes(color = genotype), width = 0.2, height = 0, size = 1, alpha = 0.8) +
#   geom_violin(alpha = 0.5) +
#   geom_errorbar(data = summary_stats,
#                 aes(x = genotype, ymin = mean - se, ymax = mean + se),
#                 width = 0.2, color = "black") +
#   facet_wrap( ~ celltype.l2_ref) +
#   ggpubr::stat_compare_means(comparisons = list(c("comutant", "tet2"),
#                                                 c("comutant", "tp53"),
#                                                 c("comutant", "WT")), method = "wilcox")
#
# colData(cds_main_human_unaligned) |> glimpse()

###############################################################################
#Do the authors see the same patterns in human patients carrying TET2+TP53 mutations?
#Similar to clusters 3 and 6 observed in mice

# helpful
# bb_var_umap(cds_combined, "leiden.1",
#             facet_by = "data_set", overwrite_labels = TRUE) / bb_var_umap(obj = cds_combined, var = "genotype", facet_by = "data_set")
#
# bb_var_umap(filter_cds(cds_combined, cells = bb_cellmeta(cds_combined) |> filter(celltype.l1_ref %in% c("DC", "Mono", "other"))), "leiden.1", overwrite_labels = TRUE) +
#
# bb_var_umap(filter_cds(cds_combined, cells = bb_cellmeta(cds_combined) |> filter(celltype.l1_ref %in% c("DC", "Mono", "other"))), "genotype")
#
# bb_cellmeta(cds_combined) |> count(celltype.l1_ref)
# bb_var_umap(filter_cds(cds_combined, cells = bb_cellmeta(cds_combined) |> filter(celltype.l1_ref == "B")), "leiden.1", overwrite_labels = TRUE)
#
# bb_var_umap(cds_combined, "celltype.l1_ref") + bb_var_umap(cds_combined, "leiden.1", overwrite_labels = TRUE)
# bb_var_umap(cds_combined, "dataset")


# need to run this---------------------
leiden.1_tibble <- bind_rows(
  tibble(
    leiden.1 = c(
      "49",
      "22",
      "81",
      "67",
      "57",
      "52",
      "91",
      "19",
      "72",
      "46",
      "27",
      "17",
      "37",
      "83",
      "61",
      "3",
      "51",
      "41",
      "40",
      "62",
      "42",
      "53",
      "35",
      "9",
      "36",
      "71",
      "60",
      "32",
      "92"
    ),
    aml_leiden_enrichment = "comutant_aml"
  ),
  tibble(
    leiden.1 = c("56", "6", "77", "14", "16", "4", "105", "2", "88"),
    aml_leiden_enrichment = "WT_aml"
  ),
  tibble(
    leiden.1 = c("39","47", "21", "89", "20", "86", "95", "69", "25"),
    aml_leiden_enrichment = "tet2_aml"
  ),
  tibble(
         leiden.1 = c("48", "26", "18", "76", "93", "87", "24", "90", "33", "5", "99", "38"),
         aml_leiden_enrichment = "tp53_aml"
  ),
  tibble(
    leiden.1 = c("106", "97", "75", "82", "13", "63", "23", "44", "7", "100"),
    aml_leiden_enrichment = "B cells"
  ),

  tibble(
    leiden.1 = c(
      "65",
      "45",
      "73",
      "103",
      "59",
      "68",
      "11",
      "98",
      "58",
      "15",
      "79",
      "102",
      "55",
      "78",
      "96",
      "28",
      "50",
      "80",
      "10",
      "85",
      "94",
      "101",
      "43",
      "66"
    ),
    aml_leiden_enrichment = "TNK"
  )
)

# colData(cds_main) |> glimpse()
# unique(colData(cds_main)$barcode)
colData(cds_combined)$aml_leiden_enrichment <- NULL

cds_combined <- left_join(
  bb_cellmeta(cds_combined),
  leiden.1_tibble,
  by = join_by(leiden.1)
) |>
  mutate(aml_leiden_enrichment = replace_na(aml_leiden_enrichment, "other")) |> select(cell_id, aml_leiden_enrichment) |>
  bb_tbl_to_coldata(obj = cds_combined, min_tbl = _)

# ms_aml_mapped_cds <- filter_cds(
#   cds_combined,
#   cells = bb_cellmeta(cds_combined) |>
#     filter(data_set == "mouse") |>
#     filter(partition %in% c(3,6)) |>
#     filter(aml_leiden_enrichment != "other")
# )

ms_mapped_aml_dat <- bb_cellmeta(filter_cds(
  cds_combined,
  cells = bb_cellmeta(cds_combined) |>
    filter(data_set == "mouse") |>
    filter(partition %in% c(3,6)) |>
    filter(aml_leiden_enrichment != "other")
)) |>
  select(cell_id, aml_leiden_enrichment) |> mutate(aml_leiden_enrichment = "Murine dKO TET2/TP53 AML")

all_cells_enrich <- bb_cellmeta(cds_combined) |> select(cell_id, aml_leiden_enrichment)

mapped_enrich <- all_cells_enrich |>
  left_join(ms_mapped_aml_dat, by = "cell_id") |>
  mutate(aml_leiden_enrichment2 = coalesce(aml_leiden_enrichment.y, aml_leiden_enrichment.x)) |>
  select(cell_id, aml_leiden_enrichment2)

# unique(mapped_enrich$aml_leiden_enrichment2)
# bb_var_umap(cds_combined, "aml_leiden_enrichment")

colData(cds_combined)$aml_leiden_enrichment2 <- NULL

cds_combined <- left_join(
  bb_cellmeta(cds_combined),
  mapped_enrich,
  by = join_by(cell_id)
) |>
   select(cell_id, aml_leiden_enrichment2) |>
  bb_tbl_to_coldata(obj = cds_combined, min_tbl = _)

#unique(colData(cds_combined)$aml_leiden_enrichment2)

colData(cds_combined)$aml_leiden_enrichment2 <- factor(
  colData(cds_combined)$aml_leiden_enrichment2,
  levels = c("B cells", "TNK", "WT_aml", "tet2_aml", "tp53_aml", "comutant_aml", "Murine dKO TET2/TP53 AML", "other")
)

# human_ms_aml_mapped_umap<- bb_var_umap(filter_cds(
#   cds_combined,
#   cells = bb_cellmeta(cds_combined) |>
#     filter(!aml_leiden_enrichment2 %in% c("other"))),
#   "aml_leiden_enrichment2", palette = experimental_group_palette_2)
#
# save_plot(
#   #filename = "temp.pdf",
#   filename = fs::path(figs_out, "human_ms_aml_mapped_umap.pdf"),
#   plot = human_ms_aml_mapped_umap,
#   base_width = 6,
#   base_height = 4
# )

human_ms_aml_mapped_umap_NoTB <- bb_var_umap(filter_cds(
  cds_combined,
  cells = bb_cellmeta(cds_combined) |>
    filter(!aml_leiden_enrichment2 %in% c("other", "B cells", "TNK"))),
  "aml_leiden_enrichment2",
    palette = experimental_group_palette_2,
  rasterize = T)

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "human_ms_aml_mapped_umap_NoTB.pdf"),
  plot = human_ms_aml_mapped_umap_NoTB,
  base_width = 6,
  base_height = 4
)

# good -----------------
colData(cds_combined)$aml_leiden_enrichment <- factor(
  colData(cds_combined)$aml_leiden_enrichment,
  levels = c("B cells", "TNK", "WT_aml", "tet2_aml", "tp53_aml", "comutant_aml")
)

# human_ms_geno_aml_umap <- bb_var_umap(
#   filter_cds(
#     cds_combined,
#     cells = bb_cellmeta(cds_combined) |> filter(aml_leiden_enrichment != "other")
#   ),
#   "aml_leiden_enrichment",
#   foreground_alpha = 0.5
# )

# save_plot(
#   #filename = "temp.pdf",
#   filename = fs::path(figs_out, "human_ms_geno_aml_umap.pdf"),
#   plot = human_ms_geno_aml_umap,
#   base_width = 6,
#   base_height = 4
# )

# good - filtered out the "other" cells---------------------
ms_all_aml_dat <- bb_cellmeta(filter_cds(
  cds_combined,
  cells = bb_cellmeta(cds_combined) |>
    filter(data_set == "mouse") |>
    filter(partition %in% c(3,6)) #|>
    #filter(aml_leiden_enrichment != "other")
)) |>
  select(cell_id) |> mutate(dko_aml = "Murine dKO TET2/TP53 AML")

colData(cds_combined)$dko_aml <- NULL

cds_combined <- left_join(
  bb_cellmeta(cds_combined),
  ms_all_aml_dat,
  by = join_by(cell_id)
) |>
  select(cell_id, dko_aml) |>
  bb_tbl_to_coldata(obj = cds_combined, min_tbl = _)

human_ms_umap2 <- bb_var_umap(obj = filter_cds(
  cds_combined,
  cells = bb_cellmeta(cds_combined) |> filter(data_set == "human") #|> filter(aml_leiden_enrichment != "other")
),
var = "genotype",
foreground_alpha = 0.5,
palette = experimental_group_palette_1,
rasterize = T) +labs(title = "Human") +
  bb_var_umap(
  filter_cds(
    cds_combined,
    cells = bb_cellmeta(cds_combined) |> filter(data_set == "mouse") #|> filter(aml_leiden_enrichment != "other")
  ),
  "dko_aml",
  value_to_highlight = "Murine dKO TET2/TP53 AML",
  #pallete = experimental_group_palette_2,
  rasterize = T
) +labs(title = "Mouse")

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "human_ms_umap2.pdf"),
  plot = human_ms_umap2,
  base_width = 10,
  base_height = 4.0
)

# good
# human_ms_mapping_QC_barplot <- bb_cellmeta(filter_cds(
#   cds_combined,
#   cells = bb_cellmeta(cds_combined) |> filter(aml_leiden_enrichment != "other"))) |>
#   # filter(aml_leiden_enrichment %in% c("comutant_aml", "other_aml")) |>
#   filter(data_set == "human") |>
#   count(aml_leiden_enrichment, genotype) |>
#   ggplot(aes(x = aml_leiden_enrichment, y = n, fill = genotype)) +
#   geom_bar(stat = "identity", position = "fill") + labs(y = "Proportion")
#
# save_plot(
#   #filename = "temp.pdf",
#   filename = fs::path(figs_out, "human_ms_mapping_QC_barplot.pdf"),
#   plot = human_ms_mapping_QC_barplot,
#   base_width = 6,
#   base_height = 4.0
# )

# good
human_ms3.6_region_enrich_bp <- bb_cellmeta(cds_combined) |>
  filter(aml_leiden_enrichment %in% c("comutant_aml", "tet2_aml", "WT_aml", "tp53_aml")) |>
  filter(data_set == "mouse") |>
  filter(partition %in% c("3", "6")) |>
  count(aml_leiden_enrichment, genotype) |>
  ggplot(aes(x = genotype, y = n, fill = aml_leiden_enrichment)) +
  geom_bar(stat = "identity", position = "fill") + labs(y = "Proportion") +
  scale_fill_manual(values = experimental_group_palette_2)

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "human_ms3.6_region_enrich_bp.pdf"),
  plot = human_ms3.6_region_enrich_bp,
  base_width = 4,
  base_height = 4.0
)

data <- bb_cellmeta(cds_combined) |>
  filter(aml_leiden_enrichment %in% c("comutant_aml", "tet2_aml", "WT_aml", "tp53_aml")) |>
  filter(data_set == "mouse") |>
  filter(partition %in% c("3", "6")) |>
  mutate(new_group = ifelse(aml_leiden_enrichment == "comutant_aml", "comutant_aml", "other_aml")) |>
  count(new_group) |>
  pivot_wider(names_from = new_group, values_from = n)

binom.test(x = data$comutant_aml, data$comutant_aml + data$other_aml, p = 0.25)
