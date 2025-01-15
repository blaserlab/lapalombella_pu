suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRtemplates"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
#suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("todor"))


# run this to update the data package
blaseRtemplates::project_data("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/datapkg")

# graphical parameters####
the_font_size <-10
theme_set(theme_cowplot(font_size = the_font_size))

# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "AML" = "#DC0000",
  "WT" = "#3C5488",
  "MDSC" = brewer.pal(n = 8, name = "Dark2")[1],
  "other" = brewer.pal(n = 8, name = "Dark2")[2],
  "comutant" = brewer.pal(n = 8, name = "Dark2")[3],
  "tet2" = brewer.pal(n = 8, name = "Dark2")[4],
  "tp53" = brewer.pal(n = 8, name = "Dark2")[5],
  "WT" = "orange3"


)

experimental_group_palette_1 <- c(
  "MDSC" = brewer.pal(n = 8, name = "Dark2")[1],
  "other" = brewer.pal(n = 8, name = "Dark2")[2],
  "comutant" = brewer.pal(n = 8, name = "Dark2")[3],
  "tet2" = brewer.pal(n = 8, name = "Dark2")[4],
  "tp53" = brewer.pal(n = 8, name = "Dark2")[5],
  "WT" = brewer.pal(n = 8, name = "Dark2")[6]


)


jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#041dff","white","#ff0505")

# conflicts ---------------------------------------------------------------
# resolve conflicting function names here
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("exprs", "monocle3")
conflicts_prefer(base::as.data.frame)
conflicts_prefer(base::intersect)

# output directories

figs_out <- fs::path("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/figs/revisions_final")
tables_out <- "/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/tables"


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
    rasterize = TRUE
  ) + labs(title = "MDSC Clusters") +
  bb_var_umap(
    cds_main_human_unaligned,
    "genotype",
    cell_size = 0.1,
    foreground_alpha = 0.05,
    overwrite_labels = F,
    palette = experimental_group_palette_1,
    rasterize = TRUE
  ) + labs(title = "Genotype")


save_plot(filename = fs::path(figs_out, "mdsc_umap.pdf"), mdsc_leiden_clusts, base_width = 7, base_height = 2.5)

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
  # expression_threshold = 0.6,
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
#   ggplot(aes(x = genotype, y = mean_exhaustion_score)) +
#   geom_jitter() +
#   # geom_violin() +
#   facet_wrap( ~ celltype.l2_ref) +
#   ggpubr::stat_compare_means(comparisons = list(c("comutant", "tet2"),
#                                                 c("comutant", "tp53"),
#                                                 c("comutant", "WT")), method = "wilcox")



###############################################################################
#Do the authors see the same patterns in human patients carrying TET2+TP53 mutations?
#Similar to clusters 3 and 6 observed in mice

bb_var_umap(cds_combined, "partition", facet_by = "data_set", value_to_highlight = c("3", "6"))
bb_var_umap(cds_combined, "leiden.1", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "partition", facet_by = "data_set", overwrite_labels = TRUE)
bb_var_umap(cds_combined, "genotype", facet_by = "data_set")
bb_var_umap(cds_combined, "genotype") + bb_var_umap(cds_combined, "leiden.1")


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
