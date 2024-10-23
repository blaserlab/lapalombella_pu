#Revisions
########################################################################################################################
aml <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12')))
########################################################################################################################

#HSC/MPP marker gene bubble plot

#Lineage Markers gene bubble plot
bb_genebubbles(
  aml,
  genes = c(
    "Ptprc", #Cd45
    "Cd3e", #T cell marker
    "Cd19", #B cell marker
    "Cd11b", #myeloid/granulocyte marker
    "Ly6g", #(Gr-1) neutrophil/granulocyte marker
    "Ter119", #erythroid lineage marker
    "Itga2b", #(Cd41) megakaryocyte marker
    "Cd8a", #T cell marker
    "Cd4", #T cell marker
    "Ncr1" #(Nk1.1) NK cell marker
  ),
  cell_grouping = c("leiden_assignment2","partition"),
  return_value = "data"
) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "Ptprc" ~ "Ptprc (Cd45)",
    gene_short_name == "Itga2b" ~ "Itga2b (Cd41)",
    gene_short_name == "Ly6g" ~ "Ly6g (Gr-1)",
    gene_short_name == "Ncr1" ~ "Ncr1 (Nk1.1)",
    TRUE ~ gene_short_name
  )) |>
  ggplot(mapping = aes(x = partition,
                       y = gene_short_name,
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~leiden_assignment2, scales = "free_x") +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(title = "AML Leiden/Partition Clusters - Lineage Markers", x = NULL, y = NULL, size = "Proportion", color = "Expression")

bb_genebubbles(
  aml,
  genes = c(
    "Flt3", #(Flk2)
    "Slamf1", #(Cd150)
    "Cd48",
    "Ly6a", #(Sca-1)
    "Kit"
  ),
  cell_grouping = c("leiden_assignment2","partition"),
  return_value = "data"
) |>
  dplyr::mutate(gene_short_name = dplyr::case_when(
    gene_short_name == "Slamf1" ~ "Slamf1 (Cd150)",
    gene_short_name == "Ly6a" ~ "Ly6a (Sca-1)",
    TRUE ~ gene_short_name
  )) |>
  ggplot(mapping = aes(x = partition,
                       y = gene_short_name,
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~leiden_assignment2, scales = "free_x") +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(title = "AML Leiden/Partition Clusters - HSC/MPP Markers", x = NULL, y = NULL, size = "Proportion", color = "Expression")
##############################################################################################################################################
#Sample summary
sample_numbers <- bb_cellmeta(cds_main) |>
  group_by(genotype, leukemia_phenotype, primary_or_engraftment, tissue, sample) |>
  summarise(n = n())
kable(sample_numbers)

F6D_by.sample <- bb_var_umap(blaseRtools::filter_cds(cds = cds_main,
                                                     cells = bb_cellmeta(cds_main) |>
                                                       filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12'))), "leiden_assignment2", overwrite_labels = F, group_label_size = 35, facet_by = "sample")
F6D_by.sample
######################################################################################################################################################################################################
#Lineage MArkers
bb_genebubbles(
  filter_cds(
    cds = cds_p568,
    cells = bb_cellmeta(cds_p568) |>
      filter(louvain %in% c(53, 48, 4, 21, 33, 38, 10, 49, 24, 9, 34, 19, 26, 20, 51, 15, 31)
    ) |> filter(leukemia_phenotype %in% c("AML", "No leukemia"))
  ),
  genes = c(
    "Ptprc",
    #Cd45
    "Cd3e",
    #T cell marker
    "Cd19",
    #B cell marker
    "Cd11b",
    #myeloid/granulocyte marker
    "Ly6g",
    #(Gr-1) neutrophil/granulocyte marker
    "Ter119",
    #erythroid lineage marker
    "Itga2b",
    #(Cd41) megakaryocyte marker
    "Cd8a",
    #T cell marker
    "Cd4",
    #T cell marker
    "Ncr1" #(Nk1.1) NK cell marker
  ),
  cell_grouping = c("louvain", "leukemia_phenotype"),
  return_value = "data"
) |>
  dplyr::mutate(
    gene_short_name = dplyr::case_when(
      gene_short_name == "Ptprc" ~ "Ptprc (Cd45)",
      gene_short_name == "Itga2b" ~ "Itga2b (Cd41)",
      gene_short_name == "Ly6g" ~ "Ly6g (Gr-1)",
      gene_short_name == "Ncr1" ~ "Ncr1 (Nk1.1)",
      TRUE ~ gene_short_name
    )
  ) |>
  ggplot(mapping = aes(
    x = louvain,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~ leukemia_phenotype, scales = "free_x") +
  theme_minimal_grid(font_size = the_font_size) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(
    title = "Lineage Markers - Louvain Clusters",
    x = NULL,
    y = NULL,
    size = "Proportion",
    color = "Expression"
  )

#hCD16 = FCGR3A (mouse: Fcgr4 & Fcgr3) & FCGR3B
#hCD32 = FCGR2A (mouse: Fcgr2b) & FCGR2B

#GMP Marker Genes - Louvain Clusters
bb_genebubbles(
  filter_cds(
    cds = cds_p568,
    cells = bb_cellmeta(cds_p568) |> filter(
      louvain %in% c(53, 48, 4, 21, 33, 38, 10, 49, 24, 9, 34, 19, 26, 20, 51, 15, 31)
    ) |> filter(leukemia_phenotype %in% c("AML", "No leukemia"))
  ),
  genes = c("Fcgr4", #CD16
            "Fcgr3", #CD16
            "Fcgr2b", #CD32
            "Ly6a", #(Sca-1)
            "Kit",
            "Cd34"),
  cell_grouping = c("louvain", "leukemia_phenotype"),
  return_value = "data"
) |>
  dplyr::mutate(
    gene_short_name = dplyr::case_when(
      gene_short_name == "Fcgr4" ~ "Fcgr4 (Cd16)",
      gene_short_name == "Fcgr3" ~ "Fcgr3 (Cd16)",
      gene_short_name == "Fcgr2b" ~ "Fcgr2b (Cd32)",
      gene_short_name == "Ly6a" ~ "Ly6a (Sca-1)",
      TRUE ~ gene_short_name
    )
  ) |>
  ggplot(mapping = aes(
    x = louvain,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~ leukemia_phenotype, scales = "free_x") +
  theme_minimal_grid(font_size = the_font_size) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(
    title = "GMP Markers - Louvain Clusters",
    x = NULL,
    y = NULL,
    size = "Proportion",
    color = "Expression"
  )

#lineage negative subset of louvain cds_p568 clusters:
##GMP Marker Genes - Louvain Clusters
bb_genebubbles(
  filter_cds(
    cds = cds_p568,
    cells = bb_cellmeta(cds_p568) |> filter(
      louvain %in% c(9, 10, 15, 19, 20, 24, 26, 31, 49, 51, 53))
   |> filter(leukemia_phenotype %in% c("AML", "No leukemia"))
  ),
genes = c("Fcgr4", #CD16
          "Fcgr3", #CD16
          "Fcgr2b", #CD32
          "Ly6a", #(Sca-1)
          "Kit",
          "Cd34"),
cell_grouping = c("louvain", "leukemia_phenotype"),
return_value = "data") |>
  dplyr::mutate(
    gene_short_name = dplyr::case_when(
      gene_short_name == "Fcgr4" ~ "Fcgr4 (Cd16)",
      gene_short_name == "Fcgr3" ~ "Fcgr3 (Cd16)",
      gene_short_name == "Fcgr2b" ~ "Fcgr2b (Cd32)",
      gene_short_name == "Ly6a" ~ "Ly6a (Sca-1)",
      TRUE ~ gene_short_name
    )
  ) |>
  ggplot(mapping = aes(
    x = louvain,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~ leukemia_phenotype, scales = "free_x") +
  theme_minimal_grid(font_size = the_font_size) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(
    title = "Lineage negative Louvain Clusters - GMP Markers",
    x = NULL,
    y = NULL,
    size = "Proportion",
    color = "Expression"
  )
################################################################################################################################################

#Highlight GMP cells found in clusters 9 & 19
head(colData(cds_p568))

GMP <- filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |> filter(
    louvain %in% c(9, 19)) |> filter(leukemia_phenotype %in% c("AML", "No leukemia")))

unique(colData(GMP)$louvain)

gmp_clust_bcode <- colData(GMP)$barcode

#Add column to cds_main to identify cells GMP cells identified in clust 9 & 19 from cds_p568
colData(cds_main)$GMP <- colData(cds_main)$barcode %in% gmp_clust_bcode

#Plot locations of GMP cells (louvain clust 9 & 19) from cds_p568 on cds_main dimensions
bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |>
               filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
  "GMP",
  value_to_highlight = TRUE,
  legend_pos = "none",
  plot_title = "GMP cells identified in louvain clusters 9 & 19 from Supp 4A",
  palette = "#EF8A62"
) +
  facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(color = "black"),
                          legend.position = "none")

#Louvain cluster 19 from cds_p568 on cds_main dimensions
clust19 <- filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |> filter(
    louvain %in% c(19)))

unique(colData(GMP)$louvain)

clust19_bcode <- colData(clust19)$barcode

#Add column to cds_main to identify cells GMP cells identified in clust 9 & 19 from cds_p568
colData(cds_main)$clust19 <- colData(cds_main)$barcode %in% clust19_bcode

bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |>
               filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
  "clust19",
  value_to_highlight = TRUE,
  legend_pos = "none",
  plot_title = "Louvain cluster 19 from Supp 4A",
  palette = "#EF8A62"
) +
  facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(color = "black"),
                          legend.position = "none")

######################################################################################################################################
#TODO perform pseudobulk analysis on GMP clusters between AML dKO and p53 KO/WT

#Pseudobulk Volcano

# unique(colData(GMP)$tissue)
#
# exp_design <-
#   bb_cellmeta(GMP) |>
#   group_by(leukemia_phenotype) |>
#   summarise()
#
# pseudobulk_res <-
#   bb_pseudobulk_mf(cds = GMP,
#                    pseudosample_table = exp_design,
#                    design_formula = "~ leukemia_phenotype",
#                    result_recipe = c("leukemia_phennotype", "AML", "No leukemia"))
#
# # Differential expression results.
# genes_to_highlight <- unique(c("Jtga4", "Klf4"))
# genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(pseudobulk_res$Result, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]
#
# volcano_data_RTvCLL <- pseudobulk_res$Result %>%
#   mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
#   mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))
#
# library(ggtext)
# F1H <-
#   ggplot(
#     volcano_data_RTvCLL,
#     aes(
#       x = log2FoldChange,
#       y = -log10(padj),
#       colour = threshold,
#       fill = threshold,
#       label = text_label
#     )
#   ) +
#   geom_point(shape = 21,
#              size = 0.5,
#              alpha = 0.4) +
#   geom_text_repel(color = "black",
#                   fontface = "italic",
#                   box.padding = 0.22, #0.5
#                   point.padding = 0.1, #0.25
#                   min.segment.length = 0,
#                   max.overlaps = 20000,
#                   size = 3,
#                   segment.size = 0.25,
#                   force = 2,
#                   seed = 1234,
#                   segment.curvature = -0.1,
#                   segment.square = TRUE,
#                   segment.inflect = TRUE) +
#   xlab("log<sub>2</sub> fold change") +
#   ylab("-log<sub>10</sub> adjusted p-value") +
#   theme(axis.title.x =  element_markdown()) +
#   theme(axis.title.y = element_markdown()) +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c("grey80", "#DC0000")) +
#   scale_fill_manual(values = c("transparent", "#DC0000")) +
#   labs(caption = "\U21D0 Up in CLL\nUp in RT \U21D2",title = "Pseudobulk:RT clusters 3 & 11 vs CLL clusters 1 & 9")+
#   theme(plot.caption.position = "panel") +
#   theme(plot.caption = element_text(hjust = 0.5)) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_cartesian(xlim = c(-1.0 * max(abs(
#     range(
#       volcano_data_RTvCLL %>% dplyr::filter(!is.na(padj)) %>% pull(log2FoldChange)
#     )
#   )), 1.0 * max(abs(
#     range(
#       volcano_data_RTvCLL %>% filter(!is.na(padj)) %>% pull(log2FoldChange)
#     )
#   ))))
