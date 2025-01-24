Supp_Fig_5A1 <- bb_var_umap(filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |>
    filter(louvain %in% c(53, 48, 4, 21, 33, 38, 10, 49, 24, 9, 34, 19, 26, 20, 51, 15, 31)
    ) |> filter(leukemia_phenotype %in% c("AML", "No leukemia"))
), "louvain", overwrite_labels = T)

Supp_Fig_5A1

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_5A1.pdf"),
  plot = Supp_Fig_5A1,
  base_width = 5,
  base_height = 4
)

Supp_Fig_5B <- bb_var_umap(filter_cds(
  cds = cds_p568,
  cells = bb_cellmeta(cds_p568) |>
    filter(louvain %in% c(53, 48, 4, 21, 33, 38, 10, 49, 24, 9, 34, 19, 26, 20, 51, 15, 31)
    ) |> filter(leukemia_phenotype %in% c("AML", "No leukemia"))
),"density", facet_by = "leukemia_phenotype")

Supp_Fig_5B

save_plot(
  filename = "temp.pdf",
  #filename = fs::path(figs_out, "Supp_Fig_5B.pdf"),
  plot = Supp_Fig_5B,
  base_width = 10,
  base_height = 4
)

#Supp 5 GMP Cluster Identification
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
               filter(leukemia_phenotype %in% c("AML"))),
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

Supp_Fig_5A2 <- bb_var_umap(
  filter_cds(cds_main, cells = bb_cellmeta(cds_main) |>
               filter(leukemia_phenotype %in% c("AML"))),
  "clust19",
  value_to_highlight = TRUE,
  legend_pos = "none",
  plot_title = "Louvain cluster 19",
  palette = "#EF8A62",
  rasterize = T
) +
  facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(color = "black"),
                          legend.position = "none")
Supp_Fig_5A2

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_5A2.pdf"),
  plot = Supp_Fig_5A2,
  base_width = 4,
  base_height = 4
)

#Supplemental Figure 5C
colData(cds_p568)$louvain <- factor(
  colData(cds_p568)$louvain,
  levels = c("10", "38", "48", "53", "21","33", "24", "4", "49", "19", "34", "20","26", "9","51","15","31")
)
Supp_Fig_5C <- bb_genebubbles(
  filter_cds(cds_p568, cells = bb_cellmeta(cds_p568) |>
               filter(louvain %in% c("10", "38", "48", "53", "21","33", "24", "4", "49", "19", "34", "20","26", "9","51","15","31"))),
               #filter(leukemia_phenotype %in% c("AML"))),
  genes = c("Klf4", "Flt3", "Fcgr4", "Itga4"),
  cell_grouping = "louvain",
  experiment_type = "Gene Expression",
  scale_expr = F) + geom_point(pch = 21) +
  scale_size_area() +
  scale_fill_viridis_c(option = "A") +
  theme_minimal_grid()

Supp_Fig_5C

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_5C.pdf"),
  plot = Supp_Fig_5C,
  base_width = 6.5,
  base_height = 4.0
)

