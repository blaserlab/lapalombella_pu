#Supplemental Figure 8

#Supp Fig 8A - microenv umap
microenv_populations <- bb_cellmeta(cds_main) |>
  group_by(leiden_assignment2) |>
  summarise() |>
  filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
  filter(str_detect(leiden_assignment2, "ALL", negate = TRUE)) |>
  pull(leiden_assignment2)

# cds_main$leukemia_phenotype <- factor(cds_main$leukemia_phenotype,
#                                   levels = c("No leukemia","AML","pre-B ALL"))
Supp_Fig_8A <- bb_var_umap(
  filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |>
      filter(geno_pheno %in% c("dKO: AML", "Wildtype", "P53 KO: T ALL"))
  ),
  var = "leiden_assignment2", cell_size = 0.5,
  value_to_highlight = microenv_populations,
  facet_by = "genotype", overwrite_labels = F, rasterize = T
) #+theme(legend.text = element_text(size = 8))
Supp_Fig_8A

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_8A.pdf"),
  plot = Supp_Fig_8A,
  base_width = 12,
  base_height = 3
)
#cell composition stacked bar chart
# pheno_sums <- bb_cellmeta(cds_main) |>
#   count(leukemia_phenotype, name = "sum")
#
# SF8_microenv_bp<- bb_cellmeta(cds_main) |>
#   count(leiden_assignment2, leukemia_phenotype) |>
#   left_join(pheno_sums) |>
#   mutate(total = sum(pheno_sums$sum)) |>
#   mutate(ratio = sum/total) |>
#   mutate(normal = n/ratio/4) |>
#   filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
#   filter(str_detect(leiden_assignment2, "ALL", negate = TRUE)) |>
#   ggplot(aes(x = leiden_assignment2, y = normal, fill = leukemia_phenotype)) +
#   geom_bar(position="fill", stat="identity") +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
#   labs(x = "ScType Assignment")

#Supp Fig 8B
Supp_Fig_8B <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(geno_pheno %in% c("dKO: AML", "Wildtype", "P53 KO: T ALL")) |>
        filter(str_detect(leiden_assignment2, "T")) |>
        # filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
        filter(str_detect(leiden_assignment2, "ALL", negate = TRUE))
    ),
    genes = c(
      "Tox",
      "Eomes",
      "Havcr2",
      "Pdcd1",
      "Cd3e",
      "Tigit",
      "Sell",
      "Slamf6"),
    cell_grouping = c("leiden_assignment2", "genotype"),
    return_value = "data"
  ) |>
  mutate(across(genotype, factor, levels=c("TP53-/-/TET2-/-", "TP53-/-", "WT")))|>
  ggplot(mapping = aes(
    x = leiden_assignment2,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~genotype, scales = "free_x",) +
  theme_minimal_grid(font_size = 8) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic", size = 9)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9)) +
  labs(x = NULL,
       y = NULL,
       size = "Proportion",
       color = "Expression") +
  theme(plot.title = element_text(size = 140))
Supp_Fig_8B

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_8B.pdf"),
  plot = Supp_Fig_8B,
  base_width = 8,
  base_height = 4
)

#Supp Fig 8C
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
#Supp Fig 8C
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
#Supp Fig 8D
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
