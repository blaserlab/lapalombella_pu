#Supplemental Figure 3
########################################################################################################################
aml <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12')))
########################################################################################################################

#Murine Lineage/HSC/MPP marker gene bubble plots
#Lineage Markers gene bubble plot
Supp_Fig_3A<- bb_genebubbles(
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
  labs(title = "Lineage Markers", x = NULL, y = NULL, size = "Proportion", color = "Expression")

Supp_Fig_3A

Supp_Fig_3B <- bb_genebubbles(
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
  labs(title = "HSC/MPP Markers", x = NULL, y = NULL, size = "Proportion", color = "Expression")

Supp_Fig_3B

#Human CITEseq/Mouse scRNAseq combined data
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

colData(cds_combined)$aml_leiden_enrichment <- NULL

cds_combined <- left_join(
  bb_cellmeta(cds_combined),
  leiden.1_tibble,
  by = join_by(leiden.1)
) |>
  mutate(aml_leiden_enrichment = replace_na(aml_leiden_enrichment, "other")) |> select(cell_id, aml_leiden_enrichment) |>
  bb_tbl_to_coldata(obj = cds_combined, min_tbl = _)

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

colData(cds_combined)$aml_leiden_enrichment2 <- NULL

cds_combined <- left_join(
  bb_cellmeta(cds_combined),
  mapped_enrich,
  by = join_by(cell_id)
) |>
  select(cell_id, aml_leiden_enrichment2) |>
  bb_tbl_to_coldata(obj = cds_combined, min_tbl = _)

colData(cds_combined)$aml_leiden_enrichment2 <- factor(
  colData(cds_combined)$aml_leiden_enrichment2,
  levels = c("B cells", "TNK", "WT_aml", "tet2_aml", "tp53_aml", "comutant_aml", "Murine dKO TET2/TP53 AML", "other")
)

colData(cds_combined)$aml_leiden_enrichment <- factor(
  colData(cds_combined)$aml_leiden_enrichment,
  levels = c("B cells", "TNK", "WT_aml", "tet2_aml", "tp53_aml", "comutant_aml")
)

# filtered out the unmapped "other" cells---------------------------------------
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

#Supplemental Figure 3C & 3D
Supp_Fig_3CD <- bb_var_umap(obj = filter_cds(
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

Supp_Fig_3CD

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_3CD.pdf"),
  plot = Supp_Fig_3CD,
  base_width = 10,
  base_height = 4.0
)

#Supplemental Figure 3E
Supp_Fig_3E <- bb_var_umap(filter_cds(
  cds_combined,
  cells = bb_cellmeta(cds_combined) |>
    filter(!aml_leiden_enrichment2 %in% c("other", "B cells", "TNK"))),
  "aml_leiden_enrichment2",
  palette = experimental_group_palette_2,
  rasterize = T)

Supp_Fig_3E

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_3E.pdf"),
  plot = Supp_Fig_3E,
  base_width = 6,
  base_height = 4
)

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

#Supplemental Figure 3F
Supp_Fig_3F <- bb_cellmeta(cds_combined) |>
  filter(aml_leiden_enrichment %in% c("comutant_aml", "tet2_aml", "WT_aml", "tp53_aml")) |>
  filter(data_set == "mouse") |>
  filter(partition %in% c("3", "6")) |>
  count(aml_leiden_enrichment, genotype) |>
  ggplot(aes(x = genotype, y = n, fill = aml_leiden_enrichment)) +
  geom_bar(stat = "identity", position = "fill") + labs(y = "Proportion") +
  scale_fill_manual(values = experimental_group_palette_2)

Supp_Fig_3F

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_3F.pdf"),
  plot = Supp_Fig_3F,
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
