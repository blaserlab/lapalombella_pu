#Figure 7
# cell composition stacked bar chart----------------------
pheno_sums <- bb_cellmeta(cds_main) |>
  count(leukemia_phenotype, name = "sum")

F7_microenv_bp<- bb_cellmeta(cds_main) |>
  count(leiden_assignment2, leukemia_phenotype) |>
  left_join(pheno_sums) |>
  mutate(total = sum(pheno_sums$sum)) |>
  mutate(ratio = sum/total) |>
  mutate(normal = n/ratio/4) |>
  filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
  filter(str_detect(leiden_assignment2, "ALL", negate = TRUE)) |>
  ggplot(aes(x = leiden_assignment2, y = normal, fill = leukemia_phenotype)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "ScType Assignment")

microenv_populations <- bb_cellmeta(cds_main) |>
  group_by(leiden_assignment2) |>
  summarise() |>
  filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
  filter(str_detect(leiden_assignment2, "ALL", negate = TRUE)) |>
  pull(leiden_assignment2)

cds_main$leukemia_phenotype <- factor(cds_main$leukemia_phenotype,
                                  levels = c("No leukemia","AML","pre-B ALL"))
F7_microenv_map <- bb_var_umap(
  filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |>
      filter(leukemia_phenotype %in% c("No leukemia", "AML", "pre-B ALL"))
  ),
  var = "leiden_assignment2", cell_size = 0.1,
  #value_to_highlight = microenv_populations,
  facet_by = "leiden", overwrite_labels = F
) #+theme(legend.text = element_text(size = 8))

marrow_gb_plot <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(tissue %in% c("marrow")) |>
        #filter(leukemia_phenotype %in% c("AML", "pre-B ALL")) |>
        filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
        filter(str_detect(leiden_assignment2, "ALL", negate = TRUE))
    ),
    genes = c(
      "Cd37",
      "Cd79a",
      "Cd3e",
      "Cd4",
      "Gzma",
      "Icos",
      "Izumo1r",
      "Trbc1"
    ),
    #"Itgam", "Ly6g", "Gzma"),
    cell_grouping = c("leiden_assignment2", "leukemia_phenotype"),
    return_value = "data"
  ) |>
  ggplot(mapping = aes(
    x = leiden_assignment2,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~leukemia_phenotype, scales = "free_x",) +
  theme_minimal_grid(font_size = 6) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL,
       y = NULL,
       size = "Proportion",
       color = "Expression")
marrow_gb_plot

gb_plot2 <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(leukemia_phenotype %in% c("AML", "pre-B ALL", "No leukemia")) |>
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
    cell_grouping = c("leiden_assignment2", "leukemia_phenotype"),
    return_value = "data"
  ) |>
  mutate(across(leukemia_phenotype, factor, levels=c("No leukemia","AML","pre-B ALL")))|>
  ggplot(mapping = aes(
    x = leiden_assignment2,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~leukemia_phenotype, scales = "free_x",) +
  theme_minimal_grid(font_size = 8) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL,
       y = NULL,
       size = "Proportion",
       color = "Expression")

gb_plot2

#Figure 7 scratch:
# aml_T_map <- bb_gene_umap(
#   filter_cds(
#     cds_main,
#     cells = bb_cellmeta(cds_main) |>
#       filter(leukemia_phenotype %in% c("AML")) |>
#       filter(str_detect(leiden_assignment2, "T")) |>
#       filter(str_detect(leiden_assignment2, "ALL", negate = TRUE))
#   ),
#   gene_or_genes = c("Cd3e",
#                     "Tigit",
#                     "Pdcd1", #Pd1
#                     "Havcr2" #Tim3)
#                     )) +
#     labs(title = "dKO AML: T cells") +
#     panel_border()
#
# ball_T_map <- bb_gene_umap(
#       filter_cds(
#         cds_main,
#         cells = bb_cellmeta(cds_main) |>
#           filter(leukemia_phenotype %in% c("pre-B ALL")) |>
#           filter(str_detect(leiden_assignment2, "T")) |>
#           filter(str_detect(leiden_assignment2, "ALL", negate = TRUE))
#       ),
#       gene_or_genes = c("Cd3e",
#                         "Tigit",
#                         "Pdcd1", #Pd1
#                         "Havcr2" #Tim3)
#       )) +
#         labs(title = "dKO pre-B ALL: T cells")+
#        panel_border()
#
# aml_ballT_map <- aml_T_map | ball_T_map
# aml_ballT_map
#Figure 7B
# cellcount_aml <- bb_cellmeta(blaseRtools::filter_cds(
#   cds = cds_main,
#   cells = bb_cellmeta(cds_main) |>
#     filter(leukemia_phenotype %in% c("AML"))
# )) |> #all clusters
#   group_by(leiden, leiden_assignment2) |>
#   summarise(n = n())
# cellcount_aml <- filter(cellcount_aml, n > 10)
#
# aml <-
#   blaseRtools::filter_cds(
#     cds = blaseRtools::filter_cds(
#       cds = cds_main,
#       cells = bb_cellmeta(cds_main) |>
#         filter(leukemia_phenotype %in% c("AML"))
#     ),
#     cells = bb_cellmeta(
#       blaseRtools::filter_cds(
#         cds = cds_main,
#         cells = bb_cellmeta(cds_main) |>
#           filter(leukemia_phenotype %in% c("AML"))
#       )
#     ) |>
#       filter(leiden %in% dput(as.character(
#         cellcount_aml$leiden
#       )))
#   )
#
# cellcount_ball <- bb_cellmeta(blaseRtools::filter_cds(
#   cds = cds_main,
#   cells = bb_cellmeta(cds_main) |>
#     filter(leukemia_phenotype %in% c("pre-B ALL"))
# )) |>
#   group_by(leiden, leiden_assignment2) |>
#   summarise(n = n())
# cellcount_ball <- filter(cellcount_ball, n > 10)
# allb <-
#   blaseRtools::filter_cds(
#     cds = blaseRtools::filter_cds(
#       cds = cds_main,
#       cells = bb_cellmeta(cds_main) |>
#         filter(leukemia_phenotype %in% c("pre-B ALL"))
#     ),
#     cells = bb_cellmeta(
#       blaseRtools::filter_cds(
#         cds = cds_main,
#         cells = bb_cellmeta(cds_main) |>
#           filter(leukemia_phenotype %in% c("pre-B ALL"))
#       )
#     ) |>
#       filter(leiden %in% dput(as.character(
#         cellcount_ball$leiden
#       )))
#   )
#
# a <-
#   bb_var_umap(aml, "leiden_assignment2", overwrite_labels = T) +
#   labs(title = "AML") +
#   ylim(-13.5, 11) +
#   xlim(-12, 15)
# b <-
#   bb_var_umap(allb, "leiden_assignment2", overwrite_labels = T) +
#   labs(title = "pre-B ALL") +
#   ylim(-13.5, 11) +
#   xlim(-12, 15)
# F7B1 <- a + b
#
# F7B1
# cellcount_aml_ball <- bb_cellmeta(blaseRtools::filter_cds(
#   cds = cds_main,
#   cells = bb_cellmeta(cds_main) |>
#     filter(leukemia_phenotype %in% c("AML", "pre-B ALL"))
# )) |> #all clusters
#   group_by(leiden, leiden_assignment2) |>
#   summarise(n = n())
# cellcount_aml_ball <- filter(cellcount_aml_ball, n > 10)
#
# aml_allb <-
#   blaseRtools::filter_cds(
#     cds = blaseRtools::filter_cds(
#       cds = cds_main,
#       cells = bb_cellmeta(cds_main) |>
#         filter(leukemia_phenotype %in% c("AML", "pre-B ALL"))
#     ),
#     cells = bb_cellmeta(
#       blaseRtools::filter_cds(
#         cds = cds_main,
#         cells = bb_cellmeta(cds_main) |>
#           filter(leukemia_phenotype %in% c("AML","pre-B ALL"))
#       )
#     ) |>
#       filter(leiden %in% dput(as.character(
#         cellcount_aml_ball$leiden
#       )))
#   )
#
# F7B1 <- bb_var_umap(aml_allb, "leiden_assignment2", overwrite_labels = T) +
#   labs(title = "AML and pre-B ALL")
# F7B1

# F7B2 <- bb_var_umap(aml_allb,"density", facet_by = "leukemia_phenotype")
# F7B2

####Pu 7C-7F:
#
# #Annotation
# #cds_p53tet2AMLann<-bb_cds_anno(query_cds=cds_p53tet2AML, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )
# #cds_mainann<-bb_cds_anno(query_cds=cds_main, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )
# #cds_preBALLann<-bb_cds_anno(query_cds=cds_preBALL, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )
#
#
# #subset T/NK clusters in AML and study exhaustive markers
# cds_aml_T <- cds_p53tet2AMLann[,colData(cds_p53tet2AMLann)$predicted.leiden_assignment %in% c("T/NK 1", "T/NK 2")]
#
# F7C <- bb_var_umap(cds_aml_T, var="louvain")
#
# #subset T/NK clusters in Pre-B ALL and study exhaustive markers
# # cds_BALL_T <- blaseRtools::filter_cds(cds = cds_main,
# #                                 cells = bb_cellmeta(cds_main) |>
# #                                   filter(leukemia_phenotype %in% c("pre-B ALL")) |> filter(predicted.leiden_assignment %in% c("T/NK 1", "T/NK 2")))
# cds_preBALL<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "pre-B ALL"]
# cds_BALL_T <- cds_preBALLann[,colData(cds_preBALLann)$predicted.leiden_assignment %in% c("T/NK 1", "T/NK 2")]
#
# F7F <- bb_var_umap(cds_BALL_T, var="louvain")
#
# bb_gene_umap(cds=cds_BALL_T, gene_or_gene = c("Cd3e","Cd8a","Cd4", "Cd34", "Cd19", "Pax5"), cell_size = 0.1)
#
# bb_gene_violinplot(cds=cds_BALL_T, variable="leiden", genes_to_plot = c("Pdcd1", "Lag3", "Tox", "Ctla4", "Tigit", "Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")
# bb_gene_violinplot(cds=cds_BALL_T, variable="louvain", genes_to_plot = c("Tox","Eomes", "Klrg1", "Pdcd1","Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")
