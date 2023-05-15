#Figure 7B
cellcount_aml<- bb_cellmeta(blaseRtools::filter_cds(cds = cds_main,
                                                    cells = bb_cellmeta(cds_main) |>
                                                      filter(leukemia_phenotype %in% c("AML")))) |> #all clusters
  group_by(leiden, leiden_assignment2) |>
  summarise(n = n())
cellcount_aml <- filter(cellcount_aml, n > 10)
aml <- blaseRtools::filter_cds(cds = blaseRtools::filter_cds(cds = cds_main,
                                                             cells = bb_cellmeta(cds_main) |>
                                                               filter(leukemia_phenotype %in% c("AML"))),
                               cells = bb_cellmeta(blaseRtools::filter_cds(cds = cds_main,
                                                                           cells = bb_cellmeta(cds_main) |>
                                                                             filter(leukemia_phenotype %in% c("AML")))) |>
                                 filter(leiden %in% dput(as.character(cellcount_aml$leiden))
                                          ))

cellcount_ball<- bb_cellmeta(blaseRtools::filter_cds(cds = cds_main,
                                                     cells = bb_cellmeta(cds_main) |>
                                                       filter(leukemia_phenotype %in% c("pre-B ALL")))) |>
  group_by(leiden, leiden_assignment2) |>
  summarise(n = n())
cellcount_ball <- filter(cellcount_ball, n > 10)
allb <-
  blaseRtools::filter_cds(
    cds = blaseRtools::filter_cds(
      cds = cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(leukemia_phenotype %in% c("pre-B ALL"))
    ),
    cells = bb_cellmeta(
      blaseRtools::filter_cds(
        cds = cds_main,
        cells = bb_cellmeta(cds_main) |>
          filter(leukemia_phenotype %in% c("pre-B ALL"))
      )
    ) |>
      filter(leiden %in% dput(as.character(
        cellcount_ball$leiden
      )))
  )

a <-
  bb_var_umap(aml, "leiden_assignment2", overwrite_labels = T) +
  labs(title = "AML") +
  ylim(-13.5, 11) +
  xlim(-12, 15)
b <-
  bb_var_umap(allb, "leiden_assignment2", overwrite_labels = T) +
  labs(title = "pre-B ALL") +
  ylim(-13.5, 11) +
  xlim(-12, 15)
F7B<- a+b

F7B
#ggsave("F7B.pdf", path = figs_out)

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
