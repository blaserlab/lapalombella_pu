
bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "pre-B ALL")
    )), "leiden_assignment2", overwrite_labels = F) + facet_wrap(~geno_pheno)
#ggsave("F7B.pdf", path = T_figs)#, width = , height = )

aml <- blaseRtools::filter_cds(cds = cds_main,
                                     cells = bb_cellmeta(cds_main) |>
                                       filter(leukemia_phenotype %in% c("AML")))
cellcount_aml<- bb_cellmeta(aml) |> #all clusters
  group_by(leiden, leiden_assignment2) |>
  summarise(n = n())
cellcount_aml <- filter(cellcount_aml, n > 10)
aml <- blaseRtools::filter_cds(cds = aml,
                               cells = bb_cellmeta(aml) |>
                                 filter(leiden %in% dput(as.character(cellcount_aml$leiden))
                                          ))
allb <- blaseRtools::filter_cds(cds = cds_main,
                                cells = bb_cellmeta(cds_main) |>
                                  filter(leukemia_phenotype %in% c("pre-B ALL")))
cellcount_ball<- bb_cellmeta(allb) |>
  group_by(leiden, leiden_assignment2) |>
  summarise(n = n())
cellcount_ball <- filter(cellcount_ball, n > 10)
allb <- blaseRtools::filter_cds(cds = allb,
                                cells = bb_cellmeta(allb) |>
                                  filter(leiden %in% dput(as.character(cellcount_ball$leiden))
                                  ))
aml_allb <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML", "pre-B ALL")))
cellcount_aml_ball<- bb_cellmeta(aml_allb) |>
  group_by(leiden, leiden_assignment2, leukemia_phenotype) |>
  summarise(n = n())
cellcount_aml_ball <- filter(cellcount_aml_ball, n > 10)

aml_allb <- blaseRtools::filter_cds(cds = aml_allb,
                               cells = bb_cellmeta(aml_allb) |>
                                 filter(leiden %in% dput(as.character(cellcount_aml_ball$leiden))
                                 ))
# aml_allb <- combine_cds(list(aml, allb))
# colData(aml_allb)

a<- bb_var_umap(aml, "leiden_assignment2", overwrite_labels = T) + labs(title = "AML") + ylim(-13.5, 11) +xlim(-12, 15)
b<- bb_var_umap(allb, "leiden_assignment2", overwrite_labels = T) + labs(title = "B ALL") + ylim(-13.5, 11) +xlim(-12, 15)
F7B<- a+b
# ggsave("F7B.pdf")
# p <- p1 + p2
# ggsave('figures/test.png', p)
