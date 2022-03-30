cds_wt_aml_spleen <- cds_main[,colData(cds_main)$specimen %in% c("P4", "P2")]
colData(cds_wt_aml_spleen)$main_partition <- colData(cds_wt_aml_spleen)$partition
colData(cds_wt_aml_spleen)$main_leiden <- colData(cds_wt_aml_spleen)$leiden
colData(cds_wt_aml_spleen)$main_louvain <- colData(cds_wt_aml_spleen)$louvain

colData(cds_wt_aml_spleen)$partition <- NULL
colData(cds_wt_aml_spleen)$leiden <- NULL
colData(cds_wt_aml_spleen)$louvain <- NULL

cds_wt_aml_spleen <- preprocess_cds(cds_wt_aml_spleen)
cds_wt_aml_spleen <- bb_align(cds_wt_aml_spleen, align_by = "specimen")
cds_wt_aml_spleen <- reduce_dimension(cds_wt_aml_spleen, cores = 39)
cds_wt_aml_spleen <- bb_triplecluster(cds_wt_aml_spleen, outfile = "data/tm_wt_aml_spleen.csv", n_cores = 24)
tm_wt_aml_spleen <- read_csv("data/tm_wt_aml_spleen.csv")
fastSave::save.pigz(tm_wt_aml_spleen, file = "data/tm_wt_aml_spleen.rda", n.cores = 8)


colData(cds_wt_aml_spleen)$leiden_assignment <-
  recode(colData(cds_wt_aml_spleen)$leiden,
         "1" = "pre-Neu2/3",
         "2" = "Neu",
         "3" = "immNeu",
         "4" = "pre-Neu1",
         "5" = "HSC/Prog",
         "6" = "Blast-like",
         "7" = "Mono",
         "8" = "cMoP/DC",
         "9" = "B",
         "10" = "T/NK 1",
         "11" = "PC",
         "12" = "T/NK 2"
         )
colData(cds_wt_aml_spleen)$blast_like_other <-
  recode(colData(cds_wt_aml_spleen)$leiden,
         "1" = "other",
         "2" = "other",
         "3" = "other",
         "4" = "other",
         "5" = "other",
         "6" = "Blast-like",
         "7" = "other",
         "8" = "other",
         "9" = "other",
         "10" = "other",
         "11" = "other",
         "12" = "other"
         )

tm_blast_like_other_AML_spleen <- top_markers(cds_wt_aml_spleen, group_cells_by = "blast_like_other", genes_to_test_per_group = 200, cores = 39)
tm_blast_like_hsc_prog_AML_spleen <- top_markers(cds_wt_aml_spleen[,colData(cds_wt_aml_spleen)$leiden_assignment %in% c("Blast-like", "HSC/Prog")], group_cells_by = "leiden_assignment", genes_to_test_per_group = 200, cores = 39)
tm_blast_like_neu_AML_spleen <- top_markers(cds_wt_aml_spleen[,colData(cds_wt_aml_spleen)$leiden_assignment %in% c("Blast-like", "Neu")], group_cells_by = "leiden_assignment", genes_to_test_per_group = 200, cores = 39)

fastSave::save.pigz(cds_wt_aml_spleen, file = "data/cds_wt_aml_spleen.rda", n.cores = 8)
fastSave::save.pigz(tm_blast_like_other_AML_spleen, file = "data/tm_blast_like_other_AML_spleen.rda", n.cores = 8)
fastSave::save.pigz(tm_blast_like_hsc_AML_spleen, file = "data/tm_blast_like_hsc_prg_AML_spleen.rda", n.cores = 8)
fastSave::save.pigz(tm_blast_like_neu_AML_spleen, file = "data/tm_blast_like_neu_AML_spleen.rda", n.cores = 8)

write_csv(tm_blast_like_other_AML_spleen, file = "data/tm_blast_like_other_AML_spleen.csv")
write_csv(tm_blast_like_hsc_prog_AML_spleen, file = "data/tm_blast_like_hsc_prg_AML_spleen.csv")
write_csv(tm_blast_like_neu_AML_spleen, file = "data/tm_blast_like_neu_AML_spleen.csv")


muench_gene_groups <- left_join(as_tibble(muench_cluster_tm), as_tibble(rowData(cds_wt_aml_spleen))) %>%
  select(id, gene_grouping = cell_group) %>%
  filter(!is.na(id))


agg_mat_wt_aml_spleen  <- aggregate_gene_expression(cds = cds_wt_aml_spleen,
                          gene_group_df = muench_gene_groups,
                          cell_group_df = data.frame(cell = rownames(colData(cds_wt_aml_spleen)),
                                                     cell_grouping = colData(cds_wt_aml_spleen)$leiden_assignment),scale_agg_values = T, norm_method = "log")

fastSave::save.pigz(muench_cds, file = "data/muench_cds.rda", n.cores = 8)
fastSave::save.pigz(agg_mat_wt_aml_spleen, file = "data/agg_mat_wt_aml_spleen.rda", n.cores = 8)

cds_wt_aml_spleen_p4 <- cds_main[,colData(cds_main)$specimen %in% c("P4")]
