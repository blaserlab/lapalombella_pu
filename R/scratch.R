# load the orthology data http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt; 20210930
homology <- read_tsv("/workspace/workspace_pipelines/lapalombella.pu.datapkg/inst/extdata/HOM_AllOrganism.rpt.txt")

human_mouse_homology <-
  homology %>%
  select(key = `DB Class Key`, organism = `Common Organism Name`, mouse_symbol = Symbol) %>%
  filter(organism == "mouse, laboratory") %>%
  select(-organism) %>%
  left_join(homology %>% select(key = `DB Class Key`, organism = `Common Organism Name`, human_symbol = Symbol) %>% filter(organism == "human") %>% select(-organism)) %>%
  select(-key)

# load the gene sets

human_gene_sets <- fgsea::gmtPathways(gmt.file = "/workspace/workspace_pipelines/lapalombella.pu.datapkg/inst/extdata/msigdb.v7.4.symbols.gmt")
human_gene_sets
mouse_gene_sets <-
  map(.x = human_gene_sets,
    .f = function(x, homology = human_mouse_homology) {
      tbl <- tibble(human_symbol = x)
      result <- left_join(tbl, homology) %>% filter(!is.na(mouse_symbol))
      return(result %>% pull(mouse_symbol))
    })


blast_like_p5p8_aggscore <- aggregate_gene_expression(cds = cds_wt_aml_marrow[,colData(cds_wt_aml_marrow)$specimen %in% c("P5", "P8") & colData(cds_wt_aml_marrow)$leiden_assignment == "Blast-like"],
                          cell_group_df = data.frame(cell = rownames(colData(cds_wt_aml_marrow[,colData(cds_wt_aml_marrow)$specimen %in% c("P5", "P8") & colData(cds_wt_aml_marrow)$leiden_assignment == "Blast-like"])),
                                                     cell_grouping = colData(cds_wt_aml_marrow[,colData(cds_wt_aml_marrow)$specimen %in% c("P5", "P8") & colData(cds_wt_aml_marrow)$leiden_assignment == "Blast-like"])$specimen)) %>%
  as.matrix() %>%
  as_tibble(rownames = "id") %>%
  mutate(P5_minus_P8 = P5-P8) %>% left_join(rowData(cds_wt_aml_marrow) %>% as_tibble())
blast_like_p5p8_aggscore

gsea_stats <- blast_like_p5p8_aggscore$P5
names(gsea_stats) <- blast_like_p5p8_aggscore$gene_short_name

blast_like_gsea_res <- fgsea::fgsea(stats = gsea_stats, pathways = mouse_gene_sets)

blast_like_gsea_res %>%
  as_tibble() %>%
  arrange(padj) %>%
  View()
fgsea::plotEnrichment(pathway = mouse_gene_sets[["GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"]], stats = gsea_stats)
