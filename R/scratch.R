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

pander(analysis_configs)

bb_gene_violinplot(cds=cds_main, variable="partition", genes_to_plot = c("Cd34","Cd47", "Vsir"), rows=3, palette = c("1" = "#F0E442", "2" = "#0072B2", "3" = "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "grey80")) + scale_fill_viridis_d()
