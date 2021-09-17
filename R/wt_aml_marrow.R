#' ---
#' title: "Overview and Question 1"
#' author: "Brad Blaser"
#' date: "9/17/2021"
#' output: pdf_document
#' ---
#'
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
source("R/dependencies.R")
source("R/configs.R")

#' ## Introduction and Overview
#'
#' Pu asked a number of questions:
#'
#' 1)	TP53/TET2 AML vs WT: How does double knockout alter the states of cells in terms of immune subpopulation distribution (T cell, B cell, NK cell, monocyte and DC) and transcriptome changes in progenitor cells (GMP, MEP, CMP, MPP and HSC)? Can we identify unique T/B ALL progenitors within AML phenotype?
#' 2)	Engrafted TP53/TET2 AML vs WT cell engrafted PepboyJ: How does leukemia remodel immune environment? Can we identify unique T/B cell progenitors within AML phenotype
#' 3)	TP53/TET2 AML vs TP53 T ALL (genotype effect): Although they develop leukemia of different lineages, are there any common shared progenitor and immune populations in AML and T ALL?
#' 4)	TP53/TET2 AML vs TP53/TET2 pre-B ALL (phenotype effect) : Although they develop leukemia of different lineages, are there any  common shared progenitor and immune populations in AML and pre-B ALL?
#' 5)	Primary transgenic TP53/TET2 AML vs engrafted TP53/TET2 AML vs WT cell engrafted PepboyJ (stage effect): Evolution of leukemia pathway and transcriptional programs in progenitors . How do pre-leukemia and leukemia cells interact with immune environment?
#' 6)	Engrafted TP53/TET2 AML spleen vs bone marrow (Organ effect): relative cell type abundance and cell cycle properties in different tissues

#' These are a lot of questions so we will take them in order and determine progressively if we are able to generate adequate answers.
#'
#' The challenge with the data you have aquired is that every sample has a different combination of metadata features.  Meaning:  different tissue, genotype, transplant status, leukemia phenotype, etc:
#'
#+ echo=TRUE
kable(analysis_configs[,-c(2,7)])
write_csv(analysis_configs[,-c(2,7)], file = str_glue("{tables_out}/analysis_configs.csv"))

#' There are no useable biological replicates that I can identify.  So this is more of a survey and we would have to consider anything done here to be confirmatory of other findings or hypothesis generating.
#'
#' I would say that the samples are usable from a qc standpoint.  Several still fell short even after all of the extra sequencing.  We aim for 15,000 reads per cell minimum:
#'

#+ echo=TRUE
kable(sequencing_metrics[,c(1,2,12)])
write_csv(sequencing_metrics, file = str_glue("{tables_out}/sequencing_metrics.csv"))

#' Here are all of the samples aligned into the same UMAP space:

#+ echo=FALSE
overview_umap <- bb_var_umap(cds_main, "specimen")

#+ echo=TRUE, dev="png", dpi=300, fig.height=4.5, fig.width=7.0
overview_umap

#' This figure is descriptive of the diversity of the samples you have acquired. However it doesn't help us with the individual questions asked above because when we ask questions about AML cells for example, the dimension reduction (UMAP) will be distorted by other samples which are irrelevant to the question.
#'
#' So what we will have to do in each case is isolate the samples we want to look at to answer the questions, perform dimensionality reduction and go from there.
#'
#' ## Question 1
#'
#' In order to look at the AML vs WT samples I subset the main dataset to include only samples P5 (AML marrow) and P8 (WT marrow).
#'
#' The following plot shows the cells from these samples with color encoding local cell density in UMAP space.  All cells from each sample are shown.  There are many more of the WT cells than the AML cells, however the color scale is relative to each sample so is comparable on the two plots.

#+ echo=FALSE
wt_aml_marrow_density_umap <- bb_var_umap(cds_wt_aml_marrow, "density", facet_by = "specimen")

#+ echo=TRUE, dev="png", dpi=300, fig.width=7.5, fig.height=3.5
wt_aml_marrow_density_umap

#' Clearly there are areas of overlap , but the distribution of cells is very different.
#'
#' The question is what cells are these.  To help answer we can use annotated reference data.  There are many ways to do this in theory, but in the end it depends on what is available.  I selected a WT mouse bone marrow dataset from Muench et al., Nature 2020 (Grimes lab, GEO accession:  GSE142341) because it was well-annotated, published in a high-level study, and available in a format we can use.
#'
#'  Here is the reference dataset from that article.  Of note, these data are from facs-enriched marrow cells.  Cluster identification is from the reference article.
#'

#+ echo=FALSE
ref_umap <- bb_var_umap(muench_cds, "muench_cluster1", overwrite_labels = T, group_label_size = 4)

#+ echo=TRUE, dev="png", dpi=300, fig.height=4.0, fig.width=4.5
ref_umap

#' Then we can identify the top 50 specific markers from each of the reference clusters.
#'
#+ echo=TRUE
write_csv(muench_cluster_tm, file = str_glue("{tables_out}/muench_cluster_tm.csv"))

#' Next we identify clusters in our dataset with suitable granularity, such as this:
#'

#+ echo=FALSE
wt_aml_marrow_leiden_umap <-
  bb_var_umap(
    cds_wt_aml_marrow,
    "leiden",
    overwrite_labels = T,
    group_label_size = 4)

#+ echo=TRUE, dev="png", dpi=300, fig.height=4.0, fig.width = 4.5
wt_aml_marrow_leiden_umap

#' Next we calculate the aggregate score for each of the reference gene lists across all of our cell clusters and plot as a heatmap:
#'

#+ echo=FALSE
col_fun_wt_aml_marrow <- colorRamp2(
  breaks = c(
    min(as.matrix(agg_mat_wt_aml_marrow)),
    0,
    max(as.matrix(agg_mat_wt_aml_marrow))
  ),
  colors = heatmap_3_colors
)
wt_aml_marrow_ref_heatmap <-
grid.grabExpr(draw(
  Heatmap(
    matrix = as.matrix(agg_mat_wt_aml_marrow),
    name = "Aggregate\nScore",
    col = col_fun_wt_aml_marrow,
    column_title = "Cell Cluster",
    column_title_side = "bottom",
    row_title = "Reference Gene Set"
  )
), wrap = T)

#+ echo=TRUE, dev="png", fig.width=7.5, fig.height=7.5
plot_grid(wt_aml_marrow_ref_heatmap)

#' From this we can come up with a table of cluster identities:
#'
#+ echo=FALSE
colData(cds_wt_aml_marrow) %>%
  as_tibble() %>%
  group_by(leiden, leiden_assignment) %>%
  summarise() %>%
  kable()

#' NB:  The reference data set did not include leukemia cells, mature neutrophils, B, or T cells, so we add these based on some inference from what our samples are and the top genes expressed in each cluster.
#'
#' Here are our samples with assigned identities:
#'

#+ echo=FALSE
wt_aml_marrow_leiden_assignment_umap <-
  bb_var_umap(
    cds_wt_aml_marrow,
    "leiden_assignment",
    overwrite_labels = T,
    group_label_size = 4)

#+ echo=TRUE, dev="png", dpi=300, fig.height=4.0, fig.width = 4.5
wt_aml_marrow_leiden_assignment_umap

#' Here are the top markers for each cluster:
#+ echo=TRUE
write_csv(tm_wt_aml_marrow, file = str_glue("{tables_out}/tm_wt_aml_marrow.csv"))

#' The next thing to do is to determine the relative proportion of cells in each of these defined clusters in the two samples.  Fold change indicates normalized cell count in AML over WT.
#'

#+ echo=FALSE
wt_aml_cluster_proportions <- colData(cds_wt_aml_marrow) %>%
  as_tibble() %>%
  group_by(leiden_assignment, specimen) %>%
  summarise(n = n()) %>%
  left_join(
    colData(cds_wt_aml_marrow) %>%
      as_tibble() %>%
      group_by(specimen) %>%
      summarise(specimen_total =  n())
  ) %>%
  mutate(overall_total = nrow(colData(cds_wt_aml_marrow))) %>%
  mutate(normalized_count = n*overall_total/specimen_total/2) %>%
  select(leiden_assignment, specimen, normalized_count) %>%
  pivot_wider(names_from = specimen, values_from = normalized_count, values_fill = 1) %>%# adds pseudocount of 1 cell
  mutate(fold_change = P5/P8) %>%
  mutate(l2fc = log2(fold_change)) %>%
  mutate(enriched = ifelse(fold_change>1, "AML", "WT")) %>%
  rename(P5_normalized = P5, P8_normalized = P8)


#+ echo=TRUE
kable(wt_aml_cluster_proportions)


#+ echo=FALSE
wt_aml_cluster_proportion_plot <-
  ggplot(wt_aml_cluster_proportions,
         mapping = aes(
           x = reorder(leiden_assignment, l2fc),
           y = l2fc,
           fill = enriched
         )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = experimental_group_palette, breaks = c("AML", "WT")) +
  labs(x = "Cluster", y = "Population L2FC") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#+ echo=TRUE, dev="png", dpi=300
wt_aml_cluster_proportion_plot

#' These numbers are dominated by the fact that there are so many mature neutrophils in the WT sample.  Hopefully we are safe in the assumption that this is not an artifact of processing, but rather is a feature of teh differentiation arrest in this mouse.  Because of this discrepancy in the number of neutrophils, the numbers of HCSs and lymphocyte subsets are correspondingly low in WT and higher in AML.
#'
#' More interesting are the Blast-like cells which had no clear cognate in the reference data and subset of which seem to be highly specific to AML.  I calcualted the top specific markers for the Blast cells vs all other cells, versus HSC/Prog and versus Neutrophils.
#'

#+ echo=TRUE
write_csv(tm_blast_like_other, str_glue("{tables_out}/tm_blast_like_other.csv"))
write_csv(tm_blast_like_hsc_prog, str_glue("{tables_out}/tm_blast_like_hsc_prog.csv"))
write_csv(tm_blast_like_neu, str_glue("{tables_out}/tm_blast_like_neu.csv"))

#' ## Conclusion
#'
#' You should check out the gene lists I have generated.  They are in my X drive in the tables directory for this project.
#'
#' GSEA/IPA analysis are good ways to go.
#'
#' Please let me know if you have any questions about this analysis.  I will work on question 2 and any followups to question 1 over the next several days as time allows.
