#AML vs WT all cell pseudobulk comparision

# bb_var_umap(
#   filter_cds(
#     cds = cds_main,
#     cells = bb_cellmeta(cds_main) |>
#       filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
#   "leiden_assignment2",
#   #alt_dim_x = "aggr_UMAP_1",
#   #alt_dim_y = "aggr_UMAP_2",
#   overwrite_labels = F,
#   facet_by = "geno_pheno"
# )/
#   bb_var_umap(
#     filter_cds(
#       cds = cds_main,
#       cells = bb_cellmeta(cds_main) |>
#         filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
#     var = "density",
#     facet_by = "geno_pheno",
#     #alt_dim_x = "aggr_UMAP_1",
#     #alt_dim_y = "aggr_UMAP_2"
#   )

#Supp Fig 6A
bb_cellmeta(filter_cds(
  cds = cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML", "No leukemia")))) |>
  group_by(sample, tissue, leukemia_phenotype) |>
  summarise(n = n())
unique(colData(wt_aml)$leukemia_phenotype)

exp_design <-
  bb_cellmeta(filter_cds(
    cds = cds_main,
    cells = bb_cellmeta(cds_main) |>
      filter(leukemia_phenotype %in% c("AML", "No leukemia")))) |> #wt_aml
  group_by(sample, leukemia_phenotype) |>
  summarise()
#exp_design

# pseudobulk_res <-
#   bb_pseudobulk_mf(cds = wt_aml,
#                    pseudosample_table = exp_design,
#                    design_formula = "~ leukemia_phenotype",
#                    result_recipe = c("leukemia_phenotype", "AML", "No leukemia"))
pseudobulk_res <-
  bb_pseudobulk_mf(cds = filter_cds(
    cds = cds_main,
    cells = bb_cellmeta(cds_main) |>
      filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
    pseudosample_table = exp_design,
    design_formula = "~ leukemia_phenotype",
    result_recipe = c("leukemia_phenotype", "AML", "No leukemia"))

#less conservative approach (pseudobulk is a very conservative approach)
#bb_monocle_regression(cds = wt_aml, gene_or_genes = "Mpo", form = "~genotype")

pseudobulk_res$Header

pseudobulk_res$Result |> filter(gene_short_name == "Mpo")
pseudobulk_res$Result |> filter(gene_short_name == "Cd34")
# Fig6_wt_aml_clust1_all_pseudobulk<- pseudobulk_res$Result
F6_aml_wt_pseudobulk<- pseudobulk_res$Result

#write.csv(F6_aml_wt_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F6_aml_wt_pseudobulk.csv")
F6_aml_wt_pseudobulk<- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F6_aml_wt_pseudobulk.csv")

#Volcano Plot:
# Differential expression results.  Positive L2FC indicates up in B vs T upregulated
genes_to_highlight <- unique(c("Prox1", "Ifi44l", "Tdrd5", "Etv4", "Etv5")) #"Cd34", "Mpo", "Klf4", "Il7r" #"Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5"
#genes_to_highlight <- filter(F6_aml_wt_pseudobulk, padj < 0.001 & abs(log2FoldChange) >= 5)|>pull(gene_short_name)
#genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(pseudobulk_res$Result, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]
#genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(F6_aml_wt_pseudobulk, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]


volcano_data <- F6_aml_wt_pseudobulk %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))
#write.csv(volcano_data, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/volcano_data.csv")

#Pu wants to highlight the annotated genes and expand the volcano data through removal of Slc4a8
volcano_data<-volcano_data[!(volcano_data$gene_short_name=="Slc4a8"),]

#install.packages("ggbreak")
#library(ggbreak)
library(ggtext)
volcano_pseudobulk <-
  ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21,
             size = 0.5,
             alpha = 0.4) +
  geom_text_repel(color = "black",
                  fontface = "italic",
                  box.padding = 0.5, #0.5
                  point.padding = 0.25, #0.25
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  size = 3,
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log<sub>2</sub> fold change") +
  ylab("-log<sub>10</sub> adjusted p-value") +
  theme(axis.title.x =  element_markdown()) +
  theme(axis.title.y = element_markdown()) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(title = "Pseudobulk dKO AML vs WT")+ #caption = "\U21D0 Up in WT\nUp in AML \U21D2",
  theme(plot.caption.position = "panel") +
  #theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-1.0 * max(abs(range(volcano_data |> dplyr::filter(!is.na(padj)) |> pull(log2FoldChange)))), 1.0 * max(abs(range(volcano_data |> filter(!is.na(padj)) |> pull(log2FoldChange)))))) #+
#ggbreak::scale_y_cut(breaks=c(35), which=c(1), scales=c(0.1, 10))
volcano_pseudobulk

#Supp Figure 6B
#install.packages("msigdbr")
library(msigdbr)
library(topGO)
goenrichment <-
  bb_goenrichment(
    query = dplyr::filter(pseudobulk_res$Result, padj < 0.05 &
                            log2FoldChange >= 0.58) |> pull(gene_short_name),
    reference = bb_rowmeta(cds_main),
    go_db = "org.Mm.eg.db"
  )


gosummary <- bb_gosummary(x = goenrichment,
                               reduce_threshold = 0.85,
                               go_db = "org.Mm.eg.db")
pseudobulk_AML_vs_WT_GO_PCA <-
  bb_goscatter(simMatrix = gosummary$simMatrix,
               reducedTerms = gosummary$reducedTerms)
pseudobulk_AML_vs_WT_GO_PCA

#Supp Figure 6C-F
#Need to recover from Pu's previous code

#Supp Figure 6G
