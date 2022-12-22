source("R/dependencies.R")
source("R/configs.R")
# original fig6b:  cluster umap faceted by sample ---------------------
# this is a loupe file you will have to find.

# modifications to CDS.  Should incorporate into data package -----------
aggr_umap_tbl <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/cellranger_aggr/lapalombella_pu_aggr/outs/count/analysis/umap/2_components/projection.csv", col_names = c("cell_id", "aggr_UMAP_1", "aggr_UMAP_2"), skip = 1) |>
   mutate(barcode_truncated = str_remove(cell_id, "[0-9]+|[0-9]+$")) |>
   mutate(sample_num = str_extract(cell_id,  '[0-9]+|[0-9]+$')) |>
   mutate(sample_name = recode(sample_num,
                               "1" = "P1",
                               "2" = "P2",
                               "3" = "P3",
                               "4" = "P4",
                               "5" = "P10_R0906_SPN",
                               "6" = "P12_R0909_SPN",
                               "7" = "P14_S1303_SPN_AML",
                               "8" = "P16_S1310_SPN_WT",
                               "9" = "P5",
                               "10" = "P6",
                               "11" = "P7",
                               "12" = "P8",
                               "13" = "P9_R0906_BM",
                               "14" = "P11_R0909_BM",
                               "15" = "P13_S1303_BM_AML",
                               "16" = "P15_S1310_BM_WT"
                               )) |>
   mutate(cell_id = paste0(barcode_truncated, "1_", sample_name))|> #|> count(sample_name)
   select(cell_id, aggr_UMAP_1, aggr_UMAP_2)

aggr_cluster_tbl <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/cellranger_aggr/lapalombella_pu_aggr/outs/count/analysis/clustering/kmeans_10_clusters/clusters.csv", col_names = c("cell_id", "cluster"), skip = 1) |>
  mutate(barcode_truncated = str_remove(cell_id, "[0-9]+|[0-9]+$")) |>
  mutate(sample_num = str_extract(cell_id, "[0-9]+|[0-9]+$")) |>
  mutate(sample_name = recode(sample_num,
                              "1" = "P1",
                              "2" = "P2",
                              "3" = "P3",
                              "4" = "P4",
                              "5" = "P10_R0906_SPN",
                              "6" = "P12_R0909_SPN",
                              "7" = "P14_S1303_SPN_AML",
                              "8" = "P16_S1310_SPN_WT",
                              "9" = "P5",
                              "10" = "P6",
                              "11" = "P7",
                              "12" = "P8",
                              "13" = "P9_R0906_BM",
                              "14" = "P11_R0909_BM",
                              "15" = "P13_S1303_BM_AML",
                              "16" = "P15_S1310_BM_WT"
                              )) |>
  mutate(cell_id = paste0(barcode_truncated, "1_", sample_name)) |>
  select(cell_id, kmeans10_cluster = cluster) |>
  mutate(kmeans10_cluster = as.character(kmeans10_cluster))

 cds_main <- bb_tbl_to_coldata(obj = cds_main, min_tbl = aggr_umap_tbl)
 cds_main <- bb_tbl_to_coldata(obj = cds_main, min_tbl = aggr_cluster_tbl)

#Create genotype/phenotype column in cds
unique(colData(cds_main)$leukemia_phenotype)
colData(cds_main)$leukemia_phenotype <-
  recode(
    colData(cds_main)$leukemia_phenotype,
    "PreB ALL" = "pre-B ALL",
    "T cell leukemia" = "T ALL"
  )
 colData(cds_main)$geno_pheno <-
   paste0(colData(cds_main)$genotype, " ", colData(cds_main)$leukemia_phenotype)

colData(cds_main)$geno_pheno <-
   recode(colData(cds_main)$geno_pheno,
          "WT No leukemia" = "Wildtype",
          "TP53-/-/TET2-/- AML" = "dKO: AML",
          "TP53-/-/TET2-/- pre-B ALL" = "dKO: pre-B ALL",
          "TP53-/- T ALL" = "P53 KO: T ALL"

          )
unique(colData(cds_main)$geno_pheno)
 #Order factor levels
 colData(cds_main)$kmeans10_cluster <- factor(colData(cds_main)$kmeans10_cluster,
                                              levels = 1:10)
 colData(cds_main)$geno_pheno <-
   factor(
     colData(cds_main)$geno_pheno,
     levels = c(
       "dKO: AML",
       "dKO: pre-B ALL",
       "P53 KO: T ALL",
       "Wildtype"
     )
   )
 colData(cds_main)$genotype <-
   factor(colData(cds_main)$genotype,
          levels = c("TP53-/-/TET2-/-", "TP53-/-", "WT"))
 #Make genotype/phenotype/tissue specific column:
 colData(cds_main)$pheno_tissue <-
   paste0(colData(cds_main)$geno_pheno, " ", colData(cds_main)$tissue)

 #where did these assignments come from
 #-partition 3 displays AML Blast Markers
 colData(cds_main)$partition_assignment <-
   recode(colData(cds_main)$partition,
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
          "12" = "T/NK 2",
          "13" = "Th1",
          "14" = "Th2",
          "15" = "Th3",
          "16" = "plasma cells"
   )
 ##############################################################################################

 #Figure 5A
bb_cellmeta(cds_main) |>
   group_by(pheno_tissue, specimen, primary_or_engraftment) |>
   summarise(n = n())

#By genotype
 A1 <-  bb_var_umap(
   cds_main,
   var = "kmeans10_cluster",
   alt_dim_x = "aggr_UMAP_1",
   alt_dim_y = "aggr_UMAP_2", cell_size = 0.1,
   overwrite_labels = T
 ) + facet_grid(col = vars(genotype)) +
   theme_minimal() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
   theme(panel.background = element_rect(color = "black")) +
   theme(legend.position = "none") +
   theme(axis.title.x = element_blank()) +
   theme(axis.title.y = element_blank())+
   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

 #By genotype on primary mice
 bb_var_umap(
   filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(primary_or_engraftment == "primary")),
   var = "kmeans10_cluster",
   alt_dim_x = "aggr_UMAP_1",
   alt_dim_y = "aggr_UMAP_2", cell_size = 0.1,
   overwrite_labels = T
 ) + facet_grid(col = vars(genotype)) +
   theme_minimal() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
   theme(panel.background = element_rect(color = "black")) +
   theme(legend.position = "none") +
   theme(axis.title.x = element_blank()) +
   theme(axis.title.y = element_blank())+
   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#Primary mice faceted by phenotype/tissue
 bb_var_umap(
   filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(primary_or_engraftment == "primary")),
   var = "kmeans10_cluster",
   alt_dim_x = "aggr_UMAP_1",
   alt_dim_y = "aggr_UMAP_2", cell_size = 0.1,
   overwrite_labels = T
 ) + facet_grid(col = vars(pheno_tissue)) +
   theme_minimal() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
   theme(panel.background = element_rect(color = "black")) +
   theme(legend.position = "none") +
   theme(axis.title.x = element_blank()) +
   theme(axis.title.y = element_blank())+
   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

 #Engrafted mice faceted by phenotype/tissue
   bb_var_umap(
     filter_cds(cds_main, cells = bb_cellmeta(cds_main)|> filter(primary_or_engraftment == "engraftment")),
     var = "kmeans10_cluster",
     alt_dim_x = "aggr_UMAP_1",
     alt_dim_y = "aggr_UMAP_2", cell_size = 0.1,
     overwrite_labels = T
   ) + facet_grid(col = vars(pheno_tissue)) +
     theme_minimal() + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()) +
     theme(panel.background = element_rect(color = "black")) +
     theme(legend.position = "none") +
     theme(axis.title.x = element_blank()) +
     theme(axis.title.y = element_blank())+
     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

 A2<-bb_var_umap(
   cds_main,
   var = "density",
   facet_by = "genotype",
   alt_dim_x = "aggr_UMAP_1",
   alt_dim_y = "aggr_UMAP_2", cell_size = 0.1
 ) + facet_grid(col = vars(genotype)) +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(panel.background = element_rect(color = "black"))+
   theme(axis.title.x = element_blank()) +
   theme(axis.title.y = element_blank())+
   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
   theme(legend.title=element_text(size=8)) +
   theme(legend.key.size = unit(0.3,"cm"))

F5A <-
   as_ggplot(grid.arrange(
     patchworkGrob(A1 / A2),
     left = textGrob(A1$labels$y, rot=90, vjust = 1.5, hjust=0.35),
     bottom = textGrob(A1$labels$x, hjust = 0.8, vjust = -0.5))
   )

# make the gene expression umaps
aml_plotlist <- map(.x = c("Cd34","Mpo", "Kit", "Elane", "Calr","Ctsg"),
                     .f = \(x, dat = cds_main) {
                       p <- bb_gene_umap(
                         dat,
                         gene_or_genes = x,
                         alt_dim_x = "aggr_UMAP_1",
                         alt_dim_y = "aggr_UMAP_2", cell_size = 0.1
                       ) +
                         scale_color_distiller(palette = "Oranges",
                                               direction = 1,
                                               na.value = "grey80",
                                               limits = c(0,2.5)) +
                         facet_wrap(~geno_pheno, labeller = labeller(group = label_wrap_gen(width = 5, multi_line = TRUE))) +
                         scale_y_continuous(breaks = c(-15,0, 10)) +
                         scale_x_continuous(breaks = c(-10,0, 10))+
                         theme(panel.spacing = unit(0.5, "lines"))+
                         theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())+
                         theme(panel.background = element_rect(color = "black",
                                                               fill = "white"))+
                         theme(axis.line = element_blank()) +
                         #theme(axis.ticks = element_blank()) +
                         #theme(axis.text = element_blank()) +
                         labs(x =NULL, y = NULL, subtitle = x) +
                         theme(plot.subtitle = element_text(face = "italic", hjust = 0.5))
                       # if (x != "Calr") p <- p + theme(legend.position = "none")
                       p
                     })

 aml_gexp_umap <- ggarrange(aml_plotlist[[1]],
           aml_plotlist[[2]],
           aml_plotlist[[3]],
           aml_plotlist[[4]],
           aml_plotlist[[5]],
           aml_plotlist[[6]],
           ncol = 3,
           nrow=2,
           common.legend = TRUE,
           legend="right")

aml_gexp_umap

#Heatmap
F5_topmarkers_k10 <-
  monocle3::top_markers(
    cds_main,
    group_cells_by = "kmeans10_cluster",
    genes_to_test_per_group = 50,
    cores = 12
  )
#write_csv(F5_topmarkers_k10, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data", "F5_topmarkers_k10.csv"))
F5_topmarkers_k10 <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F5_topmarkers_k10.csv")

# F6_topmarkers2 <-
#   monocle3::top_markers(
#     cds_main,
#     group_cells_by = "leukemia_phenotype",
#     genes_to_test_per_group = 20,
#     cores = 10
#   )

markers <- F5_topmarkers_k10 |> filter(
   cell_group %in% c('6','8','1','5','3'))|>pull(gene_short_name)

#Expression Matrix
  # mat <-
  #   bb_aggregate(
  #     obj = filter_cds(
  #       cds_main,
  #       cells = bb_cellmeta(cds_main) |>
  #         filter(
  #           leukemia_phenotype %in% c("AML", "B-ALL", "T-ALL", "No leukemia")
  #         ),
  #       #clusters of interest instead?
  #       genes = bb_rowmeta(cds_main) |>
  #         filter(gene_short_name %in% markers)
  #     ),
  #     cell_group_df = bb_cellmeta(cds_main) |>
  #       select(cell_id, leukemia_phenotype)
  #   ) |>
  #   t() |>
  #   scale() |>
  #   t()
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          kmeans10_cluster %in% c('3','6','8','1','5')
        ),
      genes = bb_rowmeta(cds_main) |>
        filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, kmeans10_cluster)
  ) |>
  t() |>
  scale() |>
  t()
  rownames(mat) <-
    tibble(feature_id = rownames(mat)) |>
    left_join(bb_rowmeta(cds_main) |>
                select(feature_id, gene_short_name)) |>
    pull(gene_short_name
    )

heatmap_3_colors <-
    c("#313695", "white", "#A50026")

  colfun = circlize::colorRamp2(breaks = c(min(mat),
                                           0,
                                           max(mat)),
                                colors = heatmap_3_colors)
  heatmap_highlights <- c(
    "Mpo",
    "Ctsg",
    "Elane",
    "Ybx1",
    "Vpreb1",
    "Cd3d",
    "Tnnt1"
  )

  anno <-
    ComplexHeatmap::rowAnnotation(link =  anno_mark(
      at = which(rownames(mat) %in% heatmap_highlights),
      labels = rownames(mat)[rownames(mat) %in% heatmap_highlights],
      labels_gp = gpar(fontsize = 7),
      padding = 4
    ))

 F5heatmap<-grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
    mat,
    col = colfun,
    name = "Expression",
    show_row_names = F,
    show_column_names = F,
    right_annotation = anno,
    #top_annotation = hmap_bp,
    row_dend_width = unit(4, "mm"),
    column_dend_height = unit(4, "mm"),
    heatmap_legend_param = list(
      legend_direction = "vertical",
      #legend_width = unit(1, "mm"),
      title_position = "lefttop-rot",
      title_gp = gpar(fontsize = 10)
    )
  )))

#stacked bar chart
####fraction of cells contributed to each cluster by leukemia_phenotype
cellcount<- bb_cellmeta(cds_main) |>
  group_by(kmeans10_cluster, leukemia_phenotype) |>
  summarise(n = n()) |> filter(kmeans10_cluster %in% c("5", "6", "1", "8","3"))
# mutate(cellcount, frac = n/group_by(kmeans10_cluster) |> summarise(n = sum(n)))
#create fraction column
library(data.table)
setDT(cellcount)[, frac := n / sum(n), by=kmeans10_cluster]

#factor levels
cellcount$kmeans10_cluster <- factor(cellcount$kmeans10_cluster,
                                     levels = c("5", "6","1","8","3"))
cellcount$leukemia_phenotype <- factor(cellcount$leukemia_phenotype,
                                     levels = c("AML", "B-ALL","T-ALL", "No leukemia"))
#plot
hmap_bp <-
  ggplot(cellcount,
         aes(x = kmeans10_cluster, y = frac, fill = leukemia_phenotype)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 12)) +
  coord_fixed(ratio = 0.45)+
  theme(plot.margin = unit(c(0,0.4,0,0), "cm")) +
  labs(x = "Clusters (k-means10)")+
  scale_x_discrete(expand = c(0.13,0))

#library(gtable)
#F5D<- F5heatmap / hmap_bp + plot_layout(heights = c(3,1))

#library(cowplot)
F5hm<- plot_grid(F5heatmap,hmap_bp, ncol=1, rel_heights = c(3,1))

#Pseudobulk:
##AML cluster 1 (17401 cells) vs WT cluster 1 (400 cells)

# bb_cellmeta(cds_main) |>
#   group_by(kmeans10_cluster, leukemia_phenotype) |>
#   summarise(n = n()) |> filter(kmeans10_cluster %in% c("1"))

wt_aml <- blaseRtools::filter_cds(
  cds = cds_main,
  cells = bb_cellmeta(cds_main) |> filter(kmeans10_cluster %in% c("1")) |>
    filter(leukemia_phenotype %in% c("AML", "No leukemia")))

bb_var_umap(wt_aml,
            "kmeans10_cluster",
            overwrite_labels = T,
            cell_size = 0.2, alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2") + facet_grid(
              col = vars(tissue),
              row = vars(leukemia_phenotype, primary_or_engraftment)
            ) + theme(legend.position = "right")

bb_cellmeta(wt_aml) |>
  group_by(sample, tissue, leukemia_phenotype) |>
  summarise(n = n())

exp_design <-
  bb_cellmeta(wt_aml) |>
  group_by(sample, leukemia_phenotype) |>
  summarise()
exp_design

pseudobulk_res <-
  bb_pseudobulk_mf(cds = wt_aml,
                   pseudosample_table = exp_design,
                   design_formula = "~ leukemia_phenotype",
                   result_recipe = c("leukemia_phenotype", "AML", "No leukemia"))

#less conservative approach (pseudobulk is a very conservative approach)
#bb_monocle_regression(cds = wt_aml, gene_or_genes = "Mpo", form = "~genotype")

pseudobulk_res$Header

pseudobulk_res$Result |> filter(gene_short_name == "Mpo")
pseudobulk_res$Result |> filter(gene_short_name == "Cd34")
Fig6_wt_aml_clust1_all_pseudobulk<- pseudobulk_res$Result
write.csv(Fig6_wt_aml_clust1_all_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/Fig6_wt_aml_clust1_all_pseudobulk.csv")

#Volcano Plot:
# Differential expression results.  Positive L2FC indicates up in B vs T upregulated
genes_to_highlight <- unique(c("Cd34", "Mpo", "Klf4", "Il7r"))
genes_to_highlight <- genes_to_highlight[genes_to_highlight %in% (filter(pseudobulk_res$Result, padj < 0.1 & abs(log2FoldChange) >= 0.58)|>pull(gene_short_name))]


volcano_data <- pseudobulk_res$Result %>%
  mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))
#write.csv(volcano_data, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/volcano_data.csv")

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
  labs(caption = "\U21D0 Up in WT\nUp in AML \U21D2",title = "Pseudobulk Cluster 1: dKO AML vs WT")+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-1.0*max(abs(range(volcano_data |> dplyr::filter(!is.na(padj)) |> pull(log2FoldChange)))), 1.0*max(abs(range(volcano_data |> filter(!is.na(padj)) |> pull(log2FoldChange))))))
volcano_pseudobulk
ggsave("volcano_pseudob_clust1_wt_vs_aml.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Figures")

#fgsea:
install.packages("msigdbr")
install.packages("fgsea")
library(msigdbr); library(fgsea)

#pseudbulk_data<- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/Fig6_wt_aml_clust1_all_pseudobulk.csv")
#pseudobulk_data[duplicated(GSEA_input)]
distinct(Fig6_wt_aml_clust1_all_pseudobulk$gene_short_name)
pseudobulk_data<- Fig6_wt_aml_clust1_all_pseudobulk[!duplicated(Fig6_wt_aml_clust1_all_pseudobulk$gene_short_name), ]
##Rank
GSEA_input <- mutate(pseudobulk_data, Rank = -log10(pseudobulk_data$padj)*pseudobulk_data$log2FoldChange)
GSEA_input <- GSEA_input[complete.cases(GSEA_input), ]
#GSEA_input[duplicated(GSEA_input)]
rankData <- GSEA_input$Rank
names(rankData) <- GSEA_input$gene_short_name
# rankData[duplicated(rankData)]

msigdbr_df <- rbind(msigdbr(species = "Mus musculus", category = "H"))#, #Hallmark
                 # msigdbr(species = "Mus musculus", category = "C6"), #Cancer relevant
                 # msigdbr(species = "Mus musculus", category = "C5", subcategory= "GO:BP"), #Gene ontology-biological process
                 # msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")) #immune relevant

pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

fgseaRes <- fgsea(pathwaysH,
                  rankData,
                  minSize=15,
                  maxSize = 300, scoreType = "std", eps = 0) #, nperm=1000)
#fgsea has a default lower bound eps=1e-10 for estimating P-values.
#If you need to estimate P-value more accurately, you can set the eps argument to zero
#in the fgsea function.

#Collapse redundant pathways
fgseaRes2 <- fgseaRes[match(collapsePathways(fgseaRes, pathwaysH, rankData, pval.threshold = 0.05, nperm = 200, gseaParam = 1)$mainPathways, fgseaRes$pathway), ]

# Number of Results w/significant hits after FDR correction
sum(fgseaRes[, padj < 0.01])
sum(fgseaRes[, padj < 0.05])
sum(fgseaRes2[, padj < 0.01])
sum(fgseaRes2[, padj < 0.05])

fgseaRes2 |>
  arrange(desc(abs(NES))) |>
  top_n(10, -padj) #can adjust top_n # for sig hits after FDR correction

fgseaResTidy <- fgseaRes2 |>
  as_tibble() |>
  arrange(desc(NES))
# Show in a nice table:
TT<-as.tibble(fgseaResTidy)

topPathways <- fgseaRes |>
  top_n(10, wt=-padj) |>
  arrange(-NES) |>
  pull(pathway)

plotGseaTable(pathwaysH[topPathways],
              rankData,
              fgseaRes,
              gseaParam = 0.5)

GSEA_HM<- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA Hallmark Pathways")

# Number of Results w/significant hits after FDR correction
sum(fgseaRes[, padj < 0.01])
sum(fgseaRes2[, padj < 0.05])
fgseaRes <- fgseaRes |> top_n(8, wt=-padj) #filter for # of sig hits after FDR correction
fgseaRes3 <- fgseaRes |> top_n(7, wt=-padj) #filter for # of sig hits after FDR correction
HM_PseudoBulkGSEA<- ggplot(fgseaRes3, aes(y = reorder(pathway, NES), x = NES, size = size)) +
  geom_point(aes(color = padj), alpha = 1) +
  scale_color_gradient(low = '#FDE725',high = "#414487") +
  #scale_color_viridis_b() +
  ggtitle("Pseudobulk Cluster 1: AML vs WT") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 9)) +labs(y = "Pathways (Hallmark)")
HM_PseudoBulkGSEA

plotEnrichment(pathwaysH[["HALLMARK_GLYCOLYSIS"]], rankData)  +
  labs(title =
         "AML vs WT: Hallmark Glycolysis",
       subtitle = paste0(
         "NES = ",
         round(fgseaRes$NES[19], digits = 2),
         ",  Adjusted pval = ",
         formatC(fgseaRes$padj[19], format = "e", digits = 2)
       ))

#T cell exhaustion score Zheng et al 2021 (https://www.science.org/doi/pdf/10.1126/science.abe6474):




######################################################################
#Pu partition assignment
 bb_genebubbles(
   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
   genes = c("Cd19",
             "Cd79a",
             "Cd3d",
             "Cd4",
             "Cd14",
             "Itgam",
             "Cd8a",
             "Cd177", #(Netrophils: Cd177>Itgam>Cd14, Monocytes: Cd14>Itgam>Cd177)
             "Pdcd1",# Cd4 T
             #"Foxp3", #T?
             "Cd34", #AML Blast/HSC
             #"Cd33", #AML Blast/Myeloid lineage
             "Kit", #CD117/cKit - Myeloid Blast Differentiation
             "Il3ra" #Cd123 - Leukemic Stem Cells
   ), cell_grouping = "partition") + labs(x = "Partition Clusters", y = NULL)


 # colData(cds_main)$partition_assignment <-
 #   recode(colData(cds_main)$partition,
 #          "1" = "B",
 #          "2" = "Cd4+ T",
 #          "3" = "Blast-like",
 #          "4" = "Cd8+ T",
 #          "5" = "B",
 #          "6" = "Neu/Blast-like",
 #          "7" = "Neu/Cd33+ Blast-like",
 #          "8" = "T",
 #          "9" = "Neu",
 #          "10" = "B",
 #          "11" = "Unassigned",
 #          "12" = "Unassigned",
 #          "13" = "Unassigned",
 #          "14" = "Unassigned",
 #          "15" = "Mixed?",
 #          "16" = "B"
 #   )

 #kmeans10 assignment
 bb_genebubbles(
   obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
   genes = c("Cd19",
             "Cd79a",
             "Cd3d",
             "Cd4",
             "Cd14",
             "Itgam",
             "Cd8a",
             "Cd177", #(Netrophils: Cd177>Itgam>Cd14, Monocytes: Cd14>Itgam>Cd177)
             "Pdcd1",# Cd4 T
             #"Foxp3", #T?
             "Cd34", #AML Blast/HSC
             "Cd33", #AML Blast/Myeloid lineage
             "Kit", #CD117/cKit - Myeloid Blast Differentiation
             "Il3ra" #Cd123 - Leukemic Stem Cells
   ), cell_grouping = "kmeans10_cluster") + labs(x = "kmeans10 Clusters", y = NULL)

 colData(cds_main)$k10_assignment <-
   recode(
     colData(cds_main)$kmeans10_cluster,
     "1" = "AML Blasts",
     "2" = "Unassigned",
     "3" = "B-ALL",
     "4" = "Cd33+",
     "5" = "Cd4+ T",
     "6" = "B-ALL", #Cd8+ cells here too
     "7" = "Cd8+ T",
     "8" = "B-ALL",
     "9" = "Mixed?",
     "10" = "Unassigned"
   )

#Figure 6A
 #By Leukemia Phenotype
   # bb_var_umap(
   #   cds_main,
   #   var = "kmeans10_cluster",
   #   alt_dim_x = "aggr_UMAP_1",
   #   alt_dim_y = "aggr_UMAP_2",
   #   overwrite_labels = T
   # ) + facet_grid(col = vars(leukemia_phenotype)) + labs(x = "UMAP 1", y = "UMAP 2") +
   # theme_minimal() + theme(panel.grid.major = element_blank(),
   #                         panel.grid.minor = element_blank()) +
   # theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")


A2 <- bb_var_umap(
  cds_main,
  var = "partition",
  overwrite_labels = T
) + facet_grid(col = vars(genotype)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")
A5 <- bb_var_umap(
  cds_main,
  var = "partition_assignment",
  overwrite_labels = T
) + facet_grid(col = vars(genotype)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")

  #6A Density by genotype
A3<-bb_var_umap(
     cds_main,
     var = "density",
     facet_by = "genotype",
     alt_dim_x = "aggr_UMAP_1",
     alt_dim_y = "aggr_UMAP_2"
   ) + facet_grid(col = vars(genotype)) + labs(x = "UMAP 1", y = "UMAP 2")+ theme_minimal() + theme(panel.grid.major = element_blank(),
                                                                            panel.grid.minor = element_blank()) +theme(panel.background = element_rect(color = "black"))
A4<-bb_var_umap(
  cds_main,
  var = "density",
  facet_by = "genotype"
) + facet_grid(col = vars(genotype)) + labs(x = "UMAP 1", y = "UMAP 2")+ theme_minimal() + theme(panel.grid.major = element_blank(),
                                                                                                 panel.grid.minor = element_blank()) +theme(panel.background = element_rect(color = "black"))
A<-A1/A3 #kmeans10 - alt dims
A2<- A2/A4 #partition
ggsave("figure_6A_partition.pdf", path = pu_figs)
A3<-A5/A4
#ggsave("figure_6A_partition_assignment.pdf", path = pu_figs)


#Figure 6B
# bb_gene_umap(cds_main, gene_or_genes = c("Cd34"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+ facet_wrap(~leukemia_phenotype)+labs(title = "Cd34", y = "UMAP2")+
#   theme(legend.position = "none")+
#   theme(axis.title.x =element_blank())+
#   scale_color_distiller(palette = "Oranges",
#                         direction = 1,
#                         na.value = "grey80", limits = c(0,2.5))

#Phenotypic Markers
#AML: c(Mpo,Cd34,Kit,Cd11b,Elane,Calr, Ctsg)
#T-ALL: c(Sca1 = Atxn1, Ly6a, Cd7, Cd3, Tdt, Cd34)
#B-ALL: c(Cd19, Pax5, Cd24a, Cd79a, Vpreb1, Vpreb2, Vpreb3, Tdt)

#bb_gene_umap(cds_main, gene_or_genes = c("Cd34"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")

#bb_gene_umap(cds_main, gene_or_genes = c("Cd19"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")
b_all_plotlist <- map(.x = c("Cd19", "Pax5", "Cd24a", "Cd79a", "Vpreb1", "Vpreb2", "Vpreb3", "Dntt"),
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.2
                      ) +
                        scale_color_distiller(palette = "Greens",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) +
                        facet_wrap(~geno_pheno) + scale_y_continuous(breaks = c(-15,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, -10))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      if (x != "Vpreb3") p <- p + theme(legend.position = "none")
                      p
                    })

b_all_plotlist[[1]] / b_all_plotlist[[2]] |
  b_all_plotlist[[3]] / b_all_plotlist[[4]] |
                           b_all_plotlist[[5]] / b_all_plotlist[[6]] |
                           b_all_plotlist[[7]] / b_all_plotlist[[8]]

bb_gene_umap(cds_main, gene_or_genes = c("Dntt"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+ facet_wrap(~leukemia_phenotype)

t_all_plotlist <- map(.x = c("Atxn1", "Ly6a", "Cd7", "Cd3d", "Dntt", "Cd34"), #Atxn1 = Sca1
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        alt_dim_x = "aggr_UMAP_1",
                        alt_dim_y = "aggr_UMAP_2", cell_size = 0.2
                      ) +
                        scale_color_distiller(palette = "Blues",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) + #check limits!
                        facet_wrap(~geno_pheno) + scale_y_continuous(breaks = c(-15,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, -10))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      if (x != "Tdt") p <- p + theme(legend.position = "none")
                      p
                    })
t_all_plotlist[[1]] / t_all_plotlist[[2]] |
  t_all_plotlist[[3]] / t_all_plotlist[[4]] |
t_all_plotlist[[5]] / t_all_plotlist[[6]]

bb_var_umap(
  cds_main,
  var = "k10_assignment",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  overwrite_labels = T
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")

bb_var_umap(
  cds_main,
  var = "kmeans10_cluster",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  overwrite_labels = T
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")

#AML Markers:
####

celltypes <- bb_genebubbles(
  cds_main,
  genes = c(
    "Cd19",
    "Cd79a",
    "Cd3d",
    "Cd4",
    "Cd14",
    "Itgam",
    "Cd8a",
    "Cd177", #(Netrophils: Cd177>Itgam>Cd14, Monocytes: Cd14>Itgam>Cd177)
    "Pdcd1",# Cd4 T
    #"Foxp3", #T?
    "Cd34", #AML Blast
    "Cd33", #AML Blast (Cd34+ but Cd33- is APL) (Cd34+ but Cd33- - Acute Leukemia)
    "Kit", #CD117/cKit - Myeloid Blast Differentiation
    "Il3ra"
  ),
  cell_grouping = c("kmeans10_cluser", "k10_assignment"),
  return_value = "data"
) |>
  ggplot(mapping = aes(x = kmeans10_cluster,
                       y = gene_short_name,
                       color = expression,
                       size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap(~k10_assignment, scales = "free_x", ) +
  theme_minimal_grid(font_size = the_font_size) +
  theme(strip.background = ggh4x::element_part_rect(side = "b", colour = "black", fill = "transparent")) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")
celltypes


#Pu Celltype assignment:
# subset to include only wt and aml cells with old dimensions
olddims_wt_aml <- blaseRtools::filter_cds(cds = cds_main,
                                          cells = bb_cellmeta(cds_main) |>
                                            filter(leukemia_phenotype %in% c("AML", "No leukemia")))

# plot the wt aml cells with old dimensions
bb_var_umap(
  olddims_wt_aml,
  "partition",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  facet_by = "genotype"
)

colData(olddims_wt_aml)$partition_assignment <-
  recode(colData(olddims_wt_aml)$partition,
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
         "12" = "T/NK 2",
         "13" = "Th1",
         "14" = "Th2",
         "15" = "Th3",
         "16" = "plasma cells"
  )
bb_var_umap(
  olddims_wt_aml,
  "partition_assignment",
  overwrite_labels = T,
  facet_by = "geno_pheno"
)
#with loupe dims
bb_var_umap(
  olddims_wt_aml,
  "partition_assignment",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2",
  overwrite_labels = T,
  facet_by = "geno_pheno"
)/
bb_var_umap(
  olddims_wt_aml,
  var = "density",
  facet_by = "geno_pheno",
  alt_dim_x = "aggr_UMAP_1",
  alt_dim_y = "aggr_UMAP_2"
)

#Pu Heatmap Code---------------------------------------------------------
cds_p53tet2AML<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "AML"]

cds_p53tet2AMLann<-bb_cds_anno(query_cds=cds_p53tet2AML, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )

# I added the assignment for top20
top20 <- tm_wt_aml_marrow |>
  filter(cluster_method == "partition") |>
  group_by(cell_group) |>
  slice_min(order_by = marker_test_q_value, n = 20) |>
  pull(gene_short_name)
bb_cellmeta(cds_p53tet2AMLann)
colData(cds_p53tet2AMLann)

top20

# put in a new row metadata column
# you didnt define your top20 genes here
cds_p53tet2AML<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "AML"]
cds_p53tet2AMLann<-bb_cds_anno(query_cds=cds_p53tet2AML, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )

rowData(cds_p53tet2AMLann)$top20 <- ifelse(rowData(cds_p53tet2AMLann)$gene_short_name %in% top20,
                                           "yes",
                                           "no")

# filter the cds and pipe into aggregate gene expression
agg_mat_wt_aml_olddim <-
  olddims_wt_aml |>
  filter_cds(genes = bb_rowmeta(cds_p53tet2AMLann) |>
               filter(top20 == "yes")) |>
  aggregate_gene_expression(cell_group_df = bb_cellmeta(olddims_wt_aml) |>
                              select(cell_id, partition_assignment))


max(agg_mat_wt_aml_olddim)

min(agg_mat_wt_aml_olddim)


# convert from sparse to regular matrix
agg_mat_wt_aml_olddim <- as.matrix(agg_mat_wt_aml_olddim)
max(agg_mat_wt_aml_olddim)

agg_mat_wt_aml_olddim


# fix the rownames
rownames(agg_mat_wt_aml_olddim) <-
  left_join(tibble(feature_id = rownames(agg_mat_wt_aml_olddim)),
            bb_rowmeta(olddims_wt_aml)) |>
  pull(gene_short_name)
max(agg_mat_wt_aml_olddim)
# transpose and then put all of the genes (columns) on the same scale
agg_mat_wt_aml_olddim <- scale(t(agg_mat_wt_aml_olddim))



# make a list of genes you want to point out
# put whatever genes you want from top20 here
heatmap_highlights <- c(
  "S100a8",
  "S100a9",
  "Mpo",
  "Ctsg",
  "Elane",
  "Ybx1",
  "Gzmb",
  "Tox",
  "Cd8b1",
  "Vpreb1",
  "Cd3d"
)
# make the heatmap color scale
# see configs.R for what these colors are
col_fun_heatmap <-
  colorRamp2(breaks = c(min(agg_mat_wt_aml_olddim),
                        0,
                        max(agg_mat_wt_aml_olddim)),
             colors = heatmap_3_colors)


# make the annotation object
heatmap_anno_df <-
  map(
    .x = heatmap_highlights,
    .f = function(x) {
      index <- which(colnames(agg_mat_wt_aml_olddim) == x)
      return(index)
    }
  ) %>% set_names(heatmap_highlights) %>%
  bind_cols() %>%
  pivot_longer(everything()) %>%
  as.data.frame()

heatmap_gene_anno <- HeatmapAnnotation(
  foo = anno_mark(
    at = heatmap_anno_df$value,
    labels = heatmap_anno_df$name,
    labels_gp = gpar(fontsize = 8),
    padding = 1.5,
    labels_rot = 45
  ),
  which = "column"
)


# make the heatmap finally
partition_heatmap <- grid.grabExpr(draw(
  Heatmap(
    matrix = agg_mat_wt_aml_olddim, # you had the wrong matrix here
    col = col_fun_heatmap,
    name = "Expression",
    heatmap_legend_param = list(
      title_gp = gpar(fontface = "plain", fontsize = 9),
      grid_width = unit(0.14, "in"),
      labels_gp = gpar(fontsize = 8)
    ),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    column_dend_side = "bottom",
    show_row_names = T,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = F,
    top_annotation = heatmap_gene_anno,
    row_dend_gp = gpar(lwd = 0.5),
    column_dend_gp = gpar(lwd = 0.5),
    row_title = "Partition",
    column_title = "Top 20 Genes"
  )
), wrap = T)

# plot the heatmap
plot_grid(partition_heatmap)

# save the heatmap
save_plot(
  plot_grid(partition_heatmap),
  filename = file.path(figs_out, "wt_aml_heatmap.png"),
  base_width = 5.5,
  base_height = 3.5
)



#################################################################################
#generate cds subset by pt & B cells (via clonotype_id)
cds_subset2712 <- cds_main[, colData(cds_main)$patient == "pt_2712" &
                             colData(cds_main)$clonotype_id %in% "clonotype1"]

cds_subset1245 <- cds_main[, colData(cds_main)$patient == "pt_1245" &
                             colData(cds_main)$clonotype_id %in% "clonotype1"]

#bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP"))))
#blaseRtools::bb_gene_pseudotime

cds_subset2712<- order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP")))
colData(cds_subset2712)

cds_subset2712 <- cluster_cells(cds_subset2712)
cds_subset2712 <- learn_graph(cds_subset2712)
colData(cds_subset2712)
plot_cells(cds_subset2712,
           color_cells_by = "sample",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds_subset2712 <- order_cells(cds_subset2712)

#manually selected root nodes
plot_cells(cds_subset2712,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

monocle3::graph_test
#use trace and substitute Matrix::rBind with rbind
trace('calculateLW', edit = T, where = asNamespace("monocle3"))

gtest <- monocle3::graph_test(cds_subset2712, neighbor_graph="principal_graph", cores=4)

pr_deg_ids <- row.names(subset(gtest, q_value < 0.05))
view(pr_deg_ids)
pr_deg_ids_q_0.01 <- row.names(subset(gtest, q_value < 0.01))
pr_deg_ids_q_0.01
view(gtest)
write.csv(gtest)
pr_deg_ids_q_0.0001 <- row.names(subset(gtest, q_value < 0.0001))
gtestOK <- filter(gtest, status == "OK")
#q0.gtest <- filter(gtest, q_value < 4.407045e-307)
gene_module_df <- find_gene_modules(cds_subset2712[pr_deg_ids_q_0.0001,], resolution=c(0,10^seq(-6,-1)))

q0.gtest4 <- filter(gtest, q_value < 4.407045e-307 & module == '4')
view(gtest)
write.csv(q0.gtest4)
view(pr_deg_ids_q_0)
write.table(pr_deg_ids_q_0)
#packageVersion("monocle3")
#cds_subset <- choose_cells(cds_subset)

plot_cells(cds_subset2712, genes=c("CDK1","AURKB","AURKA","NME1","PLK1","CCL3", "CCL4","KIF4A","CCNB1","UBE2C","FAM72D","MND1","MTFR2","POC1A"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#UBE2C correlated gene expression from TCGA - A Comprehensive Bioinformatics Analysis of UBE2C in Cancers
plot_cells(cds_subset2712, genes=c("UBE2C","FAM72D","MND1","MTFR2","POC1A","FOXM1","CCNB1","CCNA1","CCNA2","MKI67","CDK1","AURKB","AURKA"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
#tabula sapiens LN B cell genes
plot_cells(cds_subset2712, genes=c("FOSB","NR4A2","NR4A1","AREG","HSP90AA1","DNAJB1","HSPA8","LY9","CD83","RHOB","HSP90AB1","HSPA1B","SERTAD1","LINC01781","HSPE1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#bb_gene_dotplot(cds_subset2712, gene_or_genes)
#bb_var_umap(cds_subset1245, var = "sample")
#bb_gene_umap(cds_subset1245, gene_or_genes = c("CCL3", "CCL4"))
#bb_gene_umap(cds_main, gene_or_genes = c("CD57", "CD3"))

#nice labeling
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

#dotplot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
             colData(cds_main)$clonotype_id %in% "clonotype1"],
  markers = c("CCL3", "CCL4", "AURKB","AURKA","NME1","CDK1"),
  group_cells_by = "nice_label",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) )) |> pull(gene_short_name)



#####################################################################
#Violin Plots
bb_gene_violinplot(cds_main,
  # filter_cds(
  #   cds_main,
  #   cells = bb_cellmeta(cds_main) |>
  #     filter( == "B")
  # ),
  variable = "kmeans10_cluster",
  genes_to_plot = "S100a8",
  pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
C1<- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                        genes_to_plot = "Elane",
                        pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
C2 <- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                         genes_to_plot = "Ctsg",
                         pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
C3<- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                        genes_to_plot = "Mpo",
                        pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
C4<- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                        genes_to_plot = "Calr",
                        pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)

C5 <- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                   genes_to_plot = "S100a8",
                   pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
C6<- bb_gene_violinplot(cds_main, variable = "kmeans10_cluster",
                        genes_to_plot = "S100a9",
                        pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
)
FC <- (C1|C2|C3)/(C4|C5|C6)
FC

bb_genebubbles(
  obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
  genes = c("Elane",
            "Ctsg",
            "Mpo",
            "Calr",
            "S100a8",
            "S100a9"
  ), cell_grouping = "kmeans10_cluster") + labs(x = "kmeans10 Clusters", y = NULL)

unique(colData(cds_main)$partition)

# subset to include only wt and B ALL cells with old dimensions
olddims_wt_aml <- blaseRtools::filter_cds(cds = cds_main,
                                          cells = bb_cellmeta(cds_main) |>
                                            filter(leukemia_phenotype %in% c("AML", "No leukemia")))

bb_gene_umap(cds_main, gene_or_genes = c("Cd34"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+ facet_wrap(~leukemia_phenotype)+labs(title = "Cd34", y = "UMAP2")+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank())+
  scale_color_distiller(palette = "Oranges",
                        direction = 1,
                        na.value = "grey80", limits = c(0,2.5))
#kmeans10 assignment
bb_genebubbles(
  obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main)),
  genes = c("Cd19",
            "Cd79a",
            "Cd24",
            ""
  ), cell_grouping = "kmeans10_cluster") + labs(x = "kmeans10 Clusters", y = NULL)

#Pseudotime Code - Ethan

#for use in group_cells_by = "nice_label"
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

#generate cds subset by pt & B cells (via clonotype_id)
cds_subset2712 <- cds_main[, colData(cds_main)$patient == "pt_2712" &
                             colData(cds_main)$clonotype_id %in% "clonotype1"]

cds_subset1245 <- cds_main[, colData(cds_main)$patient == "pt_1245" &
                             colData(cds_main)$clonotype_id %in% "clonotype1"]

#bb_gene_pseudotime(order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP"))))
#blaseRtools::bb_gene_pseudotime

cds_subset2712<- order_cells(learn_graph(cluster_cells(cds_subset2712, reduction_method = "UMAP")))
colData(cds_subset2712)

cds_subset2712 <- cluster_cells(cds_subset2712)
cds_subset2712 <- learn_graph(cds_subset2712)
colData(cds_subset2712)
plot_cells(cds_subset2712,
           color_cells_by = "sample",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds_subset2712 <- order_cells(cds_subset2712)

#manually selected root nodes
plot_cells(cds_subset2712,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

monocle3::graph_test
#use trace and substitute Matrix::rBind with rbind
trace('calculateLW', edit = T, where = asNamespace("monocle3"))

gtest <- monocle3::graph_test(cds_subset2712, neighbor_graph="principal_graph", cores=4)

pr_deg_ids <- row.names(subset(gtest, q_value < 0.05))
view(pr_deg_ids)
pr_deg_ids_q_0.01 <- row.names(subset(gtest, q_value < 0.01))
pr_deg_ids_q_0.01
view(gtest)
write.csv(gtest)
pr_deg_ids_q_0.0001 <- row.names(subset(gtest, q_value < 0.0001))
gtestOK <- filter(gtest, status == "OK")
#q0.gtest <- filter(gtest, q_value < 4.407045e-307)
gene_module_df <- find_gene_modules(cds_subset2712[pr_deg_ids_q_0.0001,], resolution=c(0,10^seq(-6,-1)))

q0.gtest4 <- filter(gtest, q_value < 4.407045e-307 & module == '4')
view(gtest)
write.csv(q0.gtest4)
view(pr_deg_ids_q_0)
write.table(pr_deg_ids_q_0)
#packageVersion("monocle3")
#cds_subset <- choose_cells(cds_subset)

plot_cells(cds_subset2712, genes=c("CDK1","AURKB","AURKA","NME1","PLK1","CCL3", "CCL4","KIF4A","CCNB1","UBE2C","FAM72D","MND1","MTFR2","POC1A"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#UBE2C correlated gene expression from TCGA - A Comprehensive Bioinformatics Analysis of UBE2C in Cancers
plot_cells(cds_subset2712, genes=c("UBE2C","FAM72D","MND1","MTFR2","POC1A","FOXM1","CCNB1","CCNA1","CCNA2","MKI67","CDK1","AURKB","AURKA"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
#tabula sapiens LN B cell genes
plot_cells(cds_subset2712, genes=c("FOSB","NR4A2","NR4A1","AREG","HSP90AA1","DNAJB1","HSPA8","LY9","CD83","RHOB","HSP90AB1","HSPA1B","SERTAD1","LINC01781","HSPE1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#bb_gene_dotplot(cds_subset2712, gene_or_genes)
#bb_var_umap(cds_subset1245, var = "sample")
#bb_gene_umap(cds_subset1245, gene_or_genes = c("CCL3", "CCL4"))
#bb_gene_umap(cds_main, gene_or_genes = c("CD57", "CD3"))

#nice labeling
colData(cds_main)$nice_label <-
  recode(
    colData(cds_main)$sample,
    "L34_19972712RTPBMC" = "RT PBMC",
    "L33_19972712RTLN" = "RT LN",
    "L35_19972712CLLPBMC" = "CLL PBMC"
  )

#dotplot
bb_gene_dotplot(
  cds_main[, colData(cds_main)$patient == "pt_2712" &
             colData(cds_main)$clonotype_id %in% "clonotype1"],
  markers = c("CCL3", "CCL4", "AURKB","AURKA","NME1","CDK1"),
  group_cells_by = "nice_label",
  group_ordering = c("CLL PBMC", "RT PBMC", "RT LN"),
  colorscale_name = "Expression",
  sizescale_name = "Proportion\nExpressing",
) + labs(x = NULL, y = NULL)


