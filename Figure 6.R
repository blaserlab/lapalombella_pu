T_Figs <- "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Figures"

source("R/dependencies.R")
source("R/configs.R")
source("R/cds_mods.R")
##################################################################################################################
#Violin Plots
# bb_gene_violinplot(cds_main, variable = "partition",
#                    genes_to_plot = c("S100a8", "S100a9"),
#                    pseudocount = 0, jitter_fill = "transparent", violin_alpha = 0.55, jitter_alpha = 0.1, include_jitter = TRUE
# ) + theme(axis.title = element_text(size = 12)) +
#   theme(axis.text = element_text(size = 12))
##################################################################################################################
#Figure 6A
F6A1 <- bb_var_umap(
  cds_main,
  var = "partition",
  overwrite_labels = T
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")
F6A1
F6A2<-bb_var_umap(
  cds_main,
  var = "density",
  facet_by = "geno_pheno"
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2")+ theme_minimal() + theme(panel.grid.major = element_blank(),
                                                                                                   panel.grid.minor = element_blank()) +theme(panel.background = element_rect(color = "black"))
F6A2
#####################################################################################################################
#Figure 6B
aml_plotlist <- map(.x = c("Cd34","Mpo", "Kit", "Elane", "Calr","Ctsg"),
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        #alt_dim_x = "aggr_UMAP_1",
                        #alt_dim_y = "aggr_UMAP_2",
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80",
                                              limits = c(0,2.5)) +
                        facet_wrap(~geno_pheno, labeller = labeller(group = label_wrap_gen(width = 5, multi_line = TRUE))) +
                        scale_y_continuous(breaks = c(-10, 0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, 10))+
                        theme(panel.spacing = unit(0.25, "lines"))+
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black",
                                                              fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        #labs(x =NULL, y = NULL, subtitle = x) +
                        #theme(plot.subtitle = element_text(face = "italic", hjust = 0.5))+
                        theme(legend.position = "none")+
                        coord_fixed(ratio = 1)+ labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
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

F6B1 <-
  (aml_plotlist[[1]] |
     plot_spacer() |
     aml_plotlist[[2]] |
     plot_spacer() |
     aml_plotlist[[3]])+ theme(
       legend.position = "right",
       legend.key.size = unit(4, "mm")
     ) + plot_layout(widths = c(1, 0, 1, 0, 1))
#ggsave("F6C1.pdf", path = T_Figs, width = 8.25, height = 4.5)

F6B2 <-
  (aml_plotlist[[4]] |
     plot_spacer() |
     aml_plotlist[[5]] |
     plot_spacer() |
     aml_plotlist[[6]])+ theme(
       legend.position = "right",
       legend.key.size = unit(4, "mm")
     ) + plot_layout(widths = c(1, 0, 1, 0, 1))
#ggsave("F6C2.pdf", path = T_Figs, width = 8.25, height = 4.5)
########################################################################################################################
#Figure 6C
#Heatmap
##top markers
F6_topmarkers_part <-
  monocle3::top_markers(
    cds_main,
    group_cells_by = "partition",
    genes_to_test_per_group = 20,
    cores = 12
  )
#write_csv(F6_topmarkers_part, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data", "F6_topmarkers_part.csv"))
F6_topmarkers_part <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F6_topmarkers_part.csv")


# markers <- F5_topmarkers_k10 |> filter(
#    cell_group %in% c('6','8','1','5','3'))|>pull(gene_short_name)
# markers <- F5_topmarkers_k10 |> filter(
#   cell_group %in% c('6','8','1','5','3','2','4','7','9','10'))|>pull(gene_short_name)
markers <- F6_topmarkers_part |> pull(gene_short_name)

#Expression Matrix
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          partition %in% c(1:16)
        ),
      genes = bb_rowmeta(cds_main) |>
        filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, partition)
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

#Annotation: Top 5 marker genes per cluster
heatmap_highlights <-
  F6_topmarkers_k10 |> group_by(cell_group) |>
  slice_max(order_by = marker_score, n=5)|>
  pull(gene_short_name)

anno <-
  ComplexHeatmap::rowAnnotation(link =  anno_mark(
    at = which(rownames(mat) %in% heatmap_highlights),
    labels = rownames(mat)[rownames(mat) %in% heatmap_highlights],
    labels_gp = gpar(fontsize = 7),
    padding = 0.15
  ))

#F6heatmap<-grid.grabExpr(draw(
ComplexHeatmap::Heatmap(
  mat,
  col = colfun,
  name = "Expression",
  show_row_names = F,
  show_column_names = F, #check column clustering order for bp
  right_annotation = anno,
  #top_annotation = hmap_bp,
  #width = ncol(mat)*unit(0.1, "mm"),
  height = nrow(mat)*unit(0.8, "mm"),
  row_dend_width = unit(6, "mm"),
  column_dend_height = unit(3, "mm"),
  heatmap_legend_param = list(
    legend_direction = "vertical",
    #legend_height = unit(1, "cm"),
    legend_width = unit(0.25, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 8.5)
  )
)#))

#stacked bar chart
####fraction of cells contributed to each cluster by leukemia_phenotype
# cellcount<- bb_cellmeta(cds_main) |>
#   group_by(kmeans10_cluster, leukemia_phenotype) |>
#   summarise(n = n()) |> filter(kmeans10_cluster %in% c("5", "6", "1", "8","3")) #select clusters
cellcount<- bb_cellmeta(cds_main) |> #all clusters
  group_by(partition, leukemia_phenotype) |>
  summarise(n = n()) |> filter(partition %in% c(1:16))|> mutate(leukemia_phenotype = recode(leukemia_phenotype,
                                                                                            "No leukemia" = "WT"))
#create fraction column
library(data.table)
setDT(cellcount)[, frac := n / sum(n), by=partition]

#factor levels
cellcount$partition <- factor(cellcount$partition,
                                     levels = c("12","2","4","15","6","3","9","14", "11","1","16","10","5","8","13","7")) #for all
cellcount$leukemia_phenotype <- factor(cellcount$leukemia_phenotype,
                                       levels = c("AML", "pre-B ALL","T ALL", "WT"))
#plot
hmap_bp <-
  ggplot(cellcount,
         aes(x = partition, y = frac, fill = leukemia_phenotype)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme_minimal()+
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(legend.justification = "top")+
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 9))+
  coord_fixed(ratio = 1.5)+ #0.45 #1.2
  theme(plot.margin = unit(c(-0.1,3.23,0,.23), "cm")) + #0.4 for select clusters
  #theme(plot.margin = unit(c(1,3.31,0,0.13), "cm")) + #0.4 for select clusters
  labs(x = "Partition Clusters")+
  scale_x_discrete(expand = c(0.07,0))+ #(0.13,0) for select clusters
  theme(legend.margin=margin(t = -0.22, unit='cm'))
hmap_bp
#ggsave("hmap_bp.pdf", path = T_Figs, width = 3.56, height = 3.85)

# F6hm<- plot_grid(F6heatmap,hmap_bp, ncol=1, rel_heights = c(3,0.75))
# F6hm

########################################################################################################################
#Figure 6D & 6E

#AML leiden clusters heatmap
aml <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12')))
#Figure 6D
bb_var_umap(aml, "leiden_assignment2", overwrite_labels = T, group_label_size = 5)
ggsave("Figure6D-aml_leiden_clusters.pdf", path = T_Figs)

# AML_topmarkers <-
#   monocle3::top_markers(
#     aml,
#     group_cells_by = "leiden",
#     genes_to_test_per_group = 20,
#     cores = 12)

#write_csv(AML_leiden_clusts_topmarkers, file = file.path("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data", "AML_leiden_clusts_topmarkers.csv"))
AML_topmarkers <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/AML_leiden_clusts_topmarkers.csv")

markers <- AML_topmarkers |> pull(gene_short_name)

#Expression Matrix
mat <-
  bb_aggregate(
    obj = filter_cds(
      aml,
      cells = bb_cellmeta(aml) |>
        filter(
          leiden %in% c('4', '5', '8', '24', '12')
        ),
      genes = bb_rowmeta(aml) |>
        filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(aml) |>
      select(cell_id, leiden)
  ) |>
  t() |>
  scale() |>
  t()
rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(aml) |>
              select(feature_id, gene_short_name)) |>
  pull(gene_short_name
  )

heatmap_3_colors <-
  c("#313695", "white", "#A50026")

colfun = circlize::colorRamp2(breaks = c(min(mat),
                                         0,
                                         max(mat)),
                              colors = heatmap_3_colors)

#Annotation: Top 5 marker genes per cluster
# heatmap_highlights <-
#   AML_topmarkers |> group_by(cell_group) |>
#   slice_max(order_by = marker_score, n=5)|>
#   pull(gene_short_name)
highlights1 <- c("Ccl6","Mpo", "Ctsg", "Elane", "Cd34", "Birc5", "S100A8", "S100A9", "Lyz2", "Cd52", "Stmn1", "Cd34", "Mpo", "Klf4", "Il7r", "Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5")

highlights2 <- c("Ifi27l2a", "S100a8", "S100a9", "Wfdc17", "Trem2", "Lyz2", "Ccl6", "Fcer1g", "Ftl1", "Cd52","Lgals3", "Fn1", "Ifi213","Hbb-bt","Hbb-bs","Hba-a1","Hba-a2","Lcn2", "Hmgn2","Ngp", "Camp", "Ngp","Ltf","Lcn","Ifitm6","Wfdc21","Elane","Mpo","Top2a","Tubb5", "Tubb4b", "Birc5", "Ran", "Stmn1" ,"Cox5b","Atp5c1","Chchd2","Mdh2","Atp5j2","Mif","Npm1","Atp5g1","Uqcrb")
heatmap_highlights <- unique(c(highlights1, highlights2))
#heatmap_highlights <- unique(c(heatmap_highlights, highlights))
anno <-
  ComplexHeatmap::rowAnnotation(link =  anno_mark(
    at = which(rownames(mat) %in% heatmap_highlights),
    labels = rownames(mat)[rownames(mat) %in% heatmap_highlights],
    labels_gp = gpar(fontsize = 6.5),
    padding = 0.15
  ))

#Figure 6E
F6E<-grid.grabExpr(draw(
ComplexHeatmap::Heatmap(
  mat,
  col = colfun,
  name = "Expression",
  show_row_names = F,
  show_column_names = T, #check column clustering order for bp
  right_annotation = anno,
  #top_annotation = hmap_bp,
  #width = ncol(mat)*unit(0.1, "mm"),
  height = nrow(mat)*unit(1.25, "mm"),
  row_dend_width = unit(6, "mm"),
  column_dend_height = unit(3, "mm"),
  heatmap_legend_param = list(
    legend_direction = "vertical",
    #legend_height = unit(1, "cm"),
    legend_width = unit(0.25, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 8.5)
  )
)))
F6E<-as_ggplot(F6E)
F6E
#ggsave("F6E.pdf", path = T_Figs)
##############################################################################################################################################
#Fig 6F Volcano
# exp_design <-
#   bb_cellmeta(filter_cds(
#     cds = cds_main,
#     cells = bb_cellmeta(cds_main) |>
#       filter(leukemia_phenotype %in% c("AML", "No leukemia")))) |> #wt_aml
#   group_by(sample, leukemia_phenotype) |>
#   summarise()
#
# pseudobulk_res <-
#   bb_pseudobulk_mf(cds = filter_cds(
#     cds = cds_main,
#     cells = bb_cellmeta(cds_main) |>
#       filter(leukemia_phenotype %in% c("AML", "No leukemia"))),
#     pseudosample_table = exp_design,
#     design_formula = "~ leukemia_phenotype",
#     result_recipe = c("leukemia_phenotype", "AML", "No leukemia"))
#
# #write.csv(F6_aml_wt_pseudobulk, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F6_aml_wt_pseudobulk.csv")
# F6_aml_wt_pseudobulk<- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/F6_aml_wt_pseudobulk.csv")
#
# genes_to_highlight <- unique(c("Prox1", "Ifi44l", "Tdrd5", "Etv4", "Etv5"))
#
# volcano_data <- F6_aml_wt_pseudobulk %>%
#   mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
#   mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight, gene_short_name, ""))
# #write.csv(volcano_data, "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/volcano_data.csv")
#
# volcano_data<-volcano_data[!(volcano_data$gene_short_name=="Slc4a8"),]
#
# library(ggtext)
# volcano_pseudobulk <-
#   ggplot(
#     volcano_data,
#     aes(
#       x = log2FoldChange,
#       y = -log10(padj),
#       colour = threshold,
#       fill = threshold,
#       label = text_label
#     )
#   ) +
#   geom_point(shape = 21,
#              size = 0.5,
#              alpha = 0.4) +
#   geom_text_repel(color = "black",
#                   fontface = "italic",
#                   box.padding = 0.5, #0.5
#                   point.padding = 0.25, #0.25
#                   min.segment.length = 0,
#                   max.overlaps = 20000,
#                   size = 3,
#                   segment.size = 0.25,
#                   force = 2,
#                   seed = 1234,
#                   segment.curvature = -0.1,
#                   segment.square = TRUE,
#                   segment.inflect = TRUE) +
#   xlab("log<sub>2</sub> fold change") +
#   ylab("-log<sub>10</sub> adjusted p-value") +
#   theme(axis.title.x =  element_markdown()) +
#   theme(axis.title.y = element_markdown()) +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c("grey80", "#DC0000")) +
#   scale_fill_manual(values = c("transparent", "#DC0000")) +
#   labs(title = "Pseudobulk dKO AML vs WT")+ #caption = "\U21D0 Up in WT\nUp in AML \U21D2",
#   theme(plot.caption.position = "panel") +
#   #theme(plot.caption = element_text(hjust = 0.5)) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_cartesian(xlim = c(-1.0 * max(abs(range(volcano_data |> dplyr::filter(!is.na(padj)) |> pull(log2FoldChange)))), 1.0 * max(abs(range(volcano_data |> filter(!is.na(padj)) |> pull(log2FoldChange)))))) #+
# volcano_pseudobulk
##############################################################################################################################################
#Fig 6G
AML_topmarkers <- read.csv("~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Data/AML_leiden_clusts_topmarkers.csv")

library(topGO)
goenrichment <-
  bb_goenrichment(
    query = dplyr::filter(pseudobulk_res$Result, padj < 0.1 &
                            log2FoldChange >= 0.58) |> pull(gene_short_name),
    reference = bb_rowmeta(cds_main),
    go_db = "org.Mm.eg.db"
  )


gosummary_0.89 <- bb_gosummary(x = goenrichment,
                              reduce_threshold = 0.89,
                              go_db = "org.Mm.eg.db")
pseudobulk_AML_vs_WT_GO_PCA <-
  bb_goscatter(simMatrix = gosummary_0.89$simMatrix,
               reducedTerms = gosummary_0.89$reducedTerms)
pseudobulk_AML_vs_WT_GO_PCA

gosummary_0.9 <- bb_gosummary(x = goenrichment,
                                         reduce_threshold = 0.9,
                                         go_db = "org.Mm.eg.db")
pseudobulk_AML_vs_WT_GO_PCA_0.9 <-
  bb_goscatter(simMatrix = gosummary_0.9$simMatrix,
               reducedTerms = gosummary_0.9$reducedTerms)
pseudobulk_AML_vs_WT_GO_PCA_0.9

#pseudoGO barplot
# pseudo3n11_goenrichment$res_table$classicFisher[pseudo3n11_goenrichment$res_table$classicFisher=="< 1e-30"]<-"1.0e-30"
# pseudo3n11_goenrichment$res_table$Rank <- as.numeric(as.character(pseudo3n11_goenrichment$res_table$Rank))
# pseudob_top25 <- filter(pseudo3n11_goenrichment$res_table, as.numeric(pseudo3n11_goenrichment$res_table$Rank) <= 25) |> mutate(neg_log10_pval = -log10(as.numeric(classicFisher))) |> rename(GO_Term = Term, Genes_Mapped = Significant)
# pseudob_top25
#
# pseudob_3n11.Up_GObp<- ggplot(data=pseudob_top25, aes(reorder(x= GO_Term, y= neg_log10_pval, neg_log10_pval), y= neg_log10_pval, fill = Genes_Mapped)) +
#   geom_bar(stat="identity") +
#   coord_flip() + scale_fill_viridis_c() + labs(x = "GO Terms", y = "-log(pval)")#theme(axis.title.x = element_text("GO Terms"), axis.title.y = element_text("-log(pval)"))
# pseudob_3n11vs1n9_GObp
# #ggsave("pseudob_3n11v.Up_GObp.pdf", path = "~/network/T/Labs/EHL/Rosa/Ethan/EHL/PRMT5/Hing et al manuscript - NatComm/10X Project Update/Figs/Composed Figs/Figure1_Supp1/S1_GO_&_Pathway_Analysis/pseudobulk")


##############################################################################################################################################

#cell type assignment: scType package
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
)
source(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
)

# get cell-type-specific gene sets from ScType database
gs_list = gene_sets_prepare(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
  "Immune system"
)

#louvain
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          partition %in% c(1:16)
          #leukemia_phenotype %in% c("AML", "pre-B ALL", "T ALL", "No leukemia")
        ),
      genes = bb_rowmeta(cds_main) #|>
      #   filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, louvain)#barcode)
  ) |>
  t() |>
  scale() |>
  t()

rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(cds_main) |>
              select(feature_id, gene_short_name)) |> pull(gene_short_name)

# assign cell types
es.max = sctype_score(scRNAseqData = mat, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
#rownames(es.max)

#louvain
es.max3<-as.data.frame(es.max|> t())
es.max3 <- as.data.frame(colnames(es.max3)[max.col(es.max3)])
es.max3$louvain <- rownames(es.max3)
colnames(es.max3)[1] <- "louvain_assignment"
