#cell type assignment: ScType package
#install.packages("HGNChelper")
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
) # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

#gs_list[["gs_positive"]][["Myeloid Dendritic cells"]]

#unique(colData(cds_main)$partition)
#Expression Matrix - partition
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
      select(cell_id, partition)#barcode)
  ) |>
  t() |>
  scale() |>
  t()

rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(cds_main) |>
              select(feature_id, gene_short_name)) |> pull(gene_short_name)

# #count matrix
# mat <- monocle3::exprs(cds_main)

# assign cell types
es.max = sctype_score(scRNAseqData = mat, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
rownames(es.max)

es.max2 <- as.data.frame(es.max)
es.max2$Largest_Column <-colnames(es.max2)[apply(es.max2,1,which.max)]

colData(cds_main)$partition_assignment2 <-
  recode(colData(cds_main)$partition,
         "1" = "Pre/Pro-B cells", #close rankings
         "2" = "Naive/Effector CD4+ T cells", #close rankings
         "3" = "Blast-like Neu",
         "4" = "Effector CD8+ T cells",
         "5" = "Myeloid Dendritic cells",
         "6" = "Blast-like", #didnt change but annotation says platelets?
         "7" = "Neutrophils",
         "8" = "Naive CD4+ T cells",
         "9" = "Macrophages",
         "10" = "Naive B cells",
         "11" = "Neutrophils",
         "12" = "Erythroid-like and erythroid precursor cells",
         "13" = "Natural killer  cells",
         "14" = "Granulocytes",
         "15" = "Effector CD8+ T cells",
         "16" = "Pro-B cells" #pre? less specific than 16
  )

harmony_colfun <- circlize::colorRamp2(breaks = c(0, 1), colors = c("grey80", "red"))
heatmap_3_colors <-
  c("#313695", "white", "#A50026")

# colfun = circlize::colorRamp2(breaks = c(min(es.max),
#                                          0,
#                                          max(es.max)),
#                               colors = heatmap_3_colors)
# ComplexHeatmap::Heatmap(es.max,
#                         col = colfun)

ComplexHeatmap::Heatmap(es.max,
                        col = colfun)

bb_var_umap(cds_main,
               var = "partition",
               #alt_dim_x = "aggr_UMAP_1",
               #alt_dim_y = "aggr_UMAP_2",
               cell_size = 0.1,
               overwrite_labels = T)
bb_var_umap(cds_main,
            var = "partition_assignment2",
            #alt_dim_x = "aggr_UMAP_1",
            #alt_dim_y = "aggr_UMAP_2",
            cell_size = 0.1,
            overwrite_labels = T) + facet_wrap(~geno_pheno)
  #facet_grid(rows = vars(tissue), cols = vars(geno_pheno))

bb_var_umap(cds_main,
            var = "partition_assignment",
            #alt_dim_x = "aggr_UMAP_1",
            #alt_dim_y = "aggr_UMAP_2",
            cell_size = 0.1,
            overwrite_labels = T)

#leiden
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
      select(cell_id, leiden)#barcode)
  ) |>
  t() |>
  scale() |>
  t()

rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(cds_main) |>
              select(feature_id, gene_short_name)) |> pull(gene_short_name)

# #count matrix
# mat <- monocle3::exprs(cds_main)

# assign cell types
es.max = sctype_score(scRNAseqData = mat, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
rownames(es.max)

es.max2 <- as.data.frame(es.max)
es.max2 <- es.max2[,39:43]
es.max2$Largest_Column <-colnames(es.max2)[apply(es.max2,1,which.max)]
view(es.max2)

colData(cds_main)$leiden_assignment2 <-
  recode(colData(cds_main)$leiden,
         "1" = "Pro-B cells", #4.0
         "2" = "Pre-B cells", #3.4
         "3" = "Naive CD4+ T cells", #2.9
         "4" = "Neutrophils", #1.8
         "5" = "Neutrophils", #1.8
         "6" = "Effector CD8+ T cells", #2.3
         "7" = "Immature B cells", #4.2
         "8" = "Neutrophils", #2.65
         "9" = "Myeloid Dendritic cells", #2.4
         "10" = "Effector CD4+ T cells", #2.3
         "11" = "Myeloid Dendritic cells", #2.2
         "12" = "Granulocytes", #2.0
         "13" = "Neutrophils",#4.4
         "14" = "Neutrophils",#5.1
         "15" = "Naive CD4+ T cells",#6.0
         "16" = "Neutrophils",#3.5
         "17" = "Macrophages",#2.2
         "18" = "Neutrophils",#4.8
         "19" = "Naive CD4+ T cells",#5
         "20" = "Neutrophils",#4.2
         "21" = "Naive B cells",#3.7
         "22" = "Macrophages",#5.6
         "23" = "Naive CD4+ T cells",#2.4
         "24" = "ISG expressing immune cells",#2.9
         "25" = "Neutrophils", #0.4 ???????
         "26" = "Erythroid-like and erythroid precursor cells", #5.5
         "27" = "Macrophages", #4.7
         "28" = "Pro-B cells", #4.4
         "29" = "Effector CD8+ T cells", #3.0
         "30" = "γδ-T cells", #1.7
         "31" = "ISG expressing immune cells",#5.1
         "32" = "Natural killer  cells",#12
         "33" = "Platelets", #9.0
         "34" = "Immature B cells",#2.9
         "35" = "Naive B cells",#7.5
         "36" = "Platelets", #4.1
         "37" = "Neutrophils",#2.5
         "38" = "Megakaryocyte",#4.7
         "39" = "Immature B cells",#4.7
         "40" = "Effector CD8+ T cells",#0.46
         "41" = "Myeloid Dendritic cells",#5.0
         "42" = "CD8+ NKT-like cells",#1.9
         "43" = "Pro-B cells"#4.7
      )

colfun = circlize::colorRamp2(breaks = c(min(es.max),
                                         0,
                                         max(es.max)),
                              colors = heatmap_3_colors)
heatmap_3_colors <-
  c("#313695", "white", "#A50026")
ComplexHeatmap::Heatmap(es.max,
                        col = colfun)
a<-bb_var_umap(cds_main,
            var = "leiden_assignment2",
            #alt_dim_x = "aggr_UMAP_1",
            #alt_dim_y = "aggr_UMAP_2",
            cell_size = 0.1,
            overwrite_labels = T) + facet_wrap(~geno_pheno)
b<-bb_var_umap(cds_main,
            var = "leiden",
            #alt_dim_x = "aggr_UMAP_1",
            #alt_dim_y = "aggr_UMAP_2",
            cell_size = 0.1,
            overwrite_labels = T) + facet_wrap(~geno_pheno)
a/b
bb_var_umap(cds_main,
            var = "partition_assignment",
            #alt_dim_x = "aggr_UMAP_1",
            #alt_dim_y = "aggr_UMAP_2",
            cell_size = 0.1,
            overwrite_labels = T)

colData(cds_main)$leiden_assignment3 <-
  recode(colData(cds_main)$leiden,
         "1" = "Pro-B cells", #4.0
         "2" = "Pre-B cells", #3.4
         "3" = "Naive CD4+ T cells", #2.9
         "4" = "AML Blast-like 1", #Neutrophils", #1.8
         "5" = "AML Blast-like 2", #Neutrophils", #1.8
         "6" = "Effector CD8+ T cells", #2.3
         "7" = "Immature B cells", #4.2
         "8" = "AML Blast-like 3", #Neutrophil-like", #2.65
         "9" = "Myeloid Dendritic cells", #2.4
         "10" = "Effector CD4+ T cells", #2.3
         "11" = "Myeloid Dendritic cells", #2.2
         "12" = "AML Blast-like 4",# Granulocytes", #2.0
         "13" = "Neutrophils",#4.4
         "14" = "Neutrophils",#5.1
         "15" = "Naive CD4+ T cells",#6.0
         "16" = "Neutrophils",#3.5
         "17" = "Macrophages",#2.2
         "18" = "Neutrophils",#4.8
         "19" = "Naive CD4+ T cells",#5
         "20" = "Neutrophils",#4.2
         "21" = "Naive B cells",#3.7
         "22" = "Macrophages",#5.6
         "23" = "Naive CD4+ T cells",#2.4
         "24" = "AML Blast-like 5", #"ISG expressing immune cells",#2.9
         "25" = "Neutrophils", #0.4 ???????
         "26" = "Erythroid-like and erythroid precursor cells", #5.5
         "27" = "Macrophages", #4.7
         "28" = "Pro-B cells", #4.4
         "29" = "Effector CD8+ T cells", #3.0
         "30" = "γδ-T cells", #1.7
         "31" = "ISG expressing immune cells",#5.1
         "32" = "Natural killer  cells",#12
         "33" = "Platelets", #9.0
         "34" = "Immature B cells",#2.9
         "35" = "Naive B cells",#7.5
         "36" = "Platelets", #4.1
         "37" = "Neutrophils",#2.5
         "38" = "Megakaryocyte",#4.7
         "39" = "Immature B cells",#4.7
         "40" = "Effector CD8+ T cells",#0.46
         "41" = "Myeloid Dendritic cells",#5.0
         "42" = "CD8+ NKT-like cells",#1.9
         "43" = "Pro-B cells"#4.7
  )

bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "No leukemia"))),
               var = "leiden_assignment3",
               #alt_dim_x = "aggr_UMAP_1",
               #alt_dim_y = "aggr_UMAP_2",
               cell_size = 0.1,
               overwrite_labels = F) + facet_wrap(~geno_pheno)
bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "No leukemia"))),
  var = "leiden",
  #alt_dim_x = "aggr_UMAP_1",
  #alt_dim_y = "aggr_UMAP_2",
  cell_size = 0.1,
  overwrite_labels = T) + facet_wrap(~geno_pheno)

bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "No leukemia"))),
  var = "density",
  #alt_dim_x = "aggr_UMAP_1",
  #alt_dim_y = "aggr_UMAP_2",
  cell_size = 0.1,
  overwrite_labels = T, facet_by = geno_pheno) + facet_wrap(~geno_pheno)


#Pu's requested gene expression
unique(colData(cds_main)$leukemia_phenotype)
bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("pre-B ALL"))), gene_or_genes = c("Cd79a", "Pax5", "Cd19"))#+facet_wrap(~geno_pheno)

bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "No leukemia"))), gene_or_genes = c("Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5"))+facet_wrap(~geno_pheno)
# make the gene expression umaps
aml_plotlist2 <- map(.x = c("Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5"),
                    .f = \(x, dat = filter_cds(
                      cds_main,
                      cells = bb_cellmeta(cds_main) |>
                        filter(
                          leukemia_phenotype %in% c("AML", "No leukemia")))) {
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

aml_gexp_umap2 <- ggarrange(aml_plotlist2[[1]],
          aml_plotlist2[[2]],
          aml_plotlist2[[3]],
          aml_plotlist2[[4]],
          aml_plotlist2[[5]],
          aml_plotlist2[[6]],
          aml_plotlist2[[7]],
          aml_plotlist2[[8]],
          aml_plotlist2[[9]],
          aml_plotlist2[[10]],
          aml_plotlist2[[11]],
          ncol = 3,
          nrow=4,
          common.legend = TRUE,
          legend="right")
aml_gexp_umap2
#Make Heatmap of AML leiden clusters
