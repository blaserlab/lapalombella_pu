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
#Expression Matrix
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

#Refine celltype annotation:

mat2 <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          partition %in% c(1, 7, 9, 5, 11)
          #leukemia_phenotype %in% c("AML", "pre-B ALL", "T ALL", "No leukemia")
        ),
      genes = bb_rowmeta(cds_main) #|>
      #   filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, leiden)
  ) |>
  t() |>
  scale() |>
  t()
