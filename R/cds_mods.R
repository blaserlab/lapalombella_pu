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
#unique(colData(cds_main)$geno_pheno)

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

#where did these assignments come from?
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

#Assignment from ScType Package
colData(cds_main)$partition_assignment2 <-
  recode(colData(cds_main)$partition,
         "1" = "Pre/Pro-B cells", #close rankings
         "2" = "Naive/Effector CD4+ T cells", #close rankings
         "3" = "Neutrophils",
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
##############################################################################################
