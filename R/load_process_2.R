library("blaseRtools")
library("blaseRtemplates")
library("tidyverse")
library("monocle3")
library("conflicted")
library("Seurat")
library("SeuratObject")

conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("exprs", "monocle3")
conflict_prefer("rename", "dplyr")

# read in the config file -------------------------------------------------
analysis_configs <- read_csv("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/analysis_configs.csv") |>
  mutate(pipestance_path = bb_fix_file_path(pipestance_path)) |>
  mutate(pipestance_path = fs::path(pipestance_path)) |>
  mutate(citeseq_path = fs::path(citeseq_path)) |>
  mutate(sample = as.character(sample))

# identify the cds files -------------------------------
cds_dirs <- fs::dir_ls("~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/cds")

cds_list <- map(.x = cds_dirs,
    .f = \(x) {
      cds_path <- fs::path(x, "cds.rda")
      load(cds_path)
      cds
    }) |>
  set_names(analysis_configs$sample)


cds_main_disk <- combine_cds(cds_list = cds_list,
                        sample_col_name = "check_and_remove")

colData(cds_main_disk)$check_and_remove <- NULL
rowData(cds_main_disk)$data_type <- rowData(cds_main_disk)$experiment_type
rowData(cds_main_disk)$experiment_type <- NULL
cds_main_disk <- bb_split_citeseq(cds_main_disk)

cds_main_disk <- convert_counts_matrix(cds = cds_main_disk, matrix_control=list(matrix_class='BPCells'))
cds_main_disk <- monocle3::estimate_size_factors(cds_main_disk)

cds_main_disk <- preprocess_cds(cds_main_disk)
cds_main_disk <- reduce_dimension(cds_main_disk)
# save_monocle_objects(cds_main_disk, directory_path = "/workspace/brad_workspace/lapalombella.pu.datapkg/inst/extdata/cds_main_disk")
# cds_main_disk <- load_monocle_objects(directory_path = "/workspace/brad_workspace/lapalombella.pu.datapkg/inst/extdata/cds_main_disk")

# bb_cellmeta(cds_main_disk)
# bb_var_umap(cds_main_disk, "uid", foreground_alpha = 0.1)
# bb_gene_umap(cds_main_disk, "CD3E")

# # align ---------------------------------------------------
cds_main_disk<- bb_align(cds_main_disk, align_by = "pool", n_cores = 20)
# bb_var_umap(cds_main_disk, "pool")

# # Identify clusters and calculate top markers ---------------------------
cds_main_disk <-
  bb_triplecluster(
    cds_main_disk,
    n_top_markers = 50,
    outfile = "data/cds_main_disk_top_markers.csv",
    n_cores = 8
  )
cds_main_disk_top_markers <- read_csv("data/cds_main_disk_top_markers.csv")
save_monocle_objects(cds_main_disk, directory_path = "/workspace/brad_workspace/lapalombella.pu.datapkg/inst/extdata/cds_main_disk")

# # Identify gene modules and add them to the gene metadata. ---------------
# not working
# cds_main_disk <- bb_gene_modules(cds_main_disk, n_cores = 4)
# pr_graph_test_res <- graph_test(cds = cds_main_disk, neighbor_graph="knn", cores=4, )
# pr_graph_test_res |> View()
#
# pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#
# save_monocle_objects(cds_main_disk, directory_path = "/workspace/brad_workspace/lapalombella.pu.datapkg/inst/extdata/cds_main_disk")

# # align to seurat reference ---------------------------------------------
# # NB:  human PBMC only
cds_main_disk <-
  bb_seurat_anno(
    cds_main_disk,
    reference = system.file("extdata", "pbmc_multimodal.h5seurat",
                            package = "blaseRextras")
  )
cds_main_disk

bb_seurat_anno

# save the objects -------------------------------------------------------
save(cds_main_disk_top_markers, file = "data/cds_main_top_markers.rda", compress = "bzip2")
save_monocle_objects(cds_main_disk, directory_path = "/workspace/brad_workspace/lapalombella.pu.datapkg/inst/extdata/cds_main_disk")

bb_rowmeta(cds_main_disk)
bb_cellmeta(cds_main_disk)
