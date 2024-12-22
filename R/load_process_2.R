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

save_dir <- fs::path("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024")

project_data("~/network/X/Labs/Blaser/share/resources/datapkg/blaseRextras/")
# read in the config file -------------------------------------------------
analysis_configs <- read_csv("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/analysis_configs.csv") |>
  mutate(pipestance_path = bb_fix_file_path(pipestance_path)) |>
  mutate(pipestance_path = fs::path(pipestance_path)) |>
  mutate(citeseq_path = fs::path(citeseq_path)) |>
  mutate(sample = as.character(sample))

# read in the updated patient sample metadata and clean it -----------------------
ps_meta <- read_csv("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/patient_sample_metadata.csv",
  col_names = c("pid", "hash_id", "hashtag_barcode", "sample", "PZ", "date")
) |>
  select(-c(hashtag_barcode, date)) |>
  mutate(sample = recode(sample,
    "COMUTANT-pool1" = "comutant-pool1",
    "COMUTANT-pool2" = "comutant-pool2",
    "COMUTANT-pool3" = "comutant-pool3",
    "COMUTANT-pool4" = "comutant-pool4",
    "TET2-pool1" = "tet2-pool1",
    "TET2-pool2" = "tet2-pool2",
    "TP53-pool1" = "tp53-pool1",
    "TP53-pool2" = "tp53-pool2"
  )) |>
  mutate(sample = str_replace(sample, "-", "_")) |>
  left_join(analysis_configs, by = join_by("sample", "PZ")) |>
  mutate(uid = paste0(sample, "_", hash_id)) |>
  select(uid, pid)



# identify the cds files -------------------------------
cds_dirs <- fs::dir_ls("~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/cds")


cds_list <- map(
  .x = cds_dirs,
  .f = \(x) {
    cds_path <- fs::path(x, "cds.rda")
    load(cds_path)
    cds
  }
) |>
  set_names(fs::path_file(cds_dirs))


cds_main_human <- combine_cds(
  cds_list = cds_list,
  sample_col_name = "check_and_remove"
)

waldo::compare(
  bb_cellmeta(cds_main_human)$sample,
  bb_cellmeta(cds_main_human)$check_and_remove
)
rm(cds_list)
gc()

colData(cds_main_human)$check_and_remove <- NULL
rowData(cds_main_human)$data_type <- rowData(cds_main_human)$experiment_type
rowData(cds_main_human)$experiment_type <- NULL
cds_main_human <- bb_split_citeseq(cds_main_human)

# cds_main_human <- convert_counts_matrix(cds = cds_main_human, matrix_control=list(matrix_class='BPCells'))
cds_main_human <- monocle3::estimate_size_factors(cds_main_human)

cds_main_human <- preprocess_cds(cds_main_human)
cds_main_human <- reduce_dimension(cds_main_human,
  build_nn_index = TRUE
)

# add in additional cell metadata --------------------------
new_cellmeta <- left_join(bb_cellmeta(cds_main_human), ps_meta, by = "uid") |>
  select(cell_id, pid)
cds_main_human <- bb_tbl_to_coldata(obj = cds_main_human, min_tbl = new_cellmeta)

# # align ---------------------------------------------------
cds_main_human <- bb_align(cds_main_human, align_by = "PZ", n_cores = 20)

# # Identify clusters and calculate top markers ---------------------------
cds_main_human <-
  bb_triplecluster(
    cds_main_human,
    n_top_markers = 50,
    outfile = fs::path(save_dir, "cds_main_human_top_markers.csv"),
    n_cores = 8
  )
cds_main_human_top_markers <- read_csv(fs::path(save_dir, "cds_main_human_top_markers.csv"))

# reference mapping

cds_main_human <- bb_monocle_anno(
  cds_qry = cds_main_human,
  cds_ref = pbmc_ref,
  labels = c("celltype.l1", "celltype.l2", "celltype.l3")
)

# save the objects -------------------------------------------------------
save_monocle_disk(cds_main_human, data_directory = "../lapalombella.pu.datapkg/data/", extdata_directory = "../lapalombella.pu.datapkg/inst/extdata/")
save(cds_main_human_top_markers, file = "../lapalombella.pu.datapkg/data/cds_main_human_top_markers.rda", compress = "bzip2")
