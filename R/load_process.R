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

save_dir <- fs::path("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024")
fs::dir_create(fs::path(save_dir, "qc_lists"))

# read in the config file -------------------------------------------------
analysis_configs <- read_csv("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/analysis_configs.csv") |>
  mutate(pipestance_path = bb_fix_file_path(pipestance_path)) |>
  mutate(pipestance_path = fs::path(pipestance_path)) |>
  mutate(citeseq_path = fs::path(citeseq_path)) |>
  mutate(sample = as.character(sample))


# generate sequencing qc table --------------------------------------------
seq_qc <- map_dfr(.x = analysis_configs$sample,
                  .f = \(x, data = analysis_configs) {
                    pipestance <- data |>
                      filter(sample == x) |>
                      pull(pipestance_path)
                    stopifnot(fs::dir_exists(pipestance))
                    read_csv(
                      list.files(
                        path = pipestance,
                        pattern = "metrics_summary.csv",
                        recursive = TRUE,
                        full.names = T
                      )
                    ) |>
                      mutate(sample = x) |>
                      relocate(sample)

                  })

# read in the hashtagging data -------------------------------------
hto_file <- fs::dir_ls("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/hto_lookup/")

# get the sample names directly from the hto files
sample <- str_remove(fs::path_file(hto_file), ".csv")

# join hto data onto analysis configs
analysis_configs <- left_join(analysis_configs, tibble(hto_file, sample), by = join_by("sample"))

hto_list <- map2(
  .x = analysis_configs$hto_file,
  .y = analysis_configs$sample,
  .f  = \(x,y) {
  read_csv(x) |>
    mutate(sample = y) |>
    mutate(cell_id = paste0(barcode_sequence, "-1")) |>
    select(cell_id, hash_id, sample)
}) |>
  set_names(analysis_configs$sample)

# Generate a list of CDS objects using purrr::map ------------------------
walk2(
  .x = analysis_configs$sample,
  .y = hto_list,
  .f = \(x, y, conf = analysis_configs, sd = save_dir) {
    conf_filtered <- conf |>
      filter(sample == x)
    h5 <- list.files(
      conf_filtered$pipestance_path,
      pattern = "filtered_feature_bc_matrix.h5",
      recursive = T,
      full.names = T
    )
    cds <- bb_load_tenx_h5(
      filename = h5,
      sample_metadata_tbl = conf_filtered |>
        select(-c(pipestance_path))
    ) |>
      bb_tbl_to_coldata(min_tbl = y)

    # filter out the doublets and cells without hashtags
    cds <- filter_cds(cds, cells = bb_cellmeta(cds) |>
      filter(hash_id != "Doublet") |>
      filter(!is.na(hash_id)))

    # standard QC
    qc_res <- bb_qc(
      cds,
      cds_name = x,
      genome = "human"
    )

    # identify doublets
    anticipated_doublet_rate <- ncol(cds) / 100000

    doubletfinder_res <- bb_doubletfinder(
      cds,
      doublet_prediction = anticipated_doublet_rate,
      qc_table = qc_res[[1]])


    cellmeta_list <- list(hto = y, qc = qc_res, df = doubletfinder_res)

    save_dir <- fs::path(sd, "qc_lists", x)
    fs::dir_create(save_dir)
    save(cellmeta_list,
         file = fs::path(save_dir, "cellmeta_list.rda"),
         compress = "bzip2")
    gc()
  }
)

save(seq_qc, file = fs::path(save_dir, "seq_qc.rda"), compress = "bzip2")
save(analysis_configs, file = fs::path(save_dir, "analysis_configs.rda"), compress = "bzip2")
