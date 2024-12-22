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

# read in the config file -------------------------------------------------
analysis_configs <- read_csv("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/analysis_configs.csv") |>
  mutate(pipestance_path = bb_fix_file_path(pipestance_path)) |>
  mutate(pipestance_path = fs::path(pipestance_path)) |>
  mutate(citeseq_path = fs::path(citeseq_path)) |>
  mutate(sample = as.character(sample))


# get the citeseq barcode list ------------------------


biolegend_adts <-
  read_csv(
    "~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/TotalSeq_C_Human_Universal_Cocktail_v2_184_Antibodies_399910_Barcodes.csv",
    skip = 1,
    col_names = c(
      "short_id",
      "description",
      "clone",
      "barcode",
      "ensembl_id",
      "new_or_old"
    )
  ) |>
  mutate(id = paste0(short_id, "-", barcode)) |>
  mutate(
    gene_short_name = str_remove(
      description,
      "anti-human |anti-human/mouse |anti-mouse/human |anti-human/mouse/rat "
    )
  ) |>
  mutate(gene_short_name = str_replace(gene_short_name, "  ", " "))


# get a list of ensembl gene ids and gene short names ----------------------
ensembl <-
  biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl_lookup <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               mart = ensembl) |>
  as_tibble() |>
  rename(id = ensembl_gene_id, gene_short_name = hgnc_symbol) |>
  filter(!is.na(gene_short_name)) |>
  filter(gene_short_name != "")



# read in the data -------------------------
walk(.x = analysis_configs$sample,
     .f = \(x,
            conf = analysis_configs,
            ens_lookup = ensembl_lookup,
            bl_adts = biolegend_adts,
            sd = save_dir) {
       conf_filtered <- conf |>
         filter(sample == x)
       mat_gex <-
         Matrix::readMM(
           file = fs::path(
             conf_filtered$pipestance_path,
             "filtered_feature_bc_matrix/matrix.mtx.gz"
           )
         )
       mat_adt <-
         Matrix::readMM(file = fs::path(conf_filtered$citeseq_path,
                                        "matrix.mtx.gz"))
       barcodes_gex <-
         readr::read_tsv(
           fs::path(
             conf_filtered$pipestance_path,
             "filtered_feature_bc_matrix/barcodes.tsv.gz"
           ),
           col_names = "cell_id"
         )
       barcodes_adt <-
         readr::read_tsv(
           file = fs::path(conf_filtered$citeseq_path,
                           "barcodes.tsv.gz"),
           col_names = "cell_id"
         ) |>
         dplyr::mutate(cell_id = paste0(cell_id, "-1"))
       features_gex <-
         readr::read_tsv(
           file = fs::path(
             conf_filtered$pipestance_path,
             "filtered_feature_bc_matrix/features.tsv.gz"
           ),
           col_names = "feature_id"
         )
       features_adt <-
         readr::read_tsv(fs::path(conf_filtered$citeseq_path,
                                  "features.tsv.gz"),
                         col_names = "feature_id")

       colnames(mat_adt) <- barcodes_adt$cell_id
       rownames(mat_adt) <- features_adt$feature_id

       colnames(mat_gex) <- barcodes_gex$cell_id
       rownames(mat_gex) <- features_gex$feature_id

       cb_intersection <-
         base::intersect(barcodes_gex$cell_id, barcodes_adt$cell_id)

       genes_to_remove <- ens_lookup |>
         filter(gene_short_name %in% blaseRdata::hg38_remove_genes) |>
         pull(id)

       genes_to_keep <-
         rownames(mat_gex)[!rownames(mat_gex) %in% genes_to_remove]

       stopifnot(sum(rownames(mat_gex) %notin% genes_to_remove) == length(genes_to_keep))
       stopifnot(length(genes_to_keep) == length(unique(genes_to_keep)))

       mat_use <- rbind(mat_gex[genes_to_keep, cb_intersection],
                        mat_adt[, cb_intersection])

       stopifnot(length(rownames(mat_use)) == length(unique(rownames(mat_use))))

       # make the basic cellmeta ---------------------
       cellmeta_use <-
         data.frame(barcode = colnames(mat_use),
                    row.names = colnames(mat_use))

       # make the rowmeta -------------------
       full_lookup <-
         bind_rows(
           ens_lookup |>
             mutate(experiment_type = "Gene Expression"),
           bl_adts |>
             select(id, gene_short_name) |>
             mutate(experiment_type = "Antibody Capture")
         )

       dupes <- full_lookup |>
         count(id) |>
         filter(n >= 2) |>
         pull(id)

       cleaned_lookup <- full_lookup |>
         filter(id %notin% dupes)
       {
         rowmeta_use <- left_join(tibble(id = rownames(mat_use)),
                                  cleaned_lookup, by = join_by(id)) |>
           as.data.frame()
         rownames(rowmeta_use) <- rowmeta_use$id
         }

       stopifnot(nrow(rowmeta_use) == length(rownames(mat_use)))

       cds <- monocle3::new_cell_data_set(
         expression_data = mat_use,
         cell_metadata = cellmeta_use,
         gene_metadata = rowmeta_use
       )

       # identify and load the correct cellmeta list ----------

       cellmeta_list_dir <-
         fs::path(
           "~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/qc_lists/"
         )
       cellmeta_list_dir <- fs::dir_ls(cellmeta_list_dir, regexp = x)
       load(fs::path(cellmeta_list_dir, "cellmeta_list.rda"))

       # make the cds --------------------
       cds <- bb_tbl_to_coldata(
         cds,
         min_tbl = cellmeta_list$hto |>
           mutate(uid = paste0(sample, "_", hash_id)) |>
           select(cell_id, sample, hash_id, uid) |>
           left_join(analysis_configs, by = join_by("sample")) |>
           select(-c(pipestance_path, citeseq_path))
       )
       cds <-
         bb_tbl_to_coldata(cds, min_tbl = cellmeta_list$qc$qc_calls, join_col = "barcode")
       cds <-
         bb_tbl_to_coldata(cds, min_tbl = cellmeta_list$df, join_col = "barcode")

       cds <- filter_cds(
         cds,
         cells = bb_cellmeta(cds) |>
           filter(
             hash_id != "Doublet",
             !is.na(hash_id),
             qc.any == FALSE,
             doubletfinder_high_conf == "Singlet"
           ),
         genes = bb_rowmeta(cds) |>
           filter(!is.na(gene_short_name),
                  !is.na(experiment_type))
       )

       # Remove the qc and doubletfinder columns from the cell metadata --------
       colData(cds)$qc.any <- NULL
       colData(cds)$doubletfinder_low_conf <- NULL
       colData(cds)$doubletfinder_high_conf <- NULL

       # save the cds
       save_dir <- fs::path(sd, "cds", x)
       fs::dir_create(save_dir)
       save(
         cds,
         file = fs::path(save_dir, "cds.rda")
       )
       gc()


     })

save(analysis_configs, file = fs::path(save_dir, "analyis_configs.rda"), compress = "bzip2")
save(biolegend_adts, file = fs::path(save_dir, "biolegend_adts.rda"), compress = "bzip2")
gc()
