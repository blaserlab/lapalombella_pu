# # make a directory in sc_working to collect all the processed data from network
#
# processed_destination <- fs::path("~/brad_workspace/sc_working/lapalombella_pu_geo/processed")
# fs::dir_create(processed_destination)
#
# # update the file paths with the new pipestance locations
#
# pipestance_locs <-
#   analysis_configs |> pull(pipestance_location) |> str_replace("Blaser/", "Blaser/staff/") |> fs::path()
#
# walk(.x = pipestance_locs,
#      .f = \(x, pd = processed_destination) {
#        sample <- fs::path_file(x)
#        final_path <- fs::path(x, "outs", "per_sample_outs", sample, "count", "sample_feature_bc_matrix")
#        fs::file_copy(path = fs::path(final_path, "barcodes.tsv.gz"),
#                      new_path = fs::path(pd, str_glue("{sample}_barcodes.tsv.gz")))
#        fs::file_copy(path = fs::path(final_path, "features.tsv.gz"),
#                      new_path = fs::path(pd, str_glue("{sample}_features.tsv.gz")))
#        fs::file_copy(path = fs::path(final_path, "matrix.mtx.gz"),
#                      new_path = fs::path(pd, str_glue("{sample}_matrix.mtx.gz")))
#      })
#
# # make a directory in sc_working to collect all the raw data from network
#
# raw_destination <- fs::path("~/brad_workspace/sc_working/lapalombella_pu_geo/raw")
# fs::dir_create(raw_destination)
#
# # read in the location of the fastqs (raw) data
#
# raw_file_locs <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/geo/raw_data_location.csv") |>
#   mutate(flowcell = fs::path_file(fastqs)) |>
#   mutate(unique = str_glue("{fastq_id}_{flowcell}")) |>
#   mutate(fastqs = str_replace(fastqs, "workspace_pipelines/sc_working", "network/X/Labs/Blaser/staff/single_cell/lapalombella_pu")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_run1", "output_lapalombella_pu_run1_20210425050827/lapalombella_pu_run1")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_run2", "output_lapalombella_pu_run2_20210425134737/lapalombella_pu_run2")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_rerun_2_8", "output_lapalombella_pu_rerun_2_8_20210826062041/lapalombella_pu_rerun_2_8")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_rerun_1_3_5_12_13", "output_lapalombella_pu_rerun_1_3_5_12_13_20210826085244/lapalombella_pu_rerun_1_3_5_12_13")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_P9_P16/", "output_lapalombella_pu_P9_P16_20210720143419/lapalombella_pu_P9_P16/")) |>
#   mutate(fastqs = str_replace(fastqs, "lapalombella_pu_P9_P16_run2", "output_lapalombella_pu_P9_P16_run2_20210804181700/lapalombella_pu_P9_P16_run2")) |>
#   mutate(fastqs = fs::path(fastqs, fastq_id)) |>
#   select(unique, fastqs)
# raw_file_locs
#
# walk2(.x = raw_file_locs$fastqs,
#      .y = raw_file_locs$unique,
#     .f = \(x, y, rd = raw_destination) {
#       new_name <- paste0(y, "_", fs::path_file(fs::dir_ls(x)))
#       new_dest <- fs::path(rd, y)
#       files <- fs::dir_ls(x)
#       walk(.x = files,
#           .f = \(x, nd = new_dest) {
#             x
#             fs::file_copy(path = x,
#                           new_path = paste0(nd, "_", fs::path_file(x)))
#           })
#     })
#
# # md5 sum
# # processed files
#
# future::plan("multisession")
# processed_md5 <-
#   furrr::future_map_dfr(
#     .x = fs::dir_ls(processed_destination),
#     .f = \(x) {
#       name <- fs::path_file(x)
#       dig <- digest::digest(object = x, algo = "md5", file = TRUE)
#       tibble::tibble(file_name = name, file_checksum = dig)
#     },
#     .progress = TRUE
#   )
# write_csv(processed_md5, fs::path(geo_out, "processed_md5.csv"))
#
#
# raw_md5 <-
#   furrr::future_map_dfr(
#     .x = fs::dir_ls(raw_destination),
#     .f = \(x) {
#       name <- fs::path_file(x)
#       dig <- digest::digest(object = x, algo = "md5", file = TRUE)
#       tibble::tibble(file_name = name, file_checksum = dig)
#     },
#     .progress = TRUE
#   )
#
# write_csv(raw_md5, fs::path(geo_out, "raw_md5.csv"))
#
# # make a wide table with samples and files
#
#
#
# # samples processed data
# processed_files <- tibble(processed_file = fs::path_file(fs::dir_ls(processed_destination))) |>
#   mutate(sample = str_remove(processed_file, "_barcodes.tsv.gz|_features.tsv.gz|_matrix.mtx.gz"))  |>
#   mutate(type = str_extract(processed_file, "barcodes|features|matrix"))  |>
#   pivot_wider(names_from = "type", values_from = "processed_file") |> arrange(sample, )
#
# # samples raw data
# raw_files <- left_join(read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/geo/raw_data_location.csv", col_types = cols()) |> select(-fastqs),
#
# tibble(raw_file = fs::path_file(fs::dir_ls(raw_destination))) |>
#   mutate(fastq_id = str_remove(raw_file, "_H.*"))) |> distinct() |> select(-fastq_id) |>
#   group_by(sample) |>
#   mutate(file_no = paste0("raw_data_file", row_number())) |>
#   pivot_wider(names_from = "file_no", values_from = "raw_file")
#
# left_join(raw_files, processed_files) |>
#   select(sample, barcodes, features, matrix, starts_with("raw")) |> write_csv(fs::path(geo_out, "sample_files.csv"))
#
# # paired end experiments
# tibble(raw_file = fs::path_file(fs::dir_ls(raw_destination))) |>
#   mutate(fastq_id = str_remove(raw_file, "_H.*")) |>
#   mutate(cell = str_extract(raw_file, "H........")) |>
#   mutate(lane = str_extract(raw_file, "L00[12]")) |>
#   mutate(id_cell_lane = paste0(fastq_id, "_", cell, "_", lane)) |>
#   select(id_cell_lane, raw_file) |> group_by(id_cell_lane) |> mutate(file_name = paste0("file_name_", row_number())) |>
#   pivot_wider(names_from = file_name, values_from = raw_file) |>
#   mutate(orderer = str_extract(id_cell_lane, "P[:digit:]{1,2}")) |>
#   relocate(orderer) |>
#   mutate(orderer = str_remove(orderer, "P")) |>
#   mutate(orderer = as.numeric(orderer)) |>
#   arrange(orderer) |>
#   ungroup() |>
#   select(-c(orderer, id_cell_lane)) |>
#   write_csv(fs::path(geo_out, "paired_end_exp.csv"))
