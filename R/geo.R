source("R/dependencies.R")
source("R/configs.R")

analysis_configs_cds_main_human
bb_cellmeta(cds_main_human_unaligned) |> glimpse()

sra <- readr::read_tsv(fs::path("data/tet2_p53_citeseq_sra.tsv"))
View(sra)

sample_metadata <- readr::read_csv("~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/patient_sample_metadata.csv",
                                   col_names = c("pid", "hash_id", "hash", "sample", "PZ", "date"),
                                   col_types = c("cccccc")) |>
  mutate(date = paste0(date, "-2024")) |>
  mutate(date = lubridate::dmy(date)) |>
  mutate(`*sample_name` = paste0(PZ, "_", sample)) |>
  group_by(`*sample_name`, date, sample) |>
  summarize()
sample_metadata

table_with_paths <- fs::dir_map(fs::path("~/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/geo/FASTQ/"), fs::path_file) |>
  unlist() |>
  tibble(filepaths = _) |>
  mutate(library_name = stringr::str_remove(filepaths, "_IGO.*")) |>
  mutate(`*sample_name` = stringr::str_remove(library_name, "_CITE")) |>
  relocate(`*sample_name`) |>
  mutate(`*organism` = "Homo sapiens") |>
  left_join(sample_metadata) |>
  mutate(`*age` = "adult") |>
  rename(`*isolate` = sample) |>
  rename(`*collection_date` = date) |>
  mutate(`*biomaterial_provider` = "Dr. Omar Abdel-Wahab, Molecular Pharmacology Program, Sloan Kettering Institute, New York, NY, USA") |>
  mutate(`*geo_loc_name` = "USA:New York") |>
  mutate(`*sex` = "pooled male and female") |>
  mutate(`*tissue` = "blood") |>
  mutate(sample_type = "PBMC") |>
  mutate(description = "CITE-seq ADT and Gene expression libraries from pooled patient samples") |>
  mutate(bioproject_accession = "PRJNA901155")

tet2_p53_biosamples <- table_with_paths |>
  select(-c(filepaths, library_name)) |>
  distinct() |>
  mutate(`*collection_date` = as.character(`*collection_date`)) |>
  pivot_longer(everything())
tet2_p53_biosamples

left_join(tibble(name = colnames(sra)), tet2_p53_biosamples) |>
  pivot_wider(names_from = "name", values_from = "value") |> unnest(cols = everything()) |>
  write_tsv("data/tet2_p53_biosamples.tsv")

# metadata --------------

table_with_paths |>
  rename(sample_name = `*sample_name`) |>
  rename(library_ID = library_name) |>
  mutate(library_ID = ifelse(str_detect(library_ID, "CITE"), library_ID, paste0(library_ID, "_GEX"))) |>
  mutate(title = ifelse(str_detect(library_ID, "CITE"),
         "CITE-seq ADT library of Homo sapiens: adult PBMC",
         "CITE-seq Gene Expression library of Homo sapiens:  adult PBMC")) |>
  mutate(library_strategy = "OTHER") |>
  mutate(design_description = ifelse(str_detect(library_ID, "CITE"), "CITE-seq ADT library", "CITE-seq gene expression library")) |>
  mutate(library_source = "TRANSCRIPTOMIC SINGLE CELL") |>
  mutate(library_selection = "Oligo-dT") |>
  mutate(library_layout = "Paired") |>
  mutate(platform = "ILLUMINA") |>
  mutate(instrument_model = "Illumina NovaSeq X Plus") |>
  mutate(filetype = "fastq") |>
  group_by(library_ID) |>
  mutate(rank_filepath = rank(filepaths, ties.method = "first")) |>
  mutate(rank_filepath = paste0("filename", rank_filepath)) |>
  mutate(rank_filepath = recode(rank_filepath, "filename1" = "filename")) |>
  select(c(sample_name, library_ID, title, library_strategy, library_source, library_selection, library_layout, platform, instrument_model, design_description, filetype, rank_filepath, filepaths)) |>
  pivot_wider(names_from = rank_filepath, values_from = filepaths, values_fill = "") |>
  mutate(assembly = NA) |>
  mutate(fasta_file = NA) |>
  write_tsv("data/tet2_p53_metadata.tsv")


count(library_ID)
  # count(lane_read)
  mutate(filepath_names = )
