# read in the analysis configs
analysis_configs <- read_excel("~/network/X/Labs/Blaser/collaborators/lapalombella_pu_network/analysis_configs.xlsx", sheet = "specimen_metadata") %>%
  mutate(pipestance_location = bb_fix_file_path(pipestance_location))

sequencing_metrics <- map_dfr(
  .x = analysis_configs$pipestance_location,
  .f = function(x) {
    read_csv(
      paste0(
        x,
        "/outs/per_sample_outs/",
        str_replace(x, pattern = ".*/", ""),
        "/metrics_summary.csv"
      )
    ) %>%
      filter(`Library or Sample` == "Sample") %>%
      mutate(specimen_fullname = str_replace(x, pattern = ".*/", "")) %>%
      select(specimen_fullname, `Metric Name`, `Metric Value`) %>%
      mutate(`Metric Value` = str_replace_all(string = `Metric Value`, pattern = ",|%", replacement = "")) %>%
      mutate(`Metric Value` = as.numeric(`Metric Value`)) %>%
      pivot_wider(names_from = `Metric Name`, values_from = `Metric Value`)

  }
)

sum(sequencing_metrics$`Number of reads assigned to the sample`)
