# graphical parameters####
theme_set(theme_cowplot(font_size = 10))

# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "AML" = "#DC0000",
  "WT" = "#3C5488"

)

jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#313695","white","#A50026")

# unmask important functions
filter <- dplyr::filter
mutate <- dplyr::mutate
group_by <- dplyr::group_by
select <- dplyr::select
rename <- dplyr::rename

# output directories

figs_out <- "~/network/X/Labs/Blaser/collaborators/lapalombella_pu_network/figs"
tables_out <- "~/network/X/Labs/Blaser/collaborators/lapalombella_pu_network/tables"
