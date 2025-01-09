# graphical parameters####
the_font_size <-10
theme_set(theme_cowplot(font_size = the_font_size))

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
heatmap_3_colors <- c("#041dff","white","#ff0505")

# conflicts ---------------------------------------------------------------
# resolve conflicting function names here
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("exprs", "monocle3")
conflicts_prefer(base::as.data.frame)
conflicts_prefer(base::intersect)

# output directories

figs_out <- "/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/figs"
tables_out <- "/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/tables"


