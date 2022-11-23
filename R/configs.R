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



# output directories

figs_out <- "~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/figs"
T_figs <- "~/network/T/Labs/EHL/Rosa/Ethan/10X/Tet2_P53/Figures"
pu_figs <- "~/network/T/Labs/EHL/Mouse Group/1.MOUSE COLONIES/p53 colonies RL/Mechanistic studies/scRNA-seq/scRNA-seq analysis/2022.11 Ethan's analysis for manuscript"
tables_out <- "~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/tables"
geo_out <- "~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/geo"

# source local configs ----------------------------------------------------
# these are sourced after main configs and will overwrite duplicate entries if
# present. The file local_configs.R is ignored by git and so is useful for user-
# specific configurations such as output directories or formatting.

#source("R/local_configs.R")
