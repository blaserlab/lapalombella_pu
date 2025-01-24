colData(cds_main)$genotype <- factor(
  colData(cds_main)$genotype,
  levels = c("WT", "TP53-/-", "TP53-/-/TET2-/-")
)

#Figure 6A
F6A1<- bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia"))
),
"partition",
overwrite_labels = T) + facet_wrap( ~ genotype)  +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = "black"),
    legend.position = "none"
  )
F6A2<-bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia"))
),
  var = "density",
  facet_by = "genotype") +
  facet_grid(col = vars(genotype)) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black"))

F6A<-F6A1/plot_spacer()/F6A2 + plot_layout(heights = c(2,-0.2, 2))

F6A

#####################################################################################################################
#Figure 6B
aml_plotlist <- map(.x = c("Cd34", "Mpo", "Kit", "Elane", "Calr", "Ctsg"),
                    .f = \(x, dat = filter_cds(
                      cds_main,
                      cells = bb_cellmeta(cds_main) |>
                        filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia"))
                    )) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80",
                                              limits = c(0, 2.5)) +
                        facet_wrap(~genotype, labeller = labeller(group = label_wrap_gen(width = 5, multi_line = TRUE))) +
                        scale_y_continuous(breaks = c(-10, 0, 10)) +
                        scale_x_continuous(breaks = c(-10, 0, 10)) +
                        theme(panel.spacing = unit(0.25, "lines")) +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
                        theme(panel.background = element_rect(color = "black",
                                                              fill = "white")) +
                        theme(axis.line = element_blank()) +
                        theme(legend.position = "none") +
                        labs(x = NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))

                      # Remove x-axis text and ticks for all plots except Ctsg
                      if (x != "Ctsg") {
                        p <- p + theme(axis.text.x = element_blank(),
                                       axis.ticks.x = element_blank())
                      }

                      # Remove facet title for all plots except Cd34
                      if (x != "Cd34") {
                        p <- p + theme(strip.text = element_blank())
                      }

                      p
                    })

F6B <- ggarrange(aml_plotlist[[1]],
                           aml_plotlist[[2]],
                           aml_plotlist[[3]],
                           aml_plotlist[[4]],
                           aml_plotlist[[5]],
                           aml_plotlist[[6]],
                           ncol = 1,
                           nrow = 6,
                           common.legend = TRUE,
                           legend = "right")

#ggsave("F6B.pdf", path = figs_out, width = 8.25, height = 4.5)
#####################################################################################################################
# #Figure 6C: Heatmap
# F6_topmarkers_part <-
#   monocle3::top_markers(
#     filter_cds(
#       cds_main,
#       cells = bb_cellmeta(cds_main) |>
#         filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia"))
#     ),
#     group_cells_by = "partition",
#     genes_to_test_per_group = 15,
#     cores = 12
#   )

F6_topmarkers_part <- read.csv("~/network/T/Labs/EHL/Senior Associate group/Rosa/Ethan/1. EHL/Tet2_P53/Data/F6_topmarkers_part.csv")

markers <- F6_topmarkers_part |> pull(gene_short_name)

#Expression Matrix
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")),
      genes = bb_rowmeta(cds_main) |>
        filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, partition)
  ) |>
  t() |>
  scale() |>
  t()
rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(cds_main) |>
              select(feature_id, gene_short_name)) |>
  pull(gene_short_name)

colfun = circlize::colorRamp2(breaks = c(min(mat),
                                         0,
                                         max(mat)),
                              colors = heatmap_3_colors)

#Annotation: Top 5 marker genes per cluster
heatmap_highlights <-
  F6_topmarkers_part |> group_by(cell_group) |>
  slice_max(order_by = marker_score, n=5)|>
  pull(gene_short_name)

anno <-
  ComplexHeatmap::rowAnnotation(link =  anno_mark(
    at = which(rownames(mat) %in% heatmap_highlights),
    labels = rownames(mat)[rownames(mat) %in% heatmap_highlights],
    labels_gp = gpar(fontsize = 7),
    padding = 0.15
  ))

F6C_heatmap<-grid.grabExpr(draw(
ComplexHeatmap::Heatmap(
  mat,
  col = colfun,
  name = "Expression",
  show_row_names = F,
  show_column_names = F, #check column clustering order for bp
  right_annotation = anno,
  #top_annotation = hmap_bp,
  #width = ncol(mat)*unit(0.1, "mm"),
  height = nrow(mat)*unit(0.8, "mm"),
  row_dend_width = unit(6, "mm"),
  column_dend_height = unit(3, "mm"),
  heatmap_legend_param = list(
    legend_direction = "vertical",
    #legend_height = unit(1, "cm"),
    legend_width = unit(0.25, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 8.5)
  )
)))

# Stacked bar chart: heatmap cell composition
#normalized to cell quantity contributed from each leukemia phenotype
cellgeno_sums <- bb_cellmeta(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")))) |>
  count(genotype, name = "geno_sum")

cellcount <- bb_cellmeta(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")))) |>
  group_by(partition, genotype, leiden_assignment2) |>
  summarise(n = n()) |>
  left_join(cellgeno_sums) |>
  #mutate(leukemia_phenotype = recode(leukemia_phenotype, "No leukemia" = "WT")) |>
  mutate(total = sum(cellgeno_sums$geno_sum)) |>
  mutate(ratio = geno_sum/total) |>
  mutate(normalized_cell_frac = n/ratio)

#write.csv(cellcount, file = "~/network/T/Labs/EHL/Rosa/Ethan/10x/Tet2_P53/cellcount.csv")

#factor levels
cellcount$partition <- factor(cellcount$partition,
                                     levels = c("12","5", "2", "4", "8", "13", "6", "3", "9","11","14", "1", "16", "10", "7")) #for all
#plot
F6C_hmap_bp <-
  ggplot(cellcount,
         aes(x = partition, y = normalized_cell_frac, fill = genotype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.9) +
  theme_minimal()+
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(legend.justification = "top")+
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 9))+
  coord_fixed(ratio = 1.5)+ #0.45 #1.2
  theme(plot.margin = unit(c(-0.1,3.23,0,.23), "cm")) + #0.4 for select clusters
  #theme(plot.margin = unit(c(1,3.31,0,0.13), "cm")) + #0.4 for select clusters
  labs(x = "Partition Clusters")+
  scale_x_discrete(expand = c(0.07,0))+ #(0.13,0) for select clusters
  theme(legend.margin=margin(t = -0.22, unit='cm'))+
  theme(plot.margin = unit(c(1,0.1,1,0.25), "cm"))
F6C_hmap_bp

#ggsave("hmap_bp.pdf", path = figs_out, width = 6.8, height = 2.15)

########################################################################################################################
#Figure 6D
aml <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12')))

F6D <- bb_var_umap(blaseRtools::filter_cds(cds = cds_main,
                                    cells = bb_cellmeta(cds_main) |>
                                      filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12'))), "leiden_assignment2", overwrite_labels = T, group_label_size = 5)
#ggsave("Figure6D-aml_leiden_clusters.pdf", path = figs_out)
F6D

#Figure 6E
markers <- AML_topmarkers |> pull(gene_short_name)

bb_cellmeta(aml) |> count(leiden_assignment2, leiden)

#Expression Matrix
mat <-
  bb_aggregate(
    obj = filter_cds(
      aml,
      cells = bb_cellmeta(aml) |>
        filter(
          leiden_assignment2 %in% c('4', '5', '8', '24', '12')
        ),
      genes = bb_rowmeta(aml) |>
        filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(aml) |>
      select(cell_id, leiden_assignment2)
  ) |>
  t() |>
  scale() |>
  t()
rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(aml) |>
              select(feature_id, gene_short_name)) |>
  pull(gene_short_name
  )

colfun = circlize::colorRamp2(breaks = c(min(mat),
                                         0,
                                         max(mat)),
                              colors = heatmap_3_colors)

#Annotation: Top 5 marker genes per cluster
# heatmap_highlights <-
#   AML_topmarkers |> group_by(cell_group) |>
#   slice_max(order_by = marker_score, n=5)|>
#   pull(gene_short_name)

highlights1 <- c("Ifitm3", "Ccl6","Mpo", "Ctsg", "Elane", "Cd34", "Birc5", "S100A8", "S100A9", "Lyz2", "Cd52", "Stmn1", "Cd34", "Mpo", "Klf4", "Il7r", "Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5")

highlights2 <- c("Ifi27l2a", "S100a8", "S100a9", "Wfdc17", "Trem2", "Lyz2", "Ccl6", "Fcer1g", "Ftl1", "Cd52","Lgals3", "Fn1", "Ifi213","Hbb-bt","Hbb-bs","Hba-a1","Hba-a2","Lcn2", "Hmgn2","Ngp", "Camp", "Ngp","Ltf","Lcn","Ifitm6","Wfdc21","Elane","Mpo","Top2a","Tubb5", "Tubb4b", "Birc5", "Ran", "Stmn1" ,"Cox5b","Atp5c1","Chchd2","Mdh2","Atp5j2","Mif","Npm1","Atp5g1","Uqcrb")
heatmap_highlights2 <- unique(c(highlights1, highlights2))
#heatmap_highlights <- unique(c(heatmap_highlights, highlights))
anno <-
  ComplexHeatmap::rowAnnotation(link =  anno_mark(
    at = which(rownames(mat) %in% heatmap_highlights2),
    labels = rownames(mat)[rownames(mat) %in% heatmap_highlights2],
    labels_gp = gpar(fontsize = 6),
    padding = 0.35
  ))

#Figure 6E
F6E<-grid.grabExpr(draw(
ComplexHeatmap::Heatmap(
  mat,
  col = colfun,
  name = "Expression",
  show_row_names = F,
  show_column_names = T,
  column_names_rot = 45,
  right_annotation = anno,
  #top_annotation = hmap_bp,
  width = ncol(mat)*unit(15, "mm"),
  height = nrow(mat)*unit(1.25, "mm"),
  row_dend_width = unit(6, "mm"),
  column_dend_height = unit(3, "mm"),
  heatmap_legend_param = list(
    legend_direction = "vertical",
    #legend_height = unit(1, "cm"),
    legend_width = unit(0.25, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 8.5)
  )
)))
F6E<-as_ggplot(F6E)
F6E


