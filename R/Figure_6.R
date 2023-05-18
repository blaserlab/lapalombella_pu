#Figure 6A
F6A1 <- bb_var_umap(
  cds_main,
  var = "partition",
  overwrite_labels = T
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(color = "black")) + theme(legend.position = "none")
F6A1
F6A2<-bb_var_umap(
  cds_main,
  var = "density",
  facet_by = "geno_pheno"
) + facet_grid(col = vars(geno_pheno)) + labs(x = "UMAP 1", y = "UMAP 2")+ theme_minimal() + theme(panel.grid.major = element_blank(),
                                                                                                   panel.grid.minor = element_blank()) +theme(panel.background = element_rect(color = "black"))
F6A2
#####################################################################################################################
#Figure 6B
aml_plotlist <- map(.x = c("Cd34","Mpo", "Kit", "Elane", "Calr","Ctsg"),
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        #alt_dim_x = "aggr_UMAP_1",
                        #alt_dim_y = "aggr_UMAP_2",
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80",
                                              limits = c(0,2.5)) +
                        facet_wrap(~geno_pheno, labeller = labeller(group = label_wrap_gen(width = 5, multi_line = TRUE))) +
                        scale_y_continuous(breaks = c(-10, 0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, 10))+
                        theme(panel.spacing = unit(0.25, "lines"))+
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black",
                                                              fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        #labs(x =NULL, y = NULL, subtitle = x) +
                        #theme(plot.subtitle = element_text(face = "italic", hjust = 0.5))+
                        theme(legend.position = "none")+
                        coord_fixed(ratio = 1)+ labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      # if (x != "Calr") p <- p + theme(legend.position = "none")
                      p
                    })

aml_gexp_umap <- ggarrange(aml_plotlist[[1]],
                           aml_plotlist[[2]],
                           aml_plotlist[[3]],
                           aml_plotlist[[4]],
                           aml_plotlist[[5]],
                           aml_plotlist[[6]],
                           ncol = 3,
                           nrow=2,
                           common.legend = TRUE,
                           legend="right")


#ggsave("F6C2.pdf", path = figs_out, width = 8.25, height = 4.5)
#####################################################################################################################
#Figure 6C: Heatmap

markers <- F6_topmarkers_part |> pull(gene_short_name)

#Expression Matrix
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(partition %in% c(1:16)),
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

heatmap_3_colors <-
  c("#313695", "white", "#A50026")

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

F6heatmap<-grid.grabExpr(draw(
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
cellpheno_sums <- bb_cellmeta(cds_main) |>
  count(leukemia_phenotype, name = "pheno_sum")

cellcount <- bb_cellmeta(cds_main) |>
  group_by(partition, leukemia_phenotype) |>
  summarise(n = n()) |> filter(partition %in% c(1:16)) |>
  left_join(cellpheno_sums) |>
  mutate(leukemia_phenotype = recode(leukemia_phenotype, "No leukemia" = "WT")) |>
  mutate(total = sum(cellpheno_sums$pheno_sum)) |>
  mutate(ratio = pheno_sum/total) |>
  mutate(normalized_cell_frac = n/ratio/4)

#factor levels
cellcount$partition <- factor(cellcount$partition,
                                     levels = c("12","2","4","15","6","3","9","14", "11","1","16","10","5","8","13","7")) #for all
cellcount$leukemia_phenotype <- factor(cellcount$leukemia_phenotype,
                                       levels = c("AML", "pre-B ALL","T ALL", "WT"))
#plot
hmap_bp <-
  ggplot(cellcount,
         aes(x = partition, y = normalized_cell_frac, fill = leukemia_phenotype)) +
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
hmap_bp
#ggsave("hmap_bp.pdf", path = figs_out, width = 6.8, height = 2.15)

#F6hm<- plot_grid(F6heatmap,hmap_bp, ncol=1, rel_heights = c(3,0.75))
#F6hm

########################################################################################################################
#Figure 6D
aml <- blaseRtools::filter_cds(cds = cds_main,
                               cells = bb_cellmeta(cds_main) |>
                                 filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12')))

F6D <- bb_var_umap(blaseRtools::filter_cds(cds = cds_main,
                                    cells = bb_cellmeta(cds_main) |>
                                      filter(leukemia_phenotype %in% c("AML")) |> filter(leiden %in% c('4', '5', '8', '24', '12'))), "leiden_assignment2", overwrite_labels = T, group_label_size = 5)
#ggsave("Figure6D-aml_leiden_clusters.pdf", path = figs_out)

#Figure 6E
markers <- AML_topmarkers |> pull(gene_short_name)

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

heatmap_3_colors <-
  c("#313695", "white", "#A50026")

colfun = circlize::colorRamp2(breaks = c(min(mat),
                                         0,
                                         max(mat)),
                              colors = heatmap_3_colors)

#Annotation: Top 5 marker genes per cluster
# heatmap_highlights <-
#   AML_topmarkers |> group_by(cell_group) |>
#   slice_max(order_by = marker_score, n=5)|>
#   pull(gene_short_name)

highlights1 <- c("Ccl6","Mpo", "Ctsg", "Elane", "Cd34", "Birc5", "S100A8", "S100A9", "Lyz2", "Cd52", "Stmn1", "Cd34", "Mpo", "Klf4", "Il7r", "Fcnb", "Nedd4", "Cebpe", "Ms4a3", "Ets1", "Mapk13", "Ifit1", "Ifit3", "Ifi47", "Il6ra", "Irf5")

highlights2 <- c("Ifi27l2a", "S100a8", "S100a9", "Wfdc17", "Trem2", "Lyz2", "Ccl6", "Fcer1g", "Ftl1", "Cd52","Lgals3", "Fn1", "Ifi213","Hbb-bt","Hbb-bs","Hba-a1","Hba-a2","Lcn2", "Hmgn2","Ngp", "Camp", "Ngp","Ltf","Lcn","Ifitm6","Wfdc21","Elane","Mpo","Top2a","Tubb5", "Tubb4b", "Birc5", "Ran", "Stmn1" ,"Cox5b","Atp5c1","Chchd2","Mdh2","Atp5j2","Mif","Npm1","Atp5g1","Uqcrb")
heatmap_highlights2 <- unique(c(highlights1, highlights2))
#heatmap_highlights <- unique(c(heatmap_highlights, highlights))
anno <-
  ComplexHeatmap::rowAnnotation(link =  anno_mark(
    at = which(rownames(mat) %in% heatmap_highlights2),
    labels = rownames(mat)[rownames(mat) %in% heatmap_highlights2],
    labels_gp = gpar(fontsize = 4),
    padding = 0.15
  ))

#Figure 6E
F6E<-grid.grabExpr(draw(
ComplexHeatmap::Heatmap(
  mat,
  col = colfun,
  name = "Expression",
  show_row_names = F,
  show_column_names = F,
  column_names_rot = 45,
  right_annotation = anno,
  #top_annotation = hmap_bp,
  #width = ncol(mat)*unit(0.1, "mm"),
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


#  scratch work below this ---------------------------
#pseudotime analysis

# bb_cellmeta(cds_p568) |> glimpse()
# bb_var_umap(cds_p568, "sample")
# bb_var_umap(cds_p568, "leukemia_phenotype", value_to_highlight = "AML") +
#   bb_var_umap(cds_p568, "leukemia_phenotype", value_to_highlight = "PreB ALL") +
#   bb_var_umap(cds_p568, "leukemia_phenotype", value_to_highlight = "No leukemia")
#
# bb_var_umap(cds_p568, "density", facet_by = "leukemia_phenotype")
# bb_var_umap(cds_p568, "leiden_assignment2", facet_by = "leukemia_phenotype")
#
# bb_cellmeta(muench_cds) |> glimpse()
#
# data(muench_cds, package = "lapalombella.pu.datapkg")

