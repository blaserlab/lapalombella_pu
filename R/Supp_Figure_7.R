#Supplemental Figure 7
Supp_Fig_7A <- bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      #partition %in% c(1:16)
      leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")
    )),
  var = "leiden",
  cell_size = 0.1,
  overwrite_labels = T)

Supp_Fig_7A

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_7A.pdf"),
  plot = Supp_Fig_7A,
  base_width = 6.0,
  base_height = 6.0
)

#ScType
library(HGNChelper)
#
source(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
)
source(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
)

# get cell-type-specific gene sets from ScType database
gs_list = gene_sets_prepare(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
  "Immune system"
) # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# #count matrix - leiden clusters
mat <-
  bb_aggregate(
    obj = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")
        ),
      genes = bb_rowmeta(cds_main) #|>
      #   filter(gene_short_name %in% markers)
    ),
    cell_group_df = bb_cellmeta(cds_main) |>
      select(cell_id, leiden)#barcode)
  ) |>
  t() |>
  scale() |>
  t()

rownames(mat) <-
  tibble(feature_id = rownames(mat)) |>
  left_join(bb_rowmeta(cds_main) |>
              select(feature_id, gene_short_name)) |> pull(gene_short_name)

# assign cell types
es.max = sctype_score(scRNAseqData = mat, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
rownames(es.max)

#es.max2 <- as.data.frame(es.max)
#es.max2$Largest_Column <-colnames(es.max2)[apply(es.max2,1,which.max)]

harmony_colfun <- circlize::colorRamp2(breaks = c(0, 1), colors = c("grey80", "red"))
heatmap_3_colors <-
  c("#313695", "white", "#A50026")

colfun = circlize::colorRamp2(breaks = c(min(es.max),
                                         0,
                                         max(es.max)),
                              colors = heatmap_3_colors)

Supp_Fig_7B<- ComplexHeatmap::Heatmap(es.max,
                                      col = colfun,
                                      heatmap_legend_param = list(
                                        title = "", # Change to your desired legend title
                                        title_position = "topcenter", # Align title in the center at the top
                                        legend_direction = "vertical" # Move the legend to the bottom
                                      )
)
Supp_Fig_7B

#Supp Fig 7C
plotlist <- map(.x = c("Cd3e", "Cd4", "Gzma", "Icos", "Izumo1r", "Trbc1"),
                    .f = \(x, dat = filter_cds(
                      cds_main,
                      cells = bb_cellmeta(cds_main) |>
                        filter(
                          leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")
                        ))) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Oranges",
                                              direction = 1,
                                              na.value = "grey80",
                                              limits = c(0,2.5)) +
                        facet_wrap(~genotype, labeller = labeller(group = label_wrap_gen(width = 5, multi_line = TRUE))) +
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

Supp_Fig_7C <- ggarrange(plotlist[[1]],
                           plotlist[[2]],
                           plotlist[[3]],
                           plotlist[[4]],
                           plotlist[[5]],
                           plotlist[[6]],
                           ncol = 2,
                           nrow= 3,
                           common.legend = TRUE,
                           legend="right")
Supp_Fig_7C

#Supp Fig 7D

Supp_Fig_7D <- bb_gene_dotplot(cds=filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(
      leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")
    )), markers=c("Cd19", "Cd34", "Kit", "Cd8a", "Cd3e", "Cd4", "Cd37"), group_cells_by = "geno_pheno" )+ scale_fill_viridis_d(option = "turbo")+
  theme_minimal_grid(font_size = 12) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  )

Supp_Fig_7D

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_7D.pdf"),
  plot = Supp_Fig_7D,
  base_width = 5,
  base_height = 5
)

#Supp Fig 7E
Supp_Fig_7E <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(tissue %in% c("marrow")) |>
        filter(leukemia_phenotype %in% c("AML", "T ALL", "No leukemia")) |>
        filter(str_detect(leiden_assignment2, "T|B|NK|Natural")) |>
        filter(str_detect(leiden_assignment2, "ALL", negate = TRUE))
    ),
    genes = c(
      "Cd37",
      "Cd79a",
      "Cd3e",
      "Cd4",
      "Gzma",
      "Icos",
      "Izumo1r",
      "Trbc1"
    ),
    #"Itgam", "Ly6g", "Gzma"),
    cell_grouping = c("leiden_assignment2", "genotype"),
    return_value = "data"
  ) |>
  ggplot(mapping = aes(
    x = leiden_assignment2,
    y = gene_short_name,
    color = expression,
    size = proportion
  )) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  facet_wrap( ~genotype, scales = "free_x",) +
  theme_minimal_grid(font_size = 6) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL,
       y = NULL,
       size = "Proportion",
       color = "Expression")

Supp_Fig_7E

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_7E.pdf"),
  plot = Supp_Fig_7E,
  base_width = 10,
  base_height = 3
)

#Supp Figure 7F - dKO AML mice
Supp_Fig_7F <- bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML"))),
  gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)

Supp_Fig_7F

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_7F.pdf"),
  plot = Supp_Fig_7F,
  base_width = 10,
  base_height = 8
)

#Supp Figure 7G - TP53 KO mice
Supp_Fig_7G <- bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("T-ALL"))),
  gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)

Supp_Fig_7G

save_plot(
  #filename = "temp.pdf",
  filename = fs::path(figs_out, "Supp_Fig_7G.pdf"),
  plot = Supp_Fig_7G,
  base_width = 10,
  base_height = 8
)
