#Supplemental Figure 7
#Supp Fig 7A
leiden_clust <- bb_var_umap(cds_main, var = "leiden", overwrite_labels = T)
leiden_clust
#scType Assignment.R is a dependency
#Supp Fig 7B
leiden_ScType_hm<- ComplexHeatmap::Heatmap(es.max,
                                           col = colfun)
leiden_ScType_hm

#Supp Fig 7C
plotlist <- map(.x = c("Cd3e", "Cd4", "Gzma", "Icos", "Izumo1r", "Trbc1"),
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

gexp_umap <- ggarrange(plotlist[[1]],
                           plotlist[[2]],
                           plotlist[[3]],
                           plotlist[[4]],
                           plotlist[[5]],
                           plotlist[[6]],
                           ncol = 3,
                           nrow=2,
                           common.legend = TRUE,
                           legend="right")
gexp_umap

#Supp Fig 7D
unique(colData(cds_main)$geno_pheno)

bb_gene_dotplot(cds=cds_main, markers=c("Cd19", "Cd34", "Kit", "Cd8a", "Cd3e", "Cd4", "Cd37"), group_cells_by = "geno_pheno" )+ scale_fill_viridis_d(option = "turbo")+
  theme_minimal_grid(font_size = 12) +
  theme(
    strip.background = ggh4x::element_part_rect(
      side = "b",
      colour = "black",
      fill = "transparent"
    )
  )

#Supp Fig 7E
marrow_gb_plot <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(tissue %in% c("marrow")) |>
        #filter(leukemia_phenotype %in% c("AML", "pre-B ALL")) |>
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
    cell_grouping = c("leiden_assignment2", "leukemia_phenotype"),
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
  facet_wrap( ~leukemia_phenotype, scales = "free_x",) +
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
marrow_gb_plot

#Supp Figure 7F - AML
bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("AML"))),
  gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)

#Supp Figure 7G - B-ALL
bb_gene_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |>
    filter(leukemia_phenotype %in% c("pre-B ALL"))),
  gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)
