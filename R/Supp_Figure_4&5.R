#Supp Figure 4
b_all_plotlist <- map(.x = c("Cd19", "Pax5", "Cd24a", "Cd79a", "Vpreb1", "Vpreb2", "Vpreb3", "Dntt"),
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        #alt_dim_x = "aggr_UMAP_1",
                        #alt_dim_y = "aggr_UMAP_2",
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Greens",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) +
                        facet_wrap(~geno_pheno) + scale_y_continuous(breaks = c(-15,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, -10))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      if (x != "Vpreb3") p <- p + theme(legend.position = "none")
                      p
                    })

b_all_plotlist[[1]] / b_all_plotlist[[2]] |
  b_all_plotlist[[3]] / b_all_plotlist[[4]]

b_all_plotlist[[5]] / b_all_plotlist[[6]] |
  b_all_plotlist[[7]] / b_all_plotlist[[8]]


bb_gene_umap(cds_main, gene_or_genes = c("Dntt"), alt_dim_x = "aggr_UMAP_1", alt_dim_y = "aggr_UMAP_2")+ facet_wrap(~leukemia_phenotype)

t_all_plotlist <- map(.x = c("Atxn1", "Ly6a", "Cd7", "Cd3d"), #Atxn1 = Sca1
                    .f = \(x, dat = cds_main) {
                      p <- bb_gene_umap(
                        dat,
                        gene_or_genes = x,
                        #alt_dim_x = "aggr_UMAP_1",
                        #alt_dim_y = "aggr_UMAP_2",
                        cell_size = 0.1
                      ) +
                        scale_color_distiller(palette = "Blues",
                                              direction = 1,
                                              na.value = "grey80", limits = c(0,2)) + #check limits!
                        facet_wrap(~geno_pheno) + scale_y_continuous(breaks = c(-15,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, -10))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      if (x != "Tdt") p <- p + theme(legend.position = "none")
                      p
                    })
t_all_plotlist[[1]] / t_all_plotlist[[2]] |
  t_all_plotlist[[3]] / t_all_plotlist[[4]]
