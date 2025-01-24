#Supplemental Figure 4
#Supplemental Figure 4A
b_plotlist <- map(.x = c("Cd19", "Pax5", "Cd24a", "Cd79a", "Vpreb1", "Vpreb2", "Vpreb3", "Dntt"),
                    .f = \(x, dat = filter_cds(
                      cds_main,
                      cells = bb_cellmeta(cds_main) |>
                        filter(geno_pheno %in% c("dKO: AML", "Wildtype", "P53 KO: T ALL"))
                    )) {
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
                        facet_wrap(~genotype) + scale_y_continuous(breaks = c(-10,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, 10))+
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

Supp_Fig_4A.1 <- b_plotlist[[1]] / b_plotlist[[2]] /
  b_plotlist[[3]] / b_plotlist[[4]]

Supp_Fig_4A.2 <- b_plotlist[[5]] / b_plotlist[[6]] /
  b_plotlist[[7]] / b_plotlist[[8]]

Supp_Fig_4A.1
Supp_Fig_4A.2

#Supplemental Figure 4B
t_plotlist <- map(.x = c("Atxn1", "Ly6a", "Cd7", "Cd3d"), #Atxn1 = Sca1
                    .f = \(x, dat = filter_cds(
                      cds_main,
                      cells = bb_cellmeta(cds_main) |>
                        filter(geno_pheno %in% c("dKO: AML", "Wildtype", "P53 KO: T ALL"))
                    )) {
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
                        facet_wrap(~genotype) + scale_y_continuous(breaks = c(-10,0, 10)) +
                        scale_x_continuous(breaks = c(-10,0, 10))+
                        theme(panel.spacing = unit(0.5, "lines"))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                        theme(panel.background = element_rect(color = "black", fill = "white"))+
                        theme(axis.line = element_blank()) +
                        #theme(axis.ticks = element_blank()) +
                        #theme(axis.text = element_blank()) +
                        labs(x =NULL, y = x) +
                        theme(axis.title.y = element_text(size = 14, face = "italic"))
                      if (x != "Cd3d") p <- p + theme(legend.position = "none")
                      p
                    })
Supp_Fig_4B <- t_plotlist[[1]] / t_plotlist[[2]] /
  t_plotlist[[3]] / t_plotlist[[4]]
Supp_Fig_4B
