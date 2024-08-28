

bb_cellmeta(cds_WT_AML_bALL) |> glimpse()
bb_cellmeta(cds_WT_AML_bALL) |> count(specimen, tissue, genotype, leukemia_phenotype)

bb_var_umap(cds_WT_AML_bALL, "density", facet_by = c("genotype", "leukemia_phenotype"), cols = (vars(genotype)), rows = vars(leukemia_phenotype))


bb_var_umap(cds_WT_AML_bALL, "density", facet_by = c("genotype"))

p0 <- test_bb_var_umap(cds_WT_AML_bALL, "density")
p0
p1 <- test_bb_var_umap(cds_WT_AML_bALL, "density", facet_by = c("genotype"))
p1
p2 <- test_bb_var_umap(cds_WT_AML_bALL, "density", facet_by = c("leukemia_phenotype"))
p2
p3 <- test_bb_var_umap(cds_WT_AML_bALL, "density", facet_by = c("genotype", "leukemia_phenotype"), cell_size = 2)
p3
p3 + scale_color_viridis_c(limits = c(0, 0.2), na.value = "red") + scale_fill_viridis_c(limits = c(0, 0.2), na.value = "red")
left_join(p1$data |> as_tibble() |> select(cell_id, p1_density = density),
          p2$data |> as_tibble() |> select(cell_id, p2_density = density)) |>
  left_join(p3$data |> as_tibble() |> select(cell_id, p3_density = density)) |>
  left_join(p0$data |> as_tibble() |> select(cell_id, p0_density = density))


test_bb_var_umap(cds_WT_AML_bALL, "density")
p3$compound_variable


bb_gene_umap(cds_WT_AML_bALL, "Pvr") +
bb_gene_umap(cds_WT_AML_bALL, "Nectin2")

