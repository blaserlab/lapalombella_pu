
#wt-aml-marrow
colData(cds_wt_aml_marrow)
cds <- reduce_dimension(cds_wt_aml_marrow)# not so sure why you want to recalculate dimensions
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden_assignment")

ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Jun",
                    "Tox",
                    "Elane",
                    "Cd3e")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
cds <- cluster_cells(cds)# also not sure why you are reclustering
plot_cells(cds, color_cells_by = "partition")


cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "leiden_assignment",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

#P53 T ALL
colData(cds_p53ALL)
cds <- reduce_dimension(cds_p53ALL)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden")
ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Jun",
                    "Tox",
                    "Elane",
                    "Cd3e")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "leiden",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

#PreB ALL
colData(cds_preBALL)
cds <- reduce_dimension(cds_preBALL)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden")
ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Jun",
                    "Tox",
                    "Elane",
                    "Cd3e")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "leiden",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

#Spleen
colData(cds_spleen)
cds <- reduce_dimension(cds_spleen)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden")
ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Jun",
                    "Tox",
                    "Elane",
                    "Cd3e")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "leiden",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)


#cds_wt-aml-spleen
colData(cds_wt_aml_spleen)
cds <- reduce_dimension(cds_wt_aml_spleen, reduction_method = c("UMAP"))
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden")
ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Jun",
                    "Tox",
                    "Elane",
                    "Cd3e")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cell_trajectory(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE)


cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

genes <- c("Ybx1")

colData(cds)
rowData(cds)
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes,
                       colData(cds)$leukemia_phenotype %in% c("AML")]

plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="leiden",
                         min_expr=0.5)
cds_subset <- choose_cells(cds)



#cds_aml_all
colData(cds_aml_all)
cds <- reduce_dimension(cds_aml_all)
colData(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "sample")

cds_aml_p5 <- cds_main[,colData(cds_main)$specimen %in% c("P5")]
cds <- reduce_dimension(cds_aml_p5)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "leiden")



ciliated_genes <- c("Cd34",
                    "Mpo",
                    "Ezh2",
                    "Tox",
                    "Elane",
                    "Cd3e")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
                     color_cells_by = "leiden",
                     label_groups_by_cluster=FALSE,
                     label_leaves=TRUE,
                     label_branch_points=FALSE)


cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

genes <- c("Ybx1", "Ezh2", "Elane")

colData(cds)
rowData(cds)
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes,
                   colData(cds)$leukemia_phenotype %in% c("AML")]

plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="leiden",
                         min_expr=0.5)

