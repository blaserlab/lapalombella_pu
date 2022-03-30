requireData("lapalombella.pu.datapkg")
library(lapalombella.pu.datapkg)
renv::install("blaserlab/blaseRtools")
library(Seurat)
library("blaseRtools")

#wt-aml-marrow
colData(cds_wt_aml_marrow)
        cds <- reduce_dimension(cds_wt_aml_marrow)
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
cds <- cluster_cells(cds)
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
