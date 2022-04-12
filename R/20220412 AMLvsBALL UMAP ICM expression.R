requireData("lapalombella.pu.datapkg")
library(lapalombella.pu.datapkg)
renv::install("blaserlab/blaseRtools")
library(Seurat)
library("blaseRtools")
theme_set(theme_cowplot(font_size = 10))
colData(cds_main)
bb_var_umap(cds_main, var = "partition")
bb_var_umap(cds_p53ALL, var = "leiden")
cds_p53tet2AML<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "AML"]

cds_p53ALL<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "T cell leukemia"]
cds_preBALL<-cds_main[,colData(cds_main)$leukemia_phenotype %in% "PreB ALL"]
cds_spleen<-cds_main[,colData(cds_main)$tissue %in% "spleen"]

cds_bonemarrow<-cds_main[,colData(cds_main)$tissue %in% "bone marrow"]

bb_gene_umap(cds= cds_p53tet2AML, gene_or_gene = c("Kit", "Cd34", "Cd19", "Cd3e", "Cd4", "Cd8a", "Itgam", "S100a8", "S100a9"), cell_size = 0.1)

bb_gene_dotplot(cds=cds_p53tet2AML, markers=c("Cd34","Pvr","Kit","Cd47", "Il6","Cd19", "Vsir"), group_cells_by = "partition" )

bb_gene_dotplot(cds=cds_main, markers=c("Cd34","Pvr","Kit","Cd47", "Il6","Cd19", "Vsir"), group_cells_by = "leukemia_phenotype" )

bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Cd34","Cd47", "Vsir"), rows=3, palette =
                     c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"))

bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Cd34","Cd47", "Vsir"), rows=3)+scale_fill_viridis_d(option = "turbo")

bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Cd34","Cd47", "Vsir"), rows=3)+geom_point(aes(colour = partition))+scale_fill_viridis_d(option = "turbo")

display.brewer.pal(12,"Paired")

#plot selected genes in different partitions in AML and T ALL
bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6) + scale_fill_viridis_d(option = "turbo")
bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6, palette =
                     c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"))

bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6) + scale_fill_viridis_d(option = "turbo")
bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), rows=6, palette =
                     c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"))


bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Hif1a", "Nr4a1", "Exomes", "Tbx21", "Batf"), rows=8) + scale_fill_viridis_d(option = "turbo")

bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Hif1a", "Nr4a1", "Exomes", "Tbx21", "Batf"), rows=8) + scale_fill_viridis_d(option = "turbo")

bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Lag3", "Havcr2", "Pdcd1"), rows=8) + scale_fill_viridis_d(option = "turbo")


#plot selected genes in different leukemia_phenotypes

bb_gene_dotplot(cds=cds_main, markers=c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), group_cells_by = "leukemia_phenotype" )+ scale_fill_viridis_d(option = "turbo")

bb_gene_dotplot(cds=cds_main, markers=c("Cd19", "Cd34", "Kit", "Cd8a", "Cd3e", "Cd4", "Cd37", "Cd25"), group_cells_by = "leukemia_phenotype" )+ scale_fill_viridis_d(option = "turbo")

cowplot::save_plot(AMlumap, file = "Amlumap.png", base_height = 11.29, base_width = 8.27)

bb_gene_dotplot(cds=cds_bonemarrow, markers=c("S100a9","Tmem173"), group_cells_by = "leukemia_phenotype" )

bb_gene_dotplot(cds=cds_bonemarrow, markers=c("Elane", "Ctsg", "Mpo", "Calr", "S100a8", "S100a9"), group_cells_by = "leukemia_phenotype" )+ scale_fill_viridis_d(option = "turbo")



bb_gene_violinplot(cds=cds_p53tet2AML, variable="partition", genes_to_plot = c("Lag3", "Tigit", "Pdcd1", "Cd8a", "Cd3e"), rows=8) + scale_fill_viridis_d(option = "turbo") + scale_fill_brewer(palette = "Paired")

bb_gene_violinplot(cds=cds_p53ALL, variable="partition", genes_to_plot = c("Lag3", "Tigit", "Pdcd1", "Cd8a", "Cd3e"), rows=8) + scale_fill_viridis_d(option = "turbo") + scale_fill_brewer(palette = "Paired")

bb_gene_umap(cds= cds_p53tet2AML, gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)

bb_gene_umap(cds= cds_p53ALL, gene_or_gene = c("Cd3e","Cd8a","Pdcd1","Lag3","Tigit", "Kit", "Cd34", "Cd19",  "Cd4"), cell_size = 0.1)

#heatmap

bb_cds_heatmap(cds=cds_p53tet2AML)

exprs(cds_main)

#annotation
requireData("blaseRextras")
renv::install("/workspace/zhang_workspace/Lapalombella_Pu-1/Git clone/lapalombella_pu/R/pu_work/renv/library/R-4.1/x86_64-pc-linux-gnu/lapalombella.pu.datapkg/renv/library/R-4.1/x86_64-pc-linux-gnu/lapalombella.pu.datapkg/blaseRextras_0.0.0.9001.tar.gz")
cds_p53tet2AMLann<-bb_cds_anno(query_cds=cds_p53tet2AML, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )
cds_mainann<-bb_cds_anno(query_cds=cds_main, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )
cds_preBALLann<-bb_cds_anno(query_cds=cds_preBALL, ref=cds_wt_aml_marrow,transfer_col="leiden_assignment", unique_id = NULL )

colData(muench_cds)
colData(cds_main)
colData(cds_wt_aml_spleen)
bb_var_umap(muench_cds, var = "muench_cluster2")
colData(cds_p53tet2AMLann)
bb_var_umap(cds_p53tet2AMLann, var="predicted.leiden_assignment")
colData(cds_main)
bb_var_umap(cds_main, var="predicted.muench_cluster")
bb_var_umap(cds_preBALLann, var="predicted.leiden_assignment")

bb_var_umap(cds_p53tet2AMLann, var="predicted.leiden_assignment")+ bb_var_umap(cds_preBALLann, var="predicted.leiden_assignment")

bb_gene_umap(cds= cds_p53tet2AMLann, gene_or_gene = c("Cd3e","Cd8a","Cd4", "Cd34", "Cd19", "Pax5"), cell_size = 0.1)

bb_gene_umap(cds= cds_preBALLann, gene_or_gene = c("Cd3e","Cd8a","Cd4", "Cd34", "Cd19", "Pax5"), cell_size = 0.1)

#Violin plot for selected genes in AML and preB ALL
bb_gene_violinplot(cds=cds_p53tet2AMLann, variable="predicted.leiden_assignment", genes_to_plot = c("Tox", "Ikzf2","Ikzf3"), rows= 7) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_p53tet2AMLann, variable="predicted.leiden_assignment", genes_to_plot = c( "Slamf6", "Sell","Cd8a", "Mki67", "Cd44"), rows= 9) + scale_fill_brewer(palette = "Paired")

bb_gene_violinplot(cds=cds_preBALLann, variable="predicted.leiden_assignment", genes_to_plot = c("Tox", "Ikzf2","Ikzf3"), rows= 7) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_preBALLann, variable="predicted.leiden_assignment", genes_to_plot = c("Slamf6", "Sell","Cd8a", "Mki67", "Cd44"), rows= 9) + scale_fill_brewer(palette = "Paired")

#subset T/NK clusters in AML and study exhaustive markers
cds_aml_T <- cds_p53tet2AMLann[,colData(cds_p53tet2AMLann)$predicted.leiden_assignment %in% c("T/NK 1", "T/NK 2")]

colData(cds_aml_T)

bb_var_umap(cds_aml_T, var="louvain")

bb_gene_umap(cds=cds_aml_T, gene_or_gene = c("Cd3e","Cd8a","Cd4", "Cd34", "Cd19", "Pax5"), cell_size = 0.1)

bb_gene_violinplot(cds=cds_aml_T, variable="leiden", genes_to_plot = c("Pdcd1", "Lag3", "Tox", "Ctla4", "Tigit", "Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_aml_T, variable="leiden", genes_to_plot = c("Tox","Eomes", "Klrg1", "Pdcd1","Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")

#subset T/NK clusters in Pre-B ALL and study exhaustive markers

cds_BALL_T <- cds_preBALLann[,colData(cds_preBALLann)$predicted.leiden_assignment %in% c("T/NK 1", "T/NK 2")]

colData(cds_BALL_T)
bb_var_umap(cds_BALL_T, var="louvain")

bb_gene_umap(cds=cds_BALL_T, gene_or_gene = c("Cd3e","Cd8a","Cd4", "Cd34", "Cd19", "Pax5"), cell_size = 0.1)

bb_gene_violinplot(cds=cds_BALL_T, variable="leiden", genes_to_plot = c("Pdcd1", "Lag3", "Tox", "Ctla4", "Tigit", "Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")
bb_gene_violinplot(cds=cds_BALL_T, variable="louvain", genes_to_plot = c("Tox","Eomes", "Klrg1", "Pdcd1","Havcr2"), rows= 7) + scale_fill_brewer(palette = "Paired")
