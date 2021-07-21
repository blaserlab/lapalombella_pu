# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# install into your home directory if you haven't already
# devtools::install_github("git@github.com:blaserlab/blaseRtools.git")


# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/blaseRtools")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DESeq2")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/genefilter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/annotate")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/AnnotationDbi")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/KEGGREST")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Biostrings")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/geneplotter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DoubletFinder")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Seurat")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/SeuratDisk")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/rrvgo")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GO.db")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GOSemSim")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/scater")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/topGO")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/fastSave")

# load core packages for the analysis
library("blaseRtools")
library("tidyverse")
library("monocle3")
library("circlize")
library("ComplexHeatmap")
library("lazyData")
library("cowplot")
library("RColorBrewer")
library("ggrepel")
library("ggpubr")
library("rstatix")
library("readxl")


# run this to update the data package in renv
# bb_renv_datapkg("~/network/X/Labs/Blaser/collaborators/lapalombella_pu_network/datapkg")


# load the data set into a hidden environment
# requireData("lapalombella.pu.datapkg")

