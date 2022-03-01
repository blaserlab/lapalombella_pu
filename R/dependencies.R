# libraries-------------------------------------
# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("blaserlab/blaseRtools")

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
library("knitr")
library("pander")
library("conflicted")

# run this to update the data package in renv
bb_renv_datapkg("~/network/X/Labs/Blaser/collaborators/lapalombella_pu_network/datapkg")


# load the data set into a hidden environment
requireData("lapalombella.pu.datapkg")

