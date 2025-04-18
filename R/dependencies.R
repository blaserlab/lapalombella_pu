# blaseRtemplates::install_one_package("blaseRtools", "link")
# load core packages for the analysis ------------------------------------
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRtemplates"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
#suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("todor"))

# run this to update the data package
blaseRtemplates::project_data("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/datapkg")

#blaseRtemplates::install_one_package("todor")
# todor::todor(c("TODO", "NOTE"))

