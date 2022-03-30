# renv --------------------------------------------------------------------

# set up the renv from scratch

# renv::init(bioconductor = TRUE)

# restore the renv from the lockfile

# renv::restore()



# package installation ----------------------------------------------------

# # Try this first...it's faster:
# blaseRtemplates::easy_install("<package name>", how = "link_from_cache")

# # If you need a new package or an update, try this:
# blaseRtemplates::easy_install("<package name>", how = "new_or_update")
# blaseRtemplates::easy_install("blaserlab/blaseRtools", how = "new_or_update")

# # If you are installing from a "tarball", use this:
# blaseRtemplates::easy_install("/path/to/tarball.tar.gz")

# # use "bioc::<package name>" for bioconductor packages
# # use "<repo/package name>" for github source packages

# load core packages for the analysis ------------------------------------
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
bb_renv_datapkg("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/datapkg")


# load the data set into a hidden environment
requireData("lapalombella.pu.datapkg")
