## Analysis Code for "TP53 mutations and TET2 deficiency cooperate to drive leukemogenesis and establish an immunosuppressive environment"

Citation: J Clin Invest. 2025;135(10):e184021. <https://doi.org/10.1172/JCI184021>.

### Instructions for use

This code was developed on R v4.4.0 for Linux. It can be used to reproduce figures showing scRNA-seq and CITE-seq data from the article referenced above.

This code is meant to be used together with an R data package containing processed data. The data package is available on figshare at the following URL: <https://figshare.com/projects/TP53_mutations_and_TET2_deficiency_cooperate_to_drive_leukemogenesis_and_establish_an_immunosuppressive_environment/167030>

#### Steps to reproduce selected figures:

1.  System Requirements

-   R v4.4

-   Rstudio

-   This software has been tested on Linux Ubuntu 22.04.5

-   Loading the complete dataset occupies approximately 8 GB memory.

2.  Installation

-   download this object in a convenient location on your system.

-   clone the analysis project to your computer using git clone [https://github.com/blaserlab/lapalombella_pu.git](https://github.com/blaserlab/pkc_cxcl8.git)

-   open the R project

-   a list of the packages required for the project can be found in library_catalogs/blas02_lapalombella_pu.tsv. Filter for packages with status == "active". Install these packages and their dependencies.

-   install custom packages from our R Universe repository using these commands:

```         
install.packages('blaseRtools', repos = c('https://blaserlab.r-universe.dev/','https://cloud.r-project.org/'))

install.packages('blaseRtemplates', repos = c('https://blaserlab.r-universe.dev/,'https://cloud.r-project.org/'))

install.packages('blaseRdata', repos = c('https://blaserlab.r-universe.dev/','https://cloud.r-project.org/'))
```

-   source R/dependencies.R (the final line in that file must be edited to point to the directory containing the data package)

-   source R/configs.R (the file paths defining the output variables should be customized for your system)

-   see the named files in R/ to reproduce specific figures from the manuscript

-   typical time required for the first installation and data loading is approximately 15 minutes. This excludes the time required to download the data package.
