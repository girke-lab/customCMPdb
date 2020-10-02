# customCMPdb: Integrating Community and Custom Compound Collections

## Introduction
This package serves as a query interface for important community collections of
small molecules, while also allowing users to include custom compound
collections. At the time of writing, the following community databases are included: 

+ [DrugAge](https://genomics.senescence.info/drugs/) [@Barardo2017-xk]
+ [DrugBank](https://www.drugbank.ca/) [@Wishart2018-ap]
+ [CMAP02](https://portals.broadinstitute.org/cmap/) [@Lamb2006-du]
+ [LINCS](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) [@Subramanian2017-fu]


## Installation and Loading
`customCMPdb` is an R/Bioconductor package and can be installed using 
`BiocManager::install()`.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("customCMPdb")
```

To obtain the most recent updates immediately, one can install it directly from 
GitHub as follows.
```r
devtools::install_github("yduan004/customCMPdb", build_vignettes=TRUE)
```

After the package is installed, it can be loaded into an R session as follows.
```r
library(customCMPdb)
```
For detailed description of the package, please refer to the vignette by running
```r
browseVignettes("customCMPdb")
```

## Vignette
The vignette of this package is also available at [here](https://www.bioconductor.org/packages/release/bioc/vignettes/customCMPdb/inst/doc/customCMPdb.html)
