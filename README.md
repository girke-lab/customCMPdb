# compoundCollectionData

## Introduction
This package contains annotation and structure datasets for compounds in 
[DrugAge](https://genomics.senescence.info/drugs/), 
[DrugBank](https://www.drugbank.ca/), [CMAP02](https://portals.broadinstitute.org/cmap/) 
and [LINCS](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) databases.
The actual datasets are stored in `AnnotationHub`. 

## Installation and Loading
`compoundCollectionData` is a R/Bioconductor package and can be installed using 
`BiocManager::install()`.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("compoundCollectionData")
```

To obtain the most recent updates immediately, one can install it directly from 
GitHub as follows.
```r
devtools::install_github("yduan004/compoundCollectionData", build_vignettes=TRUE)
```

After the package is installed, it can be loaded into an R session as follows.
```r
library(compoundCollectionData)
```
For detailed description of the package, please refer to the vignette by running
```r
browseVignettes("compoundCollectionData")
```

## Vignette
The vignette of this package is also available at [here](https://www.bioconductor.org/packages/release/bioc/vignettes/compoundCollectionData/inst/doc/compoundCollectionData.html)
