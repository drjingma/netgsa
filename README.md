
# About

This package carries out Network-based Gene Set Analysis by incorporating external information about interactions among genes, as well as novel interactions learned from data.

## How to install?

You can install it directly from CRAN:
```
install.packages("netgsa",build_vignettes=T)
```

Reference manual is available on [CRAN](https://cran.r-project.org/web/packages/netgsa/index.html). A [vignette](https://cran.r-project.org/web/packages/netgsa/vignettes/netgsa.html) on how to use NetGSA is also available. A development version is available on [GitHub](https://github.com/drjingma/netgsa) and can be installed via the following:
```r
library(devtools)
devtools::install_github("drjingma/netgsa", build_vignettes=T)
```
More details about the method implemented can be found in the original paper [here](http://drjingma.com/papers/ma-netgsa) and a follow-up review paper [here](http://drjingma.com/papers/review). 

## Why should someone use `netgsa`?

 - NetGSA incorporates the rich network information curated in public databases (e.g. KEGG, reactome, etc.) and/or learned from high-throughput sequencing data, thereby gaining power in detecting active genetic/metabolic pathways.

## How does it compare to other methods?

 - NetGSA tests the self-contained null hypothesis and compares the set of genes in a given pathway with itself. 
 - NetGSA allows users to complement potentially misspecified information in public databases with high-throughput sequencing data.
 - NetGSA is particularly powerful in detecting active metabolic pathways from metabolomic data where the pathway size is relatively small and available metabolic network information is sparse.  
 - See more details in our review paper [here](http://drjingma.com/papers/review).

## Notes

 - NetGSA is based on a linear mixed model. Estimation of the variance components in this model can be done via restricted maximum likelihood or via the restricted Haseman-Elston (REHE) regression. See our recent paper [here](http://drjingma.com/papers/REHE) for fast variance components estimation with REHE. 
 - NetGSA is seamlessly integrated with external databases on gene-gene interactions and utilizes Cytoscape for interactive visualization of enriched pathways. Moreover, NetGSA can handle thousands of genes within minutes. See [here](http://drjingma.com/papers/netgsa) for details on how we improved the computation and visualization of NetGSA.  


  [![Travis-CI Build Status](https://travis-ci.org/drjingma/netgsa.svg?branch=master)](https://travis-ci.org/drjingma/netgsa)
