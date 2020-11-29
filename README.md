# netgsa
Network-based Gene Set Analysis

This package carries out Network-based Gene Set Analysis by incorporating external information about interactions among genes, as well as novel interactions learned from data.


## **Package Installation**
You can install it directly from GitHub through `devtools`:

```
library(devtools)
devtools::install_github("mikehellstern/netgsa", build_vignettes=T)
```

## Updates

The most recent implementation has optimized the NetGSA computation in the following aspects:

* Variance component estimation: Residuals are needed to estimate the variance components. This is done directly without evaluating the fixed effect coefficients. 
* Contrast vector/matrix: This is done more efficiently by leveraging the `lapply` function; products of the contrast vectors are first computed and reused to calculate the degrees of freedom and test statistics. 
* In the main function `NetGSA`, the default input `A` is a list of adjacency matrices across the tested groups. For each group, we assume that its adjacency matrix is again coded as a list of smaller matrices (or in the extreme case, one matrix of size `p`). We do not assume the adjancency matrices across groups to have the same block diagonal structure. For this particular structure, I removed the check on variable compatibility between the adjacency matrices and the input data, but this should be added later.  
* Note the use of `adj2inf` should not change if we have block diagonal adj matrix, because the list of eigenvalues remain the same. 
* The fixed effect coefficients `beta` is currently an output from the main function `NetGSA`, but do we need it? Is there a better way of estimating beta given the block diagonal structure of D?
* Although it was mentioned in the notes that we may not need to assemble the entire adjacency matrix to get the test statistic for a given pathway, in practice we may be working with more pathways. Would it be better if we assemble the entire adjacency matrix anyway to avoid repeatedly subsetting variables? In either case, before running `NetGSA`, we should make sure to filter variables in the input data matrix to keep only those that belong to at least one tested pathway.


## **References**
**Ma, Jing, Shojaie, Ali and Michailidis, George.** (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics https://doi.org/10.1093/bioinformatics/btw410


  [![Travis-CI Build Status](https://travis-ci.org/drjingma/netgsa.svg?branch=master)](https://travis-ci.org/drjingma/netgsa)
