\name{formatPathways}
\alias{formatPathways}
\title{Format cytoscape nested networks}

\description{
Format cytoscape nested networks using preset NetGSA format}

\usage{
formatPathways(x, pways, graph_layout = NULL)
}

\arguments{
  \item{x}{A NetGSA object returned from calling \code{NetGSA()}}
  \item{pways}{Character vector of pathways to format}
  \item{graph_layout}{(Optional) Layout to pass to plots. Must be a string for Cytoscape which will be passed to \code{RCy3::layoutNetwork}.}
}
  
\details{
Loads gene testing data into each pathway. Genes are tested using an F-test if there are 2 or more conditions or a two-sided one-class t-test against the null hypothesis of mean = 0 if there is only one condition. FDR corrected q-values are mapped to the color of the node. The scale ranges from 0 to 1 with red represents q-values of 0 and white representing q-values of 1. Loaded data includes: p-value from the F-test/t-test (pval), FDR corrected q-value (pFdr), test statistic from the F-test/t-test (teststat).

Custom formatting can be applied using the cytoscape GUI or the RCy3 pacakge.
}



\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Michael Hellstern
}

\seealso{
\code{\link{plot.NetGSA}}
}

\examples{
\donttest{
library(igraph)

## load the data
data("breastcancer2012")

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames    <- unique(c(pathways[[24]], pathways[[52]]))
sx           <- x[match(rownames(x), genenames, nomatch = 0L) > 0L,]

db_edges       <- obtainEdgeList(rownames(sx), databases = c("kegg", "reactome", "biocarta"))
adj_cluster    <- prepareAdjMat(sx, group, databases = db_edges, cluster = TRUE)
out_cluster    <- NetGSA(adj_cluster[["Adj"]], sx, group, pathways_mat[c(24, 52), rownames(sx)], lklMethod = "REHE", sampling = FALSE, sample_n = 0.1, sample_p = 0.1)

### Cytoscape closed or open
plot(out_cluster) 
my_layout <- function(graph) layout_with_graphopt(graph = graph, spring.length = 1000, spring.constant = 0.00004)
formatPathways(out_cluster, "ErbB signaling pathway")
}
}
