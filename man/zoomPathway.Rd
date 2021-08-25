\name{zoomPathway}
\alias{zoomPathway}
\title{Zoom in on pathway in igraph}

\description{
Plots the gene to gene interactions for a given pathway in igraph.}

\usage{
zoomPathway(x, pway, graph_layout = NULL)
}

\arguments{
  \item{x}{A NetGSA object returned from calling \code{NetGSA()}}
  \item{pway}{Name of pathway to plot}
  \item{graph_layout}{(Optional) Layout function to pass to igraph plots. This function should only take one parameter (an igraph object). For example one might create a custom layout by setting the spring.length and spring.constant with: \code{my_layout <- function(graph) layout_with_graphopt(graph = graph, spring.length = 1000, spring.constant = 0.00004)}}
}
  
\details{
Generates igraph plot for gene to gene interactions for a given pathway
}

\value{No return value, called to zoom to pathway}



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
zoomPathway(out_cluster, "ErbB signaling pathway", my_layout)
}
}
