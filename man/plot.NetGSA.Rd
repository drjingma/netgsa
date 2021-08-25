\name{plot.NetGSA}
\alias{plot.NetGSA}
\title{Generates NetGSA plots}

\description{
Generates network plots in Cytoscape and igraph}

\usage{
\method{plot}{NetGSA}(x, graph_layout = NULL, rescale_node = c(2,10), rescale_label = c(0.5,0.6), ...)
}

\arguments{
  \item{x}{An object of class "NetGSA" returned from calling \code{NetGSA()}}
  \item{graph_layout}{(Optional) Layout to pass to plots. Either a function for igraph plots (when Cytoscape not open) or a string for Cytoscape. The igraph function should only take one parameter (an igraph object). See \code{igraph::layout_} for more details. For example one might create a custom layout by setting the spring.length and spring.constant with: \code{my_layout <- function(graph) layout_with_graphopt(graph = graph, spring.length = 1000, spring.constant = 0.00004)}. The string for Cytoscape will be passed to \code{RCy3::layoutNetwork}.}
  \item{rescale_node}{(Optional) Node size rescaling to pass to igraph plots. Must be a vector of length 2 with the first element being the minimum node size and the second being the maximum.}
  \item{rescale_label}{(Optional) Label size rescaling to pass to igraph plots. Must be a vector of length 2 with the first element being the minimum node size and the second being the maximum.}
  \item{...}{Other arguments not used by this method}
  }
  
\details{
One of two options can occur.

(1) If Cytoscape is open on the user's computer, a nested network will be created. The main network is the interactions between pathways. In this graph, there is one node for each pathway. An edge is drawn between pathways if there is at least one edge between genes of each pathway. That is if gene A is in pathway 1 and gene B is in pathway 2, pathway 1 and pathway 2 will have an edge if gene A and gene B have an edge. Note self-edges are not drawn. The value of the test statistic is mapped to node color. Large negative values of the test statistic are orange, values around 0 are white and large positive values are blue. FDR corrected q-values are mapped to the border color of the node. The scale ranges from 0 to 1 with red representing q-values of 0 and white representing q-values of 1. Pathway size is mapped to node size so pathways with more genes are larger. Each pathway node is also linked to its network of genes so the user can see individual gene interactions within a pathway. These can be accessed by right clicking the node -> Nested Networks -> Go To Nested Network. Alternatively, the corresponding nested network has the same name as the pathway so the user can click on the network directly in the Control Panel/Network menu. It is important to note that \code{plot.NetGSA} generates default plots and loads in data into Cytoscape, but the user can customize the plots however they like using RCy3 or the Cytoscape GUI directly.

To save time, the nested networks are not formatted. One can apply NetGSA's formatting using \code{formatPathways}

For custom formatting, the node data that is loaded into Cytoscape includes the pathway results from NetGSA: Pathway size (pSize), p-value (pval), FDR corrected q-value (pFDR), test statistic (teststat) and pathway name. The edge data loaded into Cytoscape is: total number of edges between two pathways (weight). For example weight of 10 between pathway 1 and pathway 2 means there are 10 edges between the genes of pathway 1 and the genes of pathway 2.

There are two R plots also generated. The first is the legend for Cytoscape. The legend shows the mapping for node color (test statistic) and node border color (FDR corrected q-value). This is generated in R because there does not seem to be a reliable way to plot the legend for the main network (interactions between pathways). The second plot is a plot of the main network created in igraph. It mimics the Cytoscape plot as closely as possible. NetGSA exports the x and y coordinates of the nodes in the Cytoscape layout and uses them in the igraph layout. Custom layouts can be passed to this using the \code{graph_layout} argument. The user can also zoom-in on individual pathways in igraph using the \code{zoomPathway} function.

(2) If Cytoscape is not open, the \code{igraph::rglplot} function is used to plot the main network (interactions between pathways). The default layout used is \code{layout_on_sphere}, but custom layouts can be specified with the \code{graph_layout} argument. The other plot generated is the legend since it is difficult to plot on \code{rglplot}.
}

\value{
No return value, plots NetGSA object
}

\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Michael Hellstern
}

\seealso{
\code{\link{NetGSA}}
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
}
}
