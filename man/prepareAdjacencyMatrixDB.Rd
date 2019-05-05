\name{prepareAdjacencyMatrixDB}
\alias{prepareAdjacencyMatrixDB}
\title{Construct adjacency matrices from graphite databases and/or user provided network information}
\description{
Read the network information from any of the graphite databases specified by the user and construct the adjacency matrices needed for NetGSA.
}
\usage{
prepareAdjacencyMatrixDB(x, genes, group, databases = NULL, 
       file_e=NULL, file_ne=NULL, lambda_c=1, eta=0.5)
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix with rows referring to genes and columns to samples. Rownames should have gene ID types appended to them. The id and gene number should be separated by a colon. E.g. "ENTREZID:127550".}
  \item{genes}{Named character vector of genes. Names are ID types. Must contain same genes are rownames of \code{x}. E.g. c("ENTREZID" = "127550").}
  \item{group}{Vector of class indicators of length \eqn{n}.}
  \item{databases}{Character vector of graphite databases you wish to search for edges.}
  \item{file_e}{The name of the file which the list of edges is to read from. This should be a .csv file, with one edge in a line, the two vertices separated by a delimiter. Each observation is assumed to be a directed edge. Thus if you want an undirected edge input two separate edges.}
  \item{file_ne}{The name of the file which the list of negative edges is to read from. The edges in this file are negative in the sense that the corresponding vertices are not connected. This should be a .csv file, with one negative edge in a line, the two vertices separated by a delimiter. Each observation is assumed to be a directed edge. Thus if you want a negative undirected edge input two separate negative edges.}
  \item{lambda_c}{(Non-negative) a constant multiplied to the tuning parameter \code{lambda} needed for estimating the edge weights. By default, \code{lambda_c} is set to 1. See \code{\link{netEst.undir}} and \code{\link{netEst.dir}} for more details.}
  \item{eta}{(Non-negative) a small constant needed for estimating the edge weights. By default, \code{eta} is set to 0.5. See \code{\link{netEst.undir}} for more details.}
}
\details{
The function \code{prepareAdjacencyMatrixDB} accepts both network information from user specified sources as well as a list of graphite databases to search for edges in. \code{prepareAdjacencyMatrixDB} calculates the 0-1 adjacency matrices and then runs \code{\link{netEst.undir}} or \code{\link{netEst.dir}} if the graph is undirected or directed. Note \code{prepareAdjacencyMatrixDB} automatically detects directed graphs and rearranges both the 0-1 adjacency matrix and data matrix to the correct order. Note if the dimension of the problem, or equivalently the total number of unique genes across all pathways, is large, \code{prepareAdjacencyMatrixDB} may be a slow.

\code{prepareAdjacencyMatrixDB} returns the list of weighted adjacency matrices to be used for input in \code{NetGSA}.

}
\value{
A list with components
\item{Adj}{A list of weighted adjacency matrices estimated from either \code{netEst.undir} or \code{netEst.dir}. One weighted adjacency matrix will be returned for each condition in group. That is \code{length(Adj) = length(unique(group))}}
\item{invcov}{A list of inverse covariance matrices estimated from either \code{netEst.undir} or \code{netEst.dir}. One inverse covariance matrix will be returned for each condition in group. That is \code{length(invcov) = length(unique(group))}}
\item{lambda}{A list of values of tuning parametes used for each condition in \code{group}}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}
\author{
Mike Hellstern
}
\seealso{
\code{\link{NetGSA}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
library(glassoFast)
library(igraph)

set.seed(1)

## load the data
data("breastcancer2012")

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames <- unique(c(pathways[[24]], pathways[[52]]))
p <- length(genenames)
sx <- x[rownames(x) %in% genenames,]
genes <- setNames(rownames(sx), rep("ENTREZID", length(rownames(sx))))
rownames(sx) <- paste0("ENTREZID:", rownames(sx))

out1 <- prepareAdjacencyMatrixDB(sx, genes, group, databases = c("kegg", "reactome", "biocarta"))

}
