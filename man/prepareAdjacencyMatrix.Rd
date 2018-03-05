\name{prepareAdjacencyMatrix}
\alias{prepareAdjacencyMatrix}
\title{Construct adjacency matrices from existing databases or user provided network information}
\description{
Read the network information from KEGG or specified by the user and construct the adjacency matrices needed for NetGSA.
}
\usage{
prepareAdjacencyMatrix(x, group, pathways, import_from_kegg=FALSE, 
       file_e=c(NA, file_e), file_ne=c(NULL, file_ne, NA), 
       estimate_network=FALSE, lambda_c=1, eta=0.5, minsize=5, fileEncoding="")
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix with rows referring to genes and columns to samples.}
  \item{group}{Vector of class indicators of length \eqn{n}.}
  \item{pathways}{A list of pathways from \code{preparePathways}.}
  \item{import_from_kegg}{(Logical) whether to import network information from KEGG. }
  \item{file_e}{The name of the file which the list of edges is to read from. This should be a .txt file, with one edge in a line, the two vertices separated by a delimiter. The third column indicates the direction of the edge: directed or undirected.}
  \item{file_ne}{The name of the file which the list of negative edges is to read from. The edges in this file are negative in the sense that the corresponding vertices are not connected. This should be a .txt file, with one negative edge in a line, the two vertices separated by a delimiter. The third column indicates the direction of the negative edge: directed or undirected. If \code{NA}, the input edge list from \code{file_e} will be treated as complete network information and edge weights will be estimated only for the input edge list. If \code{NULL}, the input edge list from \code{file_e} will be treated as partial information and a network estimation procedure will be used to refit the network topology and edge weights. See details.}
  \item{estimate_network}{(Logical) whether to estimate the weighted adjacency matrices. Default is \code{FALSE}. Users are recommended to use the returned 0-1 adjacnecy matrix and the function \code{\link{netEst.undir}} and \code{\link{netEst.dir}} to estimate the weighted adjacency matrices.}
  \item{lambda_c}{(Non-negative) a constant multiplied to the tuning parameter \code{lambda} needed for estimating the edge weights. By default, \code{lambda_c} is set to 1. See \code{\link{netEst.undir}} and \code{\link{netEst.dir}} for more details.}
  \item{eta}{(Non-negative) a small constant needed for estimating the edge weights. By default, \code{eta} is set to 0.5. See \code{\link{netEst.undir}} for more details.}
  \item{minsize}{Minimum number of genes in pathways to be considered.}
  \item{fileEncoding}{Character string: if non-empty declares the encoding used on a file (not a connection) so the character data can be re-encoded. See the ‘Encoding’ section of the help for \href{http://127.0.0.1:25551/help/library/utils/help/file}{file}, the ‘R Data Import/Export Manual’ and ‘Note’.}
}
\details{
The function \code{prepareAdjacencyMatrix} accepts both network information from user specified sources and KEGG, and return 0-1 adjacency matrices needed for \code{\link{netEst.undir}} and \code{\link{netEst.dir}}. If \code{estimate_network=TRUE}, the function also returns a list of weighted adjacency matrices that can be directly used in \code{NetGSA}. Note if the dimension of the problem, or equivalently the total number of unique genes across all pathways, is large, \code{prepareAdjacencyMatrix} with \code{estimate_network=TRUE} may be a bit slow.

If \code{file_ne=NA}, the input edge list from \code{file_e} will be treated as complete network information and edge weights will be estimated only for the input edge list. If \code{file_ne=NULL}, the input edge list from \code{file_e} will be treated as partial information and a network estimation procedure will be used to refit the network topology and edge weights. When importing network information from KEGG, the network information can also be treated as complete if \code{file_ne=NA} or incomplete if \code{file_ne=NULL}. The information in \code{file_e} and \code{file_ne} should be compatible in the sense that the same edge should not appear in both files. An error will be reported if incompatibility occurs. 

When \code{estimate_network=FALSE}, the returned 0-1 adjacency matrices (for edges or non-edges) can be used as input in \code{netEst.undir} or \code{netEst.dir} to estimate the complete (and weighted) adjacency matrices under a variety of tuning parameters. The user can choose the set of tuning parameters that achieve the best balance in model fit and model complexity. 
When \code{estimate_network=TRUE}, \code{prepareAdjacencyMatrix} also returns the weighted adjacency matrices suitable for input in \code{NetGSA} using a fixed set of tuning parameters, as specified via \code{lambda_c} and \code{eta}. 

}
\value{
A list with components
\item{Amat}{A list of weighted adjacency matrices.}
\item{Adj}{A list of 0-1 adjacency matrices corresponding to the edges.}
\item{Zero_Adj}{A list of 0-1 adjacency matrices corresponding to the negative edges.}
\item{B}{The npath by \eqn{p} indicator matrix for pathways.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. \url{https://doi.org/10.1093/bioinformatics/btw410}
}
\author{
Jing Ma
}
\seealso{
\code{\link{NetGSA}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
library(glassoFast)
library(graph)
library(graphite)
library(igraph)

set.seed(1)

## load the data
data("breastcancer2012")

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames <- unique(c(pathways[[24]], pathways[[52]]))
p <- length(genenames)
sx <- x[match(genenames, rownames(x)),]
if (sum(is.na(rownames(sx)))>0){
  sx <- sx[-which(is.na(rownames(sx))),]
}
# compare the resulting matrices when file_ne=NULL vs file_ne=NA
file_e <- paste0(path.package("netgsa"), "/doc/edgelist.txt")
out1 <- prepareAdjacencyMatrix(sx, group, pathways,FALSE,file_e, NULL)
out2 <- prepareAdjacencyMatrix(sx, group, pathways,FALSE,file_e, NA)

}
