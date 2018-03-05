\name{netEst.dir}
\alias{netEst.dir}
\title{
Constrained estimation of directed networks
}
\description{
Estimates a directed network using a lasso (L1) penalty.
}
\usage{
netEst.dir(x, zero = NULL, one = NULL, lambda, verbose = FALSE, eps = 1e-08)
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix.}
  \item{zero}{(Optional) indices of entries of the matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere.}
  \item{one}{(Optional) indices of entries of the matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero}.}
  \item{lambda}{(Non-negative) numeric scalar or a vector of length \eqn{p-1} representing the regularization parameters for nodewise lasso. If \code{lambda} is a scalar, the same penalty will be used for all \eqn{p-1} lasso regressions. By default (\code{lambda=NULL}), the vector of \code{lambda} is defined as
  \deqn{\lambda_j(\alpha) = 2 n^{-1/2} Z^*_{\frac{\alpha}{2p(j-1)}}, \quad j=2,\ldots,p.}
  Here \eqn{Z^*_q} represents the \eqn{(1-q)}-th quantile of the standard normal distribution and \eqn{\alpha} is a positive constant between 0 and 1. See Shojaie and Michailidis (2010a) for details on the choice of tuning parameters.
  }
  \item{verbose}{Whether to print out information as estimation proceeds. Default = \code{FALSE}.}
  \item{eps}{(Non-negative) numeric scalar indicating the tolerance level for differentiating zero and non-zero edges: entries with magnitude \eqn{<} \code{eps} will be set to 0. }
}
\details{
The function \code{netEst.dir} performs constrained estimation of a directed network using a lasso (L1) penalty, as described in Shojaie and Michailidis (2010a). Two sets of constraints determine subsets of entries of the weighted adjacency matrix that should be exactly zero (the option \code{zero} argument), or should take non-zero values (option \code{one} argument). The remaining entries will be estimated from data.

The arguments \code{one} and/or \code{zero} can come from external knowledge on the 0-1 structure of underlying network, such as a list of edges and/or non-edges learned frm available databases. Then the function \code{prepareAdjacencyMatrix} can be used to first construct \code{one} and/or \code{zero}.

In this function, it is assumed that the columns of \eqn{x} are ordered accoring to a correct (Wald) causal order, such that no \eqn{x_j} is a parent of \eqn{x_k} (\eqn{k \le j}). Given the causal ordering of nodes, the resulting adjacency matrix is lower triangular (see Shojaie & Michailidis, 2010b). Thus, only lower triangular parts of \code{zero} and \code{one} are used in this function. For this reason, it is important that both of these matrices are also ordered according to the causal order of the nodes in \eqn{x}. To estimate the network, first each node is regressed on the known edges (\code{one}). The reisdual obtained from this regression is then used to find the additional edges, among the nodes that could potentially interact with the given node (those not in \code{zero}).
		
This function is closely related to \code{NetGSA}, which requires the weighted adjacency matrix as input. When the user does not have complete information on the weighted adjacency matrix, but has data (not necessarily the same as the \code{x} in \code{NetGSA}) and external information (\code{one} and/or \code{zero}) on the adjacency matrix, then \code{netEst.dir} can be used to estimate the remaining interactions in the adjacency matrix using the data.
Further, when it is anticipated that the adjacency matrices under different conditions are different, and data from different conditions are available, the user needs to run \code{netEst.dir} separately to obtain estimates of the adjacency matrices under each condition.

The algorithm used in \code{netEst.undir} is based on \code{glmnet}. Please refer to \code{glmnet} for computational details.
}
\value{
A list with components
\item{Adj}{The weighted adjacency matrix of dimension \eqn{p \times p}{p x p}. This is the matrix that will be used in \code{NetGSA}.}
\item{infmat}{The influence matrix of dimension \eqn{p \times p}{p x p}.}
\item{lambda}{The values of tuning parameters used.}
}
\references{
Shojaie, A., & Michailidis, G. (2010a). Penalized likelihood methods for estimation of sparse high-dimensional directed acyclic graphs. Biometrika 97(3), 519-538. \url{http://biomet.oxfordjournals.org/content/97/3/519.short}

Shojaie, A., & Michailidis, G. (2010b). Network enrichment analysis in complex experiments. Statistical applications in genetics and molecular biology, 9(1), Article 22. \url{http://www.ncbi.nlm.nih.gov/pubmed/20597848}.

Shojaie, A., & Michailidis, G. (2009). Analysis of gene sets based on the underlying regulatory network. Journal of Computational Biology, 16(3), 407-426. \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131840/}
}
\author{
Ali Shojaie
}

\seealso{
\code{\link{prepareAdjacencyMatrix}}, \code{\link{glmnet}}
}
\examples{
library(glmnet)
library(graph)
library(graphite)
library(igraph)

set.seed(1)

## load the data
data(breastcancer2012)

print(is_dag(g))

genenames <- V(g)$name
p <- length(genenames)

# reorder the variables and get the adjacency matrix
reOrder <- topo_sort(g,"in")
Adj <- as.matrix(get.adjacency(g))
Adj <- Adj[reOrder,reOrder]

B <- matrix(rep(1,p),nrow=1)
rownames(B) <- "Adrenergic signaling in cardiomyocytes"
colnames(B) <- rownames(Adj)
gx <- x[match(rownames(Adj), rownames(x)),]

Amat <- vector("list", 2)
for (k in 1:2){
  data_c <- gx[,which(group==k)]
  Amat[[k]] <- netEst.dir(data_c, one = Adj)$Adj
}
}
