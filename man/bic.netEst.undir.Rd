\name{bic.netEst.undir}
\alias{bic.netEst.undir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian information criterion to select the tuning parameters for \code{netEst.undir}
}
\description{
This function uses the Bayesian information criterion to select the optimal tuning parameters needed in \code{netEst.undir}.
}
\usage{
bic.netEst.undir(x, zero = NULL, one = NULL, lambda, rho = NULL, weight = NULL, 
                 eta = 0, verbose = FALSE, eps = 1e-08)
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix as in \code{netEst.undir}.}
  \item{zero}{(Optional) indices of entries of the matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere. The matrix must be symmetric.}
  \item{one}{(Optional) indices of entries of the matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero} and needs to be symmetric.}
  \item{lambda}{(Non-negative) user-supplied lambda sequence.}
  \item{rho}{(Non-negative) numeric scalar representing the regularization parameter for estimating the weights in the inverse covariance matrix. This is the same as \code{rho} in the graphial lasso algorithm \code{glassoFast}.}
  \item{weight}{(Optional) whether to add penalty to known edges. If NULL (default), then the known edges are assumed to be true. If nonzero, then a penalty equal to \code{lambda * weight} is added to penalize the known edges to account for possible uncertainty. Only non-negative values are accepted for the weight parameter. }
  \item{eta}{(Non-negative) a small constant added to the diagonal of the empirical covariance matrix of \code{X} to ensure it is well conditioned. By default, \code{eta} is set to 0.}
  \item{verbose}{Whether to print out information as estimation proceeds. Default=\code{FALSE}.}
  \item{eps}{Numeric scalar \eqn{>=0}, indicating the tolerance level for differentiating zero and non-zero edges: entries \eqn{<} \code{eps} will be set to 0. }
}
\details{
Let \eqn{\hat\Sigma}{\hat\Sigma} represent the empirical covariance matrix of data \code{x}. For a given \eqn{\lambda}, denote the estimated inverse covariance matrix by \eqn{\hat\Omega_{\lambda}}. the Bayesian information criterion (BIC) is defined as
\deqn{trace(\hat\Sigma \hat\Omega_{\lambda}) - \log \det (\hat\Omega_{\lambda}) + \frac{\log n}{n} \cdot df,}
where \eqn{df} represents the degrees of freedom in the selected model and can be estimated via the number of edges in \eqn{\hat\Omega_{\lambda}}. The optimal tuning parameter is selected as the one that minimizes the BIC over the range of \code{lambda}.

Note when the penalty parameter \code{lambda} is too large, the estimated adjacency matrix may be zero. The function will thus return a warning message. 
}
\value{
\item{lambda}{The values of \code{lambda} used.}
\item{weight}{The values of \code{weight} used.}
\item{BIC}{If \code{weight=NULL}, then a numeric vector of the same length as \code{lambda} with the corresponding BIC. If \code{weight} is a vector, then a matrix of size \code{length(lambda)} by \code{length(weight)} with the corresponding BIC.}
\item{df}{The degrees of freedom corresponding to each BIC.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. \url{https://doi.org/10.1093/bioinformatics/btw410}
}
\author{
Jing Ma
}

\seealso{
\code{\link{netEst.undir}}
}
\examples{
library(glassoFast)
library(graph)
library(graphite)
library(igraph)

set.seed(1)

## load the data
data(breastcancer2012)

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames <- unique(c(pathways[[24]], pathways[[52]]))
p <- length(genenames)
sx <- x[match(genenames, rownames(x)),]
if (sum(is.na(rownames(sx)))>0){
  sx <- sx[-which(is.na(rownames(sx))),]
}
file_e <- paste0(path.package("netgsa"), "/doc/edgelist.txt")
out <- prepareAdjacencyMatrix(sx, group, pathways,FALSE,file_e, NULL)
sx <- sx[match(colnames(out$B), rownames(sx)),]

ncond <- length(unique(group))
Amat <- vector("list",ncond)
for (k in 1:ncond){
  data_c <- sx[,(group==k)]
  fitBIC <- bic.netEst.undir(data_c,one=out$Adj,
                             lambda=seq(1,10)*sqrt(log(p)/ncol(data_c)),eta=0.1)
  fit <- netEst.undir(data_c,one=out$Adj,
                      lambda=which.min(fitBIC$BIC)*sqrt(log(p)/ncol(data_c)),eta=0.1)
  Amat[[k]] <- fit$Adj
}

}