\name{netEst.undir}
\alias{netEst.undir}
\title{
Constrained estimation of undirected networks
}
\description{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty.
}
\usage{
netEst.undir(x, zero = NULL, one = NULL, lambda, rho=NULL, penalize_diag = TRUE, weight = NULL, 
             eta = 0, verbose = FALSE, eps = 1e-08)
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix with rows referring to genes and columns to samples.}
  \item{zero}{(Optional) indices of entries of the weighted adjacency matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere. The matrix must be symmetric.}
  \item{one}{(Optional) indices of entries of the weighted adjacency matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero} and needs to be symmetric.}
  \item{lambda}{(Non-negative) numeric vector representing the regularization parameters for lasso. Can choose best based on BIC using \code{bic.netEst.undir}}
  \item{rho}{(Non-negative) numeric scalar or symmetric \eqn{p \times p}{p x p} matrix representing the regularization parameter for estimating the weights in the inverse covariance matrix. This is the same as \code{rho} in the graphical lasso algorithm \code{glassoFast}.}
  \item{penalize_diag}{Logical. Whether or not to penalize diagonal entries when estimating weighted adjacency matrix. If TRUE a small penalty is used, otherwise no penalty is used.}
  \item{weight}{(Optional) whether to add penalty to known edges. If NULL (default), then the known edges are assumed to be true. If nonzero, then a penalty equal to lambda * weight is added to penalize the known edges to account for possible uncertainty. Only non-negative values are accepted for the weight parameter.  }
  \item{eta}{(Non-negative) a small constant added to the diagonal of the empirical covariance matrix of \code{X} to ensure it is well conditioned. By default, \code{eta} is set to 0.}
  \item{verbose}{Whether to print out information as estimation proceeds. Default = \code{FALSE}.}
  \item{eps}{(Non-negative) numeric scalar indicating the tolerance level for differentiating zero and non-zero edges: entries with magnitude \eqn{<} \code{eps} will be set to 0. }
}
\details{
The function \code{netEst.undir} performs constrained estimation of sparse inverse covariance (concentration) matrices using a lasso (L1) penalty, as described in Ma, Shojaie and Michailidis (2016). Two sets of constraints determine subsets of entries of the inverse covariance matrix that should be exactly zero (the option \code{zero} argument), or should take non-zero values (option \code{one} argument). The remaining entries will be estimated from data.

The arguments \code{one} and/or \code{zero} can come from external knowledge on the 0-1 structure of underlying concentration matrix, such as a list of edges and/or non-edges learned from available databases.

\code{netEst.undir} estimates both the support (0-1 structure) of the concentration matrix, or equivalently, the adjacency matrix of the corresponding Gaussian graphical model, for a given tuning parameter, \code{lambda}; and the concentration matrix with diagonal entries set to 0, or equivalently, the weighted adjacency matrix.
The weighted adjacency matrix is estimated using maximum likelihood based on the estimated support. The parameter \code{rho} controls the amount of regularization used in the maximum likelihood step. A small \code{rho} is recommended, as a large value of \code{rho} may result in too much regularization in the maximum likelihood estimation, thus further penalizing the support of the weighted adjacency matrix.  
Note this function is suitable only for estimating the adjacency matrix of a undirected graph. The \code{weight} parameter allows one to specify whether to penalize the known edges. If known edges obtained from external information contain uncertainty such that some of them are spurious, then it is recommended to use a small positive \code{weight} parameter to select the most probable edges from the collection of known ones. 

This function is closely related to \code{NetGSA}, which requires the weighted adjacency matrix as input. When the user does not have complete information on the weighted adjacency matrix, but has data (\code{x}, not necessarily the same as the \code{x} in \code{NetGSA}) and external information (\code{one} and/or \code{zero}) on the adjacency matrix, then \code{netEst.undir} can be used to estimate the remaining interactions in the adjacency matrix using the data.
Further, when it is anticipated that the adjacency matrices under different conditions are different, and data from different conditions are available, the user needs to run \code{netEst.undir} separately to obtain estimates of the adjacency matrices under each condition.

The algorithm used in \code{netEst.undir} is based on \code{glmnet} and \code{glasso}. Please refer to \code{glmnet} and \code{glasso} for computational details.
}
\value{
A list with components
\item{Adj}{List of weighted adjacency matrices (partial correlations) of dimension \eqn{p \times p}{p x p}, with diagonal entries set to 0. Each element in the list is the weighted adjacency matric corresponding to each value in \code{lambda}. Each element is a matrix that will be used in \code{NetGSA}.}
\item{invcov}{List of estimated inverse covariance matrix of dimension \eqn{p \times p}{p x p}.}
\item{lambda}{List of values of tuning parameters used.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. \url{https://doi.org/10.1093/bioinformatics/btw410}
}
\author{
Jing Ma & Michael Hellstern
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prepareAdjMat}}, \code{\link{bic.netEst.undir}}, \code{\link{glmnet}}
}
\examples{
library(glassoFast)
library(graphite)
library(igraph)

set.seed(1)

## load the data
data(breastcancer2012)

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames <- unique(c(pathways[[24]], pathways[[52]]))
sx <- x[match(genenames, rownames(x)),]
if (sum(is.na(rownames(sx)))>0){
  sx <- sx[-which(is.na(rownames(sx))),]
}
p <- length(genenames)

## zero/one matrices should be based on known non-edges/known edges. Random used as an example
one <- matrix(sample(c(0,1), length(rownames(sx))**2, replace = TRUE, prob = c(0.9, 0.1)), length(rownames(sx)), dimnames = list(rownames(sx), rownames(sx)))

ncond <- length(unique(group))
Amat <- vector("list",ncond)
for (k in 1:ncond){
    data_c <- sx[,(group==k)]
    fitBIC <- bic.netEst.undir(data_c,one=one,
                               lambda=seq(1,10)*sqrt(log(p)/ncol(data_c)),eta=0.1)
    fit <- netEst.undir(data_c,one=one,
                        lambda=which.min(fitBIC$BIC)*sqrt(log(p)/ncol(data_c)),eta=0.1)
    Amat[[k]] <- fit$Adj
}

}
