\name{NetGSA}
\alias{NetGSA}
\title{Network-based Gene Set Analysis}
\description{Tests the significance of pre-defined sets of genes (pathways) with respect to an outcome variable, such as the condition indicator (e.g. cancer vs. normal, etc.), based on the underlying biological networks.}
\usage{
NetGSA(A, x, group, pathways, lklMethod = c("REML", "ML", "HE", "REHE"), 
       sampling=FALSE, sample_n = NULL, sample_p = NULL, minsize=5, 
       eta = 0.1, lim4kappa = 500)
}
\arguments{
  \item{A}{A list of weighted adjacency matrices.}
  \item{x}{The \eqn{p \times n}{p x n} data matrix with rows referring to genes and columns to samples. It is very important that the adjacency matrices \code{A} share the same rownames as the data matrix \code{x}.}
  \item{group}{Vector of class indicators of length \eqn{n}.}
  \item{pathways}{The npath by \eqn{p} indicator matrix for pathways.}
  \item{lklMethod}{Method used for variance component calculation: options are \code{ML} (maximum likelihood), \code{REML} (restricted maximum likelihood), \code{HE} (Haseman-Elston regression) or \code{REHE} (restricted Haseman-Elston regression). See details.}
  \item{sampling}{(Logical) whether to subsample the observations and/or variables. See details.}
  \item{sample_n}{The ratio for subsampling the observations if \code{sampling=TRUE}.}
  \item{sample_p}{The ratio for subsampling the variables if \code{sampling=TRUE}.}
  \item{minsize}{Minimum number of genes in pathways to be considered.}
  \item{eta}{Approximation limit for the Influence matrix. See 'Details'.}
  \item{lim4kappa}{Limit for condition number (used to adjust \code{eta}). See 'Details'.}
}
\details{
The function \code{NetGSA} carries out a Network-based Gene Set Analysis, using the method described in Shojaie and Michailidis (2009) and Shojaie and Michailidis (2010). It can be used for gene set (pathway) enrichment analysis where the data come from \eqn{K} heterogeneous conditions, where \eqn{K}{K=1,2}, or more. NetGSA differs from Gene Set Analysis (Efron and Tibshirani, 2007) in that it incorporates the underlying biological networks. Therefore, when the networks encoded in \code{A} are empty, one should instead consider alternative approaches such as Gene Set Analysis (Efron and Tibshirani, 2007).

The NetGSA method is formulated in terms of a mixed linear model. Let \eqn{X} represent the rearrangement of data \code{x} into an \eqn{np \times 1}{np x 1} column vector.  
\deqn{X=\Psi \beta + \Pi \gamma + \epsilon}
where \eqn{\beta} is the vector of fixed effects, \eqn{\gamma} and \eqn{\epsilon} are random effects and random errors, respectively. The underlying biological networks are encoded in the weighted adjacency matrices, which determine the influence matrix under each condition. The influence matrices further determine the design matrices \eqn{\Psi} and \eqn{\Pi} in the mixed linear model. Formally, the influence matrix under each condition represents the effect of each gene on all the other genes in the network and is calculated from the adjacency matrix (\code{A[[k]]} for the \eqn{k}-th condition). A small value of \code{eta} is used to make sure that the influence matrices are well-conditioned (i.e. their condition numbers are bounded by \code{lim4kappa}.)

The problem is then to test the null hypothesis \eqn{\ell\beta = 0}{l\beta = 0} against the alternative \eqn{\ell\beta \neq 0}{l\beta != 0}, where \eqn{\ell}{l} is a contrast vector, optimally defined through the underlying networks. 
For a one-sample or two-sample test, the test statistic \eqn{T} for each gene set has approximately a t-distribution under the null, whose degrees of freedom are estimated using the Satterthwaite approximation method. When analyzing complex experiments involving multiple conditions, often multiple contrast vectors of interest are considered for a specific subnetwork. Alternatively, one can combine the contrast vectors into a contrast matrix \eqn{L}. A different test statistic \eqn{F} will be used. Under the null, \eqn{F} has an F-distribution, whose degrees of freedom are calculated based on the contrast matrix \eqn{L} as well as variances of \eqn{\gamma} and \eqn{\epsilon}. The fixed effects \eqn{\beta} are estimated by generalized least squares, and the estimate depends on estimated variance components of \eqn{\gamma} and \eqn{\epsilon}. 

Estimation of the variance components (\eqn{\sigma^2_{\epsilon}} and \eqn{\sigma^2_{\gamma}}) can be done in several different ways after profiling out \eqn{\sigma^2_{\epsilon}}, including \code{REML/ML} which uses Newton's method or \code{HE/REHE} which is based on the Haseman-Elston regression method. The latter notes the fact that \eqn{Var(X)=\sigma^2_{\gamma}\Pi*\Pi' + \sigma^2_{\epsilon}I}, and uses an ordinary least squares to solve for the unknown coefficients after vectorizing both sides. In particular, \code{REHE} uses nonnegative least squares for the regression and therefore ensures nonnegative estimate of the variance components. Due to the simple formulation, \code{HE/REHE} also allows subsampling with respect to both the samples and the variables, and is recommended especially when the problem is large (i.e. large \eqn{p} and/or large \eqn{n}).

The pathway membership information is stored in \code{pathways}, which should be a matrix of \eqn{npath} x \eqn{p}. See \code{\link{prepareAdjacencyMatrix}} for details on how to prepare a suitable pathway membership object.

This function can deal with both directed and undirected networks, which are specified via the option \code{directed}. Note \code{NetGSA} uses slightly different procedures to calculate the influence matrices for directed and undirected networks. 
In either case, the user can still apply \code{NetGSA} if only partial information on the adjacency matrices is available. The functions \code{netEst.undir} and \code{netEst.dir} provide details on how to estimate the weighted adjacency matrices from data based on available network information. 
}
\value{
A list with components
\item{results}{A data frame with pathway names, pathway sizes, p-values and false discovery rate corrected q-values, and test statistic for all pathways.}
\item{beta}{Vector of fixed effects of length \eqn{2p}, of which the first half is for condition 1 and the second half for condition 2.}
\item{s2.epsilon}{Variance of the random errors \eqn{\epsilon}.}
\item{s2.gamma}{Variance of the random effects \eqn{\gamma}.}
\item{graph}{List of components needed in \code{plot.NetGSA}.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. \url{https://doi.org/10.1093/bioinformatics/btw410}

Shojaie, A., & Michailidis, G. (2010). Network enrichment analysis in complex experiments. Statistical applications in genetics and molecular biology, 9(1), Article 22. \url{http://www.ncbi.nlm.nih.gov/pubmed/20597848}.

Shojaie, A., & Michailidis, G. (2009). Analysis of gene sets based on the underlying regulatory network. Journal of Computational Biology, 16(3), 407-426. \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131840/}

}
\author{
Ali Shojaie and Jing Ma
}
\seealso{
\code{\link{prepareAdjacencyMatrix}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
library(glassoFast)
library(glmnet)
library(igraph)

set.seed(1)
data(breastcancer2012)


## -----------------Undirected networks-----------------
## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames <- unique(c(pathways[[24]], pathways[[52]]))
p <- length(genenames)
sx <- x[match(genenames, rownames(x)),]
if (sum(is.na(rownames(sx)))>0){
  sx <- sx[-which(is.na(rownames(sx))),]
}
file_e <- system.file("extdata", "edgelist.txt", package = "netgsa")
out <- prepareAdjacencyMatrix(sx, group, pathways, FALSE, file_e, NULL)
sx <- sx[match(colnames(out$B), rownames(sx)),]

ncond <- length(unique(group))
Amat <- vector("list",ncond)

## -- Not run --
# for (k in 1:ncond){
#  data_c <- sx[,(group==k)]
#  fitBIC <- bic.netEst.undir(data_c,one=out$Adj,
#                             lambda=seq(1,10)*sqrt(log(p)/ncol(data_c)),eta=0.1)
#  fit <- netEst.undir(data_c,one=out$Adj,
#                      lambda=which.min(fitBIC$BIC)*sqrt(log(p)/ncol(data_c)),eta=0.1)
#  Amat[[k]] <- fit$Adj
# }
# test <- NetGSA(Amat, sx, group, pathways = out$B, lklMethod = 'REHE')

## -------------------Directed networks------------------
## NetGSA also works for directed acyclic graphs (DAGs).
# e.g. the "Adrenergic signaling in cardiomyocytes" pathway from KEGG is a DAG.
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
# for (k in 1:2){
#   data_c <- gx[,which(group==k)]
#   Amat[[k]] <- netEst.dir(data_c, one = Adj)$Adj
# }
# test <- NetGSA(Amat, gx, group, pathways = B, lklMethod = 'REHE')

}
