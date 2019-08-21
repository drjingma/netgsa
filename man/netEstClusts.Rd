\name{netEstClusts}

\alias{netEstClusts}

\title{Estimates network for each cluster}

\description{
Estimates network using \code{netEst.dir} or \code{netEst.undir} for each cluster. This is a helper function in \code{prepareAdjacencyMatrixDB} and should not be called by the user.
}


\usage{
netEstClusts(grp, X, group, net_info, n, p, lambda_c, eta, net_clusters)
}


\arguments{
  \item{grp}{Specific group to analyze e.g. condition 1, 2, etc. Same type as \code{group}}
  \item{X}{\eqn{p \times n}{p x n} data matrix}
  \item{group}{Vector specifying which columns in the data matrix correspond to which condition or group}
  \item{net_info}{List with named elements "zero" and "one" corresponding to the zero and one information matrices used in \code{netEst.dir} and \code{netEst.undir}}
  \item{n}{Vector of the total number of observations for \code{grp}. Used to determine \code{lambda} in \code{netEst.dir} and \code{netEst.undir}}
  \item{p}{Total number of genes. Used to determine \code{lambda} in \code{netEst.dir} and \code{netEst.undir}}
  \item{lambda_c}{lambda constant. Used to determine \code{lambda} in \code{netEst.dir} and \code{netEst.undir}}
  \item{eta}{Value of eta passed to \code{netEst.undir}}
  \item{net_clusters}{Named numeric vector specifying which genes correspond to which clusters. Names are genes and the values are their corresponding clusters}
}

\details{
This function loops through each cluster and calls \code{netEst.undir} or \code{netEst.dir} with the relevant parameters.
}


\value{
A list with components
\item{Adj}{List of the weighted adjacency matrices (partial correlations) for each cluster. The structure is Adj[[cluster]]. If cluster = FALSE, Adj[[cluster]] only has one element. Note this is used in a lapply where we loop over the groups giving us the final Adj[[condition]][[cluster]] structure.}

\item{invcov}{List of estimated inverse covariance matrices for each cluster. The structure is invcov[[cluster]]. If cluster = FALSE, invcov[[cluster]] only has one element}

\item{lambda}{List of values of tuning parameters used for each cluster.}

}

\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Mike Hellstern
}

\seealso{
\code{\link{prepareAdjacencyMatrixDB}}
}
