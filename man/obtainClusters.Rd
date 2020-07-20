\name{obtainClusters}

\alias{obtainClusters}

\title{Estimate optimal gene clustering structure}

\description{
Tries six different clustering methods and chooses the one with the best results. This is a helper function in \code{prepareAdjMat} and should not be called by the user.
}


\usage{
obtainClusters(A, order, cluster)
}


\arguments{
  \item{A}{A 0-1 adjacency matrix}
}

\details{
This function tries the six different clustering methods in igraph and chooses the best one. As stated in \code{prepareAdjMat} the six methods evaluated are: \code{cluster_walktrap}, \code{cluster_leading_eigen}, \code{cluster_fast_greedy}, \code{cluster_label_prop}, \code{cluster_infomap}, and \code{cluster_louvain}. See \code{prepareAdjMat} for how the best is chosen. Even if \code{cluster = FALSE}, connected components of the 0-1 adjacency matrix are used as clusters.

It is essential that the order of the returned named numeric vector must be in the same order as the rows of the data matrix.
}


\value{
Named numeric vector of membership. The name of each element is the corresponding gene and the value is the cluster it belongs to.
}

\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Michael Hellstern
}

\seealso{
\code{\link{prepareAdjMat}}
}
