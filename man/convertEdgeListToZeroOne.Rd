\name{convertEdgeListToZeroOne}

\alias{convertEdgeListToZeroOne}

\title{Convert edgelist and list of non-edges to zero/one adjacency matrices}

\description{
Convert edgelist and list of non-edges to zero/one matrices used in \code{netEst.dir} and \code{netEst.undir}. This is a helper function in \code{prepareAdjMat} and should not be called by the user.
}


\usage{
convertEdgeListToZeroOne(edgelist, non_edges, genes, genes_not_in_dbs)
}


\arguments{
  \item{edgelist}{Data.table of all edges (user specified and those found in graphite)}
  \item{non_edges}{Data.table of all non-edges. Graphite doesn't give non-edges so this should only be user specified}
  \item{genes}{Character vector of genes}
  \item{genes_not_in_dbs}{Character vector of genes not found in the graphite databases}
}

\details{
This function converts edges and non-edges into both the zero and one matrices used in \code{netEst.dir} and \code{netEst.undir}. If two genes are in the graphite databases, but we don't see an edge between them they are classified as having a non-edge. However, if they are not in the graphite database we classify them as having neither an edge between them nor a non-edge between them. We are essentially assuming we have no information on these. The entries in both the zero and one matrices will be 0.
}


\value{
A list with components
\item{ones_freq}{A pxp matrix with the number of graphite databases a specific edge occurs in as the value}
\item{ones}{A pxp matrix. Values are 1 if the edge occurs at least once, 0 otherwise}
\item{zeros}{A pxp matrix. Values are 1 if we know this is not an edge, 0 otherwise}
\item{directed}{Logical indicating if our edgelist is a directed acyclic graph (DAG). Checked using \code{igraph::is_dag}}
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
