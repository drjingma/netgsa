\name{addUserEdges}

\alias{addUserEdges}

\title{Add user edgelist to edgelist found in graphite}

\description{
Combine user edges with those identified in graphite. This is a helper function in \code{prepareAdjacencyMatrixDB} and should not be called by the user.}

\usage{
addUserEdges(non_user_edges, user_edges)
}

\arguments{
  \item{non_user_edges}{Data.table of user edges found in graphite databases}
  \item{user_edges}{Data.table of user edges}
}  

\details{
This function reconciles conflicting information between user edges and edges found in graphite. This function gives preference to user information. For example, if the user specifies a frequency for an edge that will be used rather than the frequency calculated from the graphite databases.
}


\value{
A data.table of edges including both user specified and graphite identified edges.
}


\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Mike Hellstern
}

\seealso{
\code{\link{prepareAdjMat}}
}
