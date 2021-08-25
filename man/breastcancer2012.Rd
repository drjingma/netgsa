\name{breastcancer2012}
\alias{breastcancer2012}
\docType{data}
\title{
Breast cancer data from TCGA (2012).
}
\description{
An example data set consisting of RNA-seq gene expression data, KEGG pathways, edge list and non-edge list.
}
\usage{data(breastcancer2012)}
\format{
  A list with components
  \describe{
    \item{\code{x}}{The \eqn{p \times n}{p x n} data matrix.}
    \item{\code{group}}{The vector of class indicators of length \eqn{n}.}
    \item{\code{pathways}}{A list of KEGG pathways.}
    \item{\code{edgelist}}{A data frame of edges, each row corresponding to one edge.}
    \item{\code{nonedgelist}}{A data frame of nonedges, each row corresponding to one negative edge.}
    \item{\code{pathways_mat}}{Matrix with pathway indicators}
  }
}
\references{
Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. Nature, 490(7418), 61.
}
\examples{
data("breastcancer2012")
}
\keyword{datasets}
