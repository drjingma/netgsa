\name{obtainEdgeList}
\alias{obtainEdgeList}
\title{Obtain edgelist from graphite databases. To be used within prepareAdjacencyMatrixDB}
\description{
Find all edges between genes in the specified graphite databases.
}
\usage{
obtainEdgeList(genes, databases)
}
\arguments{
  \item{genes}{Named character vector of genes. Names are ID types. E.g. \code{"ENTREZID" = "25"}}
  \item{databases}{Character vector of graphite databases you wish to search for edges.}
}
\details{
\code{obtainEdgeList} searchs through the specified databases to find edges between genes in \code{genes}. These edges can be used to create the 0-1 adjacency matrices used in \code{netEst.undir} and \code{netEst.dir}
}
\value{
A list with components
\item{edgelist}{A \code{data.table} listing the edges. One row per edge. Edges are assumed to be directed. So if an edge is undirected there will be two rows.}
\item{genes_not_in_dbs}{A vector of genes specified, but were not found in the databases searched}
}
\author{
Mike Hellstern
}
\seealso{
\code{\link{prepareAdjacencyMatrixDB}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
genes <- setNames(c("10000", "10298", "106821730", "10718", "1398", "1399", "145957", "1839", "1950", "1956"), rep("ENTREZID", 10))

out <- obtainEdgeList(genes, c("kegg", "reactome"))
}
