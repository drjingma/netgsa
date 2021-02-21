\name{obtainEdgeList}
\alias{obtainEdgeList}
\title{Obtain edgelist from graphite databases. To be used within prepareAdjMat}
\description{
Find all edges between genes in the specified graphite databases.
}
\usage{
obtainEdgeList(genes, databases)
}
\arguments{
  \item{genes}{Character vector of gene ID and gene value. The ID and gene value should be separated by a colon. E.g. "ENTREZID:127550". It is very important to have these separated by a colon since \code{obtainEdgeList} uses regular expressions to split this into gene value and gene ID.}
  \item{databases}{Character vector of graphite databases you wish to search for edges. Options are: biocarta, kegg, nci, panther, pathbank, pharmgkb, reactome, smpdb, ndex. Note NDEx is recommended for expert users, see details.}
}
\details{
\code{obtainEdgeList} searches through the specified databases to find edges between genes in the \code{genes} argument. Since one can search in multiple databases with different identifiers, genes are converted using \code{AnnotationDbi::select} and metabolites are converted using \code{graphite:::metabolites()}. Databases are also used to specify non-edges. This function searches through \code{graphite} databases and also has the option to search NDEx (public databases only). However, since NDEx is open-source and does not contain curated edge information like \code{graphite}, NDEx database search is a beta function and is only recommended for expert users. When searching through NDEx, gene identifiers are not converted. Only, the gene identifiers passed to the \code{genes} argument are used to search through NDEx. NDEx contains some very large networks with millions of edges and extracting those of interest can be slow.

This function is particularly useful if the user wants to create an edgelist outside of \code{prepareAdjMat}. \code{graphite} and it's databases are constantly updated. Creating and storing an edgelist outside of \code{prepareAdjMat} may help reproducibility as this guarantees the same external information is used. It can also speed up computation since if only a character vector of databases is passed to \code{prepareAdjMat}, it calls \code{obtainEdgeList} each time and each call can take several minutes. The edges from \code{obtainEdgeList} are used to create the 0-1 adjacency matrices used in \code{netEst.undir} and \code{netEst.dir}.

Using \code{obtainEdgeList} to generate edge information is highly recommended as this performs all the searching and conversion of genes to common identifiers. Inclusion of additional edges, removal of edges, or other user modifications to edgelists should be through the \code{file_e} and \code{file_ne} arguments in \code{prepareAdjMat}.
}
\value{
A list of class \code{obtainedEdgeList} with components
\item{edgelist}{A \code{data.table} listing the edges. One row per edge. Edges are assumed to be directed. So if an edge is undirected there will be two rows.}
\item{genes_not_in_dbs}{A vector of genes specified, but were not found in the databases searched}
}
\author{
Michael Hellstern
}
\seealso{
\code{\link{prepareAdjMat}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
genes <- paste0("ENTREZID:", c("10000", "10298", "106821730", "10718", "1398", "1399", "145957", "1839", "1950", "1956"))

out <- obtainEdgeList(genes, c("kegg", "reactome"))
}
