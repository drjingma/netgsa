\name{prepareAdjacencyMatrixDB}
\alias{prepareAdjacencyMatrixDB}
\title{Construct adjacency matrices from graphite databases and/or user provided network information}
\description{
Read the network information from any of the graphite databases specified by the user and construct the adjacency matrices needed for NetGSA. This function also allows for clustering. See details for more information}
\usage{
prepareAdjacencyMatrixDB(x, group, databases = NULL, cluster = NULL,
       file_e=NULL, file_ne=NULL, lambda_c=1, eta=0.5)
}
\arguments{
  \item{x}{The \eqn{p \times n}{p x n} data matrix with rows referring to genes and columns to samples. Rownames should be unique and have gene ID types appended to them. The id and gene number should be separated by a colon. E.g. "ENTREZID:127550". It is very important to have these separated by a colon since \code{NetGSA} uses regular expressions to split the row names of x into gene value and gene ID.}
  
  \item{group}{Vector of class indicators of length \eqn{n}.}
  
  \item{databases}{(Optional) Character vector of graphite databases you wish to search for edges. Since one can search in multiple databases with different identifiers, converts genes using \code{AnnotationDbi::select} and convert metabolites using \code{graphite:::metabolites()}. Databases are also used to specify non-edges. See Details for more information. }
  
  \item{cluster}{(Optional) Logical indicating whether or not to cluster genes to estimate adjacency matrix. If not specified, set to TRUE if there are > 2,500 genes (p > 2,500). The main use of clustering is to speed up calcuation time. If the dimension of the problem, or equivalently the total number of unique genes across all pathways, is large, \code{prepareAdjacencyMatrixDB} may be a slow.
  
  If clustering is set to TRUE, the 0-1 adjacency matrices are used to detect clusters of genes. Once gene clusterings are chosen, the adjacency matrices are estimated for each cluster separately using \code{netEst.undir} or \code{netEst.dir}. Thus, the adjacency matrix for the full network is block diagonal with the blocks being the adjacency matrices from the clustering algorithms. Any edges between clusters are set to 0, so this can be thought of as an approximate adjacency matrix. Six clustering algorithms from the \code{igraph} package are considered: \code{cluster_walktrap}, \code{cluster_leading_eigen}, \code{cluster_fast_greedy}, \code{cluster_label_prop}, \code{cluster_infomap}, and \code{cluster_louvain}. To ensure increases in speed, the algorithms are first subset to those having a maximum cluster size of < 1,000 genes. Among those, the algorithm with the smallest edge loss is chosen. If all algorithms have a maximum cluster size > 1,000 genes the one with the smallest maximum cluster size is chosen. Edge loss is defined as the number of edges between genes of different clusters. These edges are "lost" since they are set to 0 in the block diagonal adjacency matrix.
  
  Singleton clusters are combined into one big cluster. This should not affect performance much since the gene in a singleton cluster likely does not have any edges to other genes.
  }
  
  \item{file_e}{(Optional) The name of the file which the list of edges is to read from. This file is read in with \code{data.table::fread}. Must have 4 columns in the following order. The columns do not necessarily need to be named, but they must be in this specific order:
  \itemize{
   \item 1st column - Source gene (base_gene_src), e.g. "7534""
   \item 2nd column - Gene identifier of the source gene (base_id_src), e.g. "ENTREZID"
   \item 3rd column - Destination gene (base_gene_dest), e.g. "8607"
   \item 4th column - Gene identifier of the destination gene (base_id_dest) e.g. "UNIPROT" 
  }
  
  This information cannot conflict with the user specified non-edges. That is, one cannot have the same edge in \code{file_e} and \code{file_ne}. In the case where the graph is undirected everything will be converted to an undirected edge or non-edge. Thus if the user specifies A->B as a directed non-edge it will be changed to an undirected non-edge if the graph is undirected. See Details for more information.
  }
  
  \item{file_ne}{(Optional) The name of the file which the list of non-edges is to read from. This file is read in with \code{data.table::fread}. The edges in this file are negative in the sense that the corresponding vertices are not connected. Format of the file must be the same as \code{file_e}. Again, each observation is assumed to be a directed edge. Thus if you want a negative undirected edge input two separate negative edges.
  
  In the case of conflicting information between \code{file_ne} and edges identified in a database, user non-edges are used. That is if the user specifies A->B in \code{file_ne}, but there is an edge between A->B in KEGG, the information in KEGG will be ignored and A->B will be treated as a non-edge. In the case where the graph is undirected everything will be converted to an undirected edge or non-edge. Thus if the user specifies A->B as a directed non-edge it will be changed to an undirected non-edge if the graph is undirected. See Details for more information.
  }
  
  \item{lambda_c}{(Non-negative) a constant multiplied to the tuning parameter \code{lambda} needed for estimating the edge weights. By default, \code{lambda_c} is set to 1. See \code{\link{netEst.undir}} and \code{\link{netEst.dir}} for more details.}
  
  \item{eta}{(Non-negative) a small constant needed for estimating the edge weights. By default, \code{eta} is set to 0.5. See \code{\link{netEst.undir}} for more details.}
}
\details{

The function \code{prepareAdjacencyMatrixDB} accepts both network information from user specified sources as well as a list of graphite databases to search for edges in. \code{prepareAdjacencyMatrixDB} calculates the 0-1 adjacency matrices and runs \code{\link{netEst.undir}} or \code{\link{netEst.dir}} if the graph is undirected or directed. 

When searching for network information, \code{prepareAdjacencyMatrixDB} makes some important assumptions about edges and non-edges. As already stated, the first is that user specified non-edges are given precedence. 

\code{prepareAdjacencyMatrixDB} uses \code{databases} to identify non-edges. If two genes are identified in the \code{databases} edges but there is no edge between them this will be coded as a non-edge. The rationale is that if there was an edge between these two genes it would be present.

No information about genes not identified in \code{databases} edgelists is assumed. That is, if the user passes gene A, but gene A is not found in any of the edges in \code{databases} no information about Gene A is assumed. Gene A will have neither edges nor non-edges.

Once all the network and clustering information has been compiled, \code{prepareAdjacencyMatrixDB} estimates the network. \code{prepareAdjacencyMatrixDB} will automatically detect directed graphs, rearrange them to the correct order and use \code{netEst.dir} to estimate the network. When the graph is undirected \code{netEst.undir} will be used. For more information on these methods see \code{?netEst.dir} and \code{?netEst.undir}.

Importantly, \code{prepareAdjacencyMatrixDB} returns the list of weighted adjacency matrices to be used as an input in \code{NetGSA}.

}
\value{
A list with components
\item{Adj}{A list of weighted adjacency matrices estimated from either \code{netEst.undir} or \code{netEst.dir}. That is \code{length(Adj) = length(unique(group))}. One list of weighted adjacency matrix will be returned for each condition in group. If cluster = TRUE is specified, the length of the list of adjacency matrices for each condition will be the same length as the number of clusters. If cluster = FALSE is sepcified there will only be one list element for each condition. The structure of Adj is Adj[[condition_number]][[cluster_adj_matrix]]}
\item{invcov}{A list of inverse covariance matrices estimated from either \code{netEst.undir} or \code{netEst.dir}. That is \code{length(invcov) = length(unique(group))}. One list of inverse covariance matrix will be returned for each condition in group. If cluster = TRUE is specified, the length of the list of inverse covariance matrices for each condition will be the same length as the number of clusters. If cluster = FALSE is sepcified there will only be one list element for each condition. The structure of invcov is invcov[[condition_number]][[cluster_adj_matrix]]}
\item{lambda}{A list of values of tuning parametes used for each condition in \code{group}. If cluster = TRUE is specified, the length of the list of tuning parameters for each condition will be the same length as the number of clusters. If cluster = FALSE is specified there will only be one list element for each condition.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}
\author{
Mike Hellstern
}
\seealso{
\code{\link{NetGSA}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
\examples{
## load the data
data("breastcancer2012")

## consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"
genenames    <- unique(c(pathways[[24]], pathways[[52]]))
sx           <- x[match(rownames(x), genenames, nomatch = 0L) > 0L,]

adj_cluster    <- prepareAdjacencyMatrixDB(sx, group, databases = c("kegg", "reactome", "biocarta"), cluster = TRUE)
adj_no_cluster <- prepareAdjacencyMatrixDB(sx, group, databases = c("kegg", "reactome", "biocarta"), cluster = FALSE)

}
