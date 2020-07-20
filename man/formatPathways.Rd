\name{formatPathways}
\alias{formatPathways}
\title{Format cytoscape nested networks}

\description{
Format cytoscape nested networks using preset NetGSA format}

\usage{
formatPathways(x, pways, graph_layout = NULL)
}

\arguments{
  \item{x}{A NetGSA object returned from calling \code{NetGSA()}}
  \item{pways}{Character vector of pathways to format}
  \item{graph_layout}{(Optional) Layout to pass to plots. Must be a string for Cytoscape which will be passed to \code{RCy3::layoutNetwork}.}
}
  
\details{
Loads gene testing data into each pathway. Genes are tested using an F-test if there are 2 or more conditions or a two-sided one-class t-test against the null hypothesis of mean = 0 if there is only one condition. FDR corrected q-values are mapped to the color of the node. The scale ranges from 0 to 1 with red represents q-values of 0 and white representing q-values of 1. Loaded data includes: p-value from the F-test/t-test (pval), FDR corrected q-value (pFdr), test statistic from the F-test/t-test (teststat).

Custom formatting can be applied using the cytoscape GUI or the RCy3 pacakge.
}

\value{
}


\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Michael Hellstern
}

\seealso{
\code{\link{plot.NetGSA}}
}

\examples{

#Cytoscape open
plot(out_cluster)
formatPathways(out_cluster, "ErbB signaling pathway")
}
