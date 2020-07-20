\name{checkUserEdges}
\alias{checkUserEdges}
\title{Check user edgelist for errors}

\description{
Read in user edgelist as data.table and check the column formatting and check for conflicting information. This is a helper function in \code{prepareAdjMat} and should not be called by the user.}

\usage{
checkUserEdges(edgelist, non_edges)
}

\arguments{
  \item{edgelist}{File path to user specified edgelist. Read in with \code{data.table::fread}}
  \item{non_edges}{File path to user specified non-edgelist. Read in with \code{data.table::fread}}}
  
\details{
This function checks to make sure there is the correct number of columns specified in the edgeslist (4 or 5). It also checks for non-numeric values and to see if there is inconsistent information such as an edge specified in the non-edgelist and the edgelist.
}

\value{
A list with components
\item{user_edges}{A data.table of the user specified edges}
\item{user_non_edges}{A data.table of the user specified non-edges}
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
