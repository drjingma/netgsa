\name{stackDatabases}

\alias{stackDatabases}

\title{Combine edges from databases into a data.table}

\description{
Retrieves edges from specified databases and stacks them into one data.table
}


\usage{
stackDatabases(databases)
}


\arguments{
  \item{databases}{Character vector of databases to compile. Should be one of the options from hspaiens in \code{graphite::pathwayDatabases()}}
}

\details{
This function compiles all the edges from all databases specified into one data.table
}


\value{
A data.table with columns:
\item{database}{Which database the edge comes from}
\item{src}{Source gene}
\item{src_type}{Source gene identifier type}
\item{dest}{Destination gene}
\item{dest_type}{Destination gene identifier type}
\item{direction}{Direction of edge. Either Directed or Undirected}
}

\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Mike Hellstern
}

\seealso{
\code{\link{obtainEdgeList}}
}
