\name{preparePathways}
\alias{preparePathways}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare pathway dataset needed by NetGSA
}
\description{
Prepare pathway dataset needed by NetGSA. See \code{NetGSA} for more details.
}
\usage{
preparePathways(db=c("kegg", "MSigDB"), 
           type=c("H","C1","C2","C3","C4","C5","C6","C7"),
           genename= c("EntrezID", "symbol"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{Database to build pathway from. Could be either \code{'kegg'} or \code{'MSigDB'}.}
  \item{type}{The type of pathways to choose from if \code{db=='MSigDB'}.}
  \item{genename}{Whether gene symbol or EntrezID is used.}
}
\value{
%%  ~Describe the value returned
A list of pathways.
%%  If it is a LIST, use
%% ...
}
\author{
Jing Ma (jingma@fredhutch.org)
}

\seealso{
\code{\link{NetGSA}}}

\examples{
#library(graphite)

#pathwayList <- preparePathways('kegg')
#pathwayList[[1]]
}
