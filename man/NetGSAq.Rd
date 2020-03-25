\name{NetGSAq}
\alias{NetGSAq}
\title{"Quick" Network-based Gene Set Analysis}
\description{Quick version of \code{NetGSA}}
\usage{
NetGSAq(x, group, pathways, lambda_c = 1, file_e = NULL, file_ne = NULL,
    lklMethod="REHE", cluster = TRUE, sampling = FALSE, sample_n = NULL,
    sample_p = NULL, sample_seeds = NULL, minsize=5, eta=0.1, lim4kappa=500)
}
\arguments{
  \item{x}{See \code{x} argument in \code{NetGSA}}
  \item{group}{See \code{group} argument in \code{NetGSA}}
  \item{pathways}{See \code{pathways} argument in \code{NetGSA}}
  \item{lambda_c}{See \code{lambda_c} argument in \code{prepareAdjMat}}
  \item{file_e}{See \code{file_e} argument in \code{prepareAdjMat}}
  \item{file_ne}{See \code{file_ne} argument in \code{prepareAdjMat}}
  \item{lklMethod}{See \code{lklMethod} argument in \code{NetGSA}}
  \item{sampling}{See \code{sampling} argument in \code{NetGSA}}
  \item{sample_n}{See \code{sample_n} argument in \code{NetGSA}}
  \item{sample_p}{See \code{sample_p} argument in \code{NetGSA}}
  \item{minsize}{See \code{minsize} argument in \code{NetGSA}}
  \item{eta}{See \code{eta} argument in \code{NetGSA}}
  \item{lim4kappa}{See \code{lim4kappa} argument in \code{NetGSA}}
}
\details{
This is a wrapper function to perform weighted adjacency matrix estimation and pathway enrichment in one step. For more details see \code{?prepareAdjMat} and \code{?NetGSA}.
}
\value{
A list with components
\item{results}{A data frame with pathway names, pathway sizes, p-values and false discovery rate corrected q-values, and test statistic for all pathways.}
\item{beta}{Vector of fixed effects of length \eqn{2p}, of which the first half is for condition 1 and the second half for condition 2.}
\item{s2.epsilon}{Variance of the random errors \eqn{\epsilon}.}
\item{s2.gamma}{Variance of the random effects \eqn{\gamma}.}
\item{graph}{List of components needed in \code{plot.NetGSA}.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. \url{https://doi.org/10.1093/bioinformatics/btw410}

Shojaie, A., & Michailidis, G. (2010). Network enrichment analysis in complex experiments. Statistical applications in genetics and molecular biology, 9(1), Article 22. \url{http://www.ncbi.nlm.nih.gov/pubmed/20597848}.

Shojaie, A., & Michailidis, G. (2009). Analysis of gene sets based on the underlying regulatory network. Journal of Computational Biology, 16(3), 407-426. \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131840/}

}
\author{
Michael Hellstern
}
\seealso{
\code{\link{prepareAdjMat}}, \code{\link{netEst.dir}}, \code{\link{netEst.undir}}
}
