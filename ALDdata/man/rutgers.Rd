\name{rutgers}

\alias{rutgers36}
\alias{rutgers37}

\docType{data}

\title{Rutgers genetic maps for builds 36 and 37}

\description{Location maps of the human genome with genetic and physical
  positions.}

\usage{data(rutgers36)}

\format{
  A data frame with 23186 - 24168 observations on the following 10 variables.
  \describe{
    \item{\code{marker}}{Marker name}
    \item{\code{type}}{Marker type (ie \code{RFLP} \code{SNP} \code{STS} \code{other})}
    \item{\code{n.meioses}}{Number of informative meioses}
    \item{\code{heterozygosity}}{Heterozygosity}
    \item{\code{phys.pos}}{Physical position in the human genome}
    \item{\code{cM}}{Genetic position (averaged over both sexes)}
    \item{\code{cM.female}}{Female specific genetic position}
    \item{\code{cM.male}}{Male specific genetic position}
    \item{\code{chr}}{Chromosome number}
  }
}

\source{http://compgen.rutgers.edu/RutgersMap/DownloadMap.aspx}

\references{
  Nato, AQ Jr, Buyske, S, and Matise, TC. 2013. The Rutgers
  Map: A third-generation combined linkage-physical map of the human
  genome. [in preparation as of September 1, 2013]

  Matise, TC et al. 2007. A second-generation combined linkage
  physical map of the human genome. Genome Res 17:1783-6.

  Matise et. al. {American Journal of Human Genetics}, 2003,
  73(2):271-84}

\examples{
data(rutgers36)
data(rutgers37)
}

\keyword{datasets}
