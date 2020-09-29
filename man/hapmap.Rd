\name{hapmap}

\alias{hapmap}
\alias{hapmap36}

\docType{data}

\title{Hapmap data}

\description{Reference data from the hapmap project.}

\usage{data(hapmap)}

\format{
  A data.frame with one row per SNP
  \describe{
    \item{rs}{RS number}
    \item{chr}{Chromosome number}
    \item{pos}{Chromosome position}
    \item{ref}{Reference base}
    \item{var}{Variant}
    \item{f.yri}{Variant frequency in YRI}
    \item{f.ceu}{Variant frequency in CEU}
    \item{f.chb}{Variant frequency in CHB}
    \item{f.jpt}{Variant frequency in JPT}
    \item{cM}{Approxiamte location in cM}
  }
}

\examples{
data(hapmap)
data(hapmap36)
}

\keyword{datasets}
