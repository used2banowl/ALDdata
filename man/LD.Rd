\name{LD}

\alias{LD}

\docType{data}

\title{LD data from Hapmap}

\description{Reference data from the hapmap project.}

\usage{data(ceu1)}

\format{
  A data.frame with one row per SNP
  \describe{
   \item{pos1}{Position of marker 1}
   \item{pos2}{Position of marker 2}
   \item{pop}{Hapmap population}
   \item{rs1}{Marker 1}
   \item{rs2}{Marker 2}
   \item{Dprime}{D' measure for LD between rs1 and rs2}
   \item{r2}{r-squared measure for LD between rs1 and rs2}
   \item{log}{log score for LD between rs1 and rs2}
   }
}

\details{
  There are separate data file for each chromosome for each of the following Hapmap populations: ASW, CEU, CHB, JPT, YRI.
}

\examples{
data(asw1)
data(jpt12)
}

\keyword{datasets}
