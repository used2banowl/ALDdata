#' Hapmap variant data
#'
#' Phased hapmap variant data
#'
#' @format A Large matrix with one row per haplotype and one column per SNP
# @source 
"phased"

#' Hapmap LD data
#' 
#' Linkage disequilibrium data derived from Hapmap samples
#' 
#' @format A data.frame with one row per SNP
#' #' \describe{
#'   \item{pos1}{Position of marker 1}
#'   \item{pos2}{Position of marker 2}
#'   \item{pop}{Hapmap population}
#'   \item{rs1}{Marker 1}
#'   \item{rs2}{Marker 2}
#'   \item{Dprime}{D' measure for LD between rs1 and rs2}
#'   \item{r2}{r-squared measure for LD between rs1 and rs2}
#'   \item{log}{log score for LD between rs1 and rs2}
#' }
# @source
"LD"