#' Simulated dataset of a 2 traits network with replicates
#'
#' Simulated data, for two replicates of 500 genotypes g1,...,g500. Two traits
#' were simulated (Y1 and Y2), such that G -> Y1 -> Y2
#'
#' @format A data frame of dimension \eqn{ 1000 \times 3 }. The first column is
#'   the factor G (genotype); the subsequent columns contain \eqn{Y_1} and
#'   \eqn{Y_2}.
#'
#'
#' @usage data(d2treps)
#' @keywords datasets
#' @family pcgen datasets

"d2treps"


#' Simulated dataset of a 2 traits network with replicates means
#'
#' The means of 100 markers and two traits corresponding to the dataset
#' \code{d2treps} of 500 genotypes g1,...,g500, where G -> Y1 -> Y2
#'
#' @format A data frame of dimension \eqn{ 500 \times 103}. The first column is
#'   the factor G (genotype); then \eqn{Y_1} and \eqn{Y_2}, followed by the 100
#'   markers \code{mxxxx} where \code{xxx} is the marker number
#'
#' @usage data(d2treps)
#' @keywords datasets
#' @family pcgen datasets

"d2tsnpsmeans"
