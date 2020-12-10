#' simdata: Simulated data
#'
#' Simulated data, for two replicates of 200 genotypes g1,...,g200. Three traits
#' were simulated (Y1, Y2 and Y3), using a structural equation model defined by
#' Y1 -> Y2 -> Y3, and direct genetic effects on Y1 and Y3.
#'
#' @format A data frame of dimension \eqn{ 4 \times 400 }. The first column is
#'   the factor G (genotype); the subsequent columns contain \eqn{Y_1, Y_2} and
#'   \eqn{Y_3}.
#'
#'
#' @usage data(simdata)
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#'  data(simdata)
#'  out <- pcgen(simdata)
#'  out2 <- pcRes(suffStat = simdata, alpha = 0.01, verbose= FALSE)
#' }

"simdata"


