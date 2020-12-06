#' Estimate genetic covariances between all pairs of traits, and test their
#' significance
#'
#' This will perform the conditional independence test used in pcgen algorithm,
#' assuming there are replicates, and K = Z Z^t
#'
#'
#' @inheritParams pcgen
#'
#' @param out.cor If \code{TRUE}, the output will contain estimates of genetic
#'   correlations; otherwise covariances. The pvalues are always for genetic
#'   covariance.
#'
#' @return A list with elements pvalues and out.cor, which are both p x p
#'   matrices
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @references Kruijer, W., Behrouzi, P., Bustos-Korts, D., Rodríguez-Álvarez,
#' M. X., Mahmoudi, S. M., Yandell, B., ... & van Eeuwijk, F. A. (2020).
#' Reconstruction of networks with direct and indirect genetic effects.
#' \emph{Genetics}, 214(4), 781-807.
#'
#' @examples
#' \dontrun{
#' data(simdata)
#' test <- gencovTest(suffStat= simdata, max.iter = 200, out.cor= TRUE )
#' }
#'
#' @export
#'
gencovTest <- function(suffStat, max.iter = 200, out.cor = TRUE)
{

  #QTLs <- integer()
  #covariates <- NULL

  # stop if the first column in suffStat does not have the name G
  stopifnot(names(suffStat)[1] == 'G')

  X <- as.data.frame(matrix(0,nrow(suffStat),0))

  ###################################################################
  # The case where x and y both represent real traits (no QTLs), and 1
  # (direct genetic effects: G) IS contained in S. S may also contain QTLs.
  # We do not skip any tests here.
  # Note that G will account for QTLs that may not be contained in S.

  p       <- ncol(suffStat)
  pvalues <- matrix(NA, p-1, p-1)
  gencor  <- matrix(NA, p-1, p-1)

  for (j in 1:(p-2)) {
    for (k in (j+1):(p-1)) {
      out <- gen.covar.test(x=suffStat[,j+1], y=suffStat[,k+1],
                            G=suffStat[,1], max.iter = max.iter)
      pvalues[j,k] <- pvalues[k,j] <- out$pvalue

      if (out.cor == TRUE) {
        if (out$Vg[1] > 0 & out$Vg[2] > 0) {
          gencor[j,k]  <- gencor[k,j] <- out$Vg[3] / sqrt(out$Vg[1] * out$Vg[2])
        }
      } else {
        gencor[j,k]  <- gencor[k,j] <- out$Vg[3]
      }

    }
  }

  colnames(pvalues) <- rownames(pvalues) <- colnames(suffStat)[-1]
  colnames(gencor)  <- rownames(gencor)  <- colnames(suffStat)[-1]

  return(list(pvalues = pvalues, gencor = gencor))

}

