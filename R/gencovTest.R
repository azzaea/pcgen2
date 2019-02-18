#' The conditional independence test in pcgen
#'
#' This will perform the conditional independence test used in pcgen algorithm,
#' assuming there are replicates, and K = Z Z^t
#'
#' Here we
#' put a lot
#' of details
#' COMMENT on Vg and Ve, which should be either both NULL or positive definite matrices
#'
#' @inheritParams pcgen
#'
#' @param x,y column numbers in suffStat that should be tested for conditional
#'            independence given the variables in S
#'
#' @param S vector of integers defining the conditioning set,
#'          where the integers refer to column numbers in suffStat.
#'          May be numeric(), i.e. the empty set.
#'
#' @param alpha (Default 0.01) The significance level used in the test. The test
#'              itself of course does not depend on this, but it is used in the
#'              EM-algorithm to speed up calculations. More precisely, the
#'              EM-algorithm is stopped once the p-value is below the
#'              significance level. This can be done because the PC-algorithm
#'              only needs an accept/reject decision.
#'
#' @param stop.if.significant If TRUE, the EM-algorithm used in some
#'              of the conditional independence tests will be stopped
#'              whenever the p-value becomes significant, i.e. below
#'              alpha. This will speed up calculations, and can be done
#'              because (1) the PC algorithm only needs an accept/reject
#'              decision (2) In EM the likelihood is nondecreasing. Should
#'              be put to FALSE if the precise p-values are of interest.
#'
#' @return a p-value
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references A paper on arxiv
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

