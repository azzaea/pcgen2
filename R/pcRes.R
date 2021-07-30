#' The pc algorithm applied to residuals
#'
#' The pc algorithm applied to residuals of a linear mixed model, assuming there
#' are replicates, and K = Z Z^t.
#'
#' If `use.GBLUP = FALSE`, GBLUP residuals are used as input for the pc-stable
#' algorithm of Colombo and Maathuis (2014).
#'
#' This closely resembles the residual networks of Valente et al., (2010) and
#' Topner et al., (2017) (who used different ways to predict the genetic
#' effects, and applied other causal inference algorithms to the residuals).
#'
#' When `use.GBLUP = TRUE`, pc-stable is applied to the GBLUP itself (i.e the
#' fitted or predicted genetic values, not the residuals), which resembles the
#' genomic networks of Topner et al., (2017).
#'
#' If `cov.method = "uni"\`, the GBLUP and the residuals are computed separately
#' for each trait in `suffStat`. The covariance of each trait is assumed to be
#' \deqn{\sigma_G^2 Z K Z^t + \sigma_E^2 I_n} where \eqn{Z} is a binary
#' incidence matrix, assigning plants or plots to genotypes. \eqn{Z} is based on
#' the first column in suffStat. If there is a single observation per genotype
#' (typically a genotypic mean), \eqn{Z} is the identity matrix, and the
#' relatedness matrix \eqn{K} should be specified. If there are replicates for
#' at least some of the genotypes, and no \eqn{K} is provided, independent
#' genetic effects are assumed (\eqn{K} will be the identity matrix). It is also
#' possible to have replicates and specify a non-diagonal \eqn{K}. Whenever
#' \eqn{K} is specified, sommer (mmer2) will be used; otherwise lmer (lme4).
#' mmer2 is also used when `cov.method = "us"`, in which case the multivariate
#' GBLUP is computed, for all traits in suffStat simultaneously. This is only
#' possible for a limited number of traits.
#'
#' Finally, also note that some \code{pc} function parameters are fixed here:
#' `skel.method = "stable"`, `conservative = FALSE`, `maj.rule = TRUE` and
#' `solve.confl = TRUE`. As a results, the resulting graph is fully
#' order-independent (see Colombo and Maathuis 2014). Additionally, we fixed
#' `u2pd = "relaxed"`, so that that an invalid CPDAG is output in case of
#' conflicting orientation information. Labels (defining the names of the nodes
#' of the graph) are derived from the data-frame `suffStat`, containing the
#' data.
#'
#' @references * Colombo, D., & Maathuis, M. H. (2014). Order-independent
#'   constraint-based causal structure learning. \emph{The Journal of Machine
#'   Learning Research}, 15(1), 3741-3782.
#' @references * Kruijer, W., Behrouzi, P., Bustos-Korts, D., Rodríguez-Álvarez,
#'   M. X., Mahmoudi, S. M., Yandell, B., ... & van Eeuwijk, F. A. (2020).
#'   Reconstruction of networks with direct and indirect genetic effects.
#'   \emph{Genetics}, 214(4), 781-807.
#' @references * Töpner, K., Rosa, G. J., Gianola, D., & Schön, C. C. (2017).
#'   Bayesian networks illustrate genomic and residual trait connections in
#'   maize (Zea mays L.). \emph{G3: Genes, Genomes, Genetics}, 7(8), 2779-2789.
#' @references * Valente, B. D., Rosa, G. J., de los Campos, G., Gianola, D., &
#'   Silva, M. A. (2010). Searching for recursive causal structures in
#'   multivariate quantitative genetics mixed models. \emph{Genetics}, 185(2),
#'   633-644.
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @inheritParams pcgen
#' @inheritParams getResiduals
#'
#' @param m.max maximum size of the conditioning set, in the pc-algorithm on the
#'   residuals
#'
#'
#' @param use.GBLUP Use the GBLUP itself (ie the fitted value from the linear
#'   mixed model), instead of the residuals (as in Topner et al). #These are the
#'   fitted values from the lmm (Y^)
#'
#' @return If `return.pvalues = FALSE`, the output is a graph (an object with S3
#'   class \code{"pcgen"}). If `return.pvalues = TRUE`, the output is a list
#'   with elements `gr` (the graph) and `pMax` (a matrix with the p-values).
#'
#' @examples
#' \dontrun{
#' data(simdata)
#' out <- pcRes(suffStat = simdata, alpha = 0.01, verbose= FALSE)
#' }
#'

#'
#' @export
#' @importFrom pcalg pc gaussCItest
#' @importFrom stats cor

pcRes <- function (suffStat, K = NULL, covariates = NULL, QTLs = integer(),
            cov.method = 'uni', use.GBLUP = FALSE, alpha = 0.01, m.max = Inf,
            return.pvalues = FALSE, verbose = FALSE)
  {

    stopifnot(cov.method %in% c('uni','us'))

    if (max(as.numeric(unlist(lapply(suffStat, function(x){sum(is.na(x))})))) > 0.001)
      if (cov.method == 'us')
        stop('Missing values detected. Choose cov.method = uni')


    if (length(QTLs) > 0) {
      if (is.null(covariates)) {
        res.covariates <- suffStat[, QTLs]
      } else {
        res.covariates <- cbind(covariates, suffStat[, QTLs])
      }
    } else {
      res.covariates <- covariates
    }

    res <- getResiduals(suffStat = suffStat, covariates = res.covariates,
                        K = K, cov.method = cov.method, verbose = verbose)

    if (use.GBLUP == TRUE) {res <- suffStat[, -1] - res}

    pc.res <- pc(suffStat =  list(C = cor(res), n = nrow(res)),
                 indepTest = gaussCItest, alpha = alpha,
                 labels = colnames(suffStat)[-1], verbose = verbose,
                 m.max = m.max, skel.method = "stable", u2pd = "relaxed",
                 conservative = FALSE, maj.rule = TRUE, solve.confl = TRUE)

    if (return.pvalues == TRUE) {
      return(list(gr = pc.res, pMax = pc.res@pMax))
    } else {
      return(pc.res)
    }

}
