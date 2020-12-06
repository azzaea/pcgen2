#' The pc algorithm applied to residuals
#'
#' The pc algorithm applied to residuals, assuming there are replicates, and K =
#' Z Z^t.
#'
#' If \code{use.GBLUP = FALSE}, GBLUP residuals are used as input for the
#' pc-stable algorithm of Colombo and Maathuis (2014). This closely resembles
#' the residual networks of Valente et al., (2010) and Topner et al., (2017)
#' (who used different ways to predict the genetic effects, and applied other
#' causal inference algorithms to the residuals). When \code{use.GBLUP = TRUE},
#' pc-stable is applied to the GBLUP itself, which resembles the genomic
#' networks of Topner et al., (2017). If \code{cov.method = "uni"}, the GBLUP
#' and the residuals are computed separately for each trait in suffStat. The
#' covariance of each trait is assumed to be \deqn{\sigma_G^2 Z K Z^t +
#' \sigma_E^2 I_n} where \eqn{Z} is a binary incidence matrix, assigning plants
#' or plots to genotypes. \eqn{Z} is based on the first column in suffStat. If
#' there is a single observation per genotype (typically a genotypic mean),
#' \eqn{Z} is the identity matrix, and the relatedness matrix \eqn{K} should be
#' specified. If there are replicates for at least some of the genotypes, and no
#' \eqn{K} is provided, independent genetic effects are assumed (\eqn{K} will be
#' the identity matrix). It is also possible to have replicates and specify a
#' non-diagonal \eqn{K}. Whenever  \eqn{K} is specified, sommer (mmer2) will be
#' used; otherwise lmer (lme4). mmer2 is also used when \code{cov.method =
#' "us"}, in which case the multivariate GBLUP is computed, for all traits in
#' suffStat simultaneously. This is only possible for a limited number of
#' traits.
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
#' @references Valente, B. D., Rosa, G. J., de los Campos, G., Gianola, D., &
#'   Silva, M. A. (2010). Searching for recursive causal structures in
#'   multivariate quantitative genetics mixed models. \emph{Genetics}, 185(2), 633-644.
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @inheritParams pcgen
#'
#' @param m.max maximum size of the conditioning set, in the pc-algorithm on the
#'   residuals
#'
#' @param cov.method A string, specifying which method should be used to compute
#'   the G-BLUP. Options are \code{"us"} (unstructured multi-trait model fitted
#'   using \code{sommer}) and \code{"uni"} (based on univariate GBLUPs).
#'   (Default \code{"uni"})
#'
#' @param use.GBLUP Use the GBLUP itself, instead of the residuals (as in Topner
#'   et al)
#'
#' @return If \code{return.pvalues = FALSE}, the output is a graph (an object
#'   with S3 class \code{"pcgen"}). If \code{return.pvalues = TRUE}, the output
#'   is a list with elements \code{gr} (the graph) and \code{pMax} (a matrix
#'   with the p-values).
#'
#' @examples
#' \dontrun{
#' data(simdata)
#' out <- pcRes(suffStat = simdata, alpha = 0.01, verbose= FALSE)
#' }
#'
#' @references A paper on arxiv and Topner et al (2017)
#'
#' @export
#' @import pcalg

pcRes <-
  function (suffStat, alpha= 0.01, K = NULL, m.max = Inf,
            verbose = FALSE, covariates=NULL, QTLs=integer(),
            cov.method = 'uni', use.GBLUP = FALSE, return.pvalues = FALSE)
  {
    NAdelete <- TRUE
    if (is.null(alpha)) alpha = 0.01

    stopifnot(cov.method %in% c('uni','us'))

    if (max(as.numeric(unlist(lapply(suffStat, function(x){sum(is.na(x))})))) > 0.001) {
    #if (max(as.numeric(unlist(lapply(cbind(suffStat, covariates), function(x){sum(is.na(x))})))) > 0.001) {
      if (cov.method == 'us') {stop('Missing values detected. Choose cov.method = uni')}
    }

    if (length(QTLs) > 0) {
      if (is.null(covariates)) {
        res.covariates <- suffStat[,QTLs]
      } else {
        res.covariates <- cbind(covariates, suffStat[,QTLs])
      }
    } else {
      res.covariates <- covariates
    }

    res <- getResiduals(suffStat = suffStat, covariates = res.covariates,
                        K = K, cov.method = cov.method)
    ## estimate CPDAG

    if (use.GBLUP == TRUE) {res <- suffStat[,-1] - res}

    pc.res <- pc(suffStat =  list(C = cor(res), n = nrow(res)),
                         indepTest = gaussCItest,
                         alpha = alpha, labels = colnames(suffStat)[-1],
                         verbose = verbose, m.max = m.max, skel.method = "stable",
                         u2pd = "relaxed", conservative = FALSE, maj.rule = TRUE,
                         solve.confl = TRUE)

    if (return.pvalues == TRUE) {
      return(list(gr = pc.res, pMax = pc.res@pMax))
    } else {
      return(pc.res)
    }

}
