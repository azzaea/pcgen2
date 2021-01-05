#' The pcgenFast algorithm
#'
#' pcgen with residual-based screening
#'
#' The pcgen algorithm starting with a skeleton estimated using the standard
#' pc-algorithm, based on residuals from the GBLUP.
#'
#' Here we put a lot of details COMMENT on cov.method (the same option is used
#' in the computation of the residuals and the conditional means) fixedGaps and
#' fixedEdges are NULL m.max, alpha,.... are also used in the first step
#'
#' @references * Kruijer, W., Behrouzi, P., Bustos-Korts, D., Rodríguez-Álvarez,
#'   M. X., Mahmoudi, S. M., Yandell, B., ... & van Eeuwijk, F. A. (2020).
#'   Reconstruction of networks with direct and indirect genetic effects.
#'   \emph{Genetics}, 214(4), 781-807.
#' @references * Colombo, D., & Maathuis, M. H. (2014). Order-independent
#'   constraint-based causal structure learning. \emph{The Journal of Machine
#'   Learning Research}, 15(1), 3741-3782.
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @inheritParams pcgen
#' @inheritParams getResiduals
#'
#' @param res.m.max Maximum size of the conditioning set, in the pc-algorithm on
#'   the residuals (used for prior screening).
#'
#' @param use.res If \code{FALSE}, residuals from GBLUP are only used for screening with the
#'   standard pc algotihm. After that, the standard pcgen algorithm is run on
#'   the remaining edges; the test for conditional independence of 2 traits
#'   given a set of other traits and G is based on bivariate mixed models. If
#'   \code{TRUE}, this test for conditional independence of 2 traits given a set
#'   of other traits and G is based on the residuals. In this case, no further
#'   edges between traits are removed after screening and \code{pcgen} will only
#'   infer the orientation, and the direct genetic effects.
#'
#'
#' @return If \code{return.pvalues = FALSE}, the output is a graph (an object
#'   with S3 class \code{"pcgen"}). If \code{return.pvalues = TRUE}, the output
#'   is a list with elements \code{gr} (the graph) and \code{pMax} (a matrix
#'   with the p-values).
#'
#' @examples
#' \dontrun{
#' data(simdata)
#' out <- pcgenFast(suffStat = simdata, alpha = 0.01, verbose= FALSE, use.res = TRUE)
#' }
#'
#' @seealso{\code{\link{getResiduals}}}
#'
#' @export
#' @importFrom pcalg skeleton gaussCItest
#' @importFrom methods as
#' @importFrom stats cor
#'
pcgenFast <-
function (suffStat, alpha = 0.01, m.max = Inf,
          res.m.max = Inf, verbose = FALSE, covariates = NULL,
          fixedEdges = NULL, QTLs = integer(), max.iter = 50,
          stop.if.significant = TRUE,
          cov.method = 'uni', use.res = FALSE, return.pvalues = FALSE)
{

  # suffStat=d; alpha=.01; K = NULL; m.max = 3;
  # res.m.max = Inf; verbose = FALSE; covariates = NULL;
  # fixedEdges = NULL; QTLs = integer(); cov.method = 'uni';
  # mean.adj = 'none'; max.iter = 50; stop.if.significant = TRUE; fix.skel = FALSE
  #

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
                      cov.method = cov.method)
  ## estimate CPDAG

  skel.res <- skeleton(suffStat =  list(C = cor(res), n = nrow(res)),
                       indepTest = gaussCItest,
                       alpha=alpha, labels = colnames(suffStat)[-1],
                       verbose = FALSE, m.max = res.m.max,
                       method = "stable")

  gapMatrix <- (as.matrix(as(skel.res, "amat")) == 0)
  gapMatrix <- rbind(rep(FALSE, ncol(res) + 1), cbind(rep(FALSE, ncol(res)), gapMatrix))
  rownames(gapMatrix)[1] <- colnames(gapMatrix)[1] <- 'G'


  fixMatrix <- fixedEdges

  res.cor <- cor(res)


  if (return.pvalues == TRUE) {
    ## estimate CPDAG
    pcgen.fit <- pcgen(suffStat = suffStat, alpha=alpha, verbose = verbose,
                       fixedGaps = gapMatrix, fixedEdges = fixMatrix,
                       covariates = covariates, QTLs = QTLs, m.max = m.max,
                       max.iter = max.iter,
                       stop.if.significant = stop.if.significant,
                       use.res = use.res, res.cor = res.cor,
                       return.pvalues = TRUE)
    pMax <- pcgen.fit[[2]]

    #pMax[-1, -1] <- pmax(pMax[-1, -1], skel.res@pMax)
    inf.select <- which(pMax[-1, -1] == -Inf, arr.ind = T)

    pMax[-1, -1][inf.select] <- skel.res@pMax[inf.select]

    return(list(gr = pcgen.fit[[1]], pMax = pMax))

  } else {

    pcgen.fit <- pcgen(suffStat = suffStat, alpha=alpha, verbose = verbose,
                       fixedGaps = gapMatrix, fixedEdges = fixMatrix,
                       covariates = covariates, QTLs = QTLs, m.max = m.max,
                       max.iter = max.iter,
                       stop.if.significant = stop.if.significant,
                       use.res = use.res, res.cor = res.cor,
                       return.pvalues = FALSE)
    return(pcgen.fit)
  }

}
