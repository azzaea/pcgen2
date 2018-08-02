#' The pcgenFast algorithm
#'
#' This will run the pcgen algorithm, assuming there are replicates,
#' and K = Z Z^t. Fast version, with residuals-based pre-selection.
#'
#' Here we
#' put a lot
#' of details
#' COMMENT on cov.method (the same option is
#' used in the computation of the residuals and the conditional means)
#' fixedGaps and fixedEdges are NULL
#' m.max, alpha,.... are also used in the first step
#'
#' @inheritParams pcgen
#'
#' @param m.max maximum size of the conditioning set, in pcgen
#'
#' @param res.m.max maximum size of the conditioning set, in the pc-algorithm on the residuals
#'
#' @param use.res If FALSE, residuals are only used for screening. If TRUE,
#'                residuals are also used the
#'                test for conditional independence of 2 traits given a set of
#'                other traits and G
#'
#' @param cov.method (Default 'uni') A string, specifying which method should be
#'                  used to estimate the covariance structure required for the
#'                  computation of conditional means, as well as for the
#'                  computation of residuals in the first step. Options are 'us' (unstructured multi-trait model fitted
#'                  using sommer), and uni' (based on univariate GBLUPs).
#'
#' @return A graph (an object with S3 class \code{"pcgen"})
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references A paper on arxiv
#'
#' @export

pcgenFast <-
function (suffStat, alpha = 0.01, m.max = Inf,
          res.m.max = Inf, verbose = FALSE, covariates = NULL,
          fixedEdges = NULL, QTLs = integer(), max.iter = 50,
          stop.if.significant = TRUE,
          cov.method = 'uni', use.res = FALSE)
{

  # suffStat=d; alpha=.01; K = NULL; NAdelete = TRUE; m.max = 3;
  # res.m.max = Inf; verbose = FALSE; covariates = NULL;
  # fixedEdges = NULL; QTLs = integer(); cov.method = 'uni';
  # mean.adj = 'none'; max.iter = 50; stop.if.significant = TRUE; fix.skel = FALSE
  #

  NAdelete <- TRUE

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
  rownames(gapMatrix)[1] <- colnames(gapMatrix)[1] <- 'genotype'


  fixMatrix <- fixedEdges

  res.cor <- cor(res)

  ## estimate CPDAG
  pcgen.fit <- pcgen(suffStat = suffStat, alpha=alpha, verbose = verbose,
                    fixedGaps = gapMatrix, fixedEdges = fixMatrix,
                    covariates = covariates, QTLs = QTLs, m.max = m.max,
                    max.iter = max.iter,
                    stop.if.significant = stop.if.significant,
                    use.res = use.res, res.cor = res.cor)


  return(pcgen.fit)
}
