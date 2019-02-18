#' The pc algorithm applied to residuals
#'
#' The pc algorithm applied to residuals, assuming there are replicates,
#' and K = Z Z^t.
#'
#' Here we
#' put a lot
#' of details
#' COMMENT
#'
#' @inheritParams pcgen
#'
#' @param m.max maximum size of the conditioning set, in the pc-algorithm on the residuals
#'
#' @param K kinship
#'
#' @param cov.method (Default 'fe') A string, specifying which method should be
#'                  used to compute the G-BLUP. Options are 'us'
#'                  (unstructured multi-trait model fitted
#'                  using sommer) and uni' (based on univariate GBLUPs).
#'
#' @param use.GBLUP Use the GBLUP itself, instead of the residuals (as in Topner et al)
#'
#'                
#' @param return.pvalues   
#' 
#' @return A graph (an object with S3 class \code{"pcgen"})
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references A paper on arxiv and Topner et al (2017)
#'
#' @export

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
