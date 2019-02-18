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
pcgenTest <- function(x,y,S,suffStat, QTLs=integer(), covariates=NULL,
                          alpha = 0.01, max.iter = 50,
                          stop.if.significant = TRUE,
                          use.res = FALSE, res.cor = NULL)
{
  # suffStat = d; x = 2; y = 8; S = c(1,11)
  # QTLs=integer(); covariates=NULL;
  # skip.nongenetic=1; genVar=rep(TRUE, ncol(suffStat));
  # qtlVar=rep(list(rep(TRUE, ncol(suffStat))); length(QTLs));
  # gen.var.exact=FALSE; gen.covar.exact=TRUE; alpha = 0.01;
  # mean.adj = 'none'; Vg = NULL; Ve = NULL; dec = NULL;
  # max.iter = 50; stop.if.significant = TRUE;
  # return.cond.mean = FALSE

  #suffStat = dr[, c(1, 26:36)]; covariates = dr[, 37:38]; x = 6; y = 7; S = 8:9; QTLs=integer(); covariates=NULL; alpha = 0.01; mean.adj = 'none'; Vg = NULL; Ve = NULL; dec = NULL; max.iter = 50;stop.if.significant = TRUE; return.cond.mean = FALSE

  #### First, some checking...

  if (class(QTLs)!='integer') {stop('QTLs should be a vector of integers')}
  if (1 %in% QTLs) {stop('QTLs should not contain the genotype column (G)')}

  # stop if the first column in suffStat does not have the name G
  stopifnot(names(suffStat)[1] == 'G')
  # stop if the same variables occurs more than once
  stopifnot(length(unique(c(x,y,S)))==length(c(x,y,S)))
  # 1 should always be the genotype factor (G); stop if it is in the QTL vector
  if (1 %in% QTLs) {stop()}
  # stop if both x and y are QTLs (cond. indep. of two QTLs is not tested)
  if (x %in% QTLs & y %in% QTLs) {stop()}
  # Neither we test cond. indep. of a QTL and genotype
  if (x %in% QTLs & y==1) {stop()}
  if (y %in% QTLs & x==1) {stop()}

  # put the QTL column numbers in ascending order
  QTLs <- sort(QTLs)

  ###################################################################
  # Define a data.frame X, containing the covariates. It is also defined when
  # there are no covariates, in which case it will have zero columns.

  if (is.null(covariates)) {
    X <- as.data.frame(matrix(0,nrow(suffStat),0))
  } else {
    X <- as.data.frame(covariates)
    names(X) <- paste0('covariate_', 1:ncol(covariates))
  }

  ######################################


  ###################################################################
  # The case x and y both represent real traits (no QTLs), and 1
  # (direct genetic effects) is not in S. S may contain QTLs

  if (!(1 %in% c(x,y,S)) & !(any(c(x,y) %in% QTLs))) {S <- c(1, S)}

  ###################################################################
  # The case where one of x and y is the number 1 (direct genetic effects).
  # S may also contain QTLs

  if (x==1 | y==1) {

    if (x==1) {

      if (length(S)==0) {

        return(gen.var.test(y=suffStat[,y],X=X,G=suffStat[,1]))

      } else {

        X <- as.data.frame(cbind(X, suffStat[,S]))

        return(gen.var.test(y=suffStat[,y], X = X, G=suffStat[,1]))

      }
    }

    if (y==1) {

      if (length(S)==0) {

        return(gen.var.test(y=suffStat[,x],X=X,G=suffStat[,1]))

      } else {

        X <- as.data.frame(cbind(X, suffStat[,S]))

        return(gen.var.test(y=suffStat[,x], X = X, G=suffStat[,1]))

      }

    }

  }

  ########################################################################
  # The case where one of x and y is a QTL;
  # the number 1 (direct genetic effects) may or may not be in S.
  # S may also contain additional QTLs

  if (any(QTLs %in% c(x,y))) {

    if (y %in% QTLs) {
      # In this case, make sure that x is the QTL, not y (just to simplify the code)
      temp.x <- x
      x <- y
      y <- temp.x
    }

    sns  <- names(suffStat)

    if (1 %in% S) {

      # Because of the 'hierarchy' (G would 'eat up' all QTL variance), drop 1 from S (in this case only)
      if (is.null(covariates)) {
        pval <- summary(lm(as.formula(paste(sns[y],'~',paste(sns[c(setdiff(S,1),x)], collapse='+'))), data=suffStat))$coefficients[length(S)+1,4]
      } else {
        pval <- summary(lm(as.formula(paste(sns[y],'~',paste(names(X), collapse='+'),'+',paste(sns[c(setdiff(S,1),x)], collapse='+'))), data=cbind(suffStat,X)))$coefficients[ncol(X) + length(S)+1,4]
      }

    } else {

      if (is.null(covariates)) {
        pval <- summary(lm(as.formula(paste(sns[y],'~',paste(sns[c(S,x)], collapse='+'))), data=suffStat))$coefficients[length(S)+2,4]
      } else {
        pval <- summary(lm(as.formula(paste(sns[y],'~',paste(names(X), collapse='+'),'+',paste(sns[c(S,x)], collapse='+'))), data=cbind(suffStat,X)))$coefficients[ncol(X) + length(S)+2,4]
      }

    }

    return(pval)
  }

  ###################################################################
  # The case where x and y both represent real traits (no QTLs), and 1
  # (direct genetic effects: G) IS contained in S. S may also contain QTLs.
  # We do not skip any tests here.
  # Note that G will account for QTLs that may not be contained in S.


    #if (1 %in% S) {  # redundant; at this point, 1 MUST be contained in S

  if (use.res == TRUE) {

    out <- gaussCItest(x = x - 1, y = y - 1, S = (setdiff(S, 1) - 1),
                       suffStat = list(C = res.cor, n = nrow(suffStat)))

  } else {

    if (length(S)==1) {

      out <- res.covar.test(x=suffStat[,x],y=suffStat[,y], X=X,
                            G=suffStat[,1], alpha = alpha)

    } else {

      X <- as.data.frame(cbind(X,suffStat[,setdiff(S,1)]))

      out <- res.covar.test(x = suffStat[,x], y = suffStat[,y],
                            X = X, G = suffStat[,1], alpha = alpha,
                            max.iter = max.iter,
                            stop.if.significant = stop.if.significant)
    }

  }
  return(out[1])

}

