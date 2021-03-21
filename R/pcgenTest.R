#' The conditional independence test in pcgen
#'
#' This performs the conditional independence test used in the pcgen algorithm,
#' assuming there are replicates, and independent genetic effects.
#'
#' This will perform the conditional independence test used in pcgen algorithm,
#' assuming (replicates and K = Z Z^t; or no replicates with a generic kinship
#' matrix K = A) and (independent genetic effects)
#'
#' \code{pcgenTest} tests for conditional independence between \eqn{x} and
#' \eqn{y} given \eqn{S}. It distinguishes 2 situations: (i) if one of \eqn{x}
#' and \eqn{y} (say \eqn{x}) is the factor G, \code{pcgenTest} will test if the
#' genetic variance in \eqn{y} is zero, given the traits in S. (ii) if \eqn{x}
#' and \eqn{y} are both traits, \code{pcgenTest} tests if the residual
#' covariance between them is zero, given the traits in \eqn{S} and the factor
#' G. The factor G is automatically included in the conditioning set \eqn{S}
#' (\eqn{S} does not need to contain the integer 1). This test is either based
#' on a bivariate mixed model (when \code{use.res=FALSE}), or on residuals from
#' GBLUP (\code{use.res=T}), obtained with the getResiduals function. In the
#' latter case, \code{res.cor} must be provided.
#'
#' @references * Kruijer, W., Behrouzi, P., Bustos-Korts, D., Rodríguez-Álvarez,
#'   M. X., Mahmoudi, S. M., Yandell, B., ... & van Eeuwijk, F. A. (2020).
#'   Reconstruction of networks with direct and indirect genetic effects.
#'   \emph{Genetics}, 214(4), 781-807.
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @inheritParams pcgen
#'
#' @param x,y Column numbers in \code{suffStat} that should be tested for
#'   conditional independence given the variables in \code{S}.
#'
#' @param S vector of integers defining the conditioning set, where the integers
#'   refer to column numbers in \code{suffStat}. May be numeric(), i.e. the
#'   empty set.
#'
#' @param stop.if.significant If TRUE, the EM-algorithm used in some of the
#'   conditional independence tests will be stopped whenever the p-value becomes
#'   significant, i.e. below alpha. This will speed up calculations, and can be
#'   done because (1) the PC algorithm only needs an accept/reject decision (2)
#'   In EM the likelihood is nondecreasing. Should be put to FALSE if the
#'   precise p-values are of interest.
#'
#' @return a p-value
#'
#' @seealso \code{\link{getResiduals}}, \code{\link[pcalg]{gaussCItest}},
#'   \code{\link[pcalg]{disCItest}}, \code{\link[pcalg]{binCItest}}
#'
#' @examples
#' \dontrun{
#'  data(simdata)
#'  rs <- pcgen::getResiduals(suffStat = simdata)
#'  pcgenTest(suffStat = simdata, x= 2, y= 3, S= 4)
#'  pcgenTest(suffStat = simdata, x= 2, y= 3, S= c(1,4))
#'  pcgenTest(suffStat=simdata, x= 2, y= 3, S= 4, use.res=T, res.cor= cor(rs))
#'  pcgenTest(suffStat = simdata, x= 2, y= 1, S= 4)
#' }
#'
#' @export
#' @importFrom pcalg gaussCItest
#' @importFrom stats lm as.formula

pcgenTest <- function(x, y, S, suffStat, covariates = NULL, QTLs = integer(), K = NULL,
                      replicates = TRUE, use.manova = TRUE, alpha = 0.01, max.iter = 50,
                      stop.if.significant = TRUE, use.res = FALSE, res.cor = NULL) {

  ## Input parameters checks: -------------------------------------------------

  stopifnot(names(suffStat)[1] == 'G') # First col in suffStat should be G
  if (!is.null(covariates)) stopifnot(is.data.frame(covariates))
  if (class(QTLs)!='integer') {stop("QTLs should be a vector of integers")}
  if (1 %in% QTLs) {stop("QTLs should not contain the genotype column (G)")}
  stopifnot(length(unique(c(x, y, S))) == length(c(x, y, S))) # x, y and S uniq
  if (x %in% QTLs & y %in% QTLs)
    stop("Cond. indep. of two QTLs is not tested")
  if ((x %in% QTLs & y == 1) | (y %in% QTLs & x == 1))
    stop("cond. indep. of QTL and G not tested")

  QTLs <- sort(QTLs)

  if (is.null(covariates)) {
    X <- as.data.frame(matrix(0, nrow = nrow(suffStat), ncol = 0))
  } else {
    X <- covariates
  }

  ## Type A CI test: Trait \perp G | {Traits +/- QTLs}: -----------------------
  if (x == 1 | y == 1) {
    if (length(S) != 0) X <- cbind(X, suffStat[,S])

    if (x == 1){
      return(gen.var.test(y = suffStat[, y], X = X, G = suffStat[, 1], K))
    } else {
      return(gen.var.test(y = suffStat[, x], X = X, G = suffStat[, 1], K))
    }

  }
  ## Type D CI test: Trait \perp QTL | {Traits +/- G +/- QTL} ----------------

  if (any(QTLs %in% c(x,y))) {
    if (y %in% QTLs) { # make x the QTL, not y (just to simplify the code)
      temp.x <- x
      x <- y
      y <- temp.x
    }
    sns  <- names(suffStat)

    if (1 %in% S) {
      # Because of the 'hierarchy' (G would 'eat up' all QTL variance), drop 1 from S
      # (in this case only)
      if (is.null(covariates)) {
        pval <- summary(lm(as.formula(paste(sns[y], '~',
                                            paste(sns[c(setdiff(S,1), x)],  collapse='+'))),
                           data = suffStat))$coefficients[length(S)+1,4]
      } else {
        pval <- summary(lm(as.formula(paste(sns[y], '~', paste(names(X), collapse='+'),
                                            '+', paste(sns[c(setdiff(S,1), x)], collapse='+'))),
                           data = cbind(suffStat, X)))$coefficients[ncol(X) + length(S) + 1, 4]
      }
    } else { # G is not in S
      if (is.null(covariates)) {
        pval <- summary(lm(as.formula(paste(sns[y], '~', paste(sns[c(S,x)], collapse='+'))),
                           data = suffStat))$coefficients[length(S)+2,4]
      } else {
        pval <- summary(lm(as.formula(paste(sns[y], '~', paste(names(X), collapse='+'),
                                            '+', paste(sns[c(S,x)], collapse = '+'))),
                           data = cbind(suffStat, X)))$coefficients[ncol(X) + length(S)+2,4]
      }
    }

    return(pval)
  }

  ## Type C CI test: Trait \perp Trait | {Traits +/- QTLs} --------------------

  # Absorbed in type B test below


  ## Type B CI test: Trait \perp Trait | {G + Traits +/- QTLs} ----------------
  # Note that G will account for QTLs that may not be contained in S.

  if (use.res == TRUE) {
    out <- gaussCItest(x = x - 1, y = y - 1, S = (setdiff(S, 1) - 1),
                       suffStat = list(C = res.cor, n = nrow(suffStat)))
  } else { # No residuals
    X <- cbind(X, suffStat[, setdiff(S, 1)])
    out <- res.covar.test(x = suffStat[, x], y = suffStat[, y], G = suffStat[, "G"],
                          K = K, X = X, alpha = alpha, max.iter = max.iter,
                          use.manova = use.manova,
                          replicates = replicates, stop.if.significant = stop.if.significant)
  }
  return(out)
}

