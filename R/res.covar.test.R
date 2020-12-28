#' EXACT test for residual correlation, based on a bivariate mixed model
#'
#'
#' @inheritParams pcgenTest
#' @param x,y  trait vectors whose independence is to be tested, conditional on
#'   G (genotype) and possible covariates in X. Manadatory
#'
#' @param G Factor of genotypes. Manadatotry
#'
#' @param Z.t Defaults to \code{NULL}
#'
#' @param X A dataframe of possibly the covariates, traits and QTLs used as a
#'   conditioning set. If not specified, this becomes a dependence test
#'
#' @param use.manova
#'
#' @importFrom Matrix Matrix

res.covar.test <- function(x, y, G, Z.t = NULL, K = NULL,
                           X = as.data.frame(matrix(0,length(x),0)),
                           alpha = 0.01, use.manova = TRUE,
                           max.iter = 50, stop.if.significant = TRUE) {

  if (is.null(Z.t)) {Z.t <- make.Z.matrix(G)}

  X.t           <- Matrix(cbind(rep(1, length(G)), as.matrix(X)))
  em.vec        <- c(x, y)
  names(em.vec) <- rep(as.character(G), 2)

  if (!is.null(K)) use.manova <- FALSE # seems MANOVA not suitable for generic K

  if (!use.manova) {
    fit.reduced <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = FALSE, # maybe false
                         cov.gen = FALSE, max.iter = 500)
  } else {
    X  <- as.data.frame(X)
    dd <- data.frame(genotype = G, x = x, y = y, X)
    if (ncol(X) > 0) {names(dd)[-(1:3)] <- paste0('cov', 1:(ncol(X)))}
    if (ncol(X) > 0) {names(X) <- paste0('cov', 1:(ncol(X)))}

    n.rep.vector<- as.integer(table(dd$genotype))
    n.geno      <- length(n.rep.vector)
    n.rep       <- (sum(n.rep.vector) - sum(n.rep.vector^2)/sum(n.rep.vector))/(length(n.rep.vector) - 1)
    manova.out  <- manova(as.formula(paste('cbind(x,y) ~',paste(c(names(X), 'genotype'), collapse='+'))), data = dd)
    #summary.out <- summary(manova.out)

    summary.out <-  tryCatch(summary(manova.out), error = function(e) e)

    if (any(class(summary.out) == "error")) {
      Vg.manova <- Ve.manova <- matrix(0, 2, 2)
      aov.x <- anova(lm(as.formula(paste('x~',paste(c(names(X), 'genotype'), collapse='+'))), data = dd))
      aov.y <- anova(lm(as.formula(paste('y~',paste(c(names(X), 'genotype'), collapse='+'))), data = dd))

      MS.res.x <- aov.x[["Mean Sq"]][ncol(dd)-1]
      MS.res.y <- aov.y[["Mean Sq"]][ncol(dd)-1]
      MS.gen.x <- aov.x[["Mean Sq"]][ncol(dd)-2]
      MS.gen.y <- aov.y[["Mean Sq"]][ncol(dd)-2]

      Vg.manova[1,1] <- (MS.gen.x - MS.res.x) / n.rep
      Vg.manova[2,2] <- (MS.gen.y - MS.res.y) / n.rep
      Ve.manova[1,1] <- MS.res.x
      Ve.manova[2,2] <- MS.res.y
    } else {
      d.f         <- as.numeric(summary.out[['stats']][,'Df'][c('genotype','Residuals')])
      MS.gen      <- summary.out[['SS']][['genotype']] / d.f[1]
      MS.res      <- summary.out[['SS']][['Residuals']] / d.f[2]
      Vg.manova   <- (MS.gen - MS.res) / n.rep
      Ve.manova   <- MS.res
    }

    # Negative estimates of genetic variance are set to zero
    diag(Vg.manova)[diag(Vg.manova) < 0] <- 0

    if (Vg.manova[1,1]==0 | Vg.manova[2,2]==0) {
      Vg.manova[1,2] <- Vg.manova[2,1] <- 0
      if (Vg.manova[1,1]==0) {Vg.manova[1,1] <- 0.01 * Ve.manova[1,1]}
      if (Vg.manova[2,2]==0) {Vg.manova[2,2] <- 0.01 * Ve.manova[2,2]}
    }

    fit.reduced <- fitEM(y = em.vec, X.t = X.t, Z.t = Z.t, K = K,
                         Vg.start = as.numeric(Vg.manova)[c(1,4,2)],
                         Ve.start = c(as.numeric(Ve.manova)[c(1,4)], 0),
                         cov.error = TRUE, cov.gen = FALSE, max.iter = 5)
  }

  fit.full    <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = TRUE, alpha = alpha,
                       stop.if.significant = TRUE, max.iter = max.iter,
                       null.dev = fit.reduced$deviance,
                       Vg.start = fit.reduced$variances$Vg,
                       Ve.start = fit.reduced$variances$Ve)

  #cat('\n\nIt:',fit.full$it, '\n\n')

  loglik_Full    <- -0.5*fit.full$deviance
  loglik_Reduced <- -0.5*fit.reduced$deviance

  REMLLRT <- 2 * max(loglik_Full - loglik_Reduced, 0)
  pvalue  <- (1 - pchisq(REMLLRT, df = 1))

  #cat('\n\nDiff:',pvalue - fit.full$pvalue, '\n\n')

  return(pvalue) # TRUE is never actually used. I removed it! Azza

}
