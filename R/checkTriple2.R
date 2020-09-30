checkTriple2 <- function (a, b, c, nbrsA, nbrsC, sepsetA, sepsetC, suffStat,
                          alpha, version.unf = c(NA, NA),
                          maj.rule = FALSE,verbose = FALSE,
                          covariates=NULL, QTLs = integer(), K = NULL,
                          max.iter = 50, stop.if.significant = TRUE,
                          use.res = FALSE, res.cor = NULL)
{
##16-1-18## : added Vg = Vg, Ve = Ve, dec = dec

#' param a, b, c, nbrsA, nbrsC, sepsetA, sepsetC, version.unf, maj.rule, verbose:
#                            as in the original checkTriple function
#
#' param suffStat data.frame, of which the first column is the factor genotype,
#                 and subsequent columns contain the traits. The name of the
#                 first column should be genotype
#' param alpha (Default 0.01) The significance level used in the test. The test
#              itself of course does not depend on this, but it is used in the
#              EM-algorithm to speed up calculations. More precisely, the
#              EM-algorithm is stopped once the P value is below the
#              significance level. This can be done because the PC algorithm
#              only needs an accept/ reject decision.
#' param QTLs column numbers in suffStat that correspond to QTLs
#             These may be partly in S and x and y, but not x and y both in QTLs!
#             Note: the factor genotype (column number 1) may occur in S, as well as x and y
#' param K (Default NULL) The kinship (i.e Genetic Relationship Matrix)
#' param covariates data.frame containing covariates, that should always be used
#        in each conditional independence test. Should be either NULL (default)
#        or a data.frame with the same number of rows as suffStat
#
#' param Vg, Ve (Default NULL) Genetic and residual covariance matrices of
#        dimension (ncol(suffStat)-1) x (ncol(suffStat)-1). Required if
#        cov.means equals 'exact'.
#
#' param dec (Default NULL) Contains a spectral decomposition of K = Z Z^t.
#            Should be a list with components Dk (a vector with the eigenvalues)
#            and Uk (a matrix with eigenvectors). Obtained as follows:
#            K <- Z %*% t(Z); w   <- eigen(K); Dk  <- diag(w$values); Uk  <- w$vectors; dec <- list(Dk = Dk, Uk = Uk)


   ##!##
   # Modified the original checkTriple function, by adding the following lines:
   non.collider.nodes <- sort(c(1,QTLs))
   if (any(non.collider.nodes %in% c(a,b,c))) {return(lapply(list(decision = 2, version = 1, SepsetA = sepsetA,
              SepsetC = sepsetC), as.integer))}
   #####

    nr.indep <- 0
    stopifnot(length(version.unf) == 2, version.unf %in% 1:2)
    tmp <- if (version.unf[2] == 2)
        (b %in% sepsetA || b %in% sepsetC)
    version <- 0
    if ((nn <- length(nbrsA)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        for (i in 1:nrow(allComb)) {
            S <- nbrsA[which(allComb[i, ] != 0)]

           ##!##
           # Modified the original checkTriple function:
           # Instead of pval <- indepTest(a, c, S, suffStat), we use:
            #cat("===> x=", a, " y=", c, " S=", S,"\n")
            pval <- pcgenTest(x=a, y=c, S=S, suffStat,
                                            covariates=covariates,
                                            QTLs = QTLs, K = K,
                                            alpha = alpha,
                                            max.iter = max.iter,
                                            stop.if.significant = stop.if.significant,
                                            use.res = use.res,
                                            res.cor = res.cor)
                      ##16-1-18## : addedd Vg = Vg, Ve = Ve, dec = dec)
            #####

            if (verbose)
                cat("a: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if ((nn <- length(nbrsC)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        for (i in 1:nrow(allComb)) {

            S <- nbrsC[which(allComb[i, ] != 0)]
           ##!##
           # Modified the original checkTriple function:
           # Instead of pval <- indepTest(a, c, S, suffStat), we use:
            #cat("===> x=", a, " y=", c, " S=", S,"\n")
            pval <- pcgenTest(x=a, y=c, S=S, suffStat,
                                            covariates=covariates,
                                            QTLs = QTLs,
                                            alpha = alpha,
                                            max.iter = max.iter,
                                            stop.if.significant = stop.if.significant,
                                            use.res = use.res,
                                            res.cor = res.cor)
                      ##16-1-18## : addedd Vg = Vg, Ve = Ve, dec = dec)
            #####

            if (verbose)
                cat("c: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if (version.unf[1] == 2 && nr.indep == 0) {
        version <- 2
    }
    if (is.null(tmp))
        tmp <- FALSE
    if (all(tmp)) {
        res <- 2
        if (b %nin% sepsetA)
            sepsetA <- c(sepsetA, b)
        if (b %nin% sepsetC)
            sepsetC <- c(sepsetC, b)
    }
    else {
        if (all(!tmp)) {
            res <- 1
            sepsetA <- setdiff(sepsetA, b)
            sepsetC <- setdiff(sepsetC, b)
        }
        else {
            if (!maj.rule) {
                res <- 3
            }
            else {
                if (sum(tmp)/length(tmp) < 0.5) {
                  res <- 1
                  sepsetA <- setdiff(sepsetA, b)
                  sepsetC <- setdiff(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) > 0.5) {
                  res <- 2
                  if (b %nin% sepsetA)
                    sepsetA <- c(sepsetA, b)
                  if (b %nin% sepsetC)
                    sepsetC <- c(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) == 0.5) {
                  res <- 3
                }
            }
        }
    }
    if (verbose && res == 3)
        cat("Triple ambiguous\n")
    lapply(list(decision = res, version = version, SepsetA = sepsetA,
        SepsetC = sepsetC), as.integer)
}
