pc.cons.intern2 <- function (sk, suffStat, alpha, version.unf = c(NA, NA),
                             maj.rule = FALSE, verbose = FALSE,
                             covariates=NULL,
                             QTLs = QTLs, K = NULL,
                             max.iter = 50, stop.if.significant = TRUE,
                             use.res = FALSE, res.cor = NULL)
{
##16-1-18## : added Vg = Vg, Ve = Ve, dec = dec

#' param sk, version.unf, maj.rule, verbose: as in the original pc.cons.intern function
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
#                   in each conditional independence test. Should be either NULL (default)
#                   or a data.frame with the same number of rows as suffStat
#
# sk = pc.fit2; version.unf = c(2, 1); qtlVar = rep(list(rep(TRUE, ncol(suffStat))), length(QTLs))

   ##!##
   # Modified the original pc.cons.intern function, by adding the following lines:
   non.collider.nodes <- sort(c(1,QTLs))
   #####

    g <- as(sk@graph, "matrix")
    stopifnot(all(g == t(g)))
    p <- as.numeric(dim(g)[1])
    unfTripl <- vers <- rep(NA, min(p * p, 1e+05))
    counter  <- 0

    if (sum(g) > 0) {

        ind <- which(g == 1, arr.ind = TRUE)
        tripleMatrix <- NULL

        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            allC <- setdiff(which(g[b, ] == 1), a)
            newC <- allC[g[a, allC] == 0]
            tmpMatrix <- cbind(rep(a, length(newC)), rep(b, length(newC)),
                newC)
            tripleMatrix <- rbind(tripleMatrix, tmpMatrix)
            colnames(tripleMatrix) <- c("", "", "")
        }

        if ((m <- nrow(tripleMatrix)) > 0) {

            deleteDupl <- logical(m)

            for (i in seq_len(m)) if (tripleMatrix[i, 1] > tripleMatrix[i,3])
                deleteDupl[i] <- TRUE

            if (any(deleteDupl))
                tripleMatrix <- tripleMatrix[!deleteDupl, , drop = FALSE]



            for (i in seq_len(nrow(tripleMatrix))) {

                if (counter + 1L == length(unfTripl)) {
                  n.xtra <- min(p * p, 1e+05)
                  new.len <- counter + 1L + n.xtra
                  length(unfTripl) <- new.len
                  length(vers) <- new.len
                }
                a <- tripleMatrix[i, 1]
                b <- tripleMatrix[i, 2]
                c <- tripleMatrix[i, 3]
                nbrsA <- which(g[, a] != 0)
                nbrsC <- which(g[, c] != 0)
                if (verbose) {
                  cat("\nTriple:", a, b, c, "and sepset by skelet:",
                    unique(sk@sepset[[a]][[c]], sk@sepset[[c]][[a]]),
                    "\n")
                }

                ##!##
                # Modified the original pc.cons.intern function:
                # 1. Instead of directly calling checkTriple, we first check if
                #    a and c are both genetic nodes (G or a QTL), in which case
                #    the decision is 1.
                # 2. If not, we call the checkTriple2 function, instead of checkTriple
                #    This requires all the additional arguments required for our
                #    independence test (pcgenTest)

                if (a %in% non.collider.nodes & c %in% non.collider.nodes) {

                  r.abc <- list()
                  r.abc$decision <- 1
                  r.abc$version <- 2
                  r.abc$SepsetA <- integer()
                  r.abc$SepsetC <- integer()

                } else {

                  r.abc <- checkTriple2(a, b, c, nbrsA, nbrsC, sk@sepset[[a]][[c]],
                                        sk@sepset[[c]][[a]],
                                        suffStat = suffStat,
                                        alpha = alpha,
                                        version.unf = version.unf,
                                        maj.rule = maj.rule,
                                        verbose = verbose,
                                        covariates=covariates,
                                        QTLs = QTLs, K = K,
                                        max.iter = max.iter,
                                        stop.if.significant = stop.if.significant,
                                        use.res = use.res,
                                        res.cor = res.cor)

                }
                #
                ####

                if (r.abc$decision == 3) {
                  counter <- counter + 1
                  unfTripl[counter] <- triple2numb(p, a, b, c)
                  vers[counter] <- r.abc$version
                }
                if ((version.unf[1] == 2) && (r.abc$version ==
                  2) && (r.abc$decision != 3)) {
                  counter <- counter + 1
                  unfTripl[counter] <- triple2numb(p, a, b, c)
                  vers[counter] <- r.abc$version
                }
                sk@sepset[[a]][[c]] <- r.abc$SepsetA
                sk@sepset[[c]][[a]] <- r.abc$SepsetC

            }
        }
    }

    length(unfTripl) <- length(vers) <- counter
    list(unfTripl = unfTripl, vers = vers, sk = sk)
}
