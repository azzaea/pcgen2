skeleton2  <-
function (suffStat, alpha, labels, p, method = c("stable",
    "original"), m.max = Inf, fixedGaps = NULL,
    fixedEdges = NULL, NAdelete = TRUE, verbose = FALSE,
    covariates=NULL, QTLs = integer(), dec = NULL,
    max.iter = 50, stop.if.significant = TRUE,
    use.res = FALSE, res.cor = NULL)
{
##16-1-18## : added Vg = Vg, Ve = Ve, dec = dec


#' param labels, p, method, m.max, fixedGaps, fixedEdges, NAdelete, verbose: as in the original skeleton function
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
#' param covariates data.frame containing covariates, that should always be used
#                   in each conditional independence test. Should be either NULL (default)
#                   or a data.frame with the same number of rows as suffStat
#
#' param Vg, Ve (Default NULL) Genetic and residual covariance matrices of
#        dimension (ncol(suffStat)-1) x (ncol(suffStat)-1). Required if
#        cov.means equals 'exact'.
#
#' param dec (Default NULL) Contains a spectral decomposition of K = Z Z^t.
#            Should be a list with components Dk (a vector with the eigenvalues)
#            and Uk (a matrix with eigenvectors). Obtained as follows:
#            K <- Z %*% t(Z); w   <- eigen(K); Dk  <- diag(w$values); Uk  <- w$vectors; dec <- list(Dk = Dk, Uk = Uk)


    #if (method == "stable.fast") {stop('The stable.fast option is not yet implemented.')}

    cl <- match.call()

    if (!is.null(covariates)) {
      covariates <- as.data.frame(covariates)
      stopifnot(nrow(covariates)==nrow(suffStat))
    }

    if (!missing(p))
        stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
            1, p >= 2)
    if (missing(labels)) {
        if (missing(p))
            stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    } else {
        stopifnot(is.character(labels))
        if (missing(p))
            p <- length(labels)
        else if (p != length(labels))
            stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    }

    seq_p    <- seq_len(p)
    method   <- match.arg(method)

    ##!##
    # Modified the original skeleton function, by replacing the chunk below (commented out)
    # by the following lines:

    if (is.null(fixedGaps)) {
      G <- matrix(TRUE, p, p)
    } else {
      stopifnot(identical(dim(fixedGaps), c(p, p)))
      stopifnot(identical(fixedGaps, t(fixedGaps)))
      G <- !fixedGaps
    }

    diag(G) <- FALSE

    # The main point of the modification: regardless of the value of fixedGaps,
    # there should never be edges between QTLs, or between QTLs and genotype :
    if (length(QTLs) > 0) {
      G[c(1,QTLs),c(1,QTLs)] <- FALSE
    }

    # The original pcalg chunk of code :
    #if (is.null(fixedGaps)) {
    #    G <- matrix(TRUE, nrow = p, ncol = p)
    #}
    #else if (!identical(dim(fixedGaps), c(p, p)))
    #    stop("Dimensions of the dataset and fixedGaps do not agree.")
    #else if (!identical(fixedGaps, t(fixedGaps)))
    #    stop("fixedGaps must be symmetric")
    #else G <- !fixedGaps
    #diag(G) <- FALSE
    #########################

    if (any(is.null(fixedEdges))) {
        fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    }
    else if (!identical(dim(fixedEdges), c(p, p)))
        stop("Dimensions of the dataset and fixedEdges do not agree.")
    else if (!identical(fixedEdges, t(fixedEdges)))
        stop("fixedEdges must be symmetric")

	# We have inactivated for this version because the stable.fast option is not yet implemented.
    #if (method == "stable.fast") {
    #    if (identical(indepTest, gaussCItest))
    #        indepTestName <- "gauss"
    #    else indepTestName <- "rfun"
    #    options <- list(verbose = as.integer(verbose), m.max = as.integer(ifelse(is.infinite(m.max),
    #        p, m.max)), NAdelete = NAdelete)
    #    res <- .Call("estimateSkeleton", G, suffStat, indepTestName,
    #        indepTest, alpha, fixedEdges, options)
    #    G <- res$amat
    #    sepset <- lapply(seq_p, function(i) c(lapply(res$sepset[[i]],
    #        function(v) if (identical(v, as.integer(-1))) NULL else v),
    #        vector("list", p - length(res$sepset[[i]]))))
    #    pMax <- res$pMax
    #    n.edgetests <- res$n.edgetests
    #    ord <- length(n.edgetests) - 1L
    #}
    #else {
        pval <- NULL
        sepset <- lapply(seq_p, function(.) vector("list", p))
        pMax <- matrix(-Inf, nrow = p, ncol = p)
        diag(pMax) <- 1
        done <- FALSE
        ord <- 0L
        n.edgetests <- numeric(1)


        #browser()

        while (!done && any(G) && ord <= m.max) {
            n.edgetests[ord1 <- ord + 1L] <- 0
            done <- TRUE
            ind <- which(G, arr.ind = TRUE)
            ind <- ind[order(ind[, 1]), ]
            remEdges <- nrow(ind)
            #print(ind)
            if (verbose)
                cat("Order=", ord, "; remaining edges:", remEdges,
                  "\n", sep = "")
            if (method == "stable") {
                G.l <- split(G, gl(p, p))
            }
            for (i in 1:remEdges) {
                if (verbose && (verbose >= 2 || i%%100 == 0))
                  cat("|i=", i, "|iMax=", remEdges, "\n")
                x <- ind[i, 1]
                y <- ind[i, 2]
                if (G[y, x] && !fixedEdges[y, x]) {
                  nbrsBool <- if (method == "stable")
                    G.l[[x]]
                  else G[, x]
                  nbrsBool[y] <- FALSE
                  nbrs <- seq_p[nbrsBool]
                  length_nbrs <- length(nbrs)
                  if (length_nbrs >= ord) {
                    if (length_nbrs > ord)
                      done <- FALSE
                    S <- seq_len(ord)
                    repeat {
                      n.edgetests[ord1] <- n.edgetests[ord1] + 1
                      #pval <- indepTest(x, y, nbrs[S], suffStat)


                      if (verbose) {
                        #cat("x=", labels[x], " y=", labels[y], " S=", labels[nbrs[S]],": pval =")
                        cat("x=", x, " y=", y, " S=", nbrs[S],": pval =")
                      }
                      ##!##
                      # Modified the original skeleton function: pcgenTest instead of indepTest:
                      pval <- pcgenTest(x, y, S=nbrs[S], suffStat,
                                            covariates=covariates,
                                            QTLs = QTLs,
                                            alpha = alpha,
                                            max.iter = max.iter,
                                            stop.if.significant = stop.if.significant,
                                            use.res = use.res,
                                            res.cor = res.cor)
                      ##16-1-18## : addedd Vg = Vg, Ve = Ve, dec = dec

                      if (verbose) {
                        ##!##
                        # Modified the original skeleton function: if an edge is removed,
                        # mark this with <<<<<<<<<< >>>>>>>>>>>>>>>>
                        if (pval >= alpha) {cat('<<<<<<< ',pval,' >>>>>>>>>', "\n")} else {cat(pval, "\n")}

                        #if (pval > 1.68 * 10^(-4) & pval < 1.7 * 10^(-4)) {cat('!!!!!!!!!!!!', "\n"); stop()} 
                                                #
                        #cat("x=", x, " y=", y, " S=", nbrs[S],": pval =", pval, "\n")
                      }

                      ################################

                      if (is.na(pval))
                        pval <- as.numeric(NAdelete)
                      if (pMax[x, y] < pval)
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        G[x, y] <- G[y, x] <- FALSE
                        sepset[[x]][[y]] <- nbrs[S]
                        break
                      }
                      else {
                        nextSet <- getNextSet(length_nbrs, ord,
                          S)
                        if (nextSet$wasLast)
                          break
                        S <- nextSet$nextSet
                      }
                    }
                  }
                }
            }
            ord <- ord + 1L
        }
        for (i in 1:(p - 1)) {
            for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i,
                j], pMax[j, i])
        }
    #}
    Gobject <- if (sum(G) == 0) {
        new("graphNEL", nodes = labels)
    }
    else {
        colnames(G) <- rownames(G) <- labels
        as(G, "graphNEL")
    }

    ##!##
    # Modified the original skeleton function: we do not only return the skeleton
    # and separating sets, but also genVar and qtlVar
    out1 <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
                max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
                sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))

    ##!##
    # Modified the original skeleton function:
    # Do something with the separating sets for the QTLs (check what !!)
    if (length(QTLs) > 0) {
      for (gg1 in c(1, QTLs)) {
        for (gg2 in c(1, QTLs)) {
          out1@sepset[[gg1]][[gg2]] <- integer()
        }
      }
    }
    #########################################

    return(skel.out = out1)
}
