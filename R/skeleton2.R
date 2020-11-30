skeleton2  <-
  function (suffStat, alpha, labels, p, method = c("stable", "original"), 
            m.max = Inf, fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, 
            verbose = FALSE, covariates=NULL, QTLs = integer(), K = NULL, dec = NULL,
            max.iter = 50, stop.if.significant = TRUE, use.res = FALSE, 
            res.cor = NULL) {

    #' param suffStat data.frame, of which the first column is the factor genotype,
    #                 and subsequent columns contain the traits. The name of the
    #                 first column should be genotype
    #
    #' param alpha (Default 0.01) The significance level used in the test. The test
    #              itself of course does not depend on this, but it is used in the
    #              EM-algorithm to speed up calculations. More precisely, the
    #              EM-algorithm is stopped once the P value is below the
    #              significance level. This can be done because the PC algorithm
    #              only needs an accept/ reject decision.
    #
    #' param labels, p, method, m.max, fixedGaps, fixedEdges, NAdelete, verbose: as
    #'                                            in the original skeleton function
    #
    #' param covariates data.frame containing covariates, that should always be used
    #                   in each conditional independence test. Should be either NULL (default)
    #                   or a data.frame with the same number of rows as suffStat
    #
    #' param QTLs column numbers in suffStat that correspond to QTLs
    #             These may be partly in S, x and y, but not x and y both in QTLs!
    #             Note: the factor genotype (column number 1) may occur in S, 
    #                   as well as x and y
    #
    #' param K (Default NULL) The kinship (i.e Genetic Relationship Matrix)

    #if (method == "stable.fast") 
    # {stop('The stable.fast option is not yet implemented.')}
    
    cl <- match.call()
    
    if (!is.null(covariates)) {
      covariates <- as.data.frame(covariates)
      stopifnot(nrow(covariates)==nrow(suffStat))
    }
    
    # Not REALLY needed. skeleton2 is only called from pcgen with appropriate
    # parameters ******** Azza
    
    if (!missing(p))
      stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
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
    # Modified the original skeleton function, so that regardless of the value of fixedGaps,
    # there should never be edges between QTLs, or between QTLs and genotype :
    
    if (is.null(fixedGaps)) {
      G <- matrix(TRUE, p, p)
    } else {
      stopifnot(identical(dim(fixedGaps), c(p, p)))
      stopifnot(identical(fixedGaps, t(fixedGaps)))
      G <- !fixedGaps
    }
    diag(G) <- FALSE
    if (length(QTLs) > 0) {
      G[c(1,QTLs),c(1,QTLs)] <- FALSE
    }

    #########################
    
    if (any(is.null(fixedEdges))) {
      fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    } else 
      if (!identical(dim(fixedEdges), c(p, p))) {
        stop("Dimensions of the dataset and fixedEdges do not agree.")
      } else if (!identical(fixedEdges, t(fixedEdges)))
        stop("fixedEdges must be symmetric")
    
    ##!##
    # We have inactivated for this version because the stable.fast option is not yet implemented.
    # stopifnot((is.integer(numCores) || is.numeric(numCores)) && 
    #             numCores > 0)
    # if (numCores > 1 && method != "stable.fast") {
    #   warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
    # }
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
    
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      # ind <- which((G*upper.tri(G))>0, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,
            "\n", sep = "")
      if (method == "stable") G.l <- split(G, gl(p, p))
      
      for (i in 1:remEdges) {
        if (verbose && (verbose >= 2 || i%%100 == 0))
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]; y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable")  G.l[[x]] else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              ##!##: pcgenTest instead of indepTest:
              pval <- pcgenTest(x, y, S = nbrs[S], suffStat, covariates = covariates,
                                QTLs = QTLs, K = K, alpha = alpha, max.iter = max.iter, 
                                stop.if.significant = stop.if.significant,
                                use.res = use.res, res.cor = res.cor)
              if (verbose) {
                ##!## 1. Print labels of variables (not their indices in the suffStat dataframe)
                cat("x=", labels[x], " y=", labels[y], " S=", labels[nbrs[S]],": pval =")
                ##!## 2. if an edge is removed, mark this with <<<<<<<<<< >>>>>>>>>>>>>>>>
                if (!(is.na(pval)) & pval >= alpha) {
                  cat('<<<<<<< ', pval, ' >>>>>>>>>', "\n")
                } else 
                  cat(pval, "\n")
              } # End if (verbose)
              
              if (is.na(pval))
                pval <- as.numeric(NAdelete)
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if (pval >= alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              } else { # if (pval < alpha)
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              } # End pval (?) alpha comparison
            } #End repeat
          } # End if (length_nbrs >= ord) 
        } # End if (G[y, x] && !fixedEdges[y, x])
      } # End for (i in 1:remEdges) loop
      ord <- ord + 1L
    } # End while (!done && any(G) && ord <= m.max)
    for (i in 1:(p - 1)) 
      for (j in 2:p) 
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, i])
    
    Gobject <- if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G, "graphNEL")
    }
    
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
