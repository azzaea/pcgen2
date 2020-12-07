#' Infer the DAG Of a Phenotypic Traits' network with Genetic Effects
#'
#' @import pcalg
#' @keywords internal
#' @inheritParams pcgen
#' @inheritParams pcalg::skeleton
#'
skeleton2  <- function (suffStat, alpha, labels, p, method = c("stable", "original"),
                        m.max = Inf, fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                        verbose = FALSE, covariates=NULL, QTLs = integer(), K = NULL,
                        max.iter = 50, stop.if.significant = TRUE, use.res = FALSE,
                        res.cor = NULL) {

  cl <- match.call()

  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    stopifnot(nrow(covariates)==nrow(suffStat))
  }

  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)

  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { #labels provided
    stopifnot(is.character(labels))
    if (missing(p)) p <- length(labels)
    if (p != length(labels))
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
    methods::new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    methods::as(G, "graphNEL")
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
