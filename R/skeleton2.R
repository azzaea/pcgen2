#' Infer the skeleton of a DAG of Phenotypes with Genetic Effects.
#'
#' The inferred skeleton is order-independent (equivalent to pcalg::skeleton
#' with skel.method = "stable") following the “PC-stable” approach proposed by
#' Colombo and Maathuis (2014).
#'
#' @references * Colombo, D., & Maathuis, M. H. (2014). Order-independent
#'   constraint-based causal structure learning. \emph{The Journal of Machine
#'   Learning Research}, 15(1), 3741-3782.
#'
#' @keywords internal
#' @inheritParams pcgen
#' @inheritParams pcalg::skeleton
#' @importFrom pcalg getNextSet
#'
skeleton2  <- function(suffStat, QTLs = integer(), K = NULL, replicates = TRUE,
                       use.manova = TRUE, alpha, labels, p, m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                       NAdelete = TRUE, covariates=NULL, max.iter = 50,
                       use.res = FALSE, res.cor = NULL, verbose = FALSE,
                       stop.if.significant = TRUE){

  cl <- match.call()

  # if (is.null(K))
  #   stopifnot(replicates)
  # if (!is.null(K))
  #   stopifnot(!replicates)

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

  seq_p <- seq_len(p)

  if (is.null(fixedGaps)) { # Graph G is complete
    G <- matrix(TRUE, p, p)
  } else {                  # Graph G has specific fixedGaps
    stopifnot(identical(dim(fixedGaps), c(p, p)))
    stopifnot(identical(fixedGaps, t(fixedGaps)))
    G <- !fixedGaps
  }
  diag(G) <- FALSE    # self-loops removed from graph G
  if (length(QTLs) > 0) {
    G[c(1,QTLs),c(1,QTLs)] <- FALSE # G<->QTL edges removed
  }

  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(FALSE, nrow = p, ncol = p)
  } else
    if (!identical(dim(fixedEdges), c(p, p))) {
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    } else if (!identical(fixedEdges, t(fixedEdges)))
      stop("fixedEdges must be symmetric")

  pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list", p))
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  ord <- 0L
  n.edgetests <- numeric(1)
  done <- FALSE

  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord + 1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    # ind <- which((G*upper.tri(G))>0, arr.ind = TRUE)
    ind <- ind[order(ind[, 1]), ]
    remEdges <- nrow(ind)
    if (verbose) cat("Order=", ord, "; remaining edges:", remEdges, "\n", sep = "")
    G.l <- split(G, gl(p, p))

    for (i in 1:remEdges) {
      if (verbose) cat("|i=", i, "|iMax=", remEdges, "\n")
      x <- ind[i, 1]; y <- ind[i, 2]
      if (G[y, x] && !fixedEdges[y, x]) {
        nbrsBool <- G.l[[x]]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)
          repeat {
            n.edgetests[ord1] <- n.edgetests[ord1] + 1
            pval <- pcgenTest(x, y, S = nbrs[S], suffStat, covariates = covariates,
                              QTLs = QTLs, K = K, replicates = replicates, use.manova = use.manova,
                              alpha = alpha, max.iter = max.iter,
                              stop.if.significant = stop.if.significant,
                              use.res = use.res, res.cor = res.cor)
            if (verbose) {
              cat("x=", labels[x], " y=", labels[y], " S=", labels[nbrs[S]],": pval =")
              if (!(is.na(pval)) & pval >= alpha) {
                cat('<<<<<<< ', pval, ' >>>>>>>>>', "\n") # Highlight removed edges
              } else
                cat(pval, "\n")
            } # End if (verbose)

            if (is.na(pval)) pval <- as.numeric(NAdelete)
            if (pMax[x, y] < pval) pMax[x, y] <- pval

            if (pval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            } else { # if (pval < alpha)
              nextSet <- getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
                break
              S <- nextSet$nextSet
            } # End pval<=> alpha?

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
