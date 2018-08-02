#' Check for consistency in genetic effects
#'
#' Given output from pcgen or pcgenFast, this function checks whether the
#' estimated graph is consistent with the set of traits having significant genetic variance.
#' The function detects traits that have significant genetic variance but for which there is no
#' partially directed path from G.
#'
#' Details
#'
#' @inheritParams pcgen
#'
#' @param pcgen.output a graph with nodes genotype and a number of
#'                     traits. Typically output from pcgen or pcgenFast
#'
#' @return A logical matrix of dimension (p+1) x (p+1), p being the number of traits.
#'         Most entries are FALSE, except those in the first row and column for which there
#'         are conflicts.
#'         Entries [1,j] and [j,1] are TRUE if the jth trait has
#'         significant genetic variance, but there is no
#'         partially directed path from G towards that trait. The matrix
#'         can then be used in a subsequent run of pcgen or pcgenFast, in the
#'         fixedEdges argument. The arguments suffStat, alpha and covariates
#'         should stay the same throughout
#'         (first run of pcgen, checkG, second run of pcgen).
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references A paper on arxiv
#'
#' @export

checkG <- function(pcgen.output, suffStat, alpha = 0.01,
                   covariates = NULL) {
# pcgen.output = pcgen.fit.d[[1]]; suffStat = dr[, 1:35]; covariates = dr[, 36:37]; alpha = 0.01; fixedEdges = NULL; mean.adj = 'none'; Vg = NULL; Ve = NULL

  # pcgen.output : output from pcgen or pcgenFast (only the graph!)
  # fixedEdges : matrix of TRUE/FALSE, of dimension p+1 x p+1
  #               first row/column MUST be FALSE (to do: check on this)
  #  The output will be the same matrix, but with the first row and column
  # possible modified (in case there are traits for which (1) there is marginally
  # significant genetic variance, and (2) there is no partially directed path from
  # G to that trait)

  # other to do's :
  #
  # in pcgen and pcgenFast, also return all input arguments ? ... and re-use here ?
  #
  # require that the first column is the genotypic column, although do not require the
  # name 'genotype'. Do the same check on suffStat in all functions, and document it
  #
  # checkG does not really allow for QTLs. MENTION this
  #
  # Find a suitable example for illustration. Do not forget to return Vg, Ve in the
  # initial call to pcgen, and use these Vg, Ve in the call to checkG
  #
  # document! especially, also clarify the role of fixedEdges
  #
  # test!

  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    stopifnot(nrow(covariates)==nrow(suffStat))
  }

  if (colnames(suffStat)[1]!='genotype') {stop('The first column of suffStat should be genotype')}

  labels <- colnames(suffStat)
  p      <- ncol(suffStat)

  genVar <- rep(TRUE, p)

  for (trait.number in 2:p) {

    pval <- pcgenTest(x=1, y=trait.number, S=integer(0), suffStat = suffStat,
                      covariates=covariates, alpha = alpha)
    #cat(pval,'\n')
    if (pval > alpha) {
      genVar[trait.number] <- FALSE
    }
  }

  genVarGraph <- rep(TRUE, p - 1)
  directGeffect <- rep(TRUE, p - 1)

  for (j in 2:p) {

    genVarGraph[j-1] <- !msep(alpha='genotype', beta=names(suffStat)[j],
                              a = pcgen.output@graph)

    directGeffect[j-1] <- isAdjacent(object=pcgen.output@graph, from='genotype', to=names(suffStat)[j])

  }

  # if (is.null(fixedEdges)) {
  fixedEdges.temp <- matrix(FALSE, p, p)
  colnames(fixedEdges.temp) <- rownames(fixedEdges.temp) <- names(suffStat)# [-1]
  # } else {
  #   fixedEdges.temp <- fixedEdges
  # }

  if (any(genVarGraph==FALSE & genVar[-1]==TRUE)) {

    # First, detect for which traits (column numbers in suffStat) there is no directed (actually: ancestral ?) path
    # from G to that trait
    nrs <- which(genVarGraph==FALSE & genVar[-1]==TRUE) + 1
    # Then, check for which traits there are direct genetic effects already; these should not disappear
    nrs2<- which(directGeffect==TRUE) + 1
    nrs <- sort(unique(c(nrs,nrs2)))

    fixedEdges.temp[1,nrs] <- fixedEdges.temp[nrs,1] <- TRUE

  }

  return(fixedEdges.temp)
}
