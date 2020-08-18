#' The pcgen algorithm
#'
#' The pcgen algorithm, assuming independent genetic effects.
#'
#' (1) This is an adaptation of the pc function from pcalg, with a
#' different conditional independence test (see pcgenTest), as well as
#' modified orientation rules  (2) Some parameters from the original
#' pc function are fixed here: skel.method = "stable", u2pd = "relaxed",
#' conservative = FALSE, maj.rule = TRUE, solve.confl = TRUE and
#' NAdelete = TRUE. labels (defining the names of the nodes of the
#' graph) is derived from the data-frame suffStat, containing the data
#' (3) pcgen requires phenotypic observations on multiple traits,
#' measured on the same samples. The current implementation extends the genetic
#' relatedness matrix definition from K = Z Z^t, Z being the incidence matrix
#' assigning plants to genotypes (as is common in plants data) to a 
#' a marker-based relatedness matrix definition (appropriate for human data).
#' (4) the test for conditional independence between Yj and Yk given
#' a set of traits YS can be based either on a bivariate mixed model
#' (put res.cor = NULL) or on partial correlations among the residuals.
#' In the latter case, res.cor should be the correlation matrix of the
#' residuals, obtained with getResiduals.
#'
#' @param m.max Maximum size of the conditioning sets.
#'
#' @param max.iter Maximum number of iterations in the EM-algorithm.
#'
#' @param fixedGaps  A logical matrix of dimension (p+1) * (p+1), where
#'                   p is the number of traits. The first row and column
#'                   refer to the node genotype, and subsequent rows and
#'                   columns to the traits. As in the pcalg package,
#'                   the edge i-j is removed before starting the algorithm
#'                   if the entry [i,j] or [j,i] (or both) are TRUE.
#'                   In that case, the edge is guaranteed to be absent
#'                   in the resulting graph.
#'
#' @param fixedEdges A logical matrix of dimension (p+1) * (p+1), where
#'                   p is the number of traits. The first row and column
#'                   refer to the node genotype, and subsequent rows and
#'                   columns to the traits. As in the pcalg package,
#'                   the edge i-j is never considered for removal
#'                   if the entry [i,j] or [j,i] (or both) are TRUE.
#'                   In that case, the edge is guaranteed to be present
#'                   in the resulting graph.
#'
#' @param verbose If TRUE, p-values for the conditional independence
#'                tests are printed
#'
#' @param suffStat data.frame, of which the first column is the factor genotype,
#'                  and subsequent columns contain the traits. The name of the
#'                 first column should be genotype. Should not contain covariates.
#'
#' @param alpha (Default 0.01) The significance level used in each conditional independence test.
#'
#' @param QTLs column numbers in suffStat that correspond to QTLs
#'             These may be partly in S and x and y, but x and y cannot be both QTLs.
#'
#' @param covariates data.frame containing covariates, that should always be used
#'                   in each conditional independence test. Should be either NULL (default)
#'                   or a data.frame with the same number of rows as suffStat
#'
#' @param stop.if.significant If TRUE, the EM-algorithm used in some
#'              of the conditional independence tests will be stopped
#'              whenever the p-value becomes significant, i.e. below
#'              alpha. This will speed up calculations, and can be done
#'              because (1) the PC algorithm only needs an accept/reject
#'              decision (2) In EM the likelihood is nondecreasing. Should
#'              be put to FALSE if the precise p-values are of interest.
#'
#' @param res.cor Correlation matrix of the residuals of the GBLUP. In case
#'                it is NULL (default), the test for conditional independence
#'                between Yj and Yk given a set of traits YS is based on a
#'                bivariate mixed model. If res.cor is provided, it is based
#'                on partial correlations among the residuals, tested using
#'                the gaussCItest function from pcalg. The residuals can be
#'                obtained using the function getResiduals
#'                
#' @param return.pvalues                
#'
#' @return a graph (an object with S3 class \code{"pcgen"})
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references Kruijer et al. (2018, in preparation) Reconstruction of networks with direct and indirect genetic effects
#'
#' @export

pcgen <-
function (suffStat, covariates = NULL, QTLs = integer(), alpha = 0.01, 
          m.max = Inf,  fixedEdges = NULL, fixedGaps = NULL, verbose = FALSE, 
          use.res = FALSE, res.cor = NULL, max.iter = 50, stop.if.significant = TRUE,
          return.pvalues = FALSE)
{
    NAdelete <- TRUE
    u2pd <- c("relaxed", "rand", "retry")[1]
    skel.method <- c("stable", "original", "stable.fast")[1]
    conservative <- FALSE
    maj.rule <- TRUE
    solve.confl <- TRUE

    cl <- match.call()
    if (is.null(alpha)) alpha = 0.01
    if (!is.null(covariates)) {
      covariates <- as.data.frame(covariates)
      stopifnot(nrow(covariates)==nrow(suffStat))
    }

    if (colnames(suffStat)[1]!='G') {stop('The first column of suffStat should be G (genotype)')}
    if (class(QTLs)!='integer') {stop('QTLs should be a vector of integers')}
    if (1 %in% QTLs) {stop('QTLs should not contain the genotype column (G)')}

    if (length(QTLs) > 0) {
      non.collider.nodes <- c(1,sort(QTLs))
    } else {
      non.collider.nodes <- 1
    }

    ###

    labels <- colnames(suffStat)
    p      <- ncol(suffStat)
#
#       dec  <- NULL
#       Vg   <- NULL # just in case the user accidentally specifies Vg, Ve
#       Ve   <- NULL

    #u2pd <- match.arg(u2pd)
    #skel.method <- match.arg(skel.method)

    if (u2pd != "relaxed") {
        if (conservative || maj.rule)
            stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
        if (solve.confl)
            stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
    }

    if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")

    ##################################
    # The following chunk is (i) not in the original pcalg function
    #                        (ii) ALSO in the skeleton2 function
    #                             (which now can also deal with QTLs by itself)

    #if (is.null(fixedEdges)) {
    #  stopifnot(identical(dim(fixedEdges), c(p, p)))
    #}

    #if (is.null(fixedGaps)) {
    #  gapMatrix <- matrix(FALSE, p, p)
    #} else {
    #  #stopifnot(identical(dim(fixedGaps), c(p, p)))
    #  stopifnot(identical(dim(fixedGaps), c(p, p)))
    #  gapMatrix <- as.matrix(fixedGaps)
    #}

    #if (length(QTLs) > 0) {
    #  gapMatrix[c(1,QTLs),c(1,QTLs)] <- TRUE
    #}

    #fixedGaps <- gapMatrix

    ##########################

    #skel.out <-
    skel <- skeleton2(suffStat = suffStat, alpha = alpha, labels = labels,
          method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
          NAdelete = NAdelete, m.max = m.max, verbose = verbose,
          covariates=covariates,
          QTLs = QTLs,
          max.iter = max.iter,
          stop.if.significant = stop.if.significant,
          use.res = use.res,
          res.cor = res.cor)
    ##16-1-18## : added Vg = Vg, Ve = Ve, dec = dec


    #genVar <- skel[['genVar']]
    #qtlVar <- skel[['qtlVar']]
    #skel   <- skel[['skel.out']]

    skel@call <- cl

    if (!conservative && !maj.rule) {
# this option can not occur, at least for the moment
        gr <- switch(u2pd, rand = udag2pdag(skel),
                     retry = udag2pdagSpecial(skel)$pcObj,
                     relaxed = udag2pdagRelaxed2(skel, verbose = verbose,
                                                 solve.confl = solve.confl,
                                                 non.collider.nodes = non.collider.nodes))

    } else {

        pc. <- pc.cons.intern2(skel, suffStat = suffStat, alpha = alpha,
                               version.unf = c(2, 1), maj.rule = maj.rule,
                               verbose = verbose,
                               covariates = covariates,
                               QTLs = QTLs, max.iter = max.iter,
                               stop.if.significant = stop.if.significant,
                               use.res = use.res, res.cor = res.cor)

        gr <- udag2pdagRelaxed2(pc.$sk, verbose = verbose,
                                unfVect = pc.$unfTripl,
                                solve.confl = solve.confl,
                                non.collider.nodes = non.collider.nodes)

    }
    
    if (return.pvalues == TRUE) {
      return(list(gr = gr, pMax = skel@pMax))
    } else {
      return(gr)
    }
}

