#' Create the incidence matrix
#'
#' The incidence (design) matrix links observations (replicates) to unique genotypes. In
#' expriements without replicates, this reduces to the Identity matrix
#'
#' @param genotype a vector (character or factor) with genotype labels
#'
#' @importFrom Matrix Matrix

make.Z.matrix <- function(genotype) {

    genotype <- factor(as.character(genotype), levels = as.character(unique(genotype)))

		xx <- data.frame(geno = genotype)
		mf <- model.frame(~ geno - 1, xx, drop.unused.levels = TRUE)
		mt <- terms(mf)
		f.terms <- attr(mt, "term.labels")[attr(mt,"dataClasses") == "factor"]
		Z.t <- Matrix(model.matrix(mt, data = mf,
                               contrasts.arg = lapply(mf[,f.terms, drop = FALSE],
                                                      contrasts, contrasts=FALSE)))
		# experiment which is faster: `Matrix`, `new("dtCMatrix",...)`, `new("dtTMatrix",...)`,
		# `as(x,"triangularMatrix")`, or `as(x,"dtCMatrix")`
return(Z.t)
}

