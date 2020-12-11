#'
#'
#' @param genotype a vector (character or factor) with genotype labels

make.Z.matrix <- function(genotype) {

    genotype <- factor(as.character(genotype), levels = as.character(unique(genotype)))

		xx <- data.frame(geno = genotype)
		mf <- model.frame(~ geno - 1, xx, drop.unused.levels = TRUE)
		mt <- terms(mf)
		f.terms <- attr(mt, "term.labels")[attr(mt,"dataClasses") == "factor"]
		Z.t <- Matrix(model.matrix(mt, data = mf,
                               contrasts.arg = lapply(mf[,f.terms, drop = FALSE],
                                                      contrasts, contrasts=FALSE)))
return(Z.t)
}

