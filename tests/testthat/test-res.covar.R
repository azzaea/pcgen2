test_that("res.covar.test() works", {
  dm <- d2tsnpsmeans
  x <- dm$Y1
  y <- dm$Y2
  G <- dm$G

  M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
  K <- M %*% t(M)

#
#   # # Generate the raw  0,2 dosage values (for Gaston)
#   # genos <- as.matrix(purrr::map_df(purrr::map_df(dm[, -c(1:3)], factor), as.integer))
#   # genos[genos == 1] <- 0
#   # rownames(genos) <- rownames(dm)
#   # genos.gaston <- as(genos,"bed.matrix")
#   # K2 <- gaston::GRM(genos.gaston, autosome.only = F) # calculated as: XX'/q, with X the standardized
#   # # genotype matrix and q the number of SNPs (ncol(x))
#   # all.equal(K, K2)  # mean difference of ~1
#
#
#   res.covar.test(x, y, G, K = K, use.manova = F)
})
