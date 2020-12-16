alpha <- 0.01

test_that("type B CI p-value of res.covar.test() WITH replicates", {
  #data(simdata) # G -> Y1 -> Y2 -> Y3; G -> Y3
                # 1     2     3     4
  G <- simdata$G
  Y1 <- simdata$Y1
  Y2 <- simdata$Y2
  Y3 <- simdata$Y3

  # Y1 _|_ Y2 | G ? (No, low p)
  expect_lt(res.covar.test(x = Y1, y = Y2, G = G), alpha)

  # Y1 _|_ Y3 | G  (No, low p)
  expect_lt(res.covar.test(x = Y1, y = Y3, G = G), alpha)

  # Y2 _|_ Y3 | G (No, low p)
  expect_lt(res.covar.test(x = Y2, y = Y3, G = G), alpha)


  # Y1 _|_ Y2 | (G, Y3) ? (No, low p)
  expect_lt(res.covar.test(x = Y1, y = Y2, X = Y3, G = G), alpha)

  # Y1 _|_ Y3 | (G, Y2) ? (Yes, high p)
  expect_gt(res.covar.test(x = Y1, y = Y3, X = Y2, G = G), alpha)

  # Y2 _|_ Y3 | (G, Y1)? (No, low p)
  expect_lt(res.covar.test(x = Y2, y = Y3, G = G, X = Y1), alpha)

})


test_that("type B CI p-value of res.covar.test() W/O replicates", {
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
  res.covar.test(x, y, G, K = K, use.manova = F)

})
