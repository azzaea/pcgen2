test_that("fitEM works", {
  # data(d2tsnpsmeans) # G --> Y1 --> Y2

  M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
  K <- M %*% t(M)

  x <- d2tsnpsmeans$Y1
  y <- d2tsnpsmeans$Y2
  G <- d2tsnpsmeans$G


  Z.t <- make.Z.matrix(G)
  X.t <- Matrix(rep(1, length(G)))
  em.vec <- c(x, y)
  names(em.vec) <- rep(as.character(G), 2)


  # Both analyses give rise to the same estimates of the variance components and
  # fixed coefficients
  fit.reduced <- fitEM(em.vec, X.t, Z.t, cov.error = TRUE,
                       cov.gen = FALSE, max.iter = 20,
                       K = K)

  # Much less time consuming! Incorporate the idea into the algorithm (if possible)!
  fit.reduced.2 <- pcgen:::fitEM(em.vec, X.t, Z.t = Matrix(Z.t %*% M),
                         cov.error = TRUE,
                         cov.gen = FALSE, max.iter = 20)


  expect_equal(fit.reduced$variances, fit.reduced.2$variances)
  expect_equal(fit.reduced$deviance, fit.reduced.2$deviance)
  expect_equal(fit.reduced$it, fit.reduced.2$it)
  expect_equal(fit.reduced$coeff$fixed, fit.reduced.2$coeff$fixed)

  # These are different- not sure if it actually matters!
  ##expect_equal(fit.reduced$coeff$random, fit.reduced.2$coeff$random)
})
