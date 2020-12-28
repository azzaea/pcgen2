alpha <- 0.01

test_that("type B CI p-value of res.covar.test() WITH replicates", {
  #data(simdata) # G -> Y1 -> Y2 -> Y3; G -> Y3

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
  ##expect_gt(res.covar.test(x = Y1, y = Y3, X = Y2, G = G), alpha)

  # Y2 _|_ Y3 | (G, Y1)? (No, low p)
  ##expect_lt(res.covar.test(x = Y2, y = Y3, G = G, X = Y1), alpha)

})


test_that("type B CI p-value of res.covar.test() W/O replicates", {
  #data(d2tsnpsmeans) # G -> Y1 -> Y2

  Y1 <- d2tsnpsmeans$Y1
  Y2 <- d2tsnpsmeans$Y2
  G <- d2tsnpsmeans$G

  M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
  K <- M %*% t(M)

  # Y1 _|_ Y2 | G ?  (No, low p)
  expect_lt(res.covar.test(x = Y1, y = Y2, G = G , K = K), alpha)

})
