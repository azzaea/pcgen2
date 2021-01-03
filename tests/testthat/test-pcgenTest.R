alpha <- 0.01

## Type A: Trait _|_ G | Traits -----------------------------------------------

test_that("Type A CI tests WITH replicates of pcgenTest() work", {
  #data(simdata) # G -> Y1 -> Y2 -> Y3; G -> Y3
                 # 1     2     3     4

  # Y1 _|_ G | Y2 (No, low p):
  expect_lt(pcgenTest(x = 2, y = 1, S = 3, simdata), alpha)
  # Y1 _|_ G | Y3 (No,   low p):
  expect_lt(pcgenTest(x = 2, y = 1, S = 4, simdata), alpha)
  # Y1 _|_ G | Y2, Y3 (No, low p):
  expect_lt(pcgenTest(x = 2, y = 1, S = c(3, 4), simdata), alpha)

  # Y2 _|_ G | Y1 (Yes, high p)
  expect_gt(pcgenTest(x = 3, y = 1, S = 2, simdata), alpha)
  # Y2 _|_ G | Y3 (No, low p)
  expect_lt(pcgenTest(x = 3, y = 1, S = 4, simdata), alpha)
  # Y2 _|_ G | Y1, Y3 (No, low p)
  expect_lt(pcgenTest(x = 3, y = 1, S = c(2, 4), simdata), alpha)

  # Y3 _|_ G | Y1 (No, low p)
  expect_lt(pcgenTest(x = 4, y = 1, S = 2, simdata), alpha)
  # Y3 _|_ G | Y2 (No, low p)
  expect_lt(pcgenTest(x = 4, y = 1, S = 3, simdata), alpha)
  # Y3 _|_ G | Y1, Y2 (No, low p)
  expect_lt(pcgenTest(x = 4, y = 1, S = c(2, 3), simdata), alpha)
})

test_that("Type A CI tests W/O replicates of pcgenTest() work", {
  #data("d2tsnpsmeans") # G -> Y1 -> Y2
                       # 1     2     3

  suffStat <- d2tsnpsmeans[,1:3]
  M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
  K <- M %*% t(M)

  # Y1 _|_ G | Y2 (No, low p):
  expect_lt(pcgenTest(x = 2, y = 1, S = 3, suffStat, K = K), alpha)

  # Y2 _|_ G | Y1 (Yes, high p):
  expect_gt(pcgenTest(x = 3, y = 1, S = 2, suffStat, K = K), alpha)
})

## Type B: Trait _|_ Trait | {Traits, G, Q} -----------------------------------

test_that("Type B CI tests WITH replicates of pcgenTest() work", {
  #data(simdata) # G -> Y1 -> Y2 -> Y3; G -> Y3
                 # 1     2     3     4

  # Remmeber that G is added anyways. So, both these are equivalent
  # S = {traits, G} <==> S = {traits}
  # I only test this equivalence here, since p-value calculation
  # is the job of res.covar.test
  # Y1 _|_ Y2 | (G, Y3) ? (No, low p)
  expect_equal(pcgenTest(x = 2, y = 3, S = 4, simdata),
              pcgenTest(x = 2, y = 3, S = c(1, 4), simdata))
  # # Y1 _|_ Y3 | (G, Y2) ? (Yes, high p)
  # expect_equal(pcgenTest(x = 2, y = 4, S = 3, simdata, max.iter = 5),
  #              pcgen:::pcgenTest(x = 2, y = 4, S = c(1, 3), simdata, max.iter = 5))
  # # Y2 _|_ Y3 | (G, Y1)? (No, low p)
  # expect_equal(pcgenTest(x = 3, y = 4, S = 2, simdata),
  #              pcgen:::pcgenTest(x = 3, y = 4, S = c(1, 2), simdata))

})

test_that("Type B CI tests W/O replicates of pcgenTest() work", {
  #data("d2tsnpsmeans") # G -> Y1 -> Y2
                       # 1     2     3
  suffStat <- d2tsnpsmeans[,1:3]
  M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
  K <- M %*% t(M)

  # Y1 \perp Y2 | G (No, low p)
  expect_lt(pcgenTest(x = 2, y = 3, S = 1, suffStat, K = K, use.manova = F), alpha)
})
