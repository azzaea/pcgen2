test_that("multiplication works", {
  setwd("~/github_repos/network-analysis/src/haplotypes/willem")
  load(file = 'simdata1.RData')
  str(d)
  dm[1:6, 1:6]
  dm[,4] %>% unique() # dm has 2 unique values for each snp
  # because data is simulated inbred plant

  M <- as.matrix(dm[,-(1:3)] / sqrt(ncol(dm)))
  M[1:6, 1:6]
  K <- M %*% t(M)


  x <- dm$Y1
  y <- dm$Y2
  G <- dm$G

  Z.t <- pcgen:::make.Z.matrix(G)
  Z.t[1:6,1:6]
  X.t <- Matrix(rep(1, length(G)))

  em.vec <- c(x, y)
  names(em.vec) <- rep(as.character(G), 2)

  fit.reduced <- pcgen:::fitEM(em.vec, X.t, Z.t, cov.error = TRUE,
                               cov.gen = FALSE, max.iter = 20)
  fit.reduced2 <- pcgen2:::fitEM(em.vec, X.t, Z.t, cov.error = TRUE,
                                 cov.gen = FALSE, max.iter = 20)
  all.equal(fit.reduced, fit.reduced2)

  # New code (see modification below and attached code)
  source("github_repos/network-analysis/src/haplotypes/willem/EM_Function_pcgen.R")

  # Both analyses give rise to the same estimates of the variance components and coefficients
  fit.reduced.K <- pcgen2:::fitEM(em.vec, X.t, Z.t, cov.error = TRUE,
                                  cov.gen = FALSE, max.iter = 20,
                                  K = K)



  # Changing Z defition, Much less time consuming!
  # Incorporate the idea into the algorithm (if possible)!
  fit.reduced.ZM <- pcgen2:::fitEM(em.vec, X.t, Z.t = Matrix(Z.t %*% M),
                                   cov.error = TRUE,
                                   cov.gen = FALSE, max.iter = 20)

  all.equal(fit.reduced.K$coeff$fixed,fit.reduced.ZM$coeff$fixed)

  qqplot(fit.reduced$coeff$random, fit.reduced.ZM$coeff$random)

  str(fit.reduced.K)
  str(fit.reduced.ZM)


  # COTE: For the moment, we only need to make some minor modifications in this bit
  # Need to explore in more detail the proposal by Marco.
  if (!is.null(K)) {
    #eig <- eigen(K)
    #U <- eig$vectors
    #d <- eig$values
    ##K.trans <- U %*% diag(1/sqrt(d))
    #K.trans <- U %*% diag(sqrt(d))
    #Z.t <- Matrix(Z.t %*% K.trans)

    eig <- svd(K)
    U <- eig$u
    d <- eig$d
    #K.trans <- U %*% diag(1/sqrt(d))
    K.trans <- U %*% diag(sqrt(d))
    Z.t <- Matrix(Z.t %*% K.trans)
  }


  # Illustration, to show how fitEM can work with only genotypic means:
  x <- dm$Y1
  y <- dm$Y2
  G <- dm$G
  Z.t <- pcgen:::make.Z.matrix(G)
  X.t <- Matrix(rep(1, length(G)))
  em.vec <- c(x, y)
  names(em.vec) <- rep(as.character(G), 2)

  fit.reduced <- pcgen:::fitEM(em.vec, X.t, Z.t = Matrix(Z.t %*% M),
                               cov.error = TRUE,
                               cov.gen = FALSE, max.iter = 20)


})
