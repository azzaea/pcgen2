test_that("fitEM works", {
#   # data(d2tsnpsmeans) # G --> Y1 --> Y2
#
#   M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
#   K <- M %*% t(M)
#
#   x <- d2tsnpsmeans$Y1
#   y <- d2tsnpsmeans$Y2
#   G <- d2tsnpsmeans$G
#
#
#   Z.t <- make.Z.matrix(G)
#   eig <- svd(K)
#   U <- eig$u
#   d <- eig$d
#   K.trans <- U %*% diag(sqrt(d))
#   Z <- as.matrix(Matrix(Z.t %*% K.trans)); str(Z)
#
#   X.t <- Matrix(rep(1, length(G)))
#   X <- as.matrix(X.t); str(X)
#
#   em.vec <- c(x, y)
#   names(em.vec) <- rep(as.character(G), 2)
#
#   ## pcgen code:
#   tictoc::tic()
#   fit.reduced <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = FALSE, max.iter = 200, # maybe false
#                        cov.gen = FALSE, Vg.start = c(1, 1, 0.1), Ve.start = c(1, 1, 0))
#   tictoc::toc()
#   fit.full    <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = TRUE,
#                        stop.if.significant = TRUE, max.iter = 200,
#                        null.dev = fit.reduced$deviance,
#                        Vg.start = fit.reduced$variances$Vg,
#                        Ve.start = fit.reduced$variances$Ve)
#   loglik_Full    <- -0.5*fit.full$deviance
#   loglik_Reduced <- -0.5*fit.reduced$deviance
#   REMLLRT <- 2 * max(loglik_Full - loglik_Reduced, 0)
#   pvalue  <- (1 - pchisq(REMLLRT, df = 1))
#   tictoc::toc()
#
#   ## gemma2 code:
#   pacman::p_load(gemma2)
#   e_out <- eigen2(K)
#   tictoc::tic()
#   gemma.reduced <- MphEM(eval = e_out$values, # eigen2(kinship)
#                          #X = as.matrix(g1) %*% e_out$vectors ,
#                          X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
#                          Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
#                          V_g = matrix(c(1, 0.1, 0.1, 1), nrow = 2),
#                          V_e = matrix(c(1, 0, 0, 1), nrow = 2), max_prec = 10^-4,
#                          max_iter = 300, verbose_output = T
#                          )
#   gemma.full    <- MphEM(eval = e_out$values, # eigen2(kinship)
#                          #X = as.matrix(g1) %*% e_out$vectors ,
#                          X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
#                          Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
#                          V_g = gemma.reduced[[1]]$Vg,
#                          V_e = gemma.reduced[[1]]$Ve, max_prec = 10^-4,
#                          max_iter = 300, verbose_output = F)
#                        #null.dev = fit.reduced$deviance,
#   loglik_Full    <- gemma.full[[1]]$logl_new
#   loglik_Reduced <- gemma.reduced[[1]]$logl_new
#   REMLLRT <- 2 * max(loglik_Full - loglik_Reduced, 0)
#   pvalue  <- (1 - pchisq(REMLLRT, df = 1))
#   tictoc::toc()
#
#
#
#   loglik <- sapply(FUN = function(x)x[[1]], X = gemma.reduced)
#   ggplot(data.frame(loglik = loglik, it = 1:length(gemma.reduced))) +
#     geom_point(aes(x = it, y=loglik)) + coord_cartesian(ylim = c(-1027, -1026), xlim = c(180, 200))
#   str(foo)
#   fit.reducedem$variances
#
#   ## -------- hglm code:
#
#   fit.reducedh <- hglm::hglm(y = y, X = Matrix(x),
#                              Z = (Z),
#                              calc.like = T)
#
#
#   y1 <- mdhglm::DHGLMMODELING(Model = "mean", Link = "identity",
#                               LinPred = y1 ~ y2 + (1|G), RandDist = "gaussian" ,
#                               LMatrix = Z )
#   y2 <- mdhglm::DHGLMMODELING(Model = "mean", Link = "identity",
#                               LinPred = y2 ~ (1|G), RandDist = "gaussian" ,
#                               LMatrix = Z )
#   data <- data.frame(y1 = x, y2 = y, G = as.factor(G), i = rep(1, length(G)))
#
#   h
#
#   ##
#   y
#
#   ###
#
#   Init_Corr=list(c(0))
#   ZZ1<-model.matrix(~as.factor(data$G)-1)
#   ZZCorr=list(ZZ1,ZZ1)
#   mdhglm::jointfit(RespDist = c( "gaussian", "gaussian"),
#                    DataMain = list(data, data), MeanModel = list(y1, y2),
#                    structure = "correlated",
#                    #Init_Corr = Init_Corr,
#                    convergence = 1,
#                    ZZCorr = ZZCorr
#                     )
#     # res_corr<-jointfit(RespDist=c("gaussian","binomial"),DataMain=list(eg1,eg1),
#   #                    MeanModel=list(jm1,jm2),structure="correlated",
#   #                    Init_Corr=Init_Corr,convergence=1,ZZCorr=ZZCorr)
#
#   fit.reduced.mcglm <- mcglm::mcglm()
#
#
#   # fit.full    <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = TRUE, alpha = alpha,
#   #                      stop.if.significant = TRUE, max.iter = max.iter,
#   #                      null.dev = fit.reduced$deviance,
#   #                      Vg.start = fit.reduced$variances$Vg,
#   #                      Ve.start = fit.reduced$variances$Ve)
#   fit.full    <- hglm::hglm(em.vec, X.t, Z.t, K = K, cov.error = TRUE, alpha = alpha,
#                             stop.if.significant = TRUE, max.iter = max.iter,
#                             null.dev = fit.reduced$deviance,
#                             Vg.start = fit.reduced$variances$Vg,
#                             Ve.start = fit.reduced$variances$Ve)
#
#
#
#
#
#
# ### -- Jose's test
#   # Both analyses give rise to the same estimates of the variance components and
#   # fixed coefficients
#   fit.reduced <- fitEM(em.vec, X.t, Z.t, cov.error = TRUE,
#                        cov.gen = FALSE, max.iter = 20,
#                        K = K)
#
#   # Much less time consuming! Incorporate the idea into the algorithm (if possible)!
#   fit.reduced.2 <- fitEM(em.vec, X.t, Z.t = Matrix(Z.t %*% M),
#                          cov.error = TRUE,
#                          cov.gen = FALSE, max.iter = 20)
#
#
#   expect_equal(fit.reduced$variances, fit.reduced.2$variances)
#   expect_equal(fit.reduced$deviance, fit.reduced.2$deviance)
#   expect_equal(fit.reduced$it, fit.reduced.2$it)
#   expect_equal(fit.reduced$coeff$fixed, fit.reduced.2$coeff$fixed)
#
#   # These are different- not sure if it actually matters!
#   ##expect_equal(fit.reduced$coeff$random, fit.reduced.2$coeff$random)
})
