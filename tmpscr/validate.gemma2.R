pacman::p_load(gemma2, tictoc, tidyverse)


data(d2tsnpsmeans) # G --> Y1 --> Y2

M <- as.matrix(d2tsnpsmeans[,-(1:3)] / sqrt(ncol(d2tsnpsmeans[,-(1:3)])))
K <- M %*% t(M)

x <- d2tsnpsmeans$Y1
y <- d2tsnpsmeans$Y2
G <- d2tsnpsmeans$G


Z.t <- make.Z.matrix(G)
eig <- svd(K)
U <- eig$u
d <- eig$d
K.trans <- U %*% diag(sqrt(d))
Z <- as.matrix(Matrix(Z.t %*% K.trans)); str(Z)

X.t <- Matrix(rep(1, length(G)))
X <- as.matrix(X.t); str(X)

em.vec <- c(x, y)
names(em.vec) <- rep(as.character(G), 2)

max.iter <- seq(from = 100, to = 900, by = 150)
seed <- 1223


pvaluefit <- itrsfit <- REMLLRTfit <-  vector(length = length(max.iter))
pvaluegemma <- itrsgemma <- REMLLRTgemma <-  vector(length = length(max.iter))


tic.clearlog()
for (it in 1:length(max.iter)){

  set.seed(seed)
  ## pcgen code:
  tic()

  fit.reduced <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = FALSE, max.iter = max.iter[it], # maybe false
                       cov.gen = FALSE, Vg.start = c(1, 1, 0.1), Ve.start = c(1, 1, 0), alpha = 0.001, verbose = T)
  fit.full    <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = TRUE,
                       stop.if.significant = TRUE, max.iter = max.iter[it],
                       null.dev = fit.reduced[[length(fit.reduced)]]$deviance,
                       Vg.start = fit.reduced[[length(fit.reduced)]]$variances$Vg,
                       Ve.start = fit.reduced[[length(fit.reduced)]]$variances$Ve, alpha = 0.001, verbose = T)


  loglik_Full    <- -0.5*fit.full[[length(fit.full)]]$deviance
  loglik_Reduced <- -0.5*fit.reduced[[length(fit.reduced)]]$deviance
  REMLLRTfit[it] <- 2 * max(loglik_Full - loglik_Reduced, 0)
  pvaluefit[it]  <- (1 - pchisq(REMLLRT, df = 1))
  itrsfit[it] <- mean(c(fit.reduced[[length(fit.reduced)]]$it, fit.full[[length(fit.full)]]$it))

  toc(log = TRUE, quiet = TRUE)
}

fit.log <- tic.log(format = FALSE)

fit.timings <- unlist(lapply(fit.log, function(x) x$toc - x$tic))



e_out <- eigen2(K)
tic.clearlog()
for (it in 1:length(max.iter)){

  set.seed(seed)
  tic()
  ## gemma2 code:

  gemma.reduced <- mygemma2(eval = e_out$values, # eigen2(kinship)
                            #X = as.matrix(g1) %*% e_out$vectors ,
                            X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
                            Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
                            V_g = matrix(c(1, 0.1, 0.1, 1), nrow = 2),
                            V_e = matrix(c(1, 0, 0, 1), nrow = 2), max_prec = 0.001,
                            max_iter = max.iter[it], verbose_output = T, reduced = T
  )
  gemma.full    <- mygemma2(eval = e_out$values, # eigen2(kinship)
                            #X = as.matrix(g1) %*% e_out$vectors ,
                            X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
                            Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
                            V_g = gemma.reduced[[1]]$Vg,
                            V_e = gemma.reduced[[1]]$Ve, max_prec = 0.001,
                            max_iter = max.iter[it], verbose_output = T , reduced = F)
  #null.dev = fit.reduced$deviance,
  loglik_Full    <- gemma.full[[length(gemma.full)]]$logl_new
  loglik_Reduced <- gemma.reduced[[length(gemma.reduced)]]$logl_new
  REMLLRTgemma[it] <- 2 * max(loglik_Full - loglik_Reduced, 0)
  pvaluegemma[it]  <- (1 - pchisq(REMLLRT, df = 1))
  itrsgemma[it] <-  mean(c(gemma.reduced[[length(gemma.reduced)]]$it, gemma.full[[length(gemma.full)]]$it))

  toc(log = TRUE, quiet = TRUE)
}

gemma.log <- tic.log(format = FALSE)

gemma.timings <- unlist(lapply(gemma.log, function(x) x$toc - x$tic))


data.frame(it = max.iter, gemma = REMLLRTgemma, fitEM = REMLLRTfit) %>%
  pivot_longer(!it) %>%
  ggplot(aes(x = it)) + geom_line(aes(y = value, color = name)) +
  # ggplot(data =  dataframe(it = max.iter, convergance = ))
  ggtitle("REMLLRT")

data.frame(it = max.iter, gemma = gemma.timings, fitEM = fit.timings) %>%
  pivot_longer(!it) %>%
  ggplot(aes(x = it)) + geom_line(aes(y = value, color = name)) +
  # ggplot(data =  dataframe(it = max.iter, gema.its = itrsgemma ))
  ggtitle("Time")








reps <- rep(1, 7)
for (it in 1:length(reps)){

  fit.reduced <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = FALSE, max.iter = max.iter[3], # maybe false
                       cov.gen = FALSE, Vg.start = c(1, 1, 0.1), Ve.start = c(1, 1, 0), alpha = 0.001, verbose = T)
  # fit.full   <- fitEM(em.vec, X.t, Z.t, K = K, cov.error = TRUE,
  #                      stop.if.significant = TRUE, max.iter = max.iter[3],
  #                      null.dev = fit.reduced[[length(fit.reduced)]]$deviance,
  #                      Vg.start = fit.reduced[[length(fit.reduced)]]$variances$Vg,
  #                      Ve.start = fit.reduced[[length(fit.reduced)]]$variances$Ve, alpha = 0.001, verbose = T)
  deviance <- sapply(FUN = function(x)x[[3]], X = fit.reduced)
  if (it == 1)
    plot(deviance)
  lines(deviance)
}

title(main = "fitEM- reduced model deviance")

#
# ggplot(data.frame(deviance = deviance, it = 1:length(fit.full))) +
#   geom_line(aes(x = it, y=deviance))  + ggtitle("fitEM, full model")


for (it in 1:length(reps)){
  gemma.reduced <- mygemma2(eval = e_out$values, # eigen2(kinship)
                            #X = as.matrix(g1) %*% e_out$vectors ,
                            X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
                            Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
                            V_g = matrix(c(1, 0.1, 0.1, 1), nrow = 2),
                            V_e = matrix(c(1, 0, 0, 1), nrow = 2), max_prec = 0.001,
                            max_iter = max.iter[3], verbose_output = T, reduced = T
  )
  gemma.full    <- mygemma2(eval = e_out$values, # eigen2(kinship)
                            #X = as.matrix(g1) %*% e_out$vectors ,
                            X = matrix(1, nrow = 1, ncol = nrow(d2tsnpsmeans)),
                            Y = t(as.matrix(d2tsnpsmeans[, c(2, 3)])) %*% e_out$vectors,
                            V_g = gemma.reduced[[1]]$Vg,
                            V_e = gemma.reduced[[1]]$Ve, max_prec = 0.001,
                            max_iter = max.iter[3], verbose_output = T , reduced = F)
  loglik <- sapply(FUN = function(x)x[[1]], X = gemma.reduced)

  if (it == 1)
    plot(loglik)
  lines(loglik)
}
title(main = "GEMMA- reduced model loglik")


















