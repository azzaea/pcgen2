gen.var.test <- function(y, X = data.frame(), G, K = NULL, return.fitted = FALSE) {
  # y=d[,2]; X=d[,3]; G=d[,1] # previous set up
  # y = dm[, 2]; X = dm[,3]; G=dm[,1]; K = K # new contribution

  if (is.null(K)) { # Old code: suitable for plant data with replicates (K = ZZ^t)
    X <- as.data.frame(X)
    if (ncol(X) > 0) {
      a <- data.frame(y = y, X, G = G)
      names(a) <- c('y',paste0('C',1:ncol(X)),'G')
    } else {
      a <- data.frame(y = y, G = G)
    }
    a <- a[!is.na(a$y),]
    G <- factor(as.character(G))
    lm.obj <- lm(as.formula(paste("y~", paste(names(a)[-1],collapse = "+"))), data=a)
    av <- anova(lm.obj)
    pvalue <- av[[5]][1+ncol(X)]
    if (return.fitted == TRUE)
      fitted.values <- as.numeric(as.matrix(X) %*% matrix(as.numeric(coefficients(lm.obj))[2:(ncol(X)+1)]))
  } else {  # Generic K = A
    X <- as.matrix(X)
    if (ncol(X) == 0) { # Y_S = \phi
      X <-  matrix(1, nrow = length(y))
    } else {  # Y_S \ne \phi
      X <- cbind(1, X)
    }
    lm.obj <- gaston::lmm.aireml(Y = y , X = X, K = K, verbose = F)
    lrt <- max(2*(lm.obj$logL - lm.obj$logL0), 0)
    pvalue <- pchisq(lrt, df = 1, lower.tail = F)
    if (return.fitted == TRUE)
      fitted.values <- NULL #Not sure how to calculate it here
  }

  if (return.fitted == TRUE)
    return(list(pvalue, fitted.values))
  else
    return(pvalue)


}
