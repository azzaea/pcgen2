gen.var.test <- function(y, X = data.frame(), G, K = NULL, return.fitted = FALSE) {
  # y=d[,2]; X=d[,3]; G=d[,1] # previous set up
  # y = dm[, 2]; X = dm[,3]; G=dm[,1]; K = K # new contribution

  X <- as.matrix(X)
  if (ncol(X) == 0) { # Y_S = \phi
    X <-  matrix(1, nrow = length(y))
  } else {  # Y_S \ne \phi    
    X <- cbind(1, X)
  }
  
  lm.obj <- lmm.aireml(Y = y , X = X, K = K, verbose = F)
  
  lrt <- 2*(lm.obj$logL - lm.obj$logL0) # lrt statistic
  pvalue <- pchisq(lrt, df = 1, lower.tail = F) 
  
  if (return.fitted == TRUE) {
    # fitted.values <- as.numeric(   #I don't know how to produce this from gaston output
    #   as.matrix(X) %*% matrix(as.numeric(coefficients(lm.obj))[2:(ncol(X)+1)]))
    fitted.values <- NULL
    return(list(pvalue, fitted.values))
  } else {
    return(pvalue)
  }
  
}
