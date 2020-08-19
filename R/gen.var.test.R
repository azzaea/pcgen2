gen.var.test <- function(y,X=data.frame(),G, return.fitted = FALSE) {
# y=d[,2]; X=d[,3:4]; G=d[,1]

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

  if (return.fitted == TRUE) {
    fitted.values <- as.numeric(as.matrix(X) %*% matrix(as.numeric(coefficients(lm.obj))[2:(ncol(X)+1)]))
    return(list(av[[5]][1+ncol(X)], fitted.values))
  } else {
    return(av[[5]][1+ncol(X)])
  }


}


