estimate.cond.mean <- function(suffStat, j, k=NULL, S, Vg, Ve, dec = NULL,
                               covariates = NULL) {

# estimate the conditional mean of suffStat[,j] (or suffStat[,j] and suffStat[,k]),
#    given suffStat[,S]
# j,S : refer to column numbers in suffStat (G,Y1,..,Yp), not Y1,..,Yp!
# S is a vector of integers, of length at least 1
# dec : a list with elements Dk and Uk, computed by w <- eigen(Z %*% t(Z)),
#       Dk <- diag(w$values); Uk <- w$vectors
# see Zhou and Stephens 2014, supplement.

  # To do: incorporate improvements by Bart-Jan

  ##################################################
  # Check, with undiagonalized version:
  #Z <- as.matrix(make.Z.matrix(suffStat$genotype)); K <- Z %*% t(Z)
  #Yu <- matrix(as.numeric(as.matrix(suffStat[,S])))# Yu <- matrix(c(suffStat[,j], as.numeric(as.matrix(suffStat[,S]))))
  #Xu.temp <- as.matrix(cbind(rep(1,nrow(suffStat)), covariates))
  #Xu <- as.matrix(bdiag(Xu.temp,Xu.temp)) #Xu <- as.matrix(bdiag(Xu.temp,Xu.temp,Xu.temp))#as.matrix(rbind(Xu.temp,Xu.temp,Xu.temp))
  ##V <- kronecker(Vg[c(j-1,S-1), c(j-1,S-1)], K) + kronecker(Ve[c(j-1,S-1), c(j-1,S-1)], diag(nrow(suffStat)))
  #V <- kronecker(Vg[S-1, S-1], K) + kronecker(Ve[S-1, S-1], diag(nrow(suffStat)))
  #Vinv <- solve(V)
  #betaU <- solve( t(Xu) %*% Vinv %*% Xu) %*% t(Xu) %*% Vinv %*% Yu
  ##V.S <- kronecker(Vg[S-1, S-1], K) + kronecker(Ve[S-1, S-1], diag(nrow(suffStat)))
  #V.j.S <- kronecker(t(matrix(Vg[S-1, j-1])), K) + kronecker(t(matrix(Ve[S-1, j-1])), diag(nrow(suffStat)))
  #pr.tempU <- Vinv %*% (Yu - Xu %*% betaU)
  #frU <- V.j.S %*% Vinv %*% (Yu - Xu %*% betaU) # solve(V.S) # [-(1:nrow(suffStat)), ]
  ###############################################

  # suffStat = d; j=2; k=3; S=4; Vg=q1$Vg; Ve=q1$Ve;covariates = NULL; dec = NULL
  # suffStat = d; j=2; k=NULL; S=4:5; Vg=q1$Vg; Ve=q1$Ve;covariates = NULL; dec = NULL


  if (!is.null(covariates)) {
    if (ncol(covariates)==0) {covariates <- NULL}
  }

  stopifnot(length(S) > 0)
  stopifnot(all(c(j,S) %in% 2:ncol(suffStat)))
  stopifnot(!(j %in% S))

  if (!is.null(k)) {
    stopifnot(!(k %in% S))
    stopifnot(k %in% 2:ncol(suffStat))
    stopifnot(k!=j)
  }

  if (is.null(dec)) {
    # To do: take the spectral decomposition out of this function, and put it in pc2
    Z <- as.matrix(make.Z.matrix(suffStat$genotype))
    K <- Z %*% t(Z)

    # To do : more efficient spectral decomposition for the balanced case,
    #         with n genotypes with r replicates. Then Z %*% t(Z) is of rank n
    #         with eigenvalues r (multiplicity n) and 0 (multiplicity n*(r-1))
    w <- eigen(K)
    Dk <- diag(w$values)
    Uk <- w$vectors
  } else {
    Dk <- dec$Dk
    Uk <- dec$Uk
  }

  V.inv.array <- make.V.inv.array(Vg = as.matrix(Vg[S-1,S-1]),
                                  Ve = as.matrix(Ve[S-1,S-1]), Dk = Dk)

  X <- t(matrix(rep(1, nrow(suffStat))))

  if (!is.null(covariates)) {
    stopifnot(nrow(suffStat)==nrow(covariates))
    # To do: check if covariates accidentally contains an intercept already.
    # In this case, give a warning and remove it.
    X <- rbind(X, t(as.matrix(covariates)))
  }
  #X.u <- X; Y.u <- t(as.matrix(suffStat[,S]))

  X <- X %*% Uk
  Y <- t(as.matrix(suffStat[,S])) %*% Uk # c(j,S)

  nc   <- nrow(X)
  p    <- nrow(Y)
  n    <- ncol(X)

  if (p==1 & nc==1) {
    Vbeta <- matrix(sum(sapply(1:n,function(m){kronecker((matrix(X[,m])) %*% t(matrix(X[,m])),V.inv.array[m,,])})))
  } else {
    Vbeta <- matrix(apply(sapply(1:n,function(m){kronecker((matrix(X[,m])) %*% t(matrix(X[,m])),V.inv.array[m,,])}),1,sum),ncol=p*nc)
  }

  M     <- solve(Vbeta)

  if (p==1 & nc==1) {
    v     <- sum(sapply(1:n,function(m){kronecker((matrix(X[,m])),V.inv.array[m,,] %*% matrix(Y[,m]))}))
  } else {
    v     <- apply(sapply(1:n,function(m){kronecker((matrix(X[,m])),V.inv.array[m,,] %*% matrix(Y[,m]))}),1,sum)
  }

  v <- matrix(v)

  fixed.effects <- as.numeric(M %*% v)

  res <- (Y - matrix(fixed.effects, ncol=nc) %*% X)

  pr.temp <- sapply(1:n, function(i){as.matrix(V.inv.array[i,,]) %*% matrix(res[,i])})

  if (class(pr.temp)=="numeric") {pr.temp <- t(matrix(pr.temp))}

  if (is.null(k)) {

    V.array     <- make.V.array(Vg = as.matrix(Vg[c(j-1,S-1),c(j-1,S-1)]),
                                Ve = as.matrix(Ve[c(j-1,S-1),c(j-1,S-1)]), Dk = Dk)


    fr <- sapply(1:n, function(i){t(as.matrix(V.array[i,-1,1])) %*% matrix(pr.temp[,i])})

    cond.mean <- t(matrix(fr)) %*% t(Uk)

  } else {

    V.array     <- make.V.array(Vg = as.matrix(Vg[c(j-1,k-1,S-1),c(j-1,k-1,S-1)]),
                                Ve = as.matrix(Ve[c(j-1,k-1,S-1),c(j-1,k-1,S-1)]), Dk = Dk)


    #fr <- sapply(1:n, function(i){t(as.matrix(V.array[i,-(1:2),1:2])) %*% matrix(pr.temp[,i])})
	fr <- sapply(1:n, function(i){t(matrix(V.array[i,-(1:2),1:2], ncol=2)) %*% matrix(pr.temp[,i])})

    cond.mean <- fr %*% t(Uk)

  }

return(list(cond.mean = t(cond.mean), fixed.effects = fixed.effects))
}
#
