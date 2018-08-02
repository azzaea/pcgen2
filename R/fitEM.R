# taken from EM_function_v6.R , and adapted

fitEM <- function(y, X.t, Z.t, K = NULL, Vg = NULL, Ve = NULL,
                                cov.error = TRUE, stop.if.significant = FALSE,
                                null.dev = NULL, alpha = 0.01, max.iter = 100,
                                Vg.start = NULL, Ve.start = NULL) {

#y = em.vec; K = NULL; null.dev = NULL; Vg.start = as.numeric(Vg.manova)[c(1,4,2)]; stop.if.significant= F; Vg = NULL; Ve = NULL; Ve.start = c(as.numeric(Ve.manova)[c(1,4)], 0); cov.error = FALSE; max.iter = 5

  if (stop.if.significant == TRUE & is.null(null.dev)) {stop('No null.dev given')}

  if (!is.null(Vg.start) & !is.null(Vg)) {
    warning('Vg set to NULL')
    Vg <- NULL
  }

  if (!is.null(Ve.start) & !is.null(Ve)) {
    warning('Ve set to NULL')
    Ve <- NULL
  }

	Ve.aux <- Ve
	Vg.aux <- Vg

	#start  <-  proc.time()[3]
	X.t <- Matrix(X.t)

	if(!is.null(K)) {
		eig <- eigen(K)
		U <- eig$vectors
		d <- eig$values
		K.trans <- U %*% diag(1/sqrt(d)) 
		Z.t <- Matrix(Z.t%*% K.trans)
		#Z.t <- Matrix(Z.t%*%U%*%diag(1/sqrt(d)))
	} else {
		Z.t <- Matrix(Z.t)
	}

	# Design matrices
	Z <- Diagonal(2)%x%Z.t
	X <- Diagonal(2)%x%X.t

	# Extract some information
	ngeno <- ncol(Z.t)
	N <- nrow(Z.t)

	# Needed matrices: we take the advantage of the kronecker product
	MM <- cbind(X,Z)
	XtX. <- crossprod(X.t)		#t(X.t)%*%X.t
	XtZ. <- crossprod(X.t, Z.t)	#t(X.t)%*%Z.t
	ZtZ. <- crossprod(Z.t)		#t(Z.t)%*%Z.t


	# Weights: for the moment all ones
	w <- rep(1, length(y))

	# Number of coefficients (fixed and random, per trait)
	np <- c(ncol(X), ngeno, ngeno)

	# Initial values
	devold  <- 1e10
	thr <- 1e-3

  ###########################
	if (is.null(Vg)) {
		est.Vg <- TRUE
		est.Vg.var <- TRUE
    if (is.null(Vg.start)) {
  	  Vg <- c(1,1,0.1)
    } else {
      Vg <- Vg.start
    }
	} else {
		if (length(Vg.aux) == 2) {
			est.Vg <- TRUE
			est.Vg.var <- FALSE
			Vg <- c(Vg.aux[1], Vg.aux[2], 0.1)
		} else {
      if (length(Vg.aux) == 3) {
			  est.Vg <- FALSE
		  } else {
			  stop('The specified variance/covariance matrix for the error component is not correct')
		  }
    }
	}

  ###########################
	if (is.null(Ve)) {
		est.Ve <- TRUE
		est.Ve.var <- TRUE
		if (cov.error) {
      if (is.null(Ve.start)) {
  			Ve <- c(1,1,0.1)
      } else {
        Ve <- Ve.start
      }
		} else {
      if (is.null(Ve.start)) {
  			Ve <- c(1,1,0)
      } else {
        Ve    <- Ve.start
        Ve[3] <- 0
      }
		}
	} else {
		if (length(Ve.aux) == 2) {
			est.Ve <- TRUE
			cov.error <- TRUE
			est.Ve.var <- FALSE
			Ve <- c(Ve.aux[1], Ve.aux[2], 0.1)
		} else if (length(Ve.aux) == 3) {
			est.Ve <- FALSE
		} else {
			stop('The specified variance/covariance matrix for the error component is not correct')
		}
	}

	# Precision matrices for the genetic variances (needed for SAP/Schall algorithm)
	g1 <- rep(c(Vg[1], 0), each = ngeno)
	g2 <- rep(c(0, Vg[2]), each = ngeno)

	for (it in 1:max.iter) {
		# Genotypic covariance matrix  # it=1
		Gi <- matrix(c(Vg[1], Vg[3], Vg[3], Vg[2]), ncol = 2)
		G  <- Gi%x%Diagonal(ngeno)

		Giinv <- solve(Gi)
		Ginv  <- Giinv%x%Diagonal(ngeno)

		# Error covariance matrix
		Ri <- matrix(c(Ve[1], Ve[3], Ve[3], Ve[2]), ncol = 2)
		R <- Ri%x%Diagonal(N)

		Riinv <- solve(Ri)
		Rinv  <- Riinv%x%Diagonal(N)

		# X'X X'Z
		# Z'X Z'Z

		# We take the advantage of the kronecker product
		XtRinvX. <- Riinv%x%XtX.
		XtRinvZ. <- Riinv%x%XtZ.
		ZtRinvZ. <- Riinv%x%ZtZ.

		XtRinvy. <- ((Riinv%x%t(X.t))%*%y)[,1]
		ZtRinvy. <- ((Riinv%x%t(Z.t))%*%y)[,1]

		u <- c(XtRinvy.,ZtRinvy.)

		V <- construct.block(XtRinvX., XtRinvZ., t(XtRinvZ.), ZtRinvZ.)

		#D <- Matrix:::bdiag(diag(rep(0,np[1])), Ginv)
		D <- bdiag(diag(rep(0,np[1])), Ginv)
		
		# Henderson system of equations
		H <- V + D

		Hinv <- try(solve(H))
		if(class(Hinv) == "try-error") {
			Hinv <- ginv(as.matrix(H))
		}

		# Fixed and random coefficients
		b <- Hinv%*%u

		b.fixed  <- b[1:np[1]]
		b.random <- b[-(1:np[1])]

		# Compute the deviance
		res <- (y - MM%*%b) # residuals
		dev <- deviance2(H, Gi, ngeno, Ri, N, Rinv, res, t(b.random)%*%Ginv%*%b.random)[1]
		#dev <- deviance(H, G, R, Rinv, res, t(b.random)%*%Ginv%*%b.random)[1]
		if(!est.Ve & !est.Vg) {
			break
		}
		#########################################################
		# Genotypic variance components
		#########################################################
		if (est.Vg) {
			#########################################################
			# EM algorithm
			#########################################################
			Ak <- Hinv[-(1:np[1]),-(1:np[1])]
			#if(est.Vg.var) {
				# First trait
				#A <- Ak[1:ngeno, 1:ngeno]
				#tau1 <- (1/ngeno)*(sum(diag(A)) + t(b.random[1:ngeno])%*%b.random[1:ngeno])
				#tau1 <- (1/ngeno)*(sum(diag(A)) + sum(b.random[1:ngeno]*b.random[1:ngeno]))

				# Second trait
				#A <- Ak[(ngeno+1):(2*ngeno), (ngeno+1):(2*ngeno)]
				#tau2 <- (1/ngeno)*(sum(diag(A)) + t(b.random[(ngeno+1):(2*ngeno)])%*%b.random[(ngeno+1):(2*ngeno)])
				#tau2 <- (1/ngeno)*(sum(diag(A)) + sum(b.random[(ngeno+1):(2*ngeno)]*b.random[(ngeno+1):(2*ngeno)]))
			#} else{
			#	tau1 <- Vg[1]
			#	tau2 <- Vg[2]
			#}
			# Covariance
			A <- Ak[1:ngeno, (ngeno+1):(2*ngeno)]
			#tau3 <- (1/ngeno)*(sum(diag(A)) + t(b.random[1:ngeno])%*%b.random[(ngeno+1):(2*ngeno)])
			tau3 <- (1/ngeno)*(sum(diag(A)) + sum(b.random[1:ngeno]*b.random[(ngeno+1):(2*ngeno)]))
			#########################################################
			# Schall: only for the variances (apparently faster)
			#########################################################
			if (est.Vg.var) {
				aux <- diag(G) - diag(Hinv[-(1:np[1]),-(1:np[1])])
				# First trait
				g.inv.d <- (1/Vg[1])*g1
				ed1  <- sum(aux*g.inv.d)
				ssv1 <- sum(b.random^2*g1)
				tau1 <- ssv1/ed1

				# Second trait
				g.inv.d <- (1/Vg[2])*g2
				ed2  <- sum(aux*g.inv.d)
				ssv2 <- sum(b.random^2*g2)
				tau2 <- ssv2/ed2
			} else {
				tau1 <- Vg[1]
				tau2 <- Vg[2]
			}
			Vg.new <- c(tau1, tau2, tau3)
		} else {
			Vg.new = Vg
		}
		#########################################################
		# Error variance components
		#########################################################
		if (est.Ve) {
			#########################################################
			# EM algorithm
			#########################################################
			resm <- matrix(res, ncol = 2) # One column per trait
			ress <- t(resm)%*%resm
			aux <- MM%*%Hinv #####time consuming try crossprod {base}
			if(est.Ve.var) {
				diag.var <- rowSums(aux*MM)
				# First trait
				aux1  <- sum(diag.var[1:N])
				aux2  <- ress[1,1]
				sig21 <- (1/N)*(aux1 + aux2)

				# Second trait
				aux1  <- sum(diag.var[(N+1):(2*N)])
				aux2  <- ress[2,2]
				sig22 <- (1/N)*(aux1 + aux2)

				# To be further studied
				#sig21 <- ress[1,1]/(N - ncol(X.t) - ed1)
				#sig22 <- ress[2,2]/(N - ncol(X.t) - ed2)

			} else {
				sig21 <- Ve[1]
				sig22 <- Ve[2]
			}
			if (cov.error) {
				#aux <- MM%*%Hinv
				diag.covar <- rowSums(aux[1:N,]*MM[(N+1):(2*N),])
				# Covariance
				aux1 <- sum(diag.covar)
				aux2 <- ress[1,2]
				sig212 = (1/N)*(aux1 + aux2)
			} else {
				sig212 = 0
			}
			Ve.new <- c(sig21, sig22, sig212)
		} else {
			Ve.new <- Ve
		}

		dla <- abs(devold - dev)

		#cat(sprintf("%1$3d %2$10.6f", it, dev))
		#cat(sprintf("%8.3f", c(Ve, Vg)), '\n')

		Ve <- Ve.new #!
		Vg <- Vg.new
		devold <- dev

		if (stop.if.significant) {
			loglik_Full    <- -0.5 * dev
			loglik_Reduced <- -0.5 * null.dev

			REMLLRT <- 2 * max(loglik_Full - loglik_Reduced, 0)
			pvalue  <- (1 - pchisq(REMLLRT, df = 1))

		 if (pvalue < alpha) break
		  if (it > 20 & pvalue > 0.1) break
		} 
	
		if (dla < thr) break
	}

	#end <-  proc.time()[3]
	#cat(paste('Computing time:', end - start, '\n'))
	if(!is.null(K)) {
		b.random <- (Diagonal(2)%x%K.trans)%*%b.random
	}
	
	res           <- list()
	res$coeff     <- list(fixed = b.fixed, random = b.random)
	res$variances <- list(Ve = Ve, Vg = Vg)
	res$deviance  <- dev
    res$it        <- it
    if (stop.if.significant) {res$pvalue    <- pvalue}

	res
}