# ---------------------------------------------
#
# Contact Aaron J. Molstad (amolstad@ufl.edu)
# for citation instructions 
# 
# Please raise an issues at 
# github.com/ajmolstad/BvMultinom
#
# ---------------------------------------------


# -------------------------------------------
#
# Get MLE for intercept with beta = 0
#
# -------------------------------------------
BvCat_Intercept <- function(Y){

	# ----------------------------------
	# preliminaries
	# -----------------------------------
	n <- dim(X)[1]
	p <- dim(X)[2] - 1

	# --------------------------------------
	# get number of response categories
	# --------------------------------------
	J <- length(unique(Y[[1]]))
	K <- length(unique(Y[[2]]))
	Ymat <- matrix(0, nrow=n, ncol=J*K)
	for(k in 1:n){
	Ymat[k,(Y[[2]][k]-1)*J + Y[[1]][k]] <- 1
	}
	margPi <- colSums(Ymat)/n
	beta <- log(margPi) 

	return(list("beta" = beta - mean(beta)))

}






# ----------------------------------------------
# 
# Evaluate proximal operator 
#
# ----------------------------------------------
genProx <- function(y, beta.init, gamma, lambda, L, DtDinvD, PDperp, DtD, D){
	
	# -------------------------------------------
	# create result matrix with correct dimension
	# -------------------------------------------
	beta.out <- y
	Llambda <- L*lambda

	# ------------------------------------------
	# Start with 2x2 closed form
	# ------------------------------------------
	if(dim(y)[2] == 4){

		for(kk in 1:dim(y)[1]){ 
			temp <- y[kk,]
			u <- temp[1] - temp[2] - temp[3] + temp[4]

			if(abs(u/(4*L*lambda)) < 1){
				beta <- c(temp[1] -  u/4, temp[2] + u/4, temp[3] + u/4, temp[4] - u/4)
			} else {
				if(u >= 4*L*lambda){
					 beta <- c(temp[1] - Llambda, temp[2] + Llambda, temp[3] + Llambda, temp[4] - Llambda)
				} else {
					beta <- c(temp[1] + Llambda, temp[2] - Llambda, temp[3] - Llambda, temp[4] + Llambda)
				}
			}
			beta.out[kk,] <- beta*max(1 - L*gamma/sqrt(sum(beta^2)), 0)
		}

	# ------------------------------------------
	# Solve for arbitrary J and K
	# ------------------------------------------

	} else {

		# --------------------------------------
		# exact solution screening (case i)
		# --------------------------------------
		zero.inds <- which(apply(y, 1, function(x){sqrt(sum(x^2))}) < L*gamma)
		beta.out[zero.inds,] <- 0

		# -----------------------------------------
		# First closed form solution (case ii)
		# -----------------------------------------
		exact.inds <- which(apply(tcrossprod(DtDinvD,y), 2, function(x){sqrt(sum(x^2))}) < L*lambda)
		exact.inds <- setdiff(exact.inds, zero.inds)
		temp <- t(tcrossprod(PDperp, y))
		beta.out[exact.inds,] <- temp[exact.inds,,drop=FALSE]*apply(temp[exact.inds,,drop=FALSE], 1, function(x){max(1- L*gamma/sqrt(sum(x^2)), 0)})
	
		# ----------------------------------------
		# Second closed form solution (case iii)
		# ----------------------------------------
		solve.inds <- setdiff(1:dim(y)[1], union(zero.inds, exact.inds))
		keepSVD <- svd(D)

		if(length(solve.inds) > 0){
			for(ss in solve.inds){ 
				Xtmp <- c(tcrossprod(y[ss,],t(keepSVD$u)))
				temp <- (keepSVD$d[1]*sqrt(sum(Xtmp[(abs(keepSVD$d) > 1e-8)]^2)) - L*lambda*keepSVD$d[1]^2)/(L*lambda)
				utemp <- y[ss,] - crossprod(t(D),crossprod(ginv(DtD + diag(temp, dim(D)[2])), crossprod(D, y[ss,])))
				beta.out[ss,] <- utemp*max(1 - L*gamma/sqrt(sum(utemp^2)), 0)
			}
		}

	}

	return(beta.out)

}


# ------------------------------------------------
#
# Evaluate the objective function 
#
# ------------------------------------------------
eval_obj <- function(X, beta, Ymat, D, gamma, lambda){

	n <- dim(X)[1]
	lleval <- crossprod(t(X), beta)
	inner <- -(n^(-1))*sum(log(rowSums(Ymat*exp(lleval))) - log(rowSums(exp(lleval))))
	temp <- crossprod(t(beta[-1,]),D)
	pen1 <- apply(temp, 1, function(x){sqrt(sum(x^2))})
	pen2 <- apply(beta[-1,], 1, function(x){sqrt(sum(x^2))})
	return(inner + lambda*sum(pen1) + gamma*sum(pen2))

}

# ------------------------------------------------
#
# Evaluate the negative log-likelihood
#
# ------------------------------------------------
eval_lik <- function(X, beta, Ymat){

	n <- dim(X)[1]
	lleval <- crossprod(t(X), beta)
	inner <- -(n^(-1))*sum(log(rowSums(Ymat*exp(lleval))) - log(rowSums(exp(lleval))))
	return(inner)

}




# ------------------------------------------------
#
# Evaluate the negative log-likelihood
#
# ------------------------------------------------
eval_grad <- function(X, beta, Ymat){

	n <- dim(X)[1]
	p <- dim(X)[2] - 1
	temp <- exp(crossprod(t(X), beta))
	Pmat <- temp/(tcrossprod(rowSums(temp), rep(1, dim(Ymat)[2])))
	return((n^(-1))*crossprod(X, Pmat - Ymat))

}

# ----------------------------------------------
#
# Main APG function
#
# ----------------------------------------------
BvCat_AccProxGD <- function(Y, X, D, lambda, gamma, tol = 1e-6, 
			max.iter = 1e4, beta.warm = NULL, inner.quiet = TRUE){
	
	# ------------------------------
	# preliminaries 
	# ------------------------------
	n <- dim(X)[1]
	p <- dim(X)[2] - 1
	J <- length(unique(Y[[1]]))
	K <- length(unique(Y[[2]]))
	
	# ------------------------------ 
	# create response matrix 
	# ------------------------------
	Ymat <- matrix(0, nrow=n, ncol=J*K)
	for(k in 1:n){
		Ymat[k,(Y[[2]][k]-1)*J + Y[[1]][k]] <- 1
	}
	
	# ------------------------------
	# precompute needed matrices
	# ------------------------------
	DtDinvD <- crossprod(ginv(crossprod(D)), t(D))
	PDperp <- diag(1, J*K) - D%*%DtDinvD
	DtD <- crossprod(D)

	# ------------------------------
	# get initial values 
	# ------------------------------
	if(!is.null(beta.warm)){
		beta.km2 <- matrix(beta.warm, nrow=p+1, ncol=J*K)
		beta.km1 <- beta.km2
		beta.k <- matrix(beta.warm, nrow=p+1, ncol=J*K)
	} else {
		beta.km2 <- matrix(0, nrow=p+1, ncol=J*K)
		beta.km1 <- beta.km2
		beta.k <- matrix(0, nrow=p+1, ncol=J*K)
	}
	
	L0 <- 1/(sqrt(J*K)*sum(X^2)/n)
	t <- 1
	alphakm2 <- 1
	alphakm1 <- 1
	pll <-  eval_obj(X, beta.km1, Ymat, D, gamma, lambda)
	grad <- eval_grad(X, beta.km1, Ymat)
	lik <- eval_lik(X, beta.km1, Ymat)
	pll <- rep(0, max.iter)
	pll[1] <- Inf

	for(t in 2:max.iter){
		
		# --------------------
		# Update
		# --------------------
		L <- L0*2000 # can modify to improve performance
		Atemp <- beta.km1 + ((alphakm2 - 1)/alphakm1)*(beta.km1 - beta.km2)
		grad <- eval_grad(X, Atemp, Ymat)
		lik <- eval_lik(X, Atemp, Ymat)
		linesearch <- TRUE
		beta.temp <- beta.k
		
		# ----------------------------
		# Proximal update
		# ----------------------------
		while(linesearch){
			
			beta.up <- Atemp - (L)*(grad)
			beta.temp <- matrix(0, nrow=p+1, ncol=J*K)
			beta.temp[1,] <- beta.up[1,]
			beta.temp[2:(p+1),1:(J*K)] <- genProx(beta.up[2:(p+1),1:(J*K)], 
																				 beta.km1[-1,], gamma, lambda, L, DtDinvD, PDperp, DtD, D)
			templik <- eval_lik(X, beta.temp, Ymat)     

			if(L == L0){

				linesearch <- FALSE
				pll[t] <- eval_obj(X, beta.temp, Ymat, D, gamma, lambda)
				beta.km2 <- beta.km1
				beta.k <- beta.temp

			} else {
				
				if(templik < (lik + sum(diag(crossprod(grad, beta.temp - Atemp))) + 1/(2*L)*sum((beta.temp -Atemp)^2))){
						
					linesearch <- FALSE
					beta.k <- beta.temp
					beta.km2 <- beta.km1
					pll[t] <- eval_obj(X, beta.k, Ymat, D, gamma,lambda)
						
				} else {
					L <- max(L/2, L0)
				}

			}
			
		}

		beta.km1 <- beta.k
		alphakm2 <- alphakm1
		alphakm1 <- (1 + sqrt(1 + 4*alphakm1^2))/2
		if(!inner.quiet){
			cat("# ------------- ", "\n")
			cat(pll[t], "\n")
		}

		if(t > 3){
			if(abs(pll[t] - pll[t-1]) < tol && abs(pll[t-1] - pll[t-2]) < tol && abs(pll[t-2] - pll[t-3]) < tol){
				break
			}
		}
	}

	return(list("beta" = beta.k))

}



# -----------------------------------------------------
#
# Main function for fitting solution path + CV
#
# -----------------------------------------------------
BvCat.cv <- function(X, Y, ngamma = 100, lambda.vec = seq(.01,1, length=10), nfolds = NULL, delta = .1, 
										standardize = TRUE, tol = 1e-6, quiet = TRUE, 
										inner.quiet= TRUE){
	
	 

	eval_obj <- function(X, beta, Ymat, D, gamma, lambda){
		n <- dim(X)[1]
		lleval <- crossprod(t(X), beta)
		inner <- -(n^(-1))*sum(log(rowSums(Ymat*exp(lleval))) - log(rowSums(exp(lleval))))
		temp <- crossprod(t(beta[-1,]),D)
		pen1 <- apply(temp, 1, function(x){sqrt(sum(x^2))})
		pen2 <- apply(beta[-1,], 1, function(x){sqrt(sum(x^2))})
		return(inner + lambda*sum(pen1) + gamma*sum(pen2))
	}


	eval_lik <- function(X, beta, Ymat){
		n <- dim(X)[1]
		lleval <- crossprod(t(X), beta)
		inner <- -(n^(-1))*sum(log(rowSums(Ymat*exp(lleval))) - log(rowSums(exp(lleval))))
		return(inner)
	}

	eval_grad <- function(X, beta, Ymat){
		n <- dim(X)[1]
		p <- dim(X)[2] - 1
		temp <- exp(crossprod(t(X), beta))
		Pmat <- temp/(tcrossprod(rowSums(temp), rep(1, dim(Ymat)[2])))
		return((n^(-1))*crossprod(X, Pmat - Ymat))
	}

	# ------------------------------
	# preliminaries 
	# ------------------------------
	n <- dim(X)[1]
	p <- dim(X)[2]

	# ------------------------------
	# standardize if necessary 
	# ------------------------------
	if(standardize){
		x <- (X - rep(1, n)%*%t(apply(X, 2, mean)))/(rep(1, n)%*%t(apply(X, 2, sd)))
	} else {
		x <- X
	}
	x <- cbind(1, x)
	J <- length(unique(Y[[1]]))
	K <- length(unique(Y[[2]]))
	Ymat <- matrix(0, nrow=n, ncol=J*K)
	for(k in 1:n){
		Ymat[k,(Y[[2]][k]-1)*J + Y[[1]][k]] <- 1
	}
	
	# -----------------------------------
	# Get D-matrix 
	# ----------------------------------
	D <- matrix(0, nrow=J*K, choose(J,2)*choose(K,2))
	d1Pairs <- combn(1:J, 2)
	d2Pairs <- combn(1:K, 2)
	qq <- 1
	for(k in 1:dim(d1Pairs)[2]){
		dimX <- d1Pairs[,k]
		for(j in 1:dim(d2Pairs)[2]){
			dimY <- d2Pairs[,j]
			indices <- c((dimY - 1)*J + dimX[1], (dimY - 1)*J + dimX[2])
			D[sort(indices),qq] <- c(1,-1,-1,1)
			qq <- qq + 1
		}
	}

	# -----------------------------------------
	# Get candidate tuning parameter
	# -----------------------------------------
	gamma.vec <- rep(0, length=ngamma)
	interceptTemp <- BvCat_Intercept(Y= Y)$beta
	gradTemp <- eval_grad(x, rbind(interceptTemp, matrix(0, nrow=p, ncol=J*K)), Ymat)
	matTemp <- sqrt(rowSums(gradTemp[-1,]^2))

	gamma.max <-  max(matTemp)
	gamma.min <- delta*gamma.max
	gamma.vec <- 2^seq(log2(gamma.max), log2(gamma.min), length=ngamma)
	
	beta.full <- Matrix(0, nrow = (p+1)*J*K, ncol = length(gamma.vec)*length(lambda.vec), sparse=TRUE)
	beta.old <- rbind(interceptTemp, matrix(0, nrow=p, ncol=J*K))
	beta.temp <- beta.old


	for(jj in 1:length(lambda.vec)){
		beta.old <- rbind(interceptTemp, matrix(0, nrow=p, ncol=J*K))
		for(kk in 1:length(gamma.vec)){
			temp <- BvCat_AccProxGD(
					Y = Y, X = x, D=D, lambda = lambda.vec[jj], gamma = gamma.vec[kk], tol = tol,
					max.iter = 1e4, beta.warm = beta.old, inner.quiet = inner.quiet)
			beta.old <- temp$beta
			beta.full[,(jj-1)*length(gamma.vec) + kk] <- c(beta.old)
			if(!quiet){
				cat("lambda = ", lambda.vec[jj], "; gamma = ", gamma.vec[kk],"; nonzero = ", sum(abs(beta.old) > 1e-8), "\n")
			}
		}
	}
	
	
	# ----------------------------------------------------
	# perform cross-validation
	# ----------------------------------------------------
	
	if(!is.null(nfolds)){
		

		cv.index <- list(NA)

		# ------------------------------------------
		# draw indices equally across classes
		# -----------------------------------------
		for(jj in 1:(J*K)){
			fold <- sample(rep(1:nfolds, length=colSums(Ymat)[jj]))
			cv.index[[jj]] <- split(which(Ymat[,jj]==1), fold)
		}

		tot.cv.index <- list(NA)
		for(jj in 1:nfolds){
			temp <- NULL
			for(kk in 1:(J*K)){
				temp <- c(temp, cv.index[[kk]][[jj]])
			}
			tot.cv.index[[jj]] <- temp
		}

		cv.index <- tot.cv.index
		tot.cv.index <- NULL

		# -----------------------------------------
		# Preload metric arrays
		# -----------------------------------------
		errs.joint <- array(0, dim=c(length(lambda.vec),length(gamma.vec), nfolds))
		deviance <- array(0, dim=c(length(lambda.vec),length(gamma.vec), nfolds))


		for(k in 1:nfolds){
			
			# ---------------------------------------------
			# center training X and get training Y
			# ---------------------------------------------
			ntrain <- dim(X[-cv.index[[k]],])[1]
			if(standardize){      
				x.inner <- (X[-cv.index[[k]], ] - tcrossprod(rep(1, ntrain), apply(X[-cv.index[[k]],], 2, mean)))/tcrossprod(rep(1, ntrain),apply(X[-cv.index[[k]], ], 2, sd))
			} else {
				x.innd <- X[-cv.index[[k]], ]
			}
			x.inner <- cbind(1, x.inner)
			Ytrain <- Y
			Ytrain[[1]] <- Ytrain[[1]][-cv.index[[k]]]
			Ytrain[[2]] <- Ytrain[[2]][-cv.index[[k]]]
			
			# ------------------------------------------------
			# center testing X and get testing Y
			# -----------------------------------------------
			n.test <- length(cv.index[[k]])
			if(standardize){
				x.test <- (X[cv.index[[k]], ] - tcrossprod(rep(1, n.test),apply(X[-cv.index[[k]],], 2, mean)))/tcrossprod(rep(1, n.test), apply(X[-cv.index[[k]], ], 2, sd))
			} else {
				x.test <- X[cv.index[[k]], ]
			}
			x.test <- cbind(1, x.test)
			Ytest <- Y
			Ytest[[1]] <- Ytest[[1]][cv.index[[k]]]
			Ytest[[2]] <- Ytest[[2]][cv.index[[k]]]
			YmatTest <- matrix(0, nrow=n.test, ncol=J*K)
			for(jj in 1:n.test){
				YmatTest[jj,(Ytest[[2]][jj]-1)*J + Ytest[[1]][jj]] <- 1
			}
			
			# --------------------------------------------------
			# Set initial value and compute solution path
			# --------------------------------------------------
			beta.old <- rbind(interceptTemp, matrix(0, nrow=p, ncol=J*K))
			beta.temp <- beta.old

			for(jj in 1:length(lambda.vec)){
				beta.old <- rbind(interceptTemp, matrix(0, nrow=p, ncol=J*K))
				for(kk in 1:length(gamma.vec)){
					temp <- BvCat_AccProxGD(
						Y = Ytrain, X = x.inner, D=D, lambda = lambda.vec[jj], gamma = gamma.vec[kk], tol = tol,
						max.iter = 1e4, beta.warm = beta.old, inner.quiet = inner.quiet)
					beta.old <- temp$beta
					lleval <- exp(crossprod(t(x.test), beta.old))
					Pmat <- lleval/tcrossprod(rowSums(lleval), rep(1, dim(lleval)[2]))
					errs.joint[jj,kk,k] <- sum(apply(Pmat, 1, which.max) != apply(YmatTest, 1, function(x){which(x==1)}))/dim(YmatTest)[1]
					deviance[jj,kk, k] <- -2*sum(apply(YmatTest*log(Pmat/YmatTest), 1, function(x){sum(x, na.rm=TRUE)}))
				}
			}
			if(!quiet){
				cat("Through CV fold", k, "\n")
			}
		}

	}
	
	if(!is.null(nfolds)){
	
		ind <- which(apply(errs.joint, c(1,2), mean) == min(apply(errs.joint, c(1,2), mean)), arr.ind=TRUE)
		lam.min.joint <- lambda.vec[ind[1,1]]
		gamma.min.joint <- gamma.vec[ind[1,2]]

	} else {

		errs.joint <- NULL
		deviance <- NULL
		lam.min.joint <- NULL
		gamma.min.joint <- NULL
		cv.index <- NULL

	}
		
	fit <-  list("beta" = beta.full, 
							 "J" = J, "K" = K, "p" = p,
							 "D" = D, 
							 "standardize" = standardize,
							 "errs.joint" = errs.joint,
							 "deviance" = deviance,
							 "lambda.min.joint" = lam.min.joint,
							 "gamma.min.joint" = gamma.min.joint,
							 "X.mean" = apply(X, 2, mean),
							 "X.sd" = apply(X, 2, sd),
							 "cv.index" = cv.index,
							 "lambda.vec" = lambda.vec,
							 "gamma.vec" = gamma.vec)
	
	class(fit) <- "BvCat"
	return(fit)
}


# -----------------------------------------------------
#
# Extract BvCat coefficients
#
# -----------------------------------------------------
BvCat.coef <- function(fit, lambda = NULL, gamma = NULL, type="matrix"){
	
	warnings("Note: coefficients are on the scale of standardized predictors!")
	if(is.null(lambda) & is.null(gamma)){
		lam.ind <- which(fit$lambda.vec == fit$lambda.min.joint)
		gamma.ind <- which(fit$gamma.vec == fit$gamma.min.joint)
		ind <- (lam.ind-1)*length(fit$gamma.vec) + gamma.ind
		if(is.null(fit$lambda.min.joint)){
			stop('No tuning parameters selected by CV')
		}
	}  else {
		lam.ind <- which(fit$lambda.vec == lambda)
		gamma.ind <- which(fit$gamma.vec == gamma)
		ind <- (lam.ind-1)*length(fit$gamma.vec) + gamma.ind
	}
	classref <- matrix(0, nrow=2, ncol=fit$J*fit$K)
	classref[1,] <- rep(1:fit$J, fit$K)
	classref[2,] <- rep(1:fit$K, each=fit$J)

	if(class(fit)!="BvCat"){
		stop('fit needs to be of class BvCat!')
	}

	beta.mat <- matrix(fit$beta[,ind], ncol=fit$J*fit$K)
	
	if(type=="matrix"){
		return(list("b0" = beta.mat[1,], "beta" = beta.mat[-1,], "classref" = classref))
	} else {
		beta.array <- array(0, dim=c(dim(beta.mat)[1], fit$J, fit$K))
		for(j in 1:fit$K){
			beta.array[,,j] <- beta.mat[,(1:fit$J) + ((j-1)*fit$J)]
		}
		return(list("b0" = beta.array[1,,,drop=FALSE], "beta" = beta.array[-1,,,drop=FALSE], "classref" = classref))
	}
	
}



# -----------------------------------------------------
#
# Make predictions from fitted model
#
# -----------------------------------------------------
BvCat.predict <- function(Xtest, fit, lambda = NULL, gamma=NULL, type="class"){
	
	warnings("Note: coefficients are on the scale of standardized predictors!")
	 if(is.null(lambda) & is.null(gamma)){
		lam.ind <- which(fit$lambda.vec == fit$lambda.min.joint)
		gamma.ind <- which(fit$gamma.vec == fit$gamma.min.joint)
		ind <- (lam.ind-1)*length(fit$gamma.vec) + gamma.ind
		if(is.null(fit$lambda.min.joint)){
			stop('No tuning parameters selected by CV')
		}
	}  else {
		lam.ind <- which(fit$lambda.vec == lambda)
		gamma.ind <- which(fit$gamma.vec == gamma)
		ind <- (lam.ind-1)*length(fit$gamma.vec) + gamma.ind
	}
	
	if(class(fit)!="BvCat"){
		stop('fit needs to be of class BvCat!')
	}
	
	if(type!="class" & type!="probabilities"){
		stop('type needs to be either "class" or "probabilities"')
	}
	p <- dim(Xtest)[2]
	beta.mat <- matrix(fit$beta[,ind], nrow=p+1)
	if(fit$standardize){
		x.test <- (Xtest - rep(1, dim(Xtest)[1])%*%t(fit$X.mean))/(rep(1, dim(Xtest)[1])%*%t(fit$X.sd))
	} else { 
		x.test <- Xtest
	}
	lleval <- exp(crossprod(t(cbind(1, x.test)), beta.mat))
	Pmat <- lleval/tcrossprod(rowSums(lleval), rep(1, dim(lleval)[2]))

	if(type=="class"){
		classref <- matrix(0, nrow=2, ncol=fit$J*fit$K)
		classref[1,] <- rep(1:fit$J, fit$K)
		classref[2,] <- rep(1:fit$K, each=fit$J)
		preds <- matrix(0, nrow=dim(Xtest)[1], ncol=2)
		temp <- apply(Pmat, 1, which.max)
		for(j in 1:dim(Xtest)[1]){
			preds[j,] <- classref[,temp[j]]
		}
	}
	
	if(type=="probabilities"){
		preds <- array(0, dim=c(dim(Xtest)[1], fit$J, fit$K))
		for(j in 1:K){
			preds[,,j] <- Pmat[,(1:fit$J) + ((j-1)*fit$J)]
		}
	}
	
	return(list("preds" = preds))
	
}
