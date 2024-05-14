\name{BvCat.coef}
\alias{BvCat.coef}
\title{Extract estimated regression coefficients for penalized bivariate multinomial logistic regression.}
\description{A function to extract coefficients along the solution path for the penalized bivariate multinomial logistic regression.}
\usage{
BvCat.coef(fit, lambda = NULL, gamma = NULL, type="matrix")
}
\arguments{
\item{fit}{An fitted model object of type \code{BvCat} obtained from the \code{BvCat.cv} function. }
\item{lambda}{Value of tuning parameter \eqn{\lambda} at which to extract coefficients. 
If \code{NULL}, will use values which minimized joint misclassificaiton error (assuming \code{nfolds} \eqn{ > 2}). If not \code{NULL}, must be a value from \code{fit$lambda.vec}.}
\item{gamma}{Value of tuning parameter \eqn{\gamma} at which to extract coefficients. 
If \code{NULL}, will use values which minimized joint misclassificaiton error (assuming \code{nfolds} \eqn{ > 2}). If not \code{NULL}, must be a value from \code{fit$gamma.vec}.}
\item{type}{In what form should the coefficients be returned? options are "\code{matrix}" or "\code{tensor}". 
If "\code{matrix}", will return as \eqn{p \times JK} matrix. If "tensor", will return a \eqn{p \times J \times K} array. }
}
\value{
 \item{b0}{Estimated intercept. Either \eqn{JK}-dimensional vector or \eqn{J \times K} matrix, depending on \code{type}.}
 \item{beta}{Estimated regression coefficients on the original scale of the predictors. 
 Either \eqn{p \times JK} matrix or \eqn{p \times J \times K} array, depending on \code{type}.}
}

\examples{
# ----------------------------------------
# Example for BvCategorical
# ----------------------------------------
set.seed(1)
p <- 50
n <- 300
J <- 3 
K <- 2
Results <- NULL
SigmaX <- matrix(0, nrow=p, ncol=p)
for(jj in 1:p){
  for(kk in 1:p){
    SigmaX[jj,kk] <- .5^abs(jj-kk)
  }
}
eo <- eigen(SigmaX) 
SigmaXsqrt <- tcrossprod(tcrossprod(eo$vec, diag(eo$val^.5)), eo$vec)

# generate data -------------
X <- tcrossprod(matrix(rnorm(n*p), nrow=n), SigmaXsqrt)

beta <- matrix(0, nrow=p, ncol=J*K)
temp <- sample(1:p, 10)

# predict affect joint probabilities -------------
for(j in 1:6){
  beta[temp[j],] <- runif(6, -3, 3)
}
# predictors affect marginal probabilities only ---
for(j in 7:10){
  temp1 <- runif(4, -3, 3)
  beta[temp[j],] <- c(-temp1[4] + temp1[3] + temp1[1], temp1[1], temp1[2], 
                      temp1[3], temp1[4], -temp1[1]+ temp1[4] + temp1[2])
}

lleval <- exp(crossprod(t(X), beta))
Pmat <- lleval/tcrossprod(rowSums(lleval), rep(1, J*K))
Ymat <- matrix(0, nrow=n, ncol=J*K)
Y <- list(NA)
Y[[1]] <- rep(0, n)
Y[[2]] <- rep(0, n)
refClasses <- rbind(c(1:J, 1:J), c(rep(1, J), rep(2, J)))
for(k in 1:n){
  Ymat[k,] <- rmultinom(1, 1, Pmat[k,])
  Y[[1]][k] <- refClasses[1,which(Ymat[k,]==1)]
  Y[[2]][k] <- refClasses[2,which(Ymat[k,]==1)]
}

# ---------------------------------------------------
# Running the following will take around 60 seconds
# --------------------------------------------------

\dontrun{
# ----------------------------------------
# Fit model with cross-validation
# ----------------------------------------
modfit <- BvCat.cv(X, Y, ngamma = 20, 
                  lambda.vec = 10^seq(-4, 0, length=5), 
                  nfolds = 5, delta = .01, standardize = TRUE, 
                  tol = 1e-8,  quiet = FALSE, inner.quiet = TRUE)

# ---------------------------------------
# CV errors
# ---------------------------------------
apply(modfit$errs.joint, c(1,2), mean) # default
apply(modfit$deviance, c(1,2), mean)

# -----------------------------------------------------
# Coefficients at tuning parameter pair selected by CV
# -----------------------------------------------------
temp.coef <- BvCat.coef(modfit, type="matrix")
temp.coef.tensor <- BvCat.coef(modfit, type="tensor")
temp.coef.untuned <- BvCat.coef(modfit, lambda = modfit$lambda.vec[1], 
 gamma = modfit$gamma.vec[10], type="matrix")



ntest <- 1
Xtest <- tcrossprod(matrix(rnorm(ntest*p), nrow=ntest),SigmaXsqrt)
temp.predict <- BvCat.predict(Xtest, modfit, type="class")
temp.predict.probs <- BvCat.predict(Xtest, modfit, type="probabilities") 
# first dimension of temp.predict.probs$preds length ntest, second J, and third K


# ----------------------------------------
# Fit model with no cross-validation
# ----------------------------------------
modfit.noCV <- BvCat.cv(X, Y, ngamma = 20, 
                 lambda.vec = 10^seq(-1, -4, length=5), 
                  nfolds = NULL, delta = .01, standardize = TRUE, 
                  tol = 1e-8,  quiet = TRUE, inner.quiet = TRUE)

# temp.coef <- BvCat.coef(modfit.noCV, type="matrix") # returns errror
temp.coef.noCV <- BvCat.coef(modfit.noCV, lambda = modfit$lambda.vec[1], 
 gamma = modfit$gamma.vec[10], type="matrix")
}
}
