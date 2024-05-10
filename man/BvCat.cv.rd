\name{BvCat.cv}
\alias{BvCat.cv}
\title{Function perform cross-validation and fit the penalized bivariate multinomial logistic regression.}
\description{A function for performing cross-validation to select the tuning parameter for the penalized bivariate multinomial logistic regression.}
\usage{
BvCat.cv(X, Y, ngamma = 100, lambda.vec = seq(.01,1, length=10), 
  nfolds = 5, delta = .1, 
  standardize = TRUE, tol = 1e-8, quiet = TRUE, 
  inner.quiet= TRUE)
}
\arguments{
\item{X}{An \eqn{n \times p} matrix of predictors. This should not include an intercept, nor should it be standardized.}
\item{Y}{A list of length two: each list corresponds to one of the two response variables. The first component of the list should be made up of integers from 1 to J; the second component of integers from 1 to K. }
\item{ngamma}{The number of candidate tuning parameters \eqn{\gamma} to be considered.}
\item{lambda.vec}{A set of candidate tuning parameters for \eqn{\lambda}. In general, values between \eqn{10^{-4}} and \eqn{1} tend to be best, but this may need to be explored for each particular application.}
\item{nfolds}{The number of folds for cross-validation. Default is \code{5}. If cross-validation is not desired, set equal to NULL.}
\item{delta}{The ratio of maximum to minimum \eqn{\gamma} candidate values. Smaller values of \eqn{\delta} would allow for less sparse model fits, but at the expensive of computing time.}
\item{standardize}{Should the predictors be standardized (centered and scaled to have unit standard deviation) for model fitting? Default is TRUE.}
\item{tol}{Convergence tolerance parameter. Convergence is claimed when three consecutive objective function value changes are below this tolerance. Default is \eqn{10^{-8}}, but for some applications, this may need to be adjusted.}
\item{quiet}{Should any progress be printed?}
\item{inner.quiet}{Should much more progress be printed? Use with caution. }
}

\value{
\item{beta}{A sparse matrix containing regression coefficient estimates. Note that if \code{standardize=TRUE}, these will be on the scale of the standardized predictors, so use with caution. It is highly recommended to extract coefficients with \code{BvCat.coef}.}
\item{J}{Number of response categories for the first response.}
\item{K}{Number of response categories for the second response.}
\item{D}{The constraint matrix used in the penalty. See description for more details.}
\item{standardize}{Were predictors standardized for model fitting?}
\item{errs.joint}{An array of dimension \code{ngamma} \eqn{\times} \code{length(lambda.vec)} \eqn{\times} \code{nfolds} of misclassification errors for the candidate tuning parameters in each of the \code{nfolds}.}
\item{deviance}{An array of dimension \code{ngamma} \eqn{\times} \code{length(lambda.vec)} \eqn{\times} \code{nfolds} of multinomial deviances for the candidate tuning parameters in each of the \code{nfolds}. Note that misclassification error is the default performance metric.}
\item{lambda.min.joint}{Tuning parameter value for \eqn{\lambda} which minimized average misclassification error.}
\item{gamma.min.joint}{Tuning parameter value for \eqn{\gamma} which minimized average misclassification error.}
\item{X.mean}{The columnwise averages of the input matrix \code{X}. Used for extracting coefficients and prediction when model fit with \code{standardize=TRUE}.}
\item{X.sd}{The columnwise standard deviations input matrix \code{X}. Used for extracting coefficients and prediction when model fit with \code{standardize=TRUE}.}
\item{cv.index}{Which observations belonged to which cross-validation folds?}
\item{lambda.vec}{Exactly the input \code{lambda.vec}.}
\item{gamma.vec}{The candidate values of \eqn{\gamma} used for model fitting.}
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
