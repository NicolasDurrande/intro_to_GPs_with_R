#' Prediction from a Gaussian Process Regression model
#'
#' Interpolation or approximation of a function (say \code{f}) given some observation points.
#' 
#' @param x $p \times d$ matrix of the points to predict (\code{p} is the number of prediction points and \code{d} is the input space dimension).
#' @param X $n \times d$ matrix of the design of experiments (\code{n} is the number of points and \code{d} is the input space dimension).
#' @param F $n \times 1$ matrix corresponding to the values of \code{f} at \code{X}.
#' @param kern covariance function of the GP prior for \code{f} (the prior is assumed to be centred). 
#' 
#' @return a list with the mean (\code{n*1}-matrix) and the covariance matrix (\code{n*n}) of the posterior distribution at \code{x}. 
#' 
#' @examples
#' n <- 5
#' p <- 101
#' x <- matrix(seq(from=0, to=1, length=p))
#' X <- matrix(seq(from=0.1, to=0.9, length=n))
#' F <- matrix(c(0.5, 0, 1.5, 3, 2))
#' 
#' pred <- predGPR(x,X,F,kGauss)
#'                
#' pred[[1]] # or pred$mean
#' pred[[2]] # or pred$cov
#' 
#' @export
predGPR <- function(x,X,F,kern,param=NULL,kernNoise=kWhite,paramNoise=0){
  kxX <- kern(x,X,param)
  kXX_1 <- solve(kern(X,X,param) + kernNoise(X,X,paramNoise))
  m <- kxX %*% kXX_1 %*% F
  K <- kern(x,x,param) - kxX %*% kXX_1 %*% t(kxX)
  return(list(mean=m,cov=K))
}
