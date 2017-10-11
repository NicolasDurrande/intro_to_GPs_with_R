#' Expected improvement a Gaussian Process Regression model
#'
#' Returns the expected improvement
#' 
#' @param xp $p \times d$ matrix of the points to predict (\code{p} is the number of prediction points and \code{d} is the input space dimension).
#' @param Xd $n \times d$ matrix of the design of experiments (\code{n} is the number of points and \code{d} is the input space dimension).
#' @param F $n \times 1$ matrix corresponding to the values of \code{f} at \code{Xd}.
#' @param kern covariance function of the GP prior for \code{f} (the prior is assumed to be centred). 
#' 
#' @return a \code{n*1}-matrix of the expected improvement at \code{xp}. 
#' 
#' @examples
#' n <- 5
#' p <- 101
#' x <- matrix(seq(from=0, to=1, length=p))
#' X <- matrix(seq(from=0.1, to=0.9, length=n))
#' F <- matrix(c(0.5, 0, 1.5, 3, 2))
#' 
#' ei <- EI(xp=x,Xd=X,F,kGauss)
#' 
#' @export
EI <- function(xp,Xd,F,kern,param=NULL,kernNoise=kWhite,paramNoise=0){
  kxX <- kern(xp,Xd,param)
  kXX_1 <- solve(kern(Xd,Xd,param) + kernNoise(Xd,Xd,paramNoise))
  m <- kxX %*% kXX_1 %*% F
  # var <- rep(kern(xp[1,,drop=FALSE],xp[1,,drop=FALSE],param) - rowSums((kxX %*% kXX_1)*kxX)) # rodo : this is issuing a warning, better version below 
  var <- rep(kern(xp[1,,drop=FALSE],xp[1,,drop=FALSE],param),times=length(xp)) - rowSums((kxX %*% kXX_1)*kxX)
  sd <- sqrt(pmax(1e-12,var))
  u <- (min(F) - m)/sd
  ei <- sd * (u * pnorm(u) + dnorm(u))
  return(ei)
}
