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
  xpmat <- matrix(xp,byrow=T,ncol=dim(Xd)[2])
  # kxX <- kern(xp,Xd,param) # rodo : doesn't work in more than 1d
  kxX <- kern(xpmat,Xd,param)
  kXX_1 <- solve(kern(Xd,Xd,param) + kernNoise(Xd,Xd,paramNoise))
  m <- kxX %*% kXX_1 %*% F
  # var <- rep(kern(xp[1,,drop=FALSE],xp[1,,drop=FALSE],param) - rowSums((kxX %*% kXX_1)*kxX)) # rodo : this is issuing a warning, better version below 
  var <- rep(kern(xpmat[1,,drop=FALSE],xpmat[1,,drop=FALSE],param),times=dim(xpmat)[1]) - rowSums((kxX %*% kXX_1)*kxX)
  sd <- sqrt(pmax(1e-12,var))
  u <- (min(F) - m)/sd
  ei <- sd * (u * pnorm(u) + dnorm(u))
  return(ei)
}

# Maximize EI with random restarts of the BFGS algorithm
#
# Inputs :
#   kern : kernel of the Gaussian Process
#   param : parameters of the Gaussian Process
#   Xd : n by d DoE where d is the number of variables
#   F : vector of n observations
#   xmin, xmax : lower and upper bounds for variables
#   nbtry : number of random restarts for the EI maximizations
#   maxit : maxit parameter of optim() base function
#   silent : logical (T or F) for printing or not maximization steps 
#
# Outputs :  a list best with
#   best$var : nb variables vector of parameters
#   best$EI : scalar value of EI at best$thetas
maxEI <- function(kern,Xd,F,kernNoise=kWhite,param=NULL,xmin=0,xmax=1,nbtry=50,maxit=100,silent=FALSE)
{
  nbvar <- dim(Xd)[2]
  if (length(xmin)<nbvar) {xmin <- rep(xmin,times=nbvar)}
  if (length(xmax)<nbvar) {xmax <- rep(xmax,times=nbvar)}
  best <- list()
  best$EI <- -Inf
  if (!silent) {cat("\n  MAX EI in ",nbtry," restarts of local optimization\n\n")}
  # it is important to restart this optimization as most of the time EI is too flat to allow a 
  # gradient based search to proceed. An alternative is to use the CMA-ES algorithm (TODO)
  for (i in 1:nbtry){
    xinit <- xmin+runif(nbvar)*(xmax-xmin)
    aEI <- EI(xp=xinit,Xd=Xd,F=F,kern=kern,param=param)
    if (!silent) {cat(i,"norm.xinit=",xinit," , EI=",aEI,"\n")}
    # in the optimization, it is important to remain in bounds, i.e., between 0 and 1 with the normalized variables
    tryCatch({ # optim generates errors when EI becomes too flat...
    opt_out <- optim(par=xinit,fn = EI,Xd=Xd,F=F,kern=kern,param=param, 
                     method="L-BFGS-B", lower=xmin, upper=xmax, control=list(fnscale=-1, maxit=maxit))
    }, warning = function(war) {
      print(paste("optim warning, cancelling the optim:  ",war))
      opt_out <- list()
      opt_out$value <- -Inf
      return(opt_out)
    }, error = function(err) {
      print(paste("optim error:  ",err))
      opt_out <- list()
      opt_out$value <- -Inf
      return(opt_out)
    }, finally = {}) # end tryCatch
    if (opt_out$value>best$EI){
      best$EI <- opt_out$value
      best$var <- abs(opt_out$par) # abs because some optimizers go to neg values and they are equiv to positive ones
    }
    if (!silent) {cat(" iter var=",opt_out$par," , iter EI=",opt_out$value,"\n")}
  }
  return(best)
}
