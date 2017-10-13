#' Log-likelihood of a multivariate Gaussian vector
#'
#' Log-likelihood of some parametres of a GP Given some observations
#' 
#' @param vector corresponding to the concatenation of the parameters of the kernel and the noise kernel (if not NULL).
#' @param kern covariance function of the GP prior for the signal. 
#' @param Xd $n \times d$ matrix of the design of experiments (\code{n} is the number of points and \code{d} is the input space dimension).
#' @param F $n \times 1$ matrix corresponding to the values of \code{f} at \code{X}.
#' @param kernNoise covariance function of the GP prior for the noise (can be NULL or kWhite). 
#' 
#' @return scalar corresponding to the log-likelihood. 
#' 
#' @examples
#' 
#' # Generate simple Data
#' Xd <- matrix(seq(from=0.1, to=0.9, length=5))
#' F <- matrix(c(0.5, 0, 1.5, 3, 2))
#' 
#' # optimize the model parameters
#' logLikelihood(c(1,.2), kern=kMat52, Xd=Xd, F=F)
#' opt_out <- optim(c(1,.2), logLikelihood, kern=kMat52, Xd=Xd, F=F, control=list(fnscale=-1))
#' 
#' @export
# rodo change X into Xd as function parameter name to avoid confusion when used with apply(X=)
logLikelihood <- function(params,kern,Xd,F,kernNoise=NULL){
  if(is.null(kernNoise)){
    kXX <- kern(Xd,Xd,params) + kWhite(Xd,Xd,1e-6)
  }else if(kernNoise==kWhite){
    param = params[-lenght(params)]
    paramNoise = params[lenght(params)]
    kXX <- kern(Xd,Xd,param) + kernNoise(Xd,Xd,paramNoise)
  }else{
    stop("the (log-)likelihood is only implemented for kernNoise=NULL and kernNoise=kWhite")    
  }
  ndata <- nrow(Xd)
  # rodo add something on the diagonal to avoid singularity
  #detkXX <- det(kXX)
  #eps <- 1.e-14
  #tau <- 1.e-14
  #while (detkXX<eps & tau<1) {
  #  tau <- tau*100
  #  kXX <- kXX + tau*diag(ndata)
  #  detkXX <- det(kXX)
  #}
  LL <- -1/2*ndata*log(2*pi) - 1/2*determinant(kXX)$modulus - 1/2*t(F)%*%solve(kXX)%*%F
  return(LL)
}

# Maximize logLikelihood on params with random restarts and the R optim() function
#
# Inputs :
#   tmin, tmax : lower and upper bounds for initial params in search. Note that the output can be out of bounds
#   nbtry : number of random restarts for the maximizations of the loglikelihood (LL)
#   maxit : maxit parameter of optim() base function
#   silent : logical (T or F) for printing or not maximization steps 
#
# Outputs :  a list O with
#   O$bestthetas : nb variables + 1 vector of parameters
#   O$bestLL : scalar value of LL at O$bestthetas
maxlogLikelihood <- function(kern,Xd,F,kernNoise=NULL,tmin=0.1,tmax=5,nbtry=20,maxit=500,silent=FALSE){
  npar <- dim(Xd)[2]+1
  nbvar <- dim(Xd)[2]
  if (length(tmin)<npar) {tmin <- rep(tmin,times=npar)}
  if (length(tmax)<npar) {tmax <- rep(tmax,times=npar)}
  bestLL <- -Inf
  if (!silent) cat("\n  MAX LIKELIHOOD in ",nbtry," restarts\n\n")
  for (i in 1:nbtry){
    tinit <- tmin + runif(nbvar+1)*(tmax-tmin)
    LL <- logLikelihood(params = tinit,kern=kMat52,Xd=Xd,F=F)
    if (!silent) cat(i,"theta_init=",tinit," , LL=",LL,"\n")
    # bounds on parameters in loglikelihood are actually important to EI: arg max likelihood typically generates too large
    # length scales that create really flats plateaus on EI and error in optim()
    opt_out <- optim(tinit, fn = logLikelihood, kern=kern, Xd=Xd, F=F, method="L-BFGS-B", lower=tmin, upper=tmax,control=list(fnscale=-1, maxit=maxit))
    if (opt_out$value>bestLL){
      bestLL <- opt_out$value
      bestthetas <- abs(opt_out$par) # abs because some optimizers (that do not use bounds) go to neg values and they are equiv to positive ones
    }
    cat(" iter thetas=",opt_out$par," , iter LL=",opt_out$value,"\n")
    
  }
  O <- list()
  O$bestthetas <- bestthetas
  O$bestLL <- bestLL
  return(O)
}

