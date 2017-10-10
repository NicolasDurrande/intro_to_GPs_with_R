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
    kXX <- kern(Xd,Xd,params) 
  }else if(kernNoise==kWhite){
    param = params[-lenght(params)]
    paramNoise = params[lenght(params)]
    kXX <- kern(Xd,Xd,param) + kernNoise(Xd,Xd,paramNoise)
  }else{
    stop("the (log-)likelihood is only implemented for kernNoise=NULL and kernNoise=kWhite")    
  }
  ndata <- nrow(Xd)
  # rodo add something on the diagonal to avoid singularity
  detkXX <- det(kXX)
  eps <- 1.e-14
  tau <- 1.e-14
  while (detkXX<eps & tau<1) {
    tau <- tau*100
    kXX <- kXX + tau*diag(ndata)
    detkXX <- det(kXX)
  }
  LL <- -1/2*ndata*log(2*pi) - 1/2*log(det(kXX)) - 1/2*t(F)%*%solve(kXX)%*%F
  return(LL)
}
