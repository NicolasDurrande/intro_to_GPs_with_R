#' Log-likelihood of a multivariate Gaussian vector
#'
#' Log-likelihood of some parametres of a GP Given some observations
#' 
#' @param vector corresponding to the concatenation of the parameters of the kernel and the noise kernel.
#' @param kern covariance function of the GP prior for the signal. 
#' @param kernNoise covariance function of the GP prior for the noise. 
#' @param X $n \times d$ matrix of the design of experiments (\code{n} is the number of points and \code{d} is the input space dimension).
#' @param F $n \times 1$ matrix corresponding to the values of \code{f} at \code{X}.
#' 
#' @return a list with the mean (\code{n*1}-matrix) and the covariance matrix (\code{n*n}) of the posterior distribution at \code{x}. 
#' 
#' @examples
#' 
#' # Generate simple Data
#' X <- matrix(seq(from=0.1, to=0.9, length=5))
#' F <- matrix(c(0.5, 0, 1.5, 3, 2))
#' 
#' # optimize the model parameters
#' logLikelihood(c(1,.2), kern=kMat52, X=X, F=F)
#' opt_out <- optim(c(1,.2), logLikelihood, kern=kMat52, X=X, F=F)
#' 
#' @export
logLikelihood <- function(params,kern,X,F,kernNoise=NULL){
  if(is.null(kernNoise)){
    kXX <- kern(X,X,params) 
  }else if(kernNoise==kWhite){
    param = params[-lenght(params)]
    paramNoise = params[lenght(params)]
    kXX <- kern(X,X,param) + kernNoise(X,X,paramNoise)
  }else{
    stop("the (log-)likelihood is only implemented for kernNoise=NULL and kernNoise=kWhite")    
  }
  LL <- -1/2*nrow(X)*log(2*pi) - 1/2*log(det(kXX)) - 1/2*t(F)%*%solve(kXX)%*%F
  return(LL)
}
