
##########################################
##########################################
## helpers
##########################################
##########################################

#' Distance matrix between sets of points
#'
#' Computes the euclidian distance matrix, with rescaling factors.
#' 
#' @param x $n \times d$ matrix where $n$ is the number of points in x and d is the space dimension.
#' @param y $m \times d$ matrix where $m$ is the number of points in y and d is the space dimension.
#' @param theta $d$ dimensional vecor corresponding to the input space rescaling (defaults to 1 for each dimension)
#' 
#' @return a \code{n*m} matrix with entries corresponding to the euclidian distance between x[i,] and \code{y[j,]}. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3,0,0,0,0),ncol=2)
#' y <- matrix(c(0,1,2,3,0,1,1,1),ncol=2)
#' dist(x,y)
dist <- function(x,y,theta=NULL){
  if(ncol(x) != ncol(y)) stop("x and y must have the same number of columns") 
  if(is.null(theta)) theta <- rep(1,ncol(x))
  dist2 <- matrix(0,dim(x)[1],dim(y)[1])
  for(i in 1:dim(x)[2]){
    dist2 <- dist2 + (outer(x[,i],y[,i],"-")/theta[i])^2
  }
  return(sqrt(dist2))
}

##########################################
##########################################
## One dimensional kernels
##########################################
##########################################

#' Brownian motion kernel
#'
#' This function corresponds to the Brownian motion covariance function. It can also be used to compute covariance matrices between sets of points.
#' This kernel is only defined on 1D input spaces. 
#' 
#' @param x $n \times 1$ matrix where $n$ is the number of points in x (a vector of length \code{n} is also supported).
#' @param y $m \times 1$ matrix where $m$ is the number of points in y (a vector of length \code{m} is also supported).
#' @param param positive scalar corresponding to the "variance" parameter. Default value is set to 1. 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=1)
#' kBrown(x,x)
#' 
#' @export
kBrown <- function(x,y,param=1){
  if(ncol(x)!=1 | ncol(y)!=1 ) stop("this kernel is only defined for 1-dimensional inputs: 
                      x and y must be matrices with 1 column (or simply vectors)")
  param*outer(c(x),c(y),"pmin")
}


#' Cosine kernel
#'
#' This kernel will generate periodic samples. The input space must be 1-dimensional.
#' 
#' @param x $n \times 1$ matrix where $n$ is the number of points in x (a vector of length \code{n} is also supported).
#' @param y $m \times 1$ matrix where $m$ is the number of points in y (a vector of length \code{m} is also supported).
#' @param param vector of length 2: \code{param[1]} is the process variance and \code{param[2]} is the period. Default value is set to \code{(1,2*pi)}. 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=1)
#' kCos(x,x)
#' 
#' @export
kCos <- function(x,y,param=NULL){
  if(ncol(x)!=1) stop("this kernel is only defined for 1-dimensional inputs: 
                      x and y must be matrices with 1 column")
  if(is.null(param)) param <- c(1,2*pi)
  param[1]*cos(outer(c(x),c(y),'-')/param[2]*2*pi)
}


##########################################
##########################################
## d-dimensional kernels
##########################################
##########################################

#' Exponential kernel (i.e. Matern 1/2)
#'
#' This corresponds to the Ornstein-Uhlenbek covariance function. It can also be used to compute covariance matrices between sets of points.
#' 
#' @param \code{x} \code{n * d} matrix where \code{n} is the number of points in \code{x}, and \code{d} is the input space dimension.
#' @param \code{y} \code{m * d} matrix where \code{m} is the number of points in \code{y}, and \code{d} is the input space dimension.
#' @param param \code{(d+1)}-dimensional vector of positive scalars. \code{param[1]} is the process variance and \code{param[-1]} are the length-scale parameters. Default value is \code{(1,0.2,...,0.2)} (reasonable values for [0,1]^d input spaces). 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=2)
#' kExp(x,x)
#' 
#' @export
kExp <- function(x,y,param=NULL){
  if(is.null(param)) param <- c(1,rep(.2,ncol(x)))
  param[1]*exp(-dist(x,y,param[-1]))
}


#' Gaussian kernel (i.e. RBF, squared exponential, exponentiated quadratic)
#'
#' This covariance function can also be used to compute a covariance matrice between two sets of points.
#' 
#' @param \code{x} \code{n * d} matrix where \code{n} is the number of points in \code{x}, and \code{d} is the input space dimension.
#' @param \code{y} \code{m * d} matrix where \code{m} is the number of points in \code{y}, and \code{d} is the input space dimension.
#' @param param \code{(d+1)}-dimensional vector of positive scalars. \code{param[1]} is the process variance and \code{param[-1]} are the length-scale parameters. Default value is \code{(1,0.2,...,0.2)} (reasonable values for [0,1]^d input spaces). 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=2)
#' kGauss(x,x)
#' 
#' @export
kGauss <- function(x,y,param=NULL){
  if(is.null(param)) param <- c(1,rep(.2,ncol(x)))
  param[1]*exp(-.5*dist(x,y,param[-1])^2)
}

#' Matern 3/2 kernel 
#'
#' This covariance function can also be used to compute a covariance matrice between two sets of points.
#' 
#' @param \code{x} \code{n * d} matrix where \code{n} is the number of points in \code{x}, and \code{d} is the input space dimension.
#' @param \code{y} \code{m * d} matrix where \code{m} is the number of points in \code{y}, and \code{d} is the input space dimension.
#' @param param \code{(d+1)}-dimensional vector of positive scalars. \code{param[1]} is the process variance and \code{param[-1]} are the length-scale parameters. Default value is \code{(1,0.2,...,0.2)} (reasonable values for [0,1]^d input spaces). 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=2)
#' kMat32(x,x)
#' 
#' @export
kMat32 <- function(x,y,param=NULL){
  if(is.null(param)) param <- c(1,rep(.2,ncol(x)))
  d <- sqrt(3)*dist(x,y,param[-1])
  return(param[1]*(1 + d)*exp(-d))
}

#' Matern 5/2 kernel 
#'
#' This covariance function can also be used to compute a covariance matrice between two sets of points.
#' 
#' @param \code{x} \code{n * d} matrix where \code{n} is the number of points in \code{x}, and \code{d} is the input space dimension.
#' @param \code{y} \code{m * d} matrix where \code{m} is the number of points in \code{y}, and \code{d} is the input space dimension.
#' @param param \code{(d+1)}-dimensional vector of positive scalars. \code{param[1]} is the process variance and \code{param[-1]} are the length-scale parameters. Default value is \code{(1,0.2,...,0.2)} (reasonable values for [0,1]^d input spaces). 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=2)
#' kMat52(x,x)
#' 
#' @export
kMat52 <- function(x,y,param=NULL){
  if(is.null(param)) param <- c(1,rep(.2,ncol(x)))
  d <- sqrt(5)*dist(x,y,param[-1])
  return(param[1]*(1 + d +1/3*d^2)*exp(-d))
}

#' White noise kernel 
#'
#' This covariance function can also be used to compute a covariance matrice between two sets of points.
#' 
#' @param \code{x} \code{n * d} matrix where \code{n} is the number of points in \code{x}, and \code{d} is the input space dimension.
#' @param \code{y} \code{m * d} matrix where \code{m} is the number of points in \code{y}, and \code{d} is the input space dimension.
#' @param param a positive scalar corresponding to the process variance. 
#' 
#' @return a \code{n*m} covariance matrix. 
#' 
#' @examples
#' x <- matrix(c(0,1,2,3),ncol=2)
#' kWhite(x,x)
#' 
#' @export
kWhite <- function(x,y,param=NULL){
  if(is.null(param)) param <- 1
  d <- dist(x,y,rep(1,dim(x)[2]))
  return(param*(d==0))
}
