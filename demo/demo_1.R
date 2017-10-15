library(MASS)
library(rgl)

source('../labSessions/kernels.R')
source('../labSessions/models.R')
source('../labSessions/likelihood.R')

################################
## Generate samples from a 2 dimensional 
## multivariate normal distribution

n <- 1000                   # number of samples

# linear transform of iid N(0,1) variables
X <- matrix(rnorm(2*n),ncol=2)
plot(X)

A <- matrix(c(2,5,6,3), ncol=2)
Y <- X %*% t(A)
plot(Y)

# Sampling from a N(m,K) distribution
X <- matrix(rnorm(2*n),ncol=2)

m <- c(3,2)
K <- matrix(c(4,2,2,2), nrow=2)

C <- t(chol(K))
Y <- t(C %*% t(X)) + matrix(m, nrow=n, ncol=2, byrow=TRUE)

plot(Y)
colMeans(Y)
var(Y)



#################################
## Sampling from a multivariate GP

# get data
library(R.matlab)
data <- readMat('../TestCase/data_nonoise.mat')[[1]]
X <- data[,1:2]
F <- data[,4]

xmin <- min(X[,1])
xmax <- max(X[,1])

ymin <- min(X[,2])
ymax <- max(X[,2])

xgrid <- seq(xmin, xmax, length=30)
ygrid <- seq(ymin, ymax, length=30)

Xgrid <- expand.grid(xgrid,ygrid)

K <- kMat52(Xgrid,Xgrid,param=c(1,4000,4000))

## generate samples by hand
C <- t(chol(K))
U <- matrix(rnorm(ncol(K)))
Z <- C %*% U

Z <- matrix(Z, nrow=length(xgrid))
persp3d(xgrid, ygrid, xlab='x1', ylab='x2', Z, col='wheat')

## generate a sample using build-in functions
Z <- mvrnorm(1, rep(0,nrow(K)), K)
persp3d(xgrid, ygrid, Z, col='wheat')


##################################
## Multivariate GPR

## look at the data
plot3d(X[,1], X[,2], F)

## what prior seem appropriate ?
K <- kExp(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

K <- kMat32(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

K <- kMat52(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

K <- kGauss(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

## build a GPR model
pred <- predGPR(Xgrid, X, F,
                kern=kGauss, param=c(1e-3,3000,3000))

pred_sd <- pmax(rep(0,ncol(pred$cov)), diag(pred$cov))
pred_upper95 <- pred$mean + 2*sqrt(pred_sd)
pred_lower95 <- pred$mean - 2*sqrt(pred_sd)

persp3d(xgrid, ygrid, pred$mean, col='red')
surface3d(xgrid, ygrid, pred_lower95, col='wheat', alpha=.8)
surface3d(xgrid, ygrid, pred_upper95, col='wheat', alpha=.8)
spheres3d(X[,1], X[,2], F,radius=200)


####################################
## Parameter estimation

## plot log likelihood
sig2_grid <- seq(1e-4, 1e-3, length=30)
theta_grid <- seq(500, 3000, length=30)

PARAM <- as.matrix(expand.grid(sig2_grid, theta_grid))
LL <- rep(0, nrow(PARAM))

for(i in 1:nrow(PARAM)){
    LL[i] <- logLikelihood(c(PARAM[i,],PARAM[i,2]),kGauss,X,F)
}

persp3d(sig2_grid, theta_grid, LL, col='red')

## optimisation
indmax <- which.max(LL)
paropt <- c(PARAM[indmax,], PARAM[indmax,2])
spheres3d(paropt[1], paropt[2], LL[indmax], radius=20)

## model associated with optimal parameters
pred <- predGPR(Xgrid, X, F,
                kern=kGauss, param=paropt)

pred_sd <- pmax(rep(0,ncol(pred$cov)), diag(pred$cov))
pred_upper95 <- pred$mean + 2*sqrt(pred_sd)
pred_lower95 <- pred$mean - 2*sqrt(pred_sd)

persp3d(xgrid, ygrid, pred$mean, col='red')
surface3d(xgrid, ygrid, pred_lower95, col='wheat', alpha=.8)
surface3d(xgrid, ygrid, pred_upper95, col='wheat', alpha=.8)
spheres3d(X[,1], X[,2], F, radius=200)

##################################
## Multivariate GPR with noise

Kn <- kExp(X, X, param=c(1e-4, 2000, 2000)) + kWhite(X, X, param=1e-5)
Fn <- mvrnorm(1,F,Kn)

## look at the data
plot3d(X[,1], X[,2], Fn)

## construct appropriate kernel for the signal
kSignal <- function(x, y, param){
    sig2 <- abs(param[1])
    theta <- rep(param[2], 2)
    decay <- rep(param[3], 2)

    loc <- matrix(c(366000, 7650000), ncol=2)
    locx <-exp(-.5*dist(x, loc, decay)^2)
    locy <-exp(-.5*dist(y, loc, decay)^2)
    Locx_mat <- matrix(locx, nrow=nrow(x), ncol=nrow(y))
    Locy_mat <- t(matrix(locy, nrow=nrow(y), ncol=nrow(x)))

    kgauss <- exp(-.5*dist(x, y, theta)^2)
 
    kern <- sig2 * Locx_mat * Locy_mat * kgauss

    return(kern)
}

KS <- kSignal(Xgrid, Xgrid, c(1e-3, 4000, 3000))
Z <- mvrnorm(1, rep(0,nrow(KS)), KS)
persp3d(xgrid, ygrid, Z, col='wheat')

## construct appropriate kernel for the noise
kNoise <- function(x, y, param){
    kW <- kWhite(x, y, param[1])
    kE <- kExp(x, y, param[c(2,3,3)])
    return(kW + kE)
}

KN <- kNoise(Xgrid, Xgrid, c(1e-5, 1e-4, 3000))
ZN <- mvrnorm(1, rep(0,nrow(KN)), KN)
persp3d(xgrid, ygrid, ZN, col='wheat')

## Optimize the log-likelihood
params <- c(1e-3, 4000, 3000, 1e-5, 1e-4, 3000)
logLikelihood(params, kSignal, X, Fn, kNoise, num_param_noise=3)

opt <- optim(params, logLikelihood, control=list(fnscale=-1), 
             kern=kSignal, Xd=X, F=Fn, kernNoise=kNoise, num_param_noise=3)

## build a GPR model
pred <- predGPR(Xgrid, X, Fn,
                kern=kSignal, param=opt$par[1:3],
                kernNoise=kNoise, paramNoise=opt$par[4:6])

pred_sd <- pmax(rep(0,ncol(pred$cov)), diag(pred$cov))
pred_upper95 <- pred$mean + 2*sqrt(pred_sd)
pred_lower95 <- pred$mean - 2*sqrt(pred_sd)

persp3d(xgrid, ygrid, pred$mean, col='red')
surface3d(xgrid, ygrid, pred_lower95, col='wheat', alpha=.5)
surface3d(xgrid, ygrid, pred_upper95, col='wheat', alpha=.5)
spheres3d(X[,1], X[,2], Fn, radius=200)
