library(MASS)
library(rgl)

################################
## Generate samples from a 2 dimensional 
## multivariate normal distribution

n <- 1000                   # number of samples

# linear transform of iid N(0,1) varialbles
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
source('../labSessions/kernels.R')

# get data
library(R.matlab)
data <- readMat('../TestCase/data_nonoise.mat')[[1]]
X <- data[,1:2]
F <- data[,4]

xmin <- min(X[,1])
xmax <- max(X[,1])

ymin <- min(X[,2])
ymax <- max(X[,2])

xgrid <- seq(xmin, xmax, length=50)
ygrid <- seq(ymin, ymax, length=50)

Xgrid <- expand.grid(xgrid,ygrid)

K <- kMat32(Xgrid,Xgrid,param=c(1,4000,4000))

## generate samples by hand
C <- t(chol(K))
U <- matrix(rnorm(ncol(K)))
Z <- C %*% U

Z <- matrix(Z,ncol=50)
persp3d(xgrid,ygrid,Z,col='wheat')

## generate a sample using build-in functions
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')


##################################
## Multivariate GP posterior

plot3d()
