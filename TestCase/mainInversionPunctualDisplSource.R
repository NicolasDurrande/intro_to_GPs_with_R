###########################
# MAIN CODE for
# Inversion of a punctual displacements source from 3D data
# The data used are under-sampled with the quadtree method (irregular grid)
#
# Rodolphe Le Riche, Nicolas Durrande, Valerie Cayol
#
#
# Notes: 
# * use of global variables, starting with Glb_...
###########################

rm(list=ls()) #  cleaning up

library(R.matlab)
source("./mogi_3D.R")
source("./wls_ulos.R")
source("../labSessions/kernels.R")
source("../labSessions/likelihood.R")

####### input for variables identification ###########
n2i <- list(xs=1, ys=2, zs=3, a=4, p=5) # name to index for variables
varnames <- c("xs","ys","zs","a","p")
nbvar <- 5
# optimum
xstar <- NA
xstar[n2i$xs] <- 367000 # X location of source in UTM coordinates
xstar[n2i$ys] <- 7650300 # Y location of source in UTM
xstar[n2i$zs] <- 0 # Elevation of source with respect to sea level
xstar[n2i$a] <- 500 # source radius
xstar[n2i$p] <- 20 # Source overpressure in MPa
# order of magnitude of the variables 
# (different from 0, to be used when scaling, setting bounds...)
xmag <- NA
xmag[n2i$xs] <- 367000 # X location of source in UTM coordinates
xmag[n2i$ys] <- 7650300 # Y location of source in UTM
xmag[n2i$zs] <- 1000 # Elevation of source with respect to sea level
xmag[n2i$a] <- 500 # source radius
xmag[n2i$p] <- 100 # Source overpressure in MPa
# initial point (for an identification)
xinit<-xstar+c(0.2,-0.2,-0.3,0.5,-0.4)*xmag
# bounds on variables
xmax <- NA
xmin <- NA
xmin[n2i$xs]<-364000
xmax[n2i$xs]<-366000
xmin[n2i$ys]<-7649000
xmax[n2i$ys]<-7651000
xmin[n2i$zs]<- -3000
xmax[n2i$zs]<-2000
xmin[n2i$a]<- 50
xmax[n2i$a]<- 2000
xmin[n2i$p]<- -500
xmax[n2i$p]<-500
#


Glb_var <<- list(n2i=n2i,nbvar=nbvar,xmag=xmag,xmax=xmax,xmin=xmin) # always useful stuff

####### load data ###########
data <- readMat('data_nonoise.mat') # TODO: do a read from cvs version (rodo)
Glb_xi <<- as.matrix(data$locdata[,1])
Glb_yi <<- as.matrix(data$locdata[,2])
Glb_zi <<- as.matrix(data$locdata[,3])
Glb_ulos <<- as.matrix(data$locdata[,4])

# calculate data Covariance matrix, store it in a Global variable
# covariance from exponential kernel, var = 5e-4m2, cor_length = 850 m
# and invert it
Xdata <- data$locdata[,1:2] # z's are not accounted for in Xdata
Glb_CXinv <<- solve(kExp(Xdata,Xdata,c(5e-4,850,850))) 
# # has been compared to the matlab covariance matrix below, max difference 1e-6
# data <- readMat('cov_nonoise.mat')
# CovData <- data$Cd
rm(data)
rm(Xdata)

# a sample call to the distance function would look like this:
# zvar <- list(n2i=n2i, x=xinit)
# wls <- wls_ulos(zvar$x)

######  do a design of experiments #############################

library(lhs)
nbinit <- 100 # number of points in the initial design of experiments
set.seed(0)
# do an optimal Latin Hypercube Sampling
X <- matrix(rep(xmin,times=nbinit),byrow = T,ncol=nbvar) + 
  optimumLHS(nbinit, nbvar) * matrix(rep((xmax-xmin),times=nbinit),byrow = T,ncol=nbvar)
X <- data.frame(X)
names(X) <- varnames
wls <- apply(X, 1, wls_ulos)
# a bit of plotting
par(mfrow=c(2,3))
for (i in 1:nbvar){
  plot(X[,i],log(wls),xlab=names(X)[i])
}
# The range of wls is quite large (between 0 and 1e^9), therefore I plot in log scale
# Observe that a, the source radius, is a sensitive variable. 

###### build a kriging model #######################

# normalize the output so that it is centered with a unit std dev
# (wls ranges from 0 to 10^9, might have to do a more radical scaling like log(1+wls))
mean_wls <- mean(wls)
std_wls <- sd(wls)
norm_wls <- (wls - mean_wls)/std_wls

# normalize the input so that all variables are between 0 and 1

# optimize the model parameters

# try random parameters
# kmin <- rep(0.01,times=(nbvar+1))
# kmax <- c(100,xmax)
# nbtry <- 10
# ks <- matrix(rep(kmin,times=nbtry),byrow = T,ncol=(nbvar+1)) + 
#   matrix(runif(n=(nbtry*(nbvar+1))),nrow=nbtry) * matrix(rep((kmax-kmin),times=nbtry),byrow = T,ncol=(nbvar+1))
# allLL <- apply(X=ks,MARGIN=1,FUN = logLikelihood, kern=kMat52, Xd=X, F=norm_wls)

# source(file="./cmaes.R")
# pcma <- list()
# pcma$xinit <- c(1,rep(1,times=nbvar))
# pcma$LB <- rep(0.01,times=nbvar+1)
# pcma$UB <- c(1000,xmax)
# pcma$budget <- 100000
# pcma$sigma <- 1
# opt_out <- cmaes(test_fun=logLikelihood, param=pcma,kern=kMat52, X=X, F=norm_wls)

# opt_out <- optim(c(1,rep(1,times=nbvar)), fn = logLikelihood, kern=kMat52, Xd=X, F=norm_wls, method = "Nelder-Mead",
#                  control=list(fnscale=-1, parscale=c(1,xmag)))
# opt_out <- optim(c(1,rep(20,times=nbvar)), fn = logLikelihood, kern=kMat52, Xd=X, F=norm_wls, method = "L-BFGS-B",
#                  control=list(fnscale=-1, parscale=c(1,xmag), maxit=500,factr=1))
opt_out <- optim(c(1,rep(0.2,times=nbvar)), fn = logLikelihood, kern=kMat52, Xd=X, F=norm_wls, control=list(fnscale=-1))
# param_opt <- opt_out$par


