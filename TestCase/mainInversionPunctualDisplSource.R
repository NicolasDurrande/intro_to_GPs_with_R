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

####### input for variables identification ###########
n2i <- list(xs=1, ys=2, zs=3, a=4, p=5) # name to index for variables
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
