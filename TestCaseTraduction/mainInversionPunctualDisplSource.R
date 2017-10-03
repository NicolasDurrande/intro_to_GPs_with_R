###########################
# MAIN CODE for
# Inversion of a punctual displacements source from 3D data
# Valerie Cayol, Rodolphe Le Riche, Nicolas Durrande
###########################

rm(list=ls()) #  cleaning up

library(R.matlab)
source("./mogi_3D.R")

####### data input #############
G = 2000 # Shear modulus in MPa
nu = 0.25 # Poisson's ratio
xs = 367000 # X location of source in UTM coordinates
ys = 7650300 # Y location of source in UTM
zs = 0 # Elevation of source with respect to sea level
a = 500 # source radius
p = 20 # Source overpressure in MPa
nlos = c(-0.664,-0.168,0.728) # cosines of land observing satellite

####### data input done #########

data <- readMat('dem.mat')
nrow <- nrow(data$xi)
ncol <- ncol(data$xi)
# for R plots, need to permute the order of the rows (increasing y as a vector)
# no need to permute the xi rows, they are all the same
data$yi <- data$yi[seq(nrow,1),]
data$zi <- data$zi[seq(nrow,1),]
xvec <- data$xi[1,] 
yvec <- data$yi[,1] 

# Compute surface displacements
U <- mogi_3D(G,nu,xs,ys,zs,a,p,data$xi,data$yi,data$zi)


########### PLOTS ################
# PLOTS on regular grid

# 3D rgl plot the landscape
library("rgl") # library for plots
open3d()
surface3d(xvec, yvec, t(data$zi), col= rainbow(16)[2])
# title3d(nameoffun, col="blue", font=4)
decorate3d()
# rgl.snapshot("./fileofplot.png", fmt="png", top=T)

# contour plots of displacements
par(mfrow=c(1,3))
image(xvec,yvec,t(U$x),xlab="x",ylab="y")
contour(xvec,yvec,t(U$x),add=TRUE,nlevels=20)
title("u")
#
image(xvec,yvec,t(U$y),xlab="x",ylab="y")
contour(xvec,yvec,t(U$y),add=TRUE,nlevels=20)
title("v")
#
image(xvec,yvec,t(U$z),xlab="x",ylab="y")
contour(xvec,yvec,t(U$z),add=TRUE,nlevels=20)
title("w")
