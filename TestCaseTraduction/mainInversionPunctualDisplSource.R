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

####### data input done #########

data <- readMat('dem.mat')
nrow <- 1255
ncol <- 1159


# Compute surface displacements
U <- mogi_3D(G,nu,xs,ys,zs,a,p,data$xi,data$yi,data$zi)

# contour plot on a regular grid
xvec <- data$xi[seq(0,ncol-1)*nrow+1]
yvec <- data$yi[seq(nrow,1)]
uxmat <- matrix(data=U$x,nrow=nrow,ncol=ncol)
uxmat <- uxmat[nrow(uxmat):1,]
image(xvec,yvec,t(uxmat))
contour(xvec,yvec,t(uxmat),add=TRUE,nlevels=20)
#
uzmat <- matrix(data=U$z,nrow=nrow,ncol=ncol)
uzmat <- uzmat[nrow(uxmat):1,]
image(xvec,yvec,t(uzmat))
contour(xvec,yvec,t(uzmat),add=TRUE,nlevels=20)

