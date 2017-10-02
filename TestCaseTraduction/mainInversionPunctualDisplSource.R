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

# Compute surface displacements
U <- mogi_3D(G,nu,xs,ys,zs,a,p,xi,yi,zi)
