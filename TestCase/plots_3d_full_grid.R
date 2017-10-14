###########################
# Scripts for plotting displacements
# of a punctual displacements source from 3D data on a full grid of points
#
# The root file is 'fullgrid_xyz.csv', which is a cvs file.
# The measured displacements come from the 'data_nonoise.mat' matlab format file.
# Note that the number of rows and columns used to create (and lost ;-( ) 
# in the file are nrowdata <- 1255 and ncoldata <- 1159 and 
# the line of sight cosines are 
#   nlos = c(-0.664,-0.168,0.728) 
#
# Rodolphe Le Riche, Valerie Cayol, Nicolas Durrande
#
###########################

rm(list=ls()) #  cleaning up

############ load 3D full csv data file
# dimensions associated to the file "fullgrid_xyzulos.csv"
nrowdata <- 1255
ncoldata <- 1159
# ulos (the 4th column) was generated with
nlos = c(-0.664,-0.168,0.728) # vector of direction of line of sight (satellite)

#
datacsv <- read.csv(file="fullgrid_xyzulos.csv")
if (nrow(datacsv)!=nrowdata*ncoldata) stop("nrow(data file) not equal to nrowdata*ncoldata")
if (ncol(datacsv)!=4) stop("there should be 4 columns in data file")
data<-list()
data$xi <- matrix(datacsv[,1],nrow=nrowdata)
data$yi <- matrix(datacsv[,2],nrow=nrowdata)
data$zi <- matrix(datacsv[,3],nrow=nrowdata)
ulos <- matrix(datacsv[,4],nrow=nrowdata)
xvec <- data$xi[1,]
yvec <- data$yi[,1]


###### 3D rgl plot the landscape #############
library("rgl") # library for plots
# associate ulos to a matrix of colours
nbcol<-512
uloscol <- floor((ulos-min(ulos))/(max(ulos)-min(ulos))*(nbcol-1)+1)
# uloscol <- floor((data$zi-min(data$zi))/(max(data$zi)-min(data$zi))*(nbcol-1)+1) # colours representing elevation
colorlut <- rainbow(nbcol) # ulos color lookup table
zcol <- colorlut[t(uloscol)]
# par3d(cex=2.0)
surface3d(xvec, yvec, t(data$zi), color=zcol)
# title3d(nameoffun, col="blue", font=4)
# decorate3d()
# plot vector of line of sight
#   find highest point
highcoord <- arrayInd(ind=which.max(data$zi),.dim=dim(data$zi))
m1 <- c(data$xi[highcoord],data$yi[highcoord],data$zi[highcoord])
m2 <- m1 + nlos*5000
M = t(matrix(c(m1,m2),nrow=3))
# rgl.lines(x=M)
lines3d(x=M,col=c("black"),lwd=2)
text3d(((m1+m2)/2+1000*c(1,1,1)),texts = "LOS",cex=2)
title3d(main = "target u_los projected on terrain")
# rgl.snapshot("./fileofplot.png", fmt="png", top=T)



###### contour plot with added location of measures #####
## first load data related to measured locations
# where are the points used in the residu taken
library(R.matlab)
source("../labSessions/kernels.R")
data_m <- readMat('data_nonoise.mat') # TODO: do a read from cvs version (rodo)
meas_xi <- as.matrix(data_m$locdata[,1])
meas_yi <- as.matrix(data_m$locdata[,2])
meas_zi <- as.matrix(data_m$locdata[,3])
meas_ulos <- as.matrix(data_m$locdata[,4])
Xdata <- data_m$locdata[,1:2] # z's are not accounted for in Xdata
CXinv <- solve(kExp(Xdata,Xdata,c(5e-4,850,850)))

# png(filename="./contour.png")
par(mfrow=c(1,2))
image(xvec,yvec,t(ulos),xlab="x",ylab="y",col=zcol)
zlevels <- c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.05,0.1,0.12)
zlevcol <- rainbow(length(zlevels)-1)
contour(xvec,yvec,t(ulos),add=TRUE,levels=zlevels,col=zlevcol)
points(x=meas_xi,y=meas_yi)
wls_target <- t((meas_ulos-meas_ulos))%*%CXinv%*%(meas_ulos-meas_ulos) # kinda dumb yet doing it
title(paste("Target , WLS=",wls_target))



####### compare target with another set of variables
source(file = './mogi_3D.R')
newU <- mogi_3D(G= 2000,nu=0.25,xs=365000,ys=7649800,zs=-2000,a=500,p=300,datacsv[,1],datacsv[,2],datacsv[,3])
newulos <- nlos[1]*newU$x+nlos[2]*newU$y+nlos[3]*newU$z
newulos <- matrix(newulos,nrow=nrowdata)
Ucalc <- mogi_3D(G= 2000,nu=0.25,xs=365000,ys=7649800,zs=-2000,a=500,p=300,meas_xi,meas_yi,meas_zi)
uloscal <- nlos[1]*Ucalc$x+nlos[2]*Ucalc$y+nlos[3]*Ucalc$z
wls <- t((uloscal-meas_ulos))%*%CXinv%*%(uloscal-meas_ulos) 
image(xvec,yvec,t(newulos),xlab="x",ylab="y",col=zcol)
contour(xvec,yvec,t(newulos),add=TRUE,levels=zlevels,col=zlevcol)
points(x=meas_xi,y=meas_yi)
title(paste("Trial , WLS=",format(wls,digits=5,scientific = T)))
# dev.off()


# # calculate weighted least squares
# #   and first the cova
# wls <- t((newulos-datacvs[,4]))%*%Glb_CXinv%*%(newulos-datacvs[,4])


