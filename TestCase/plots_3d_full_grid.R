###########################
# Scripts for plotting displacements
# of a punctual displacements source from 3D data on a full grid of points
#
# The root file is 'fullgrid_xyz.csv', which is a cvs file.
# Note that the number of rows and columns used to create (and lost ;-( ) 
# in the file are nrowdata <- 1255 and ncoldata <- 1159 and 
# the line of sight cosines are 
#   nlos = c(-0.664,-0.168,0.728) 
#
# Rodolphe Le Riche, Valerie Cayol, Nicolas Durrande
#
###########################

rm(list=ls()) #  cleaning up

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

########### PLOTS ################
# PLOTS on regular grid

# 3D rgl plot the landscape
library("rgl") # library for plots
#open3d()
# associate ulos to a matrix of colours
nbcol<-512
uloscol <- floor((ulos-min(ulos))/(max(ulos)-min(ulos))*(nbcol-1)+1)
# uloscol <- floor((data$zi-min(data$zi))/(max(data$zi)-min(data$zi))*(nbcol-1)+1) # colours representing elevation
colorlut <- rainbow(nbcol) # ulos color lookup table
zcol <- colorlut[t(uloscol)]
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
rgl.snapshot("./fileofplot.png", fmt="png", top=T)

# contour plots of displacements
# par(mfrow=c(1,1))
# png(filename="./contour.png")
image(xvec,yvec,t(ulos),xlab="x",ylab="y")
contour(xvec,yvec,t(ulos),add=TRUE,nlevels=20)
title("Displ. projected on LOS")
# dev.off()