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

####### load utilities ##########################

rm(list=ls()) #  cleaning up

library(R.matlab)
source("./mogi_3D.R")
source("./wls_ulos.R")
source("../labSessions/kernels.R")
source("../labSessions/likelihood.R")
source("../labSessions/models.R")
source("../labSessions/plots.R")
source("../labSessions/EI.R")


###### input for variables identification ###########
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

# bounds on variables
xmax <- NA
xmin <- NA
xmin[n2i$xs]<-364000
xmax[n2i$xs]<-366000
xmin[n2i$ys]<-7649000
xmax[n2i$ys]<-7651000
xmin[n2i$zs]<- -3000
xmax[n2i$zs]<-1000
xmin[n2i$a]<- 50
xmax[n2i$a]<- 1000
xmin[n2i$p]<- -500
xmax[n2i$p]<-500

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
Glb_CXinv <<- solve(kExp(Xdata,Xdata,c(5e-4,850,850))) # calculated once for all, used in wls_ulos
# # has been compared to the matlab covariance matrix below, max difference 1e-6
# data <- readMat('cov_nonoise.mat')
# CovData <- data$Cd
rm(data)
rm(Xdata)


#######  design of experiments #############################

library(lhs)
nbinit <- 100 # number of points in the initial design of experiments
set.seed(0)
# do an optimal Latin Hypercube Sampling
Xnorm <- optimumLHS(nbinit, nbvar)
# a function to denormalize variables
# Xnorm should either be a (n times nbvar) matrix or a numeric vector 
unnorm_var <- function(Xnorm){
  if (is.null(dim(Xnorm))) {Xn <- matrix(data = Xnorm, ncol = nbvar, byrow = T)} # numeric vector
  else {Xn <- Xnorm}
  nbrep <- dim(Xn)[1]
  Xu <- matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar) + 
    Xn * matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)  
  Xu <- data.frame(Xu)
  names(Xu) <- varnames
  return(Xu)
}
#
X <- unnorm_var(Xnorm=Xnorm)
pairs(X) # look at it

# calculate weighted least squares at X of U projected on LOS (w.r.t. target Glb_ulos)
wls <- apply(X, 1, wls_ulos)

# normalize the output so that it is centered with a unit std dev
# because wls ranges from 0 to 10^9, do a scaling in log(1+wls)
lwls <- log(1+wls)
mean_wls <- mean(lwls)
std_wls <- sd(lwls)
# wls normalization as a function as it is used at different locations in the code
normalizeWLS <- function(awls){
  lawls <- log(1+awls)
  return((lawls - mean_wls)/std_wls)
}
#
norm_wls <- normalizeWLS(wls)

# plot learning data
par(mfrow=c(2,3))
for (i in 1:nbvar){
  plot(X[,i],norm_wls,xlab=names(X)[i])
  # plot(X[,i],log(wls),xlab=names(X)[i])
  }
# empty plot + text
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("LEARNING\n","SET"), cex = 1.6, col = "black")

# Comments:
# * a, the source radius, is a sensitive variable, 
# * large zs seem correlated to large wls error. 

###### build a kriging model #######################

# optimize the model parameters by repeating local searches started from random initial points

nbtry <- 0 # make it <1 to skip max LL and take a past, a priori theta
if (nbtry<1){
  # Past results with non-logged outputs
  # oLL$bestthetas <- c(2.0680115 , 7.8302747 ,11.0337919 , 0.1464954 , 0.2639607 , 0.3328377) # LL=-99.12137
  # i.e., only zs and a are sensitive variables
  # past results with log(1+wls)
  oLL <- list()
  oLL$bestthetas <- c(1,2,2,0.3,0.3,0.5) # has LL=-123.8744 but works usually well
  # final thetas= 5.122735 6.87119 7.787931 5.655541 0.7671208 0.04670216  , final LL= -96.75596   
  oLL$bestLL <- logLikelihood(params = oLL$bestthetas,kern=kMat52,Xd=Xnorm,F=norm_wls)
} else {
  oLL <- maxlogLikelihood(kern=kMat52,Xd=Xnorm,F=norm_wls,nbtry=nbtry,maxit=500,silent=F)
}
cat("\n final thetas=",oLL$bestthetas," , final LL=",oLL$bestLL,"\n")

####### test the kriging model
# make a test set
ntest <- 110
Xtestnorm <- matrix(runif(ntest*nbvar),nrow=ntest)
Xtest <- matrix(rep(xmin,times=ntest),byrow = T,ncol=nbvar) + 
  Xtestnorm * matrix(rep((xmax-xmin),times=ntest),byrow = T,ncol=nbvar)
Xtest <- data.frame(Xtest)
names(Xtest) <- varnames
wls_test <- apply(Xtest, 1, wls_ulos)
lwls_test <- log(1+wls_test)
wls_testnorm <- (lwls_test - mean_wls)/std_wls
# wls_testnorm <- (wls_test - mean_wls)/std_wls # linear scaling version

# predict with the model and measure quality of predictions wrt true response 
pred <- predGPR(x=Xtestnorm, X=Xnorm, F=norm_wls, kern=kMat52, param=oLL$bestthetas)
# calculate RMSE and Q2
rmse <- sqrt(mean((wls_testnorm-pred$mean)^2))
q2 <- 1 - sum((wls_testnorm-pred$mean)^2)/sum((wls_testnorm-mean(wls_testnorm))^2)

# plot the test set
par(mfrow=c(2,3))
for (i in 1:nbvar){
  plot(Xtest[,i],log(wls_test),xlab=names(Xtest)[i])
}
# 
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("TEST SET\n","RMSE=",format(rmse,digits=4),
                             "\n Q2=",format(q2,digits =4), sep=""), cex = 1.2, col = "black")

par(mfrow=c(1,3))
# scatter plot
eps <- 0.2
limmin <- min(c(wls_testnorm,pred$mean))-eps
limmax <- max(c(wls_testnorm,pred$mean))+eps
par(pty="s")
plot(wls_testnorm, pred$mean, xlim=c(limmin,limmax), ylim=c(limmin,limmax), asp=1)
lines(c(limmin,limmax),c(limmin,limmax))
# residuals
plot(x=1:length(wls_testnorm),y=(pred$mean-wls_testnorm),xlab="test points no.")
# pred +/- std dev
plot(x=1:length(pred$mean),y=pred$mean,xlab="test points no.",ylab = "WLS",pch=3)
sdtest <- sqrt(diag(pred$cov))
points(x=1:length(pred$mean),y=(pred$mean-sdtest),pch="-")
points(x=1:length(pred$mean),y=(pred$mean+sdtest),pch="-")
points(x=1:length(wls_testnorm),y=wls_testnorm,pch=1,col="blue")
legend(x = "topright",legend = c("test","pred +/- std"),pch = c(1,3),col = c("blue","black"))
# TODO: other plot result with normal distribution?

# Comments : 
#   the normalization log(1+wls) helps a lot the kriging model

######### model identification ########
#         = EGO algorithm to optimize (minimize) the model-target WLS distance

EGOmaxiter <- 50
period_upd <- 5

cat("\n******** START EGO \n\n",file = "ego_listen.txt")
# x11()
# hplot_wlsn <- dev.cur()
par(mfrow=c(1,1))
uxlim <- dim(Xnorm)[1]+EGOmaxiter
ytarget <- -mean_wls/std_wls
yrange <- max(norm_wls)-min(norm_wls)
plot(x = 1:dim(Xnorm)[1],y=norm_wls,xlab="point number",ylab="norm. WLS",xlim=c(1,uxlim),
     ylim=c(ytarget-0.2*yrange,max(norm_wls)+0.2*yrange),type="l")
lines(x=c(1,uxlim),y=c(ytarget,ytarget),lty="dotted",col="red")
text(x=1,y=ytarget+0.1,labels = "ideal",cex = 0.8,col = "red",pos=4)

for (iter in 1:EGOmaxiter){

  cat("\n***** EGO iteration ",iter,"\n\n")
  # optimise EI
  oEI <- maxEI(kern=kMat52,Xd=Xnorm,F=norm_wls,param=oLL$bestthetas,xmin=0,xmax=1,nbtry=100,maxit=100,silent=F) 

  # calculate function at new point
  newX <- unnorm_var(Xnorm = oEI$var)
  newwls <- wls_ulos(as.numeric(newX))  
  newwls_norm <- normalizeWLS(newwls)
  # update data bases and save them
  X <- rbind(X,newX)
  Xnorm <- rbind(Xnorm,matrix(oEI$var,ncol=nbvar,byrow = T))
  wls <- c(wls,newwls)
  norm_wls <- c(norm_wls,newwls_norm)
  save(Xnorm,X,norm_wls,wls,file="ego_data_gathered.RData")

  # some printing and user control
  cat("\n Iteration ",iter," summary:\n")
  cat("   EI best var=",oEI$var," , has EI=",oEI$EI,"\n")
  cat("   unnormed var= ")
  for (ii in 1:nbvar) {cat(as.numeric(newX[ii]),varnames[ii],"  ")}
  cat(" , wls = ",newwls,"  (norm.wls=",newwls_norm,")\n")
  #
  cat("\n Iteration ",iter," summary:\n",file="ego_listen.txt",append = TRUE)
  cat("   EI best var=",oEI$var," , has EI=",oEI$EI,"\n",file="ego_listen.txt",append = TRUE)
  cat("   unnormed var= ",file="ego_listen.txt",append = TRUE)
  for (ii in 1:nbvar) {cat(as.numeric(newX[ii]),varnames[ii],"  ",file="ego_listen.txt",append = TRUE)}
  cat(" , wls = ",newwls,"  (norm.wls=",newwls_norm,")\n",file="ego_listen.txt",append = TRUE)
  #
  newn <- dim(Xnorm)[1]
  lines(x=c(newn-1,newn),y=c(norm_wls[newn-1],norm_wls[newn]),col="blue")
  
  # update the kriging model every period_upd
  if (iter%%period_upd==0) {
    tmax <- c(5*var(norm_wls),rep(3,times=nbvar))
    tmin <- c(0.1,rep(0.08,times=nbvar))
    oLL <- maxlogLikelihood(kern=kMat52,Xd=Xnorm,F=norm_wls,tmin=tmin,tmax=tmax,nbtry=50,maxit=500,silent=F)
    cat("   max LL, new thetas=",oLL$bestthetas," with LL=",oLL$bestLL,"\n")
    cat("   max LL, new thetas=",oLL$bestthetas," with LL=",oLL$bestLL,"\n",file="ego_listen.txt",append = TRUE)
    # user.stop <- readline(prompt="Enter x to exit, anything else otherwise ")
    # if (user.stop=="x") break
  }
  
} # end EGO loop

####### print out results
# where are the new points
cols <- character(nrow(X))
cols[] <- "black"
cols[nbinit+1:dim(X)[1]] <- "blue"
pairs(X,col=cols)
# where are the q percent good points
qwlsn <- quantile(x=norm_wls,probs = 0.2)
igood <- which(norm_wls[]<qwlsn)
pairs(X[igood,])
# look at non invertibility
# all optimal solutions have  a^3*p = 500^3*20 
plot(X[igood,n2i$a],X[igood,n2i$p],xlab="a",ylab="p")
aideal <- seq(from=xmin[n2i$a],to=xmax[n2i$a],length.out = 100)
pideal <- (xstar[n2i$a]^3*xstar[n2i$p])/(aideal^3)
lines(x = aideal,y=pideal,lty="solid",col="red")
text(x = 400,y=400,labels = "a^3 * p = const.",col = "red")
legend("topright",legend = c("20% best solutions","ideal solutions"),pch = c(1,NA),lty=c(NA,1),col=c("black","red"))
