source("kernels.R")
source("models.R")
source("likelihood.R")
source("plots.R")
source("EI.R")

# Generate simple Data
X <- matrix(seq(from=0.1, to=0.9, length=5))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

# predict on a grid with p points
p <- 101
x <- matrix(seq(from=0, to=1, length=p))
pred <- predGPR(x,X,F,kMat52)
plotGPR(x,pred)

# optimize the model parameters
logLikelihood(c(1,.2), kern=kMat52, Xd=X, F=F)
opt_out <- optim(c(1,.2), logLikelihood, kern=kMat52, Xd=X, F=F, control=list(fnscale=-1))
param_opt <- opt_out$par

pred <- predGPR(x, X, F, kMat52, param_opt)
plotGPR(x,pred)

# Compute the expected improvement
ei <- EI(xp=x, Xd=X, F, kMat52, param_opt)
plot(x, ei, type='l', main="Expected Improvement", col=darkRed)
# maximize EI
oEI <- maxEI(kern=kMat52,Xd=X,F=F,param=param_opt,xmin=0,xmax=1,nbtry=100,maxit=100,silent=TRUE)
text(x = 0.4, y = 0.04, paste("arg max EI=",format(oEI$var,digits = 5),
                              "\n max EI=",format(oEI$EI,digits=5)), pos=4, cex = 1.2, col = "black")
