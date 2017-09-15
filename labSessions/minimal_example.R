source("kernels.R")
source("models.R")
source("likelihood.R")
source("plots.R")

# Generate simple Data
X <- matrix(seq(from=0.1, to=0.9, length=5))
F <- matrix(c(0.5, 0, 1.5, 3, 2))

# predict on a grid with p points
p <- 101
x <- matrix(seq(from=0, to=1, length=p))
pred <- predGPR(x,X,F,kMat52)
plotGPR(x,pred)

# optimize the model parameters
logLikelihood(c(1,.2), kern=kMat52, X=X, F=F)
opt_out <- optim(c(1,.2), logLikelihood, kern=kMat52, X=X, F=F)

pred <- predGPR(x, X, F, kMat52, param=opt_out$par)
plotGPR(x,pred)
