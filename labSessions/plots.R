##########################################
## utils
lightblue <- rgb(114/255,159/255,207/255,.3) ; darkblue <- rgb(32/255,74/255,135/255,1) ; darkbluetr <- rgb(32/255,74/255,135/255,.3)
darkPurple <- "#5c3566" ; darkBlue <- "#204a87" ; darkGreen <- "#4e9a06" ; darkChocolate <- "#8f5902" ; darkRed  <- "#a40000" ; darkOrange <- "#ce5c00" ; darkButter <- "#c4a000"

##########################################
## plot functions

#' Plot GPR
#'
#' Plot predictions for 1-dimensional GP regression model.
#' 
#' @param x $p \times d$ matrix of the points to predict (\code{p} is the number of prediction points and \code{d} is the input space dimension).
#' @param pred output or \code{predGPR} (or list with mean values (\code{n*1}-matrix) and a \code{n*n}-covariance matrix).
#' @param xlab,ylab,ylim usual graphical parameters (optional)
#' @param add Boolean, should the graph be added to the current plot (creates new graph by default).
#' @param plot_observation_points Boolean, should the observation points be plotted (default to TRUE).
#' 
#' @examples
#' n <- 5
#' p <- 101
#' x <- matrix(seq(from=0, to=1, length=p))
#' X <- matrix(seq(from=0.1, to=0.9, length=n))
#' F <- matrix(c(0.5, 0, 1.5, 3, 2))
#' 
#' pred <- predGPR(x,X,F,kGauss)
#' plotGPR(x,pred)
#' 
#' @export
plotGPR <- function(x,pred,xlab="",ylab="",ylim=NULL,add=FALSE, plot_observation_points=TRUE){
  m <- pred[[1]]
  sdDiag <- sqrt(pmax(0,diag(pred[[2]])))
  upp95 <- m + 1.96*sdDiag
  low95 <- m - 1.96*sdDiag
  if(is.null(ylim)) ylim <- range(low95-0.5,upp95+0.5)
  par(mar=c(4.5,5.1,1.5,1.5))
  if(!add){
    plot(x, m, type="n", xlab=xlab,ylab=ylab, ylim=ylim, cex.axis=1.5,cex.lab=2)
  }
  polygon(c(x,rev(x)),c(upp95,rev(low95)),border=NA,col=lightblue)
  lines(x,m,col=darkblue,lwd=3)
  lines(x,low95,col=darkbluetr)  
  lines(x,upp95,col=darkbluetr)  
  if(plot_observation_points) points(X, F, pch=4, cex=1,lwd=3)
}



