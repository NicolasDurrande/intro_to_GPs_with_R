restartEGO <- function(fn,restart=TRUE,nbinit=nbinit)
# utility functions really attached to mainInversionPunctualDisplSource
# function to restart an EGO run from stored data file fn
# reloads the data base , calculates n_restart
{
  if (file.exists(fn) & restart==TRUE){
    load(file=fn)
    n_restart <- dim(Xnorm)[1]-nbinit
    # recalculate mean_wls and std_wls
    lwls <- log(1+wls[1:nbinit])
    mean_wls <<- mean(lwls)
    std_wls <<- sd(lwls)
    Xnorm <<- Xnorm
    X <<- X
    wls <<- wls
    norm_wls <<- norm_wls
    oLL <<- oLL
    return(n_restart)
  }
  else {return(n_restart<-0)}
}