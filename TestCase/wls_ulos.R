wls_ulos <- function(xx){
  # Weighted Least Squares distance function for ulos vectors
  #
  # The covariance matrix is passed through global variable 
  # xs,ys,zs,a and p come from the variables xx
  
  G = 2000 # Shear modulus in MPa
  nu = 0.25 # Poisson's ratio
  nlos = c(-0.664,-0.168,0.728) # vector of direction of line of sight (satellite)
  #
  n2i <- Glb_var$n2i
  
  # Compute surface displacements
  U <- mogi_3D(G,nu,xx[n2i$xs],xx[n2i$ys],xx[n2i$zs],xx[n2i$a],xx[n2i$p],
               Glb_xi,Glb_yi,Glb_zi)
  # project along los
  ulos <- nlos[1]*U$x+nlos[2]*U$y+nlos[3]*U$z
  # calculate weighted least squares
  wls <- t((ulos-Glb_ulos))%*%Glb_CXinv%*%(ulos-Glb_ulos)
  return(wls)  
}