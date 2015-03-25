 
# EM algorithm
# caculates the MLE of phi for given change-point instants
# and for a fixed number of clusters
# rupt ; sequence of change points
# P    : number of clusters
# phi0 ; initial parameter
EM.algo <- function(x,rupt,P,phi0){

  K     = nrow(rupt)
  delta = 1
  empty = 0
  dv    = 0
  tau   = matrix(1,nrow = K,ncol = P)
  phi   = phi0
  iter  = 0
  np    = apply(tau,2,sum)
  eps   = 10e-10
  

  while ( (delta>=1e-4) & (min(np)>eps) & (iter<=5000) ){
    iter       = iter+1
    phi_temp   = phi
    logdensity = t( apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,P,phi)))

    Estepout   = Estep(logdensity,phi)
    tau        = Estepout[[1]]
    lvinc      = Estepout[[2]]

    phi        = Mstep(x,rupt,tau,phi)
    np         = apply(tau,2,sum)
    delta      = max(abs(phi_temp-phi)/phi)
  }

  if (min(np)<eps){
    empty = 1
    lvinc = -Inf
  }

  if (iter>5000){
    dv    = 2
    lvinc = -Inf
  }

  rm(delta,logdensity)
  
    
  invisible(list(phi = phi,tau = tau,lvinc = lvinc,empty = empty,dv = dv))


  
}
