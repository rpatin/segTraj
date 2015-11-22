# Gmixt_simultanee
#'
#' Gmixt_simultanee calculates the cost matrix for a segmentation/clustering model
#' @param Don the bivariate  signal
#' @param lmin the minimum size for a segment
#' @param phi  the  parameters of the mixture
#' @return a matrix  G(i,j), the mixture density for segment between points (i+1) to j
#'         = deqn{\sum_{p=1}^P \log (\pi_p f(y^{ij};\theta_p))}
#'          Rq: this density if factorized in order to avoid numerical zeros in the log

Gmixt_simultanee <- function(Don,lmin,phi){
  
  P = length(phi$prop)
  m    = phi$mu
  s    = phi$sigma
  prop = phi$prop
  
  n = dim(Don)[2]
  #x = x'
  G = matrix(0,ncol=n,nrow=n)
  
  
  for (signal in 1:2){
    z=Don[signal,]
    zi  = cumsum(z) 
    lg  = lmin:n
    zi  = zi[lg]
    z2  = z^2
    z2i = cumsum(z2)
    z2i = z2i[lg]
    
    wk  = repmat(t( z2i/lg-(zi/lg)^2 ),P,1)
    
    #wk=repmat(wk,P,1)
    
    dkp   = (repmat(t(zi),P,1)/repmat(t(lg),P,1)-repmat(m[signal,],1,n-1))^2
    A     = (wk+dkp)/repmat(s[signal,]^2,1,n-lmin+1)+log(2*pi*repmat(s[signal,]^2,1,n-1))
    A     = -0.5*repmat(t(lg),P,1)*A +(repmat(log(prop),1,n-1))
    A_max = apply(A,2,max)
    A     = exp(A-repmat(t (A_max) ,P,1))
    
    #rappel: on fait du plus court chemin
    #        donc on prend -LV
    G[1,lmin:n] = -log(apply(A,2,sum)) - A_max
    
    for (i in (2:(n-lmin+1))) {
      ni  = n-i-lmin+3
      z2i = z2i[2:ni]-z2[i-1]
      zi  = zi[2:ni]-z[i-1]
      lgi = lmin:(n-i+1)
      wk  = repmat(t(z2i)/(lgi)-(zi/(lgi))^2,P,1)
      dkp = (repmat(t(zi),P,1)/repmat(t(lgi),P,1)-repmat(m[signal,],1,ni-1))^2
      A   = (wk+dkp)/repmat(s[signal,]^2,1,ni-1)+log(2*pi*repmat(s[signal,]^2,1,ni-1))
      A   = -0.5*repmat(t(lgi),P,1)*A +(repmat(log(prop),1,ni-1))
      A_max = apply(A,2,max)
      A     = exp(A-repmat(t (A_max) ,P,1))
      G[i,(i+lmin-1):n] =  G[i,(i+lmin-1):n]-log(apply(A,2,sum)) - A_max
    }}
  for (i in (lmin-1):n){
    for (j in 1:i){
      G[i,j]=Inf
    }}
  invisible(G)
}
