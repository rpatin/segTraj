# Mstep_simultanee
#' Mstep_simultanee computes the MLE within the EM framework

#' @param x the bivariate signal
#' @param rupt the rupture dataframe
#' @param phi the parameters of the mixture
#' @param tau the K*P matrix containing posterior probabilities of membership to clusters
#' @return phi the updated value of the pramaeters

Mstep_simultanee <- function (x,rupt,tau,phi) {
  
  K = nrow(tau)
  P = ncol(tau)
  m = matrix(nrow=2,ncol=P)
  s = matrix(nrow=2,ncol=P)
  prop = matrix(nrow=1,ncol=P)
  Yk = apply(rupt,1,FUN=function(y) rowSums(x[,y[1]:y[2]]))
  nk = rupt[,2]-rupt[,1]+1
  n  = sum(nk)

  #
  np    = nk %*% tau
  m=Yk %*% tau/rep(np,each=2)
  for (i in 1:2){
    s[i,]=  colSums( tau*(sapply(1:P, function(p) {apply(rupt,1,FUN=function(y) sum((x[i,y[1]:y[2]]-m[i,p])^2   ))}))) 
  }
  s=sqrt(s/rep(np,each=2))

  prop = apply(tau,2,sum)/K
  b    = order(m[1,])
  m    = m[,b]
  s    = s[,b]
  prop = prop[b]
  phi  = list(mu=m,sigma=s,prop=prop)

  invisible(phi)
}
