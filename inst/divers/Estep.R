 
# computes posterior probabilities and incomplete-data 
# log-likelihood for mixture models
# logdensity : K*P matrix containing the conditinal log-densities for each segment
# phi        : parameters of the mixture
# tau        : K*P matrix containing posterior probabilities of membership to clusters
Estep <- function (logdensity,phi){

  K = nrow(logdensity)
  P = ncol(logdensity)

  tau     = repmat( t (log( phi[(2*P+1):(3*P)] )),K,1)+logdensity
  tau_max = apply(tau,1,max)
  tau     = exp(tau-repmat(tau_max,1,P))
  lvinc  = sum(log( apply(tau,1,sum)) + tau_max)
  tau     = tau / repmat( apply(tau,1,sum) ,1,P)

  invisible(list(tau,lvinc))
}

