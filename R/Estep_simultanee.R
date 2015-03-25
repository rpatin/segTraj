# Estep_simultanee
#' Estep_simultanee computes posterior probabilities and incomplete-data  log-likelihood for mixture models
#' @param logdensity is a  K*P matrix containing the conditinal log-densities for each segment
#' @param phi   a list containing the parameters of the mixture
#' @return a list with tau a K*P matrix, tau kj is the posterior probability for segment k to belong to classe j and lvinc, the incomplete log vrais P(X=x)

Estep_simultanee <- function (logdensity,phi){
  K = nrow(logdensity)
  P = ncol(logdensity)
  tau     = matrix((log( phi$prop)),K,P,byrow=TRUE)+logdensity
  tau_max = apply(tau,1,max)
  tau     = exp(tau-matrix(tau_max,K,P))
  lvinc   = sum(log( apply(tau,1,sum)) + tau_max) #sum(log P(Xk=xk))
  tau     = sweep(tau,1,STATS = apply(tau,1,sum), FUN="/")

  invisible(list(tau,lvinc))
}

