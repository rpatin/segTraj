 
# calculates Gaussian log-densities for 
# a vector xk 
# P    : number of clusters
# phi  : parameters of the mixture
logdens_simultanee <- function(xk,P,phi){
 
m  = phi$m
s  = phi$s

nk = dim(xk)[2]
tmp = matrix(ncol=P,nrow=1)

tmp=sapply(1:P, function(p){
  -nk*log(sqrt(2*pi)*s[1,p])- 0.5*sum (    (xk[1,]  - m[1,p])^2  )/s[1,p]^2-nk*log(sqrt(2*pi)*s[2,p])- 0.5*sum (    (xk[2,]  - m[2,p])^2  )/s[2,p]^2
} )

invisible(tmp)

}
