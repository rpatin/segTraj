# calculates the cost matrix for a segmentation/clustering model
# lmin : minimum size for a segment
# phi  : parameters of the mixture
# P    : number of clusters
# G(i,j) = mixture density for segment between points (i+1) to j
#        = \sum_{p=1}^P \log (\pi_p f(y^{ij};\theta_p))
# Rq: this density if factorized in order to avoid numerical zeros in the log
Gmixt <- function(x,lmin,phi,P){

m    = phi[1:P]
s    = phi[(P+1):(2*P)]
prop = phi[(2*P+1):length(phi)]

n = length(x)
#x = x'

G = repmat(Inf,n,n)

xi  = cumsum(x) 
lg  = lmin:n
xi  = xi[lg]
x2  = x^2
x2i = cumsum(x2)
x2i = x2i[lg]

wk  = repmat(t( x2i/lg-(xi/lg)^2 ),P,1)

#wk=repmat(wk,P,1)

dkp   = (repmat(t(xi),P,1)/repmat(t(lg),P,1)-repmat(m,1,n-1))^2
A     = (wk+dkp)/repmat(s^2,1,n-lmin+1)+log(2*pi*repmat(s^2,1,n-1))
A     = -0.5*repmat(t(lg),P,1)*A +(repmat(log(prop),1,n-1))
A_max = apply(A,2,max)
A     = exp(A-repmat(t (A_max) ,P,1))

#rappel: on fait du plus court chemin
#        donc on prend -LV
G[1,lmin:n] = -log(apply(A,2,sum)) - A_max

for (i in (2:(n-lmin+1))) {
   ni  = n-i-lmin+3
   x2i = x2i[2:ni]-x2[i-1]
   xi  = xi[2:ni]-x[i-1]
   lgi = lmin:(n-i+1)
   wk  = repmat(t(x2i)/(lgi)-(xi/(lgi))^2,P,1)
   dkp = (repmat(t(xi),P,1)/repmat(t(lgi),P,1)-repmat(m,1,ni-1))^2
   A   = (wk+dkp)/repmat(s^2,1,ni-1)+log(2*pi*repmat(s^2,1,ni-1))
   A   = -0.5*repmat(t(lgi),P,1)*A +(repmat(log(prop),1,ni-1))
   A_max = apply(A,2,max)
   A     = exp(A-repmat(t (A_max) ,P,1))
   G[i,(i+lmin-1):n] = -log(apply(A,2,sum)) - A_max
}
invisible(G)
}
