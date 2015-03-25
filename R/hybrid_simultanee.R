# hybrid_simultanee
#' \code{hybrid_simultanee} performs a simultaneous seg - clustering for bivariate signals.
#' 
#' It is an algorithm which combines dynamic programming
#' and the EM algorithm to calculate the MLE of phi and T, which 
#' are the mixture parameters and the change point instants.
#' this algorithm is run for a given number of clusters,
#' and estimates the parameters for a segmentation/clustering
#' model with P clusters and 1:Kmax segments
#' 
#' @param x the two-dimensionnal signal, one line per dimension
#' @param P the number of classes 
#' @param Kmax the maximal number of segments
#' 
#' @export
#' @return  a list with tau posterior probability 
#
#' @examples
#' K  <- 5; rupt <- sample(1:20, K+1, replace=TRUE); rupt <- cumsum(rupt); 
#' n <- max(rupt)
#' muSim <- matrix(rnorm(2*K+2,  mean=20, sd=5), nrow=2) 
#' muSim <- apply(muSim,1, cumsum)
#' muSim <- t(muSim)
#' sdSim <- matrix(sqrt(1/rgamma(2*K+2, shape = 10, rate = 10)), nrow=2)
#' print(muSim)
#' pos <- lapply(1:(K+1), function(d) 1*(rupt[d]<(1:n )))
#' pos <- Reduce('+', x=pos)+1
#' x <- matrix(rnorm(2*length(pos), mean=muSim[,pos], sd=sdSim[,pos]), nrow=2)
#' bisig_plot(x = x)
#' n  = dim(x)[2]
#' res <- hybrid_simultanee(x, P=2, Kmax=10)
#' Kopt=5
#' param <- res$param[[Kopt]]
#' bisig_plot(x = x, rupt = param$rupt, mu=param$phi$mu )

hybrid_simultanee <- function(x,P,Kmax){

Linc  = matrix(-Inf,nrow=Kmax,ncol=1)
n     = dim(x)[2]
param = list()
Kmin  = P
  
if (P==1){
  rupt               = c(1, n)
  phi                = list(mu=matrix(rowMeans(x),ncol=1),sigma=matrix(apply(x,1,sd),ncol=1),prop=1)
  Linc[Kmin:Kmax]    = logdens_simultanee(x,phi)
  param[[1]]         = list(phi=phi,rupt=rupt)

#phi contient les moyennes, ecart-types et proportions de chaque composante du m?lange  
  
} else {

  # Rq: lmin=2 for the initialization step because of the hierarchical clustering
  lmin   = 2
  G      = Gmean_simultanee(x,lmin)
  out    = DynProg(G,Kmax=Kmax)
  
  for (K in (Kmin:Kmax)){
    
    cat("P",P,"K", K ,"\n")
    j      = 0
    delta  = Inf
    empty  = 0
    dv     = 0
    th     = out$t.est[K,1:K]
    rupt   = matrix(ncol=2,c(c(1,th[1:(K-1)]+1),th)) 
    phi    = EM.init_simultanee(x,rupt,K,P)
    out.EM = EM.algo_simultanee(x,rupt,P,phi)
    phi    = out.EM$phi
    tau    = out.EM$tau
  bisig_plot(x, rupt = rupt)
    while ( (delta>1e-4) & (empty==0) & (dv==0) & (j<=100)){           

      j          = j+1 
      phi.temp   = phi
      G          = Gmixt_simultanee(x,lmin=lmin,phi)
  
      out.DP     = DynProg(G,K)
      t.est      = out.DP$t.est
      J.est      = out.DP$J.est
      rupt       = matrix(ncol=2,c(c(1,t.est[K,1:(K-1)]+1),t.est[K,]))
      out.EM     = EM.algo_simultanee(x,rupt,P,phi.temp)
      phi        = out.EM$phi
      tau        = out.EM$tau
      empty      = out.EM$empty
      dv         = out.EM$dv
      lvinc.mixt = out.EM$lvinc
      bisig_plot(x, rupt = rupt)
      delta      =max(unlist(lapply(names(phi),function(d) {max(abs(phi.temp[[d]]-phi[[d]])/phi[[d]])})))
      

    }#end while
    

    
    Linc[K]=lvinc.mixt
    param[[K]] = list(phi=phi,rupt=rupt, tau=tau)

  } #end K

}  
#  Ltmp= rep(-Inf,Kmax)
# cat("tracking local maxima for P =",P,"\n")

#  while (sum(Ltmp!=Linc)>=1) {
#    # find the convex hull of the likelihood
#    Ltmp     = Linc  
#    kvfinite = which(is.finite(Linc[P:Kmax]))+P-1
#    Lfinite  = Linc[kvfinite]
#    a        = chull(Lfinite)
#    a        = kvfinite[a]
#    oumin    = which(a==(Kmin))
#    oumax    = which(a==(Kmax))
#    a        = a[oumin:oumax]
#    kvfinite = sort(a)
#    # find the coordinates of points out of the convex hull
#    Kconc    = c(1:Kmax)
#    Kconc    = Kconc[-which(Kconc %in% c(kvfinite))]
#    Kconc    = Kconc[Kconc>=Kmin]   

#    for (k in Kconc){
      
#      out.neighbors  = neighbors(L=Linc,k=k,param=param,ibp=ibp,P=P,lmin=lmin)
#      param          = out.neighbors$param
#      Linc           = out.neighbors$L
 #     #plot(1:length(Linc),Linc,col=1)     lines(1:length(Ltmp),Ltmp,col=2)     lines(k,Linc[k],col=3)
#
#    } # end k

#    out.neighbors  = neighbors(L=Linc,k=Kmax,param=param,ibp=ibp,P=P,lmin=lmin)
#    param          = out.neighbors$param
 #   Linc           = out.neighbors$L

#  } # end while
  
#} # end else Pmin ==1

invisible(list(Linc=Linc,param=param))

} #end function
