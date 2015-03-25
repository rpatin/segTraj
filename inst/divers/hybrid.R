# hybrid algorithm which combines dynamic programming
# and the EM algorithm to calculate the MLE of phi and T, which 
# are the mixture parameters and the change point instants.
# this algorithm is run for a given number of clusters,
# and estimates the parameters for a segmentation/clustering
# model with P clusters and 1:Kmax segments

hybrid <- function(x,ibp,P,Kmax){

Linc  = matrix(-Inf,nrow=Kmax,ncol=1)
V     = length(ibp)
n     = length(x)
param = list()
Kmin  = max(P,V+1)
  
if (P==1){
    

  rupt               = c(1, length(x))
  phi                = c(mean(x),sd(x),1)
  Linc[Kmin:Kmax]    = logdens(x,Pmin,phi)
  param[[1]]         = list(phi=phi,rupt=rupt)

} else {

  for (K in (Kmin:Kmax)){
    
    cat("P",P,"K", K ,"\n")
    j      = 0
    delta  = Inf
    empty  = 0
    dv     = 0
    # Rq: lmin=2 for the initialization step because of the hierarchical clustering
    lmin   = 2
    G      = Gmean(x,lmin)

    # consider the imposed breakpoints 
    # avoid some paths
    if (V>0){
      G = forbidden.path(G,ibp)
    }
    out    = DynProg(G,Kmax=K)
    t.est  = out[[2]]
    J.est  = out[[1]]
    th     = t.est[K,1:K]
    rupt   = matrix(ncol=2,c(c(1,th[1:(K-1)]+1),th)) 
    phi    = EM.init(x,rupt,K,P)
    out.EM = EM.algo(x,rupt,P,phi)
    phi    = out.EM$phi
    

    while ( (delta>1e-4) & (empty==0) & (dv==0) & (j<=100)){           

      j          = j+1 
      phi.temp   = phi
      G          = Gmixt(x,lmin=lmin,phi,P)

      # considers the imposed breakpoints
      # avoids some paths
      if (V>0){
        G          =forbidden.path(G,ibp)
      }
     
      out.DP     = DynProg(G,K)
      t.est      = out.DP$t.est
      J.est      = out.DP$J.est
      rupt       = matrix(ncol=2,c(c(1,t.est[K,1:(K-1)]+1),t.est[K,]))
      out.EM     = EM.algo(x,rupt,P,phi.temp)
      phi        = out.EM$phi
      empty      = out.EM$empty
      dv         = out.EM$dv
      lvinc.mixt = out.EM$lvinc
      delta      = max(abs(phi.temp-phi)/phi)

    }#end while
    

    
    Linc[K]=lvinc.mixt
    param[[K]] = list(phi=phi,rupt=rupt)

  } #end K

  
  Ltmp= rep(-Inf,Kmax)
  cat("tracking local maxima for P =",P,"\n")

  while (sum(Ltmp!=Linc)>=1) {
    # find the convex hull of the likelihood
    Ltmp     = Linc  
    kvfinite = which(is.finite(Linc[P:Kmax]))+P-1
    Lfinite  = Linc[kvfinite]
    a        = chull(Lfinite)
    a        = kvfinite[a]
    oumin    = which(a==(Kmin))
    oumax    = which(a==(Kmax))
    a        = a[oumin:oumax]
    kvfinite = sort(a)
    # find the coordinates of points out of the convex hull
    Kconc    = c(1:Kmax)
    Kconc    = Kconc[-which(Kconc %in% c(kvfinite))]
    Kconc    = Kconc[Kconc>=Kmin]   

    for (k in Kconc){
      
      out.neighbors  = neighbors(L=Linc,k=k,param=param,ibp=ibp,P=P,lmin=lmin)
      param          = out.neighbors$param
      Linc           = out.neighbors$L
      #plot(1:length(Linc),Linc,col=1)     lines(1:length(Ltmp),Ltmp,col=2)     lines(k,Linc[k],col=3)

    } # end k

    out.neighbors  = neighbors(L=Linc,k=Kmax,param=param,ibp=ibp,P=P,lmin=lmin)
    param          = out.neighbors$param
    Linc           = out.neighbors$L

  } # end while
  
} # end else Pmin ==1

invisible(list(Linc=Linc,param=param))

} #end function
