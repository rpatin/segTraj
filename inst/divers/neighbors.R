# tests whether neighbors of point k,P can 
# be used to re-initialize the EM algorithm and
# to improve the log-likelihood.

neighbors <- function (L,k,param,ibp,P,lmin) {


  V = length(ibp)

  
    # left neighbor
    a  = max(L[1:(k-1)])
    K1 = which.max(L[1:(k-1)])
    
    
    if ( (is.null(K1)) || (a==-Inf) || (is.na(a) ) ){
      K1            =-Inf
      phi1          = rep(-Inf,3*P)
      out.EM1$lvinc = list(lvinc=- Inf)
    }    else { 
    phi1                     = param[[K1]]$phi     
    G                        = Gmixt(x,lmin=lmin,phi1,P)
    if (V>0){
      G          =forbidden.path(G,ibp)
    }
    
    out.DP     = DynProg(G,k)
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    rupt1      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EM1    = EM.algo(x,rupt1,P,phi1)
    } #end else

    
    
    # right neighbor 
    a  = max(L[(k+1):Kmax])  
    K2 = which.max(L[(k+1):Kmax])
    K2 = K2 + k
    
    if ( (K2==0) || (a==-Inf) || (is.na(a)) ){
      K2            = -Inf
      phi2          = rep(-Inf,3*P)
      out.EM2       = list(lvinc=-Inf)
    }   else {
    phi2                     = param[[K2]]$phi
    G                        = Gmixt(x,lmin=lmin,phi2,P)
    if (V>0){
      G          =forbidden.path(G,ibp)
    }

    out.DP     = DynProg(G,k)
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    rupt2      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EM2    = EM.algo(x,rupt2,P,phi2)
    
    } #end else

    
    
    
    # choice of the best likelihood
    a = which.max( c(out.EM1$lvinc, out.EM2$lvinc,  L[k]) ) 

    
    # parameter update
    if (a==1){
      param[[k]] = list(phi = out.EM1$phi, rupt = rupt1)
      L[k] = out.EM1$lvinc}
    if (a==2) {
      param[[k]] = list(phi = out.EM2$phi,rupt = rupt2)
      L[k] = out.EM2$lvinc}
      
      
      
      
    invisible(list(L=L,param=param))  
      
      } #end function
      
