



forbidden.path <- function(G,ibp){
  V = length(ibp)
  n = nrow(G)
  G[1:ibp[1],(ibp[1]+1):n]=Inf
  for (v in (1:(V-1))){
    G[ (ibp[v]+1): ibp[v+1], (ibp[v+1]+1):n]=Inf
  }
  invisible(G) 
}#end function
