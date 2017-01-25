# Gmean_simultanee
#' Gmean_simultanee  calculate the cost matrix for a segmentation model with changes in the mean
#' @param Don the bivariate signal
#' @param lmin minimum size for a segment, default value is 2
#' @return the cost matrix G(i,j) which contains the variance of the data between point (i+1) to point j
#' @export
Gmean_simultanee<-function(Don,lmin)
{
  n = dim(Don)[2] 

  
  ## every element of the list is the cost motrix for one signal
  matD_list=lapply(1:2,function(nb) {
    Res=matrix(Inf,n,n)
    z=Don[nb,] #depend de la forme des donnees    
    z2=z^2
    z2i=cumsum(z2)
    zi=cumsum(z)
    z2i=z2i[lmin:n]
    zi=zi[lmin:n]
    Res[1,lmin:n]=z2i-((zi^2)/(lmin:n))
    nl=n-lmin+1
    for (i in 2:nl)
    {
      ni=n-i-lmin+3
      z2i=z2i[2:ni]-z2[i-1] 
      zi=zi[2:ni]-z[i-1]
      deno<-((i+lmin-1):n)-i+1
      Res[i,(i+lmin-1):n]=z2i-((zi^2)/deno)
    }
    return(Res)
  })
  
  matD=Reduce("+",matD_list)
  invisible(matD)
}
## debug note : lmin is correctly taken into account
