# Gmean_simultanee
#' Gmean_simultanee  calculate the cost matrix for a segmentation model with changes in the mean
#' @param Don the bivariate signal
#' @param lmin minimum size for a segment, default value is 2
#' @return the cost matrix G(i,j) which contains the variance of the data between point (i+1) to point j
#' @export
Gmean_simultanee<-function(Don,lmin=2)
{
  n = dim(Don)[2] 

  
  ## every element of the list is the cost motrix for one signal
  matD_list=lapply(1:2,function(nb) {
    Res=matrix(Inf,n,n)
    x=Don[nb,] #depend de la forme des donnees    
    x2=x^2
    x2i=cumsum(x2)
    xi=cumsum(x)
    x2i=x2i[lmin:n]
    xi=xi[lmin:n]
    Res[1,lmin:n]=x2i-((xi^2)/(lmin:n))
    nl=n-lmin+1
    for (i in 2:nl)
    {
      ni=n-i-lmin+3
      x2i=x2i[2:ni]-x2[i-1] 
      xi=xi[2:ni]-x[i-1]
      deno<-((i+lmin-1):n)-i+1
      Res[i,(i+lmin-1):n]=x2i-((xi^2)/deno)
    }
    return(Res)
  })
  
  matD=Reduce("+",matD_list)
  invisible(matD)
}
