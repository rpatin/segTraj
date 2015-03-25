 
# calculate the cost matrix for a segmentation model with changes in the mean
# lmin : minimum size for a segment
# G(i,j) = variance of the data between point (i+1) to point j
Gmean<-function(x,lmin=1)
{
  n = length(x)
  matD=matrix(Inf,n,n)
  x2=x^2
  x2i=cumsum(x2)
  xi=cumsum(x)
  x2i=x2i[lmin:n]
  xi=xi[lmin:n]
  matD[1,lmin:n]=x2i-((xi^2)/(lmin:n))
  nl=n-lmin+1
  for (i in 2:nl)
    {
      ni=n-i-lmin+3
      x2i=x2i[2:ni]-x2[i-1] 
      xi=xi[2:ni]-x[i-1]
      deno<-((i+lmin-1):n)-i+1
      matD[i,(i+lmin-1):n]=x2i-((xi^2)/deno)
    }
  invisible(matD)
}
