# bisig_plot
#' bisig_plot draws the plots of the bivaraite signal on the same plot (scale free)
#' @param x the signal to be plotted
#' @param rupt optionnal, if given add vertical lines at change points (rupt should)
#' 
#' @export
#' @return no value

bisig_plot<- function(x, rupt=NULL, mu=NULL, pop=NULL){
  n <- dim(x)[2]
  K <- dim(rupt)[1]
  m <- rowMeans(x)
  s <- apply(x,1,sd)
  
  x.norm <- sweep(x = x, MARGIN = 1, STATS = m)
  x.norm <- sweep(x = x.norm, MARGIN = 1, STATS = s, FUN = '/')
  plot(1:n,x.norm[1,], pch=19, col=2, ylab='', xlab = '', xaxt='n', yaxt='n')
  points(1:n, x.norm[2,], pch=19, col=3)
  if(!is.null(rupt)){
    invisible( lapply(1:K, function(d){ abline(v=rupt[d,1], lwd=1.5)}))
  }
  if(!is.null(mu) & !is.null(pop)){
    invisible( lapply(1:K, function(d){ 
      lines(rupt[d,], (rep(mu[1,pop[d]],2)-m[1])/s[1], col=2)
      lines(rupt[d,], (rep(mu[2,pop[d]],2)-m[2])/s[2], col=3)
    }))
  }
  
  
}



# matrixRupt
#' matrixRupt transforms a vector of change point into a data.frame with start and end of every segment
#' @param x the 
#' @param vectorRupt the vector containing the change point
#' 
#' @export
#' @return the matrix of change point

matrixRupt <- function(x,  vectorRupt){
 n <- ncol(x)
  posZeros <-  which(vectorRupt==0)
 if(length(posZeros)>0){
  K <- posZeros[1]-1
  vectorRupt <- vectorRupt[1:K]
 } else { }
  return(matrix(c(1, vectorRupt[1:(K-1)]+1, vectorRupt[1:(K)]), ncol=2))
}