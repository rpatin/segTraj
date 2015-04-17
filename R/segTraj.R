#' segTraj: A package for segmentation and clustering of a bivariate signal
#'
#' The foo package allows to perform a segmentation and a segmentation clustering procedure
#' 
#' @section segTraj functions:
#' hybrid_simultanee
#' DynProg
#'
#' @docType package
#' @name segTraj
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

NULL