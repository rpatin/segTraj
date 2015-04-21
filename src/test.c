#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>



/* Calculates the loglikelihood of a normal variable with mean m and standard deviation s.
   x is a R vector containing the values of the variable.
 */
double logdensParVariable(SEXP x, double m, double s)
{
    int n, i;
    double su;
    n=length(x);
    su = 0.0;
    
    for (i=0; i<n; i++) {
	su += dnorm(REAL(x)[i], m, s, 1);	
    }

    return(su);
}


/* Calculates the log-likelihood of a sample drawn from a bivariate normal distribution 
   (with zero correlation between the two dimensions). The two dimensions are stored
   in two vectors x1 and x2, their means are in the 1-element R vectors mu1 and mu2 
   respectively, and their variance is stored in the 1-element R vectors sigma1 and sigma2
   respectively.
 */ 
SEXP logdensSimultaneeCore(SEXP x1, SEXP x2, SEXP mu1, SEXP mu2, SEXP sigma1, 
			     SEXP sigma2)
{
    int P, i;
    double *m1, *m2, *s1, *s2;
    SEXP tmp;
    
    P = length(mu1);
    m1 = REAL(mu1);
    s1 = REAL(sigma1);
    m2 = REAL(mu2);
    s2 = REAL(sigma2);
        
    PROTECT(tmp = allocVector(REALSXP, P));
    
    for (i = 0; i < P; i++) {
	/* The log-density is the sum of the logdensities of the two
	   variables */
	REAL(tmp)[i] = logdensParVariable(x1, m1[i], s1[i]) + 
	    logdensParVariable(x2, m2[i], s2[i]);	
    }
    
    UNPROTECT(1);
    return(tmp);
}


/* Log-density of a sample x (matrix) drawn from a bivariate
   normal distribution with 0 correlation. The mean vector
   is stored in mu, and the vector of SD is in sigma
 */
SEXP logdensSimultanee(SEXP x, SEXP mu, SEXP sigma)
{
    SEXP x1, x2, mu1, mu2, sigma1, sigma2, resu;
    int i, n;

    PROTECT(x1 = allocVector(REALSXP, length(x)/2));
    PROTECT(x2 = allocVector(REALSXP, length(x)/2));
    PROTECT(mu1 = allocVector(REALSXP, length(mu)/2));
    PROTECT(mu2 = allocVector(REALSXP, length(mu)/2));
    PROTECT(sigma1 = allocVector(REALSXP, length(mu)/2));
    PROTECT(sigma2 = allocVector(REALSXP, length(mu)/2));
    n = length(x1);
    
    for (i = 0; i<n; i++) {
	REAL(x1)[i] = REAL(x)[i];
	REAL(x2)[i] = REAL(x)[i+n];
    }

    n = length(mu1);
    for (i = 0; i<n; i++) {
	REAL(mu1)[i] = REAL(mu)[i];
	REAL(mu2)[i] = REAL(mu)[i+n];
	REAL(sigma1)[i] = REAL(sigma)[i];
	REAL(sigma2)[i] = REAL(sigma)[i+n];
    }
    
    PROTECT(resu=logdensSimultaneeCore(x1, x2, mu1, mu2, sigma1, sigma2));
    
    UNPROTECT(7);
    return(resu);
}


/* function calculating the log-density of a sample of points x (matrix with two columns) 
   drawn form a bivariate normal distribution with zero correlation.
   The list phi contains two elements: (0) the vector mu of means, and (1) the SD of the
   two normal distributions. This logdensity is calculated only for rows of x with
   indices comprised between i1 and i2
*/
SEXP logdensSimultaneeParIndice(SEXP x, SEXP phi, int i1, int i2)
{
    SEXP sousx, mu, sd, resu;
    int i, n, k, nrx;
    n = i2-i1+1;
    nrx = nrows(x);
    
    PROTECT(sousx = allocMatrix(REALSXP, n, 2));
    k=0;
    for (i=i1; i<=i2; i++) {
	REAL(sousx)[k] = REAL(x)[i];
	REAL(sousx)[k+n] = REAL(x)[i+nrx];
	k++;
    }
    mu = VECTOR_ELT(phi,0);
    sd = VECTOR_ELT(phi,1);
    
    PROTECT(resu = logdensSimultanee(sousx, mu, sd));
    
    UNPROTECT(2);
    
    return(resu);
}





/* The function GmeanSimultanee (after the following one) calculates the cost matrix 
   for a segmentation model of a bivariate series with changes in the mean. 
   This cost matrix is calculated as the sum of the costs for each signal independently.
   The function GmeanSimultaneeCore calculates these "univariate" costs.
   Don is the bivariate signal (matrix with n rows and 2 columns) which is to be 
   segmented. lminr is the R object containing the minimum size for a segment. 
   nb is an integer defining for which "variable" (column index of Don) the cost is to be 
   calculated
*/
SEXP GmeanSimultaneeCore(SEXP Don, SEXP lminr, int nb)
{
    int lmin, n, i, j;
    SEXP x2i, xi, Res;

    lmin=INTEGER(lminr)[0];
    n = length(Don)/2;

    PROTECT(xi = allocVector(REALSXP, n));
    PROTECT(x2i = allocVector(REALSXP, n));
    PROTECT(Res=allocMatrix(REALSXP, n, n));
    
    /* Cumulative sums of the series */
    REAL(xi)[0] = REAL(Don)[nb*n];
    REAL(x2i)[0] = R_pow(REAL(Don)[nb*n],2.0);
    for (i=1; i < n; i++) {
	REAL(xi)[i] = REAL(xi)[i-1] + REAL(Don)[nb*n + i];
	REAL(x2i)[i] = REAL(x2i)[i-1] + R_pow(REAL(Don)[nb*n + i],2.0);
    }

    /* Initialisation of the cost at infinity (infinite cost) */
    for (i=0; i < (n*n); i++) {
	REAL(Res)[i] = R_PosInf;
    }
    
    /* Calculation of the cost based on the MSD */
    for (i=(lmin-1); i<n; i++) {
	REAL(Res)[i*n] = REAL(x2i)[i] - (R_pow(REAL(xi)[i], 2.0))/(((double) (i+1)));
    }
    for (i = 1; i<n; i++) {
	for (j=0; j<n; j++) {
	    REAL(x2i)[j] = REAL(x2i)[j] - R_pow(REAL(Don)[nb*n + (i-1)],2.0);
	    REAL(xi)[j] = REAL(xi)[j] - REAL(Don)[nb*n + (i-1)];
	}
	for (j=(i+lmin-1); j<n; j++) {
	    REAL(Res)[i+j*n] = REAL(x2i)[j] - (R_pow(REAL(xi)[j], 2.0))/(((double) (j-i+1)));
	}
	
    }
    UNPROTECT(3);
    
    /* Returns the cost */
    return(Res);    
}


/* The function GmeanSimultanee calculates the cost matrix 
   for a segmentation model of a bivariate series with changes in the mean. 
   This cost matrix is calculated as the sum of the costs for each signal independently.
   Don is the bivariate signal (matrix with n rows and 2 columns) which is to be 
   segmented. lminr is the R object containing the minimum size for a segment.
*/
SEXP GmeanSimultanee(SEXP Don, SEXP lminr)
{
    int i;
    SEXP res0, res1, res;
    
    PROTECT(res0 = GmeanSimultaneeCore(Don, lminr, 0));
    PROTECT(res1 = GmeanSimultaneeCore(Don, lminr, 1));
    PROTECT(res = allocVector(REALSXP, length(res0)));
    
    for (i=0; i< length(res0); i++) {
	REAL(res)[i] = REAL(res0)[i] + REAL(res1)[i];
    }
    
    UNPROTECT(3);
    return(res);
}




/* 
   Gmixt_simultanee calculates the cost matrix for a segmentation/clustering model 
   (contrarily to the previous function, it does not assume a simple change in the mean, 
   but also accounts for changes in variance). The result is  a matrix  G(i,j), the 
   mixture density for segment between points (i+1) to j
   = deqn{\sum_{p=1}^P \log (\pi_p f(y^{ij};\theta_p))}
   Remark: this density if factorized in order to avoid numerical zeros in the log.
   Don is the bivariate signal, lmin is the minimum size of the segment, and
   phi is a list with three elements: (i) the matrix of means of the models (one row per model
   one column per variable), (ii) the matrix of standard deviations of the P models (one
   row per model, and one column per variable), (iii) the vector of probabilities associated
   to each model.
 */
SEXP GmixtSimultanee(SEXP Don, SEXP lminr, SEXP phi)
{
    int P, n, signal, i, j, k, lmin, ii;
    SEXP m, s, prop, G, xi, x2i, wk, dkp, A, Amax, Asu;

    /* Initialization and memory allocation */
    P = length(VECTOR_ELT(phi,2));
    m = VECTOR_ELT(phi,0);
    s = VECTOR_ELT(phi,1);
    prop = VECTOR_ELT(phi,2);
    n = nrows(Don);
    lmin = INTEGER(lminr)[0];
    k = 0;
    
    PROTECT(G = allocMatrix(REALSXP, n, n));
    PROTECT(xi = allocVector(REALSXP, n));
    PROTECT(x2i = allocVector(REALSXP, n));
    PROTECT(wk = allocMatrix(REALSXP, P, n-lmin+1));
    PROTECT(dkp = allocMatrix(REALSXP, P, n-lmin+1));
    PROTECT(A = allocMatrix(REALSXP, P, n-lmin+1));
    PROTECT(Amax = allocVector(REALSXP, n-lmin+1));
    PROTECT(Asu = allocVector(REALSXP, n-lmin+1));
    

    /* Initialization of G = 0 */
    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    REAL(G)[i+n*j] = 0.0;
	}
    }

    
    /* For each variable in the signal */
    for (signal = 0; signal < 2; signal++) {
	
	
	REAL(xi)[0] = REAL(Don)[signal*n];
	REAL(x2i)[0] = R_pow(REAL(Don)[signal*n],2.0);
	
	for (i = 1; i < n; i++) {
	    REAL(xi)[i] = REAL(xi)[i-1] + REAL(Don)[i+signal*n];
	    REAL(x2i)[i] = REAL(x2i)[i-1] + R_pow(REAL(Don)[i+signal*n],2.0);

	}

	/* Calculation of wk */
	k = 0;
	for (j=lmin-1; j < n; j++) {
	    for (i=0; i<P; i++) {
		REAL(wk)[i+k*P] = (REAL(x2i)[j]/((double) (j+1))) - 
		    R_pow((REAL(xi)[j]/((double) (j+1))),2.0);
	    }
	    k++;
	}
	
	/* Calculation of dkp */
	k = 0;
	for (j=lmin-1; j < n; j++) {
	    for (i=0; i<P; i++) {
		REAL(dkp)[i+k*P] = R_pow((REAL(xi)[j]/((double) (j+1))) - 
					 REAL(m)[i+signal*P],2.0);
	    }
	    k++;
	}
	
	/* Calculation of A */
	k = 0;
	for (j=lmin-1; j < n; j++) {
	    for (i=0; i<P; i++) {
		REAL(A)[i+k*P] = -0.5*((double) (j+1))*(((REAL(wk)[i+k*P] + 
							  REAL(dkp)[i+k*P])/
							 R_pow(REAL(s)[i+signal*P],2.0)) + 
							log(M_2PI*R_pow(REAL(s)[i+signal*P],2.0))) + 
		    log(REAL(prop)[i]);
	    }
	    k++;
	}
	
	/* Calculation of Amax */
	for (i=0; i<ncols(A); i++) {
	    REAL(Amax)[i] = REAL(A)[P*i];
	    for (j=0; j<P; j++) {
		if (REAL(Amax)[i] < REAL(A)[j+P*i]) {
		    REAL(Amax)[i] = REAL(A)[j+P*i];
		}
	    }
	}
	
	/* Subtracts Amax */
	for (i=0; i<ncols(A); i++) {
	    for (j=0; j<P; j++) {
		REAL(A)[j+P*i] = exp(REAL(A)[j+P*i] - REAL(Amax)[i]);
	    }
	}
	
	/* First row of G */
	for (i=0; i<ncols(A); i++) {
	    REAL(Asu)[i] = REAL(A)[P*i];
	    for (j=1; j<P; j++) {
		REAL(Asu)[i] += REAL(A)[j+P*i];
	    }
	}
	k = 0;
	for (i=lmin-1; i<n; i++) {
	    REAL(G)[n*i] = -log(REAL(Asu)[k]) - REAL(Amax)[k];
	    k++;
	}
	
	/* Rest of the matrix */
	for (i = 1; i<(n-lmin+1); i++) {
	    
	    /* We move one step further, we remove the (i-1)th step */
	    for (j = 0; j < n; j++) {
		REAL(x2i)[j] = REAL(x2i)[j] - R_pow(REAL(Don)[i-1+signal*n],2.0);
		REAL(xi)[j] = REAL(xi)[j] - REAL(Don)[i-1+signal*n];
	    }
	    
	    /* Calculation of wk */
	    k = 0;
	    for (j=lmin-1+i; j < n; j++) {
		for (ii=0; ii<P; ii++) {
		    REAL(wk)[ii+k*P] = (REAL(x2i)[j]/((double) (j+1-i))) - 
			R_pow((REAL(xi)[j]/((double) (j+1-i))),2.0);
		}
		k++;
	    }
	    
	    /* Calculation of dkp */
	    k = 0;
	    for (j=lmin-1+i; j < n; j++) {
		for (ii=0; ii<P; ii++) {
		    REAL(dkp)[ii+k*P] = R_pow((REAL(xi)[j]/((double) (j+1-i))) - 
					      REAL(m)[ii+signal*P],2.0);
		}
		k++;
	    }
	    
	    /* Calculation of A */
	    k = 0;
	    for (j=lmin-1+i; j < n; j++) {
		for (ii=0; ii<P; ii++) {
		    REAL(A)[ii+k*P] = -0.5*((double) (j+1-i))*(((REAL(wk)[ii+k*P] + 
								 REAL(dkp)[ii+k*P])/
								R_pow(REAL(s)[ii+signal*P],2.0))+
							       log(M_2PI*
								   R_pow(REAL(s)[ii+signal*P],
									 2.0))) + 
			log(REAL(prop)[ii]);
		}
		k++;
	    }
	    
	    /* Calculation of Amax */
	    for (ii=0; ii<k; ii++) {
		REAL(Amax)[ii] = REAL(A)[P*ii];
		for (j=0; j<P; j++) {
		    if (REAL(Amax)[ii] < REAL(A)[j+P*ii]) {
			REAL(Amax)[ii] = REAL(A)[j+P*ii];
		    }
		}
	    }
	    
	    /* Subtracts Amax */
	    for (ii=0; ii<k; ii++) {
		for (j=0; j<P; j++) {
		    REAL(A)[j+P*ii] = exp(REAL(A)[j+P*ii] - REAL(Amax)[ii]);
		}
	    }
	    
	    
	    /* Stores the results in G */
	    for (ii=0; ii<k; ii++) {
		REAL(Asu)[ii] = REAL(A)[P*ii];
		for (j=1; j<P; j++) {
		    REAL(Asu)[ii] += REAL(A)[j+P*ii];
		}
	    }
	    k = 0;
	    for (ii=lmin-1+i; ii<n; ii++) {
		REAL(G)[i+n*ii] += -log(REAL(Asu)[k]) - REAL(Amax)[k];
		k++;
	    }
	}
    }

    for (i=lmin-2; i<n; i++) {
	for (j=0; j<=i; j++) {
	    REAL(G)[i+n*j] = R_PosInf;
	}
    }
    
    /* free memory and return */
    UNPROTECT(8);
    return(G);
    
}




/* Dynamic programming: From a cost matrix matD and a maximum number of segments
   Kmaxr, this function returns:
   * a vector with Kmax values, containing the Kth is the minimum contrast 
   for a model with K segments (-J.est is the log-likelihood) 
   * a matrix, line K are the coordinates of the change points for a model with K segments 
   (rupture in the series)
 */
SEXP DynProg(SEXP matD, SEXP Kmaxr)
{
    int N, Kmax, i, k, L, imin;
    SEXP Im, tm, test, resu, Jest;
    double tmp, tmp2;
    
    /* Initialization */
    N = sqrt(length(matD));
    Kmax = INTEGER(Kmaxr)[0];
    tmp=0.0;
    tmp2=0.0;
    
    PROTECT(Im = allocMatrix(REALSXP, Kmax, N));
    PROTECT(tm = allocMatrix(INTSXP, Kmax, N));
    PROTECT(test = allocMatrix(INTSXP, Kmax, Kmax));
    PROTECT(resu = allocVector(VECSXP,2));
    PROTECT(Jest = allocVector(REALSXP,Kmax));
    
    for (i = 0; i < (N*Kmax); i++) {
	REAL(Im)[i] = R_PosInf;
	INTEGER(tm)[i] = 0;
    }
    
    for (i=0; i<N; i++) {
	REAL(Im)[Kmax*i] = REAL(matD)[i*N];
    }

    for (i=0; i<(Kmax*Kmax); i++) {
	INTEGER(test)[i] = 0;
    }
    
    if (Kmax > 2) {
	
	for (k = 1; k < (Kmax-1); k++) {
	    for (L = k; L<N; L++) {
		tmp = REAL(Im)[k-1]+REAL(matD)[(N*L)+1];
		imin = 0;
		for (i=1; i<(L-1); i++) {
		    tmp2 = REAL(Im)[k-1+i*Kmax]+REAL(matD)[(N*L)+1+i];
		    if (tmp2<tmp) {
			tmp=tmp2;
			imin=i;
		    }
		}
		REAL(Im)[k+L*Kmax] = tmp;
		if (tmp < 1e16) {
		    INTEGER(tm)[(k-1)+L*Kmax] = ((double) imin);
		} else {
		    INTEGER(tm)[(k-1)+L*Kmax] = -1;
		}
	    }
	}
    }
    

    tmp = REAL(Im)[Kmax-2]+REAL(matD)[((N-1)*N)+1];
    imin = 0;
    for (i=1; i<(N-1); i++) {
	tmp2 = REAL(Im)[Kmax-2+i*Kmax]+REAL(matD)[((N-1)*N)+1+i];
	if (tmp2<tmp) {
	    tmp=tmp2;
	    imin=i;
	}
    }
    REAL(Im)[(Kmax-1)+(N-1)*Kmax] = tmp;
    if (tmp < 1e16) {
	INTEGER(tm)[(Kmax-2)+(N-1)*Kmax] = imin;
    }
    

    for (i=0; i<Kmax; i++) {
	INTEGER(test)[i+i*Kmax] = ((double) N-1);
    }
    
    for (k = 1; k<Kmax; k++) {
	for (i=k-1; i>=0; i--) {
	    if (INTEGER(test)[k+(i+1)*Kmax] > 0) {
		INTEGER(test)[k+(i*Kmax)] = INTEGER(tm)[i+(INTEGER(test)[k+(i+1)*Kmax])*Kmax];
	    }
	}
    }
    
    for (i=0; i<Kmax; i++) {
	REAL(Jest)[i] = REAL(Im)[(N-1)*Kmax + i];
    }

    SET_VECTOR_ELT(resu, 0, Jest);
    SET_VECTOR_ELT(resu, 1, test);
    UNPROTECT(5);

    return(resu);    
}




/* Initialization of the EM algorithm with a hierarchical clustering. 
   Takes a bivariate series x, a matrix rupt containing the coordinates 
   of the ruptures between segments, a number of segments
   Kr and a number of possible model Pr,
   and returns a list phi with starting values for the mean, the standard deviation and
   the probability of each model (see the description of Gmixt_simultanee for a detailed
   description of these elements.
 */
SEXP EMInitSimultanee(SEXP x, SEXP rupt, SEXP Kr, SEXP Pr)
{
    SEXP m, v, n, Dist, Dist2, m2, v2, Dtmp, mfin, vfin, propfin, phi0;
    int K, P, i, j, k, r, nrx, indice, imin, jmin, ntmp, ic, jc, nrcur;
    double mo1, mo2, va1, va2, nt, ybar1, ybar2, varpool, Dmin, mtmp1, mtmp2, vtmp1, vtmp2;

    /* === Initialization === */
    K = INTEGER(Kr)[0];
    P = INTEGER(Pr)[0];
    nrx = length(x)/2;

    mo1 = 0.0;
    mo2 = 0.0;
    va1 = 0.0;
    va2 = 0.0;
    nt = 0;
    Dmin = 0;
    nrcur =K;
    
    /* === Memory Allocation === */
    /* Vectors of mean, variance and sizes */
    PROTECT(m = allocVector(REALSXP, K*2));
    PROTECT(v = allocVector(REALSXP, K*2));
    PROTECT(n = allocVector(INTSXP, K));

    /* idem, for temporary storage */
    PROTECT(m2 = allocVector(REALSXP, K*2));
    PROTECT(v2 = allocVector(REALSXP, K*2));
    
    /* Distance matrices */
    PROTECT(Dist = allocVector(REALSXP, K*K));
    PROTECT(Dist2 = allocVector(REALSXP, K*K));
    PROTECT(Dtmp = allocVector(REALSXP, K));
    
    /* End results */
    PROTECT(mfin = allocMatrix(REALSXP, P,2));
    PROTECT(vfin = allocMatrix(REALSXP, P,2));
    PROTECT(propfin = allocVector(REALSXP, P));
    PROTECT(phi0 = allocVector(VECSXP, 3));
    
    /* === Calculates the mean, variance and size of the two variables in segments === */
    for (i=0; i<K; i++) {

	/* The mean and segment size */
	mo1 = 0.0;
	mo2 = 0.0;
	nt = 0;
	for (j=INTEGER(rupt)[i]; j<=INTEGER(rupt)[i+K]; j++) {
	    mo1 += REAL(x)[j];
	    mo2 += REAL(x)[j+nrx];
	    nt++;
	}
	mo1 = mo1/((double) nt);
	mo2 = mo2/((double) nt);

	/* The variance */
	va1 = 0.0;
	va2 = 0.0;
	for (j=INTEGER(rupt)[i]; j<=INTEGER(rupt)[i+K]; j++) {
	    va1 += R_pow(REAL(x)[j] - mo1,2.0);
	    va2 += R_pow(REAL(x)[j+nrx] - mo2, 2.0);
	}
	va1 = va1/((double) (nt-1)); /* xxxx to be changed in the end: replace by /nt */
	va2 = va2/((double) (nt-1));
	
	/* Store the results */
	REAL(m)[i] = mo1;
	REAL(m)[i+K] = mo2;
	REAL(v)[i] = va1;
	REAL(v)[i+K] = va2;
	INTEGER(n)[i] = nt;
    }
    
    
    /* === Calculate the distance matrix === */
    for (i=0; i<(K*K); i++) {
	REAL(Dist)[i] =  R_PosInf;
    }
    /* For every pair (k,r) of segments, with r>k, 
       calculate the distance between the segments */
    for (k = 0; k<K-1; k++) {
	for (r = k+1; r<K; r++) {

	    /* Mean of the two variables per "merged" segment */
	    ybar1 = (((double) INTEGER(n)[k])*(REAL(m)[k]) + 
		     ((double) INTEGER(n)[r])*(REAL(m)[r]))/
		((double) (INTEGER(n)[k] + INTEGER(n)[r]));
	    ybar2 = (((double) INTEGER(n)[k])*(REAL(m)[k+K]) + 
		     ((double) INTEGER(n)[r])*(REAL(m)[r+K]))/
		((double) (INTEGER(n)[k] + INTEGER(n)[r]));
	    
	    /* Pooled variance for "merged" segment */
	    varpool = ((double) INTEGER(n)[r])*(REAL(v)[r]) + 
		((double) INTEGER(n)[k])*(REAL(v)[k]) + 
		((double) INTEGER(n)[r])*R_pow(REAL(m)[r] - ybar1, 2.0) + 
		((double) INTEGER(n)[k])*R_pow(REAL(m)[k] - ybar1, 2.0);
	    varpool += ((double) INTEGER(n)[r])*(REAL(v)[r+K]) + 
		((double) INTEGER(n)[k])*(REAL(v)[k+K]) + 
		((double) INTEGER(n)[r])*R_pow(REAL(m)[r+K] - ybar2, 2.0) + 
		((double) INTEGER(n)[k])*R_pow(REAL(m)[k+K] - ybar2, 2.0);
	    varpool = varpool/((double) (INTEGER(n)[k] + INTEGER(n)[r]));
	    
	    /* The distance between segments r and k */
	    REAL(Dist)[k+r*K] = log(varpool)*((double) (INTEGER(n)[k] + 
							INTEGER(n)[r])) +
		(-((double) INTEGER(n)[k])*log(REAL(v)[k]) -
		 ((double) INTEGER(n)[r])*log(REAL(v)[r])) +
		(-((double) INTEGER(n)[k])*log(REAL(v)[k+K]) -
		 ((double) INTEGER(n)[r])*log(REAL(v)[r+K]));
	}
    }
     
    /* main loop (merging the closest segments together until there remains only P groups)  */
    if (K>P) {

	for (indice=0; indice<K; indice++) {
	    
	    /* Checks that the current number of "merged segments" is still greater than P
	       (in this cases, reduces it further, otherwise break)
	     */
	    if (nrcur <= P) 
		break;
	    
	    /* Finds minimum distance between two segments (and corresponding
	       indices) */
	    imin = 0;
	    jmin = 0;
	    Dmin = REAL(Dist)[0];
	    for (i=0; i<nrcur; i++) {
		for (j = 0; j<nrcur; j++) {
		    if (REAL(Dist)[i+j*K] < Dmin) {
			Dmin = REAL(Dist)[i+j*K];
			imin=i;
			jmin=j;
		    }
		}
	    }
		
	    /* Calculates the size, mean and variance of the grouped segments */
	    ntmp = INTEGER(n)[imin] + INTEGER(n)[jmin];
	    mtmp1 = (((double) INTEGER(n)[imin])*REAL(m)[imin] + 
		     ((double) INTEGER(n)[jmin])*REAL(m)[jmin])/((double) ntmp);
	    mtmp2 = (((double) INTEGER(n)[imin])*REAL(m)[imin+K] + 
		     ((double) INTEGER(n)[jmin])*REAL(m)[jmin+K])/((double) ntmp);
	    vtmp1 = (((double) INTEGER(n)[imin])*REAL(v)[imin] + 
		     ((double) INTEGER(n)[jmin])*REAL(v)[jmin] +
		     ((double) INTEGER(n)[imin])*(R_pow(REAL(m)[imin] - mtmp1,2.0)) +
		     ((double) INTEGER(n)[jmin])*(R_pow(REAL(m)[jmin] - mtmp1,2.0)))/
		((double) ntmp);
	    vtmp2 = (((double) INTEGER(n)[imin])*REAL(v)[imin+K] + 
		     ((double) INTEGER(n)[jmin])*REAL(v)[jmin+K] +
		     ((double) INTEGER(n)[imin])*(R_pow(REAL(m)[imin+K] - mtmp2,2.0)) +
		     ((double) INTEGER(n)[jmin])*(R_pow(REAL(m)[jmin+K] - mtmp2,2.0)))/
		((double) ntmp);
	    
	    
	    /* Removes the rows and columns of the distance matrix
	       corresponding to the two segments with minimum distance.
	    */
	    ic = 0;
	    jc = 0;

	    for (i=0; i < nrcur; i++) {
		if ((i!=imin)&&(i!=jmin)) {
		    jc = 0;
		    for (j=0; j<nrcur; j++) {
			if ((j!=imin)&&(j!=jmin))  {
			    REAL(Dist2)[ic+K*jc] = REAL(Dist)[i+K*j];
			    jc++;
			}
		    }
		    ic++;
		}
	    }

	    /* removes the means and variance and size corresponding to these
	       two segments
	     */
	    ic = 0;
	    for (i=0; i<nrcur; i++) {
		if ((i!=imin)&&(i!=jmin)) {
		    for (j=0; j<2; j++) {
			REAL(m2)[ic+K*j] = REAL(m)[i+K*j];
			REAL(v2)[ic+K*j] = REAL(v)[i+K*j];
		    }
		    INTEGER(n)[ic] = INTEGER(n)[i];
		    ic++;
		}
	    }
	    
	    /* Adds the mean, variance and size for the merged group */
	    REAL(m2)[ic] = mtmp1;
	    REAL(m2)[ic+K] = mtmp2;
	    REAL(v2)[ic] = vtmp1;
	    REAL(v2)[ic+K] = vtmp2;
	    INTEGER(n)[ic] = ntmp;
	    
	    /* Recalculates the distance vector for the merged group */
	    for (k=0; k<K; k++) {
		REAL(Dtmp)[k] = R_PosInf;
	    }
	    /* For every segment k, calculates the distance with this merged
	       segment */
	    for (k=0; k< (nrcur-2); k++) {
		
		/* Mean for the overall segment (merged + segment k) */
		ybar1 = (((double) INTEGER(n)[k])*(REAL(m2)[k]) + 
			 ((double) ntmp)*(mtmp1))/
		    (((double) INTEGER(n)[k]) + ((double) ntmp));
		ybar2 = (((double) INTEGER(n)[k])*(REAL(m2)[k+K]) + 
			 ((double) ntmp)*(mtmp2))/
		    (((double) INTEGER(n)[k]) + ((double) ntmp));
		
		/* Pooled variance for the overall segment */
		varpool = ((double) ntmp)*(vtmp1) + ((double) INTEGER(n)[k])*(REAL(v2)[k]) + 
		    ((double) ntmp)*R_pow(mtmp1 - ybar1, 2.0) + 
		    ((double) INTEGER(n)[k])*R_pow(REAL(m2)[k] - ybar1, 2.0);
		varpool += ((double) ntmp)*(vtmp2) + ((double) INTEGER(n)[k])*(REAL(v2)[k+K]) + 
		    ((double) ntmp)*R_pow(mtmp2 - ybar2, 2.0) + 
		    ((double) INTEGER(n)[k])*R_pow(REAL(m2)[k+K] - ybar2, 2.0);
		varpool = varpool/((double) (INTEGER(n)[k] + ntmp));
	    
		/* Distance between the merged segment and segment k */
		REAL(Dtmp)[k] = log(varpool)*(((double) INTEGER(n)[k]) + ((double) ntmp)) +
		    (-((double) INTEGER(n)[k])*log(REAL(v2)[k]) -
		     ((double) ntmp)*log(vtmp1)) +
		    (-((double) INTEGER(n)[k])*log(REAL(v2)[k+K]) -
		     ((double) ntmp)*log(vtmp2));
	    }
	    
	    /* Resets Dist = Inf */
	    for (i=0; i<(K*K); i++) {
		REAL(Dist)[i] = R_PosInf;
	    }
	    
	    /* adds the distances in Dist2 */
	    for (i=0; i<(nrcur-2); i++) {
		for (j=0; j<(nrcur-2); j++) {
		    REAL(Dist)[i+j*K] = REAL(Dist2)[i+j*K];
		}
	    }
	    
	    /* adds the distances in the merged group */
	    for (i=0; i<(nrcur-2); i++) {
		REAL(Dist)[i+(nrcur-2)*K] = REAL(Dtmp)[i];
	    }
	    
	    /* Updates the mean and variance */
	    for (i=0; i<(nrcur-1); i++) {
		REAL(m)[i] = REAL(m2)[i];
		REAL(m)[i+K] = REAL(m2)[i+K];
		REAL(v)[i] = REAL(v2)[i];
		REAL(v)[i+K] = REAL(v2)[i+K];
	    }
	    
	    nrcur--;
	}
    }
    
    /* Prepare the end results (output) */
    j=0;
    for (i=0; i<P; i++) {
	j+=INTEGER(n)[i];
    }
    for (i=0; i<P; i++) {
	for (k=0; k<2; k++) {
	    REAL(mfin)[i+P*k] = REAL(m)[i+K*k];
	    REAL(vfin)[i+P*k] = sqrt(REAL(v)[i+K*k]);
	}
	REAL(propfin)[i] = ((double) INTEGER(n)[i])/((double) j);
    }
    
    /* Put everything in a list */
    SET_VECTOR_ELT(phi0, 0, mfin);
    SET_VECTOR_ELT(phi0, 1, vfin);
    SET_VECTOR_ELT(phi0, 2, propfin);
    
    /* Free memory */
    UNPROTECT(12);

    /* Return results */
    return(phi0);

}




/* The E step of the EM algorithm: computes posterior probabilities and 
   incomplete-data log-likelihood for mixture models.
   Takes as argument a K*P matrix logdensity containing the conditional log-densities for
   each segment and each model
   and phi, a list containing the parameter of the mixture (see Gmixt_simultanee for a 
   description of the three elements of the list: mean, sigma and prop).
   
   a list with tau a K*P matrix, tau kj is the posterior probability for 
   segment k to belong to classe j and lvinc, the incomplete log vrais P(X=x)
*/
SEXP EstepSimultanee(SEXP logdensity, SEXP phi)
{
    int K, P, i, j;
    SEXP tau, tauMax, lvinc, resu, prop;
    double tmp;
    
    /* Initialization and memory allocation */
    K = nrows(logdensity);
    P = ncols(logdensity);
    PROTECT(tau = allocMatrix(REALSXP, K, P));
    PROTECT(tauMax = allocVector(REALSXP, K));
    PROTECT(lvinc = allocVector(REALSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 2));
    PROTECT(prop = VECTOR_ELT(phi, 2));
    
    /* Fill tau with log( pi*f ) */
    for (i = 0; i<K; i++) {
	for (j = 0; j<P; j++) {
	    REAL(tau)[i+j*K] = log(REAL(prop)[j]) + REAL(logdensity)[i+j*K];
	}
    }
    
    /* identify max value for each row (needed to avoid to small 
       probas when exponentiating) */
    for (i = 0; i<K; i++) {
	REAL(tauMax)[i] = REAL(tau)[i];
	if (P>1) {
	    for (j = 1; j<P; j++) {
		if (REAL(tau)[i+j*K] > REAL(tauMax)[i])
		    REAL(tauMax)[i] = REAL(tau)[i+j*K];
	    }
	}
    }

    /* subtract max value for each row and exponential + calculation of lvinc + 
       standardization of tau
     */
    REAL(lvinc)[0] = 0.0;
    for (i = 0; i<K; i++) {
	tmp = 0.0;
	for (j = 0; j< P; j++) {
	    REAL(tau)[i+j*K] = exp(REAL(tau)[i+j*K] - REAL(tauMax)[i]);
	    tmp+= REAL(tau)[i+j*K];
	}
	REAL(lvinc)[0] += (log(tmp) + REAL(tauMax)[i]);
	for (j = 0; j< P; j++) {
	    REAL(tau)[i+j*K] = REAL(tau)[i+j*K]/tmp;
	}
    }
    
    /* Store the results in a list */
    SET_VECTOR_ELT(resu, 0, tau);
    SET_VECTOR_ELT(resu, 1, lvinc);

    /* Free memory and return */
    UNPROTECT(5);
    return(resu);
    
}



/* M step of the EM algorithm. It takes as argument the bivariate signal x,
   a matrix rupt containing the limits (start, end) of each segment (rows),
   a matrix tau containing the posterior probabilities of membership to clusters,
   It returns the updated value of the parameters
 */
SEXP MstepSimultanee(SEXP x, SEXP rupt, SEXP tau)
{
    int K, P, i, j, g, nrx, *indx;
    SEXP m, s, prop, Yk, nk, np, resu, inds, m1, mfin, sfin, propfin;
    double tmp, tmp2;

    /* Initialization and memory allocation */    
    K = nrows(tau);
    P = ncols(tau);
    nrx = nrows(x);
    
    PROTECT(m = allocMatrix(REALSXP, P, 2));
    PROTECT(s = allocMatrix(REALSXP, P, 2));
    PROTECT(prop = allocVector(REALSXP, P));
    PROTECT(mfin = allocMatrix(REALSXP, P, 2));
    PROTECT(sfin = allocMatrix(REALSXP, P, 2));
    PROTECT(propfin = allocVector(REALSXP, P));
    PROTECT(Yk = allocVector(REALSXP, 2*K));
    PROTECT(nk = allocVector(INTSXP, K));
    PROTECT(np = allocVector(REALSXP, P));
    PROTECT(m1 = allocVector(REALSXP, P));
    PROTECT(inds = allocVector(INTSXP, P));
    PROTECT(resu = allocVector(VECSXP, 3));
    indx = INTEGER(inds);
    
    /* Calculation of Sum(x) for each variable and each segment 
       + calculation of the size of each segment
     */
    for (i=0; i<K; i++) {
	REAL(Yk)[i] = 0;
	REAL(Yk)[i+K] = 0;
	INTEGER(nk)[i] = 0;
	for (j=INTEGER(rupt)[i]; j<=INTEGER(rupt)[i+K]; j++) {
	    REAL(Yk)[i] += REAL(x)[j];
	    REAL(Yk)[i+K] += REAL(x)[j+nrx];
	    INTEGER(nk)[i]++;
	}
    }

    /* Calculation of np */
    for (i=0; i<P; i++) {
	REAL(np)[i] = 0.0;
	for (j=0; j<K; j++) {
	    REAL(np)[i] += ((double) INTEGER(nk)[j])*REAL(tau)[j+i*K];
	}
    }

    /* Calculation of m */
    for (i=0; i<P; i++) {
	REAL(m)[i] = 0.0;
	REAL(m)[i+P] = 0.0;
	for (j=0; j<K; j++) {
	    REAL(m)[i] += REAL(Yk)[j]*REAL(tau)[j+i*K]/REAL(np)[i];
	    REAL(m)[i+P] += REAL(Yk)[j+K]*REAL(tau)[j+i*K]/REAL(np)[i];
	}
    }
    
    /* Calculation of s */
    for (g=0; g<P; g++) {

	REAL(s)[g] = 0.0;
	REAL(s)[g+P] = 0.0;
	for (i=0; i<K; i++) {

	    tmp = 0.0;
	    tmp2 = 0.0;
	    for (j=INTEGER(rupt)[i]; j<=INTEGER(rupt)[i+K]; j++) {
		tmp += REAL(tau)[i+g*K] * R_pow(REAL(x)[j]-REAL(m)[g], 2.0);
		tmp2 += REAL(tau)[i+g*K] * R_pow(REAL(x)[j+nrx]-REAL(m)[g+P], 2.0);
	    }
	    REAL(s)[g] += tmp/(REAL(np)[g]);
	    REAL(s)[g+P] += tmp2/(REAL(np)[g]);
	}
	REAL(s)[g] = sqrt(REAL(s)[g]);
	REAL(s)[g+P] = sqrt(REAL(s)[g+P]);
    }

    /* Calculation of prop */
    for (i=0; i<P; i++) {
	REAL(prop)[i] = 0;
	for (j=0; j<K; j++) {
	    REAL(prop)[i] += REAL(tau)[j+i*K]/((double) K);
	}
    }

    /* Ordering as a function of the mean of the first variable */
    for (i=0; i<P; i++) {
	indx[i] = i;
	REAL(m1)[i] = REAL(m)[i];
    }
    R_orderVector(indx, length(inds), Rf_lang1(m1), TRUE, FALSE);
    for (i=0; i<P; i++) {
	REAL(mfin)[i] = REAL(m)[indx[i]];
	REAL(mfin)[i+P] = REAL(m)[indx[i]+P];
	REAL(sfin)[i] = REAL(s)[indx[i]];
	REAL(sfin)[i+P] = REAL(s)[indx[i]+P];
	REAL(propfin)[i] = REAL(prop)[indx[i]];
    }


    /* Puts results in resu and exits */
    SET_VECTOR_ELT(resu, 0, mfin);
    SET_VECTOR_ELT(resu, 1, sfin);
    SET_VECTOR_ELT(resu, 2, propfin);
    
    UNPROTECT(12);
    return(resu);    
}




/* EM algorithm: combines the E-step and M-step , and calculates the MLE of phi
   for a given change-point instant and a fixed number of clusters.
   It takes as argument the bivariate signal x,
   a matrix rupt containing the limits (start, end) of each segment (rows)
   a list phi containing the parameters of the mixture (see the help of 
   Gmixt_simultanee and previous functions for a description of phi), and
   the number of possible models Pr.
   return a list with  phi, the MLE, tau =(taukj) the probability for 
   segment k to belong to classe,lvinc = lvinc,empty = empty,dv = dv
 */
SEXP EMAlgoSimultane(SEXP x, SEXP rupt, SEXP Pr, SEXP phi)
{
    int K, i, j, P;
    double delta, eps, minp;
    SEXP resu, tau, np, phi2, phitmp, mu, sd, prop, ldt, logdensity, resEstep, lvinc, empty;
    SEXP dv, iter, tmp, tmp2;
    
    K = length(rupt)/2;
    P = INTEGER(Pr)[0];
    delta = 1;
    eps = 1e-10;
    
    PROTECT(tau = allocMatrix(REALSXP, K, P));
    PROTECT(lvinc = allocVector(REALSXP, 1));
    PROTECT(np = allocVector(REALSXP, P));
    PROTECT(logdensity = allocMatrix(REALSXP, K, P));
    PROTECT(mu = allocMatrix(REALSXP, P, 2));
    PROTECT(sd = allocMatrix(REALSXP, P, 2));
    PROTECT(prop = allocVector(REALSXP, P));
    PROTECT(phi2 = allocVector(VECSXP, 3));
    PROTECT(empty = allocVector(INTSXP, 1));
    PROTECT(dv = allocVector(INTSXP, 1));
    PROTECT(iter = allocVector(INTSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 5));


    INTEGER(empty)[0] = 0;
    INTEGER(dv)[0] = 0;
    INTEGER(iter)[0] = 0;

    /* Local copy of phi */
    for (i=0; i<P; i++) {
	REAL(mu)[i] = REAL(VECTOR_ELT(phi,0))[i];
	REAL(mu)[i+P] = REAL(VECTOR_ELT(phi,0))[i+P];
	REAL(sd)[i] = REAL(VECTOR_ELT(phi,1))[i];
	REAL(sd)[i+P] = REAL(VECTOR_ELT(phi,1))[i+P];
	REAL(prop)[i] = REAL(VECTOR_ELT(phi,2))[i];
    }
    SET_VECTOR_ELT(phi2, 0, mu);
    SET_VECTOR_ELT(phi2, 1, sd);
    SET_VECTOR_ELT(phi2, 2, prop);
    
    
    for (i = 0; i < K*P; i++) {
	REAL(tau)[i] = 1.0;
    }
    for (i = 0; i < P; i++) {
	REAL(np)[i] = K;
    }
    minp = K;

    while ( (delta >= 1e-4) && (minp > eps) && (INTEGER(iter)[0] <= 5000) ) {

	R_CheckUserInterrupt();

	INTEGER(iter)[0] = INTEGER(iter)[0]+1;

	/* Calculation of logdensity */
	for (i=0; i<K; i++) {
	    PROTECT(ldt = logdensSimultaneeParIndice(x, phi2, INTEGER(rupt)[i], 
						     INTEGER(rupt)[i+K]));
	    for (j=0; j<P; j++) {
		REAL(logdensity)[i+K*j] = REAL(ldt)[j];
	    }
	    UNPROTECT(1);
	}

	/* E-step */
	PROTECT(resEstep=EstepSimultanee(logdensity, phi2));
	tmp = VECTOR_ELT(resEstep,1);
	REAL(lvinc)[0] = REAL(tmp)[0];
	tmp = VECTOR_ELT(resEstep,0);
	for (i=0; i<K; i++) {
	    for (j=0; j<P; j++) {
		REAL(tau)[i+j*K] = REAL(tmp)[i+j*K];
	    }
	}
	
	
	/* M-step */
	PROTECT(phitmp = MstepSimultanee(x, rupt, tau));
	
	/* Recalculaltes minp */
	for (i=0; i<P; i++) {
	    REAL(np)[i] = 0.0;
	    for (j=0; j<K; j++) {
		REAL(np)[i] += REAL(tau)[j+i*K];
	    }
	}
	minp = REAL(np)[0];
	for (i=1; i<P; i++) {
	    if (REAL(np)[i] < minp)
		minp=REAL(np)[i];
	}
	
	/* stores the results in mu, sd and prop, and calculates delta at the same time */
	delta=0;
	for (i=0; i<P; i++) {
	    
	    /* Recalculates delta for mu*/
	    tmp = VECTOR_ELT(phitmp,0);
	    tmp2 = VECTOR_ELT(phi2,0);
	    if (fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]) > delta) {
		delta = fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]);
	    }
	    if (fabs((REAL(tmp)[i+P] - REAL(tmp2)[i+P])/REAL(tmp2)[i+P]) > delta) {
		delta = fabs((REAL(tmp)[i+P] - REAL(tmp2)[i+P])/
			     REAL(tmp2)[i+P]);
	    }
	    REAL(mu)[i] = REAL(tmp)[i];
	    REAL(mu)[i+P] = REAL(tmp)[i+P];


	    /* Recalculates delta for sd */
	    tmp = VECTOR_ELT(phitmp,1);
	    tmp2 = VECTOR_ELT(phi2,1);
	    if (fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]) > delta) {
		delta = fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]);
	    }
	    if (fabs((REAL(tmp)[i+P] - REAL(tmp2)[i+P])/REAL(tmp2)[i+P]) > delta) {
		delta = fabs((REAL(tmp)[i+P] - REAL(tmp2)[i+P])/
			     REAL(tmp2)[i+P]);
	    }
	    REAL(sd)[i] = REAL(tmp)[i];
	    REAL(sd)[i+P] = REAL(tmp)[i+P];

	    /* Recalculates delta for tmp */
	    tmp = VECTOR_ELT(phitmp,2);
	    tmp2 = VECTOR_ELT(phi2,2);
	    if (fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]) > delta) {
		delta = fabs((REAL(tmp)[i] - REAL(tmp2)[i])/REAL(tmp2)[i]);
	    }
	    REAL(prop)[i] = REAL(tmp)[i];

	    
	}
	
	/* Sets again the elements of phi2 */
	SET_VECTOR_ELT(phi2, 0, mu);
	SET_VECTOR_ELT(phi2, 1, sd);
	SET_VECTOR_ELT(phi2, 2, prop);
	
	/* UNPROTECTS the stuffs */
	UNPROTECT(2);
    }

    /* Checks the results of the algorithm */
    if (minp < eps) {
	INTEGER(empty)[0]=1;
	REAL(lvinc)[0] =  R_NegInf;
    }
    if (INTEGER(iter)[0]>5000) {
	INTEGER(dv)[0]=2;
	REAL(lvinc)[0] = R_NegInf;
    }
        
    /* Return the results */
    SET_VECTOR_ELT(resu, 0, phi2);
    SET_VECTOR_ELT(resu, 1, tau);
    SET_VECTOR_ELT(resu, 2, lvinc);
    SET_VECTOR_ELT(resu, 3, empty);
    SET_VECTOR_ELT(resu, 4, dv);    

    UNPROTECT(12);
    return(resu);
}




/* The function neighbors is the C equivalent of the R function neighbors. 
   It takes as arguments:
   - the bivariate series x
   - The number of classe P
   - the list entree containing the result of the segmentation (the object resu returned
   by the hybrid algorithm). This list has two elements:
   -- the first element (index 0) is a Kmax-vector containing the log-likelihood associated
      to each number of segment (Linc)
   -- the second element (index 1) is a Kmax-list containing the parameters of the fit (with
      phi, rupt and tau)
   - k is the value of k which is to be smoothed according to the neighbors. WARNING: k here is
     an index (min 0, max K-1)
   - lminr is the minimum size of a segment
   Returns:
   A list similar to entree, but with the results stored in the k-th element (i.e. the 
   k-th element of the first element is the log-likelihood, and the k-th element of the
   second element contain the parameters).
 */
SEXP neighbors(SEXP x, SEXP Pr, SEXP entree, SEXP kr, SEXP lminr) 
{
    int Kmax, i, k, K1, K2;
    double a, tmp, tmp1, tmp2;
    SEXP phi1, phi2, G, outDP, rupt1, outEM1, rupt2, outEM2, paramf, Lres, paramv,kk, resu, L, param;
    
    Kmax = length(VECTOR_ELT(entree,1));
    k = INTEGER(kr)[0];
    K1=-1;
    K2=-1;
    a=R_NegInf;
    
    /* Checks if the value of k is null */
    if (isNull(VECTOR_ELT(VECTOR_ELT(entree,1), k))) {
	return(entree);
    }
    
    PROTECT(L = VECTOR_ELT(entree, 0));
    PROTECT(param = VECTOR_ELT(entree, 1));    
    PROTECT(rupt1 = allocMatrix(INTSXP, k+1, 2));
    PROTECT(rupt2 = allocMatrix(INTSXP, k+1, 2));
    PROTECT(paramv = allocVector(VECSXP, 3));
    PROTECT(paramf = allocVector(VECSXP, 2));
    PROTECT(resu = allocVector(VECSXP, 2));
    PROTECT(kk = allocVector(INTSXP, 1));
    
    /* Left neighbor */
    if (k > 0) {
	a = REAL(L)[0];
	K1 = 0;
	for (i = 0; i<=k; i++) {
	    if (a < REAL(L)[i]) {
		a = REAL(L)[i];
		K1 = i;		
	    }
	}
    }

    if ((k==0)||(a < -(1e16))||ISNA(a)) {
	K1 = -1;
    } else {
	PROTECT(phi1 = VECTOR_ELT(VECTOR_ELT(param,K1),0));
	PROTECT(G = GmixtSimultanee(x, lminr, phi1));
	INTEGER(kk)[0] = k+1;
	PROTECT(outDP = DynProg(G,kk));/* à priori +1, car valeur et pas indice */
	
	/* Updates rupt */
	
	INTEGER(rupt1)[0] =0;
	INTEGER(rupt1)[k+1] = INTEGER(VECTOR_ELT(outDP,1))[k];

	for (i = 1; i<(k+1); i++) {
	    INTEGER(rupt1)[i] = INTEGER(VECTOR_ELT(outDP,1))[k+(i-1)*(k+1)]+1;
	    INTEGER(rupt1)[i+k+1] = INTEGER(VECTOR_ELT(outDP,1))[k+i*(k+1)];
	}
	
	
	/* EM */
	PROTECT(outEM1 = EMAlgoSimultane(x, rupt1, Pr, phi1));
    }


    /* right neighbor */
    if (k < (Kmax-1)) {
	a = REAL(L)[k+1];
	K2 = k+1;
	for (i = k+1; i<Kmax; i++) {
	    if (a < REAL(L)[i]) {
		a = REAL(L)[i];
		K2 = i;		
	    }
	}
    }

    if ((k==(Kmax-1))||(a < -(1e16))||ISNA(a)) {
	K2 = -1;
    } else {
	PROTECT(phi2 = VECTOR_ELT(VECTOR_ELT(param,K2),0));
	PROTECT(G = GmixtSimultanee(x, lminr, phi2));
	INTEGER(kk)[0] = k+1;
	PROTECT(outDP = DynProg(G,kk));/* à priori +1, car valeur et pas indice */
	
	/* Updates rupt */
	INTEGER(rupt2)[0] =0;
	INTEGER(rupt2)[k+1] = INTEGER(VECTOR_ELT(outDP,1))[k];
	for (i = 1; i<(k+1); i++) {
	    INTEGER(rupt2)[i] = INTEGER(VECTOR_ELT(outDP,1))[k+(i-1)*(k+1)]+1;
	    INTEGER(rupt2)[i+k+1] = INTEGER(VECTOR_ELT(outDP,1))[k+i*(k+1)];
	}
	
	/* EM */
	PROTECT(outEM2 = EMAlgoSimultane(x, rupt2, Pr, phi2));
    }
    

    /* Which is max: left, right or k? */
    a = 2;
    SET_VECTOR_ELT(paramv,0, VECTOR_ELT(VECTOR_ELT(param,k),0));
    SET_VECTOR_ELT(paramv,1, VECTOR_ELT(VECTOR_ELT(param,k),1));
    SET_VECTOR_ELT(paramv,2, VECTOR_ELT(VECTOR_ELT(param,k),2));
    tmp = REAL(L)[k];
    if (K1 > -1) {
	tmp1 = REAL(VECTOR_ELT(outEM1,2))[0];
    } else {
	tmp1 = R_NegInf;
    }
    if (K2 > -1) {
	tmp2 = REAL(VECTOR_ELT(outEM2,2))[0];
    } else {
	tmp2 = R_NegInf;
    }
    PROTECT(Lres=L);
        
    if (K1 > -1) {
	if ((tmp1 > tmp)&&(tmp1 > tmp2)) {
	    SET_VECTOR_ELT(paramv,0, VECTOR_ELT(outEM1, 0));
	    SET_VECTOR_ELT(paramv,1, rupt1);
	    SET_VECTOR_ELT(paramv,2, VECTOR_ELT(outEM1, 1));
	    REAL(Lres)[k] = tmp1;
	}
    }
    if (K2 > -1) {
	if ((tmp2 > tmp)&&(tmp2 > tmp1)) {
	    SET_VECTOR_ELT(paramv,0, VECTOR_ELT(outEM2, 0));
	    SET_VECTOR_ELT(paramv,1, rupt2);
	    SET_VECTOR_ELT(paramv,2, VECTOR_ELT(outEM2, 1));
	    REAL(Lres)[k] = tmp2;
	}
    }

    PROTECT(paramf=param);
    SET_VECTOR_ELT(paramf, k, paramv);    
    SET_VECTOR_ELT(resu,0,Lres);
    SET_VECTOR_ELT(resu,1,paramf);
    UNPROTECT(2);
    
    if (K2 > -1) {
	UNPROTECT(4);
    }
    if (K1 > -1) {
	UNPROTECT(4);
    }

    UNPROTECT(8);

    return(resu);
}



/* Hybrid algorithm, combining the EM algorithm and the dynamic programming for 
   segmentation.
   It is an algorithm which combines dynamic programming
   and the EM algorithm to calculate the MLE of phi and T, which 
   are the mixture parameters and the change point instants.
   this algorithm is run for a given number of clusters,
   and estimates the parameters for a segmentation/clustering
   model with P clusters and 1:Kmax segments
   Arguments:
   x the two-dimensionnal signal, one column per dimension
   Pr the number of classes 
   Kmaxr the maximal number of segments
   This is the function interacting with the R function picard. See
   this function for a description of the results (especially
   the function print.picard, which describes the content).
*/

SEXP hybridSimultanee(SEXP x, SEXP Pr, SEXP Kmaxr)
{
    int P, Kmax, Kmin, i, j, l, K, empty, dv, oukonenez;
    SEXP param, Linc, rupt, phi, mu, sigma, prop, lmin, G, out, Kr, outEM, tau;
    SEXP outDP, paramv, resu, nulle;
    double delta;
    
    P = INTEGER(Pr)[0];
    Kmax = INTEGER(Kmaxr)[0];
    Kmin = P;

    PROTECT(lmin = allocVector(INTSXP, 1));
    PROTECT(Kr = allocVector(INTSXP, 1));
    PROTECT(mu = allocMatrix(REALSXP, P, 2));
    PROTECT(sigma = allocMatrix(REALSXP, P, 2));
    PROTECT(prop = allocVector(REALSXP, P));
    PROTECT(Linc = allocVector(REALSXP, Kmax));
    PROTECT(param = allocVector(VECSXP, Kmax));
    PROTECT(resu = allocVector(VECSXP, 2));
    PROTECT(nulle = allocVector(NILSXP, 1));
    
    for (i=0;i<Kmax;i++) {
	REAL(Linc)[i]=R_NegInf;
	SET_VECTOR_ELT(param, i, nulle);
    }
    
    oukonenez = 0;
    if (P>1) { /*xxx gaffe: prévoir dans le code R la vérif du cas P=1 */

	INTEGER(lmin)[0] = 2;
	PROTECT(G = GmeanSimultanee(x, lmin));
	PROTECT(out = DynProg(G, Kmaxr));
	
	for (K = Kmin; K<=Kmax; K++) {
	    
	    /* Initialization */
	    PROTECT(paramv = allocVector(VECSXP, 3));
	    INTEGER(Kr)[0] = K;
	    j = 0;
	    delta = R_PosInf;
	    empty = 0;
	    dv = 0;
	    PROTECT(tau = allocMatrix(REALSXP, K, P));
	    
	    /* The segments limits */
	    PROTECT(rupt = allocMatrix(INTSXP, K, 2));
	    INTEGER(rupt)[0] =0;
	    INTEGER(rupt)[K] = INTEGER(VECTOR_ELT(out,1))[K-1];
	    
	    for (i = 1; i<K; i++) {
		INTEGER(rupt)[i] = INTEGER(VECTOR_ELT(out,1))[K-1+(i-1)*Kmax]+1;
		INTEGER(rupt)[i+K] = INTEGER(VECTOR_ELT(out,1))[K-1+i*Kmax];
	    }
	    
	    /* Initialization EM */
	    PROTECT(phi = EMInitSimultanee(x, rupt, Kr, Pr));
	    PROTECT(outEM = EMAlgoSimultane(x, rupt, Pr, phi));
	    
	    /* Copy in mu, sigma, prop, phi, tau */
	    for (i = 0; i < P; i++) {
		REAL(mu)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),0))[i];
		REAL(mu)[i+P] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),0))[i+P];
		REAL(sigma)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),1))[i];
		REAL(sigma)[i+P] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),1))[i+P];
		REAL(prop)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),2))[i];
		
		for (l = 0; l<K; l++) {
		    REAL(tau)[l+i*K] = REAL(VECTOR_ELT(outEM,1))[l+i*K];
		}
	    }
	    UNPROTECT(2);
	    PROTECT(phi = allocVector(VECSXP,3));
	    SET_VECTOR_ELT(phi,0,mu);
	    SET_VECTOR_ELT(phi,1,sigma);
	    SET_VECTOR_ELT(phi,2,prop);
	    

	    /* Loop for fitting */
	    while ((delta > 1e-4) && (empty==0) && (dv == 0) && (j <= 100)) {

		j++;
		
		/* Dynamic programming */
		PROTECT(G = GmixtSimultanee(x, lmin, phi));
		PROTECT(outDP = DynProg(G, Kr));

		/* Updates rupt */
		INTEGER(rupt)[0] =0;
		INTEGER(rupt)[K] = INTEGER(VECTOR_ELT(outDP,1))[K-1];
		for (i = 1; i<K; i++) {
		    INTEGER(rupt)[i] = INTEGER(VECTOR_ELT(outDP,1))[K-1+(i-1)*K]+1;
		    INTEGER(rupt)[i+K] = INTEGER(VECTOR_ELT(outDP,1))[K-1+i*K];
		}
		
		/* EM */
		PROTECT(outEM = EMAlgoSimultane(x, rupt, Pr, phi));
		
		/* Calculation of delta */
		delta=0;
		for (i=0; i<P; i++) {
		    for (l=0; l<2; l++) {

			/* Recalculates delta for mu*/
			if (fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),0))[i+l*P] - 
				  REAL(mu)[i+l*P])/REAL(mu)[i+l*P]) > delta) {
			    delta = fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),0))[i+l*P] - 
					  REAL(mu)[i+l*P])/REAL(mu)[i+l*P]);
			}
			
		    
			/* Recalculates delta for sd */
			if (fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),1))[i+l*P] - 
				  REAL(sigma)[i+l*P])/REAL(sigma)[i+l*P]) > delta) {
			    delta = fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),1))[i+l*P] - 
					  REAL(sigma)[i+l*P])/REAL(sigma)[i+l*P]);
			}
		    }
		    
		    /* Recalculates delta for prop */
		    if (fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),2))[i] - 
			      REAL(prop)[i])/REAL(prop)[i]) > delta) {
			delta = fabs((REAL(VECTOR_ELT(VECTOR_ELT(outEM, 0),2))[i] - 
				      REAL(prop)[i])/REAL(prop)[i]);
		    }
		}
		
		/* Other criteria */
		if (j==1)
		    delta=1;
		dv = INTEGER(VECTOR_ELT(outEM,4))[0];
		empty = INTEGER(VECTOR_ELT(outEM,3))[0];
		
		/* Update the elements of phi and tau */
		for (i = 0; i < P; i++) {
		    REAL(mu)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),0))[i];
		    REAL(mu)[i+P] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),0))[i+P];
		    REAL(sigma)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),1))[i];
		    REAL(sigma)[i+P] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),1))[i+P];
		    REAL(prop)[i] = REAL(VECTOR_ELT(VECTOR_ELT(outEM,0),2))[i];

		    for (l = 0; l<K; l++) {
			REAL(tau)[l+i*K] = REAL(VECTOR_ELT(outEM,1))[l+i*K];
		    }
		}
		SET_VECTOR_ELT(phi,0,mu);
		SET_VECTOR_ELT(phi,1,sigma);
		SET_VECTOR_ELT(phi,2,prop);
		
		REAL(Linc)[K-1] = REAL(VECTOR_ELT(outEM,2))[0];

		SET_VECTOR_ELT(paramv,0,VECTOR_ELT(outEM,0));
		SET_VECTOR_ELT(paramv,1,rupt);
		SET_VECTOR_ELT(paramv,2,tau);
		SET_VECTOR_ELT(param,K-1,paramv);

		UNPROTECT(3);
	    }

	    UNPROTECT(4); /* Unprotects paramv, phi, rupt and tau */
	    oukonenez++;
	}

	UNPROTECT(2); /* Unprotects G and out */
    }

    SET_VECTOR_ELT(resu,0,Linc);
    SET_VECTOR_ELT(resu,1,param);
    
    UNPROTECT(9); /* the rest */

    return(resu);
}
