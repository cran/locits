#include <R.h>

/* Parameters that control the storage for pre-computed values */


#define	MAXELL	20
#define	MAXJ	20
#define MAXK	20
#define MAXD	1800
#define DOFFSET	900 /* SHould be half MAXD */

/* A big store
#define	MAXELL	25
#define	MAXJ	25
#define MAXK	25
#define MAXD	30000
#define DOFFSET	15000 
*/
/*
#define	MAXELL	0
#define	MAXJ	0
#define MAXK	0
#define MAXD	0
#define DOFFSET	0	
*/

static double ThmStore[MAXELL][MAXJ][MAXK][MAXD];
static char ValExists[MAXELL][MAXJ][MAXK][MAXD];
double nfound, nstored, noutside;

#define TRUE	1
#define FALSE	0


#define ACCESS(array, nrows, col, row)	*(array+(col)*(nrows)+row)
#define max(a,b)        ((a) > (b) ? (a) : (b))
#define min(a,b)        ((a) > (b) ? (b) : (a))
/* Note: j is usual index (1,2,3...) and tau is usual one too */
#define ACCESSPSI(PsiJ, linPsiJ, lvPsiJ, j, tau)  *(PsiJ + *(linPsiJ+j-1) + tau  + *(lvPsiJ+j-1))


void initThmStore(int *error)
{
register int ell, j, k, d;

for(ell=0; ell < MAXELL; ++ell)
    for(j=0; j < MAXJ; ++j)
        for(k=0; k < MAXK; ++k)
            for(d=0; d < MAXD; ++d)	{
    		ValExists[ell][j][k][d] = FALSE;		    
    		ThmStore[ell][j][k][d] = 0.0;		    
	    }

nfound = nstored = noutside = 0.0;

*error = 0;
}

void StoreStatistics(double *lfound, double *lstored, double *loutside)
{

	*lfound = nfound;
	*lstored = nstored;
	*loutside = noutside;
}


void CPkPlj(double *Pk, int *lPk,
		double *psil, int *lpsil,
		double *psij, int *lpsij,
		int *d, 
		double *ans, int *verbose, int *error)
{
int mintau;
int maxtau;
int tau;
int k, lok, upk;
double Psiktau;
double Psiljdptau;
double localans=0.0;

/*
maxtau = (*lPk-1)/2;
*/
maxtau = *lPk;
mintau = -maxtau;

if (*verbose >= 2)	{
	Rprintf("lPk: %d\n", *lPk);
	Rprintf("mintau, maxtau: %d, %d\n", mintau, maxtau);
	}

for(tau=mintau; tau <= maxtau; ++tau)	{

	Psiktau = *(Pk + tau - mintau);

	lok = max(0, tau-*d);
	upk = min(*lpsil-1, *lpsij-1-*d+tau);

	Psiljdptau = (double)0.0;

	if (lok <= upk)	{
		for(k=lok; k <= upk; ++k)	{
			Psiljdptau += *(psil+k) * (*(psij+k+(*d)-tau));
		}
	}

	if (*verbose >= 3)
		Rprintf("tau: %d Psiktau=%lf Psiljdptau=%lf\n", tau, Psiktau, Psiljdptau);

	localans += Psiktau*Psiljdptau;
}

*ans = localans;
}

/* CcovIxscale: Formula in Theorem 1 of NuTestSty paper
 *
 * ell, j are the parameters in the theorem. Numbered from 1 (finest scale),
 * 		2, 3, 4, ... upwards for coarser scales. They can't be
 * 		bigger than *J
 *
 * m, n	the locations
 * II:	the raw wavelet periodogram at m (used when n=m)
 * S:  the spectral estimate at (m+n)/2
 * J:	the total number of scales (2^J = length of series)
 * PsiJ: the autocorrelation wavelet, in one big vector, finest scale first
 * lPsiJ: length of PsiJ vector
 * linPsiJ array that contains entry points into PsiJ. First entry in
 * 		linPsiJ corresponds to finest scale (actually, 0),
 * 		Second entry in linPsiJ corresponds to second scale, etc.
 * lvPsiJ: length of each scale of autocorrelation coefficients
 * psil, lpsil:	discrete wavelet at scale l, and its length
 * psij, lpsij:	discrete wavelet at scale j, and its length
 * verbose: if 1, some messages, if 2, loads of messages
 * ans:	the answer to return
 * error: the error code (which will be non zero if error occurs
 */

void CcovIxscale(int *ell, int *j, int *m, int *n,
		double *II, double *S, int *J,
		double *PsiJ, int *lPsiJ, int *linPsiJ, int *lvPsiJ,
		double *psil, int *lpsil,
		double *psij, int *lpsij,
		int *verbose,
		double *ans, int *error)
{
int k;
int d;
double bigsum;
double rh;
void CPkPlj();


*error =0;

/*
Rprintf("CcovIxscale: ell=%d, j=%d, m=%d, n=%d\n", *ell, *j, *m, *n);
*/


if (*ell > *j)	{
	*error = 1;
	return;
	}

if (*j > *J)	{
	*error = 2;
	return;
	}

if (*ell < 0)	{
	*error = 3;
	return;
	}

if (*j < 0)	{
	*error = 4;
	return;
	}


/*
 * Note: bug corrected here. Previous code did not also have the
 * *ell == *j check. The following line is only true when the
 * scales and locations are the same
 */

if (*m==*n && *ell == *j)	
	*ans = 2.0*(*(II+*ell-1))*(*(II+*j-1));

else	{

	rh = bigsum = 0.0;

	d = *n - *m;

	for(k=0; k < *J; ++k)	{


		/*
		Rprintf("*ell-1: %d MAXELL: %d\n", *ell-1, MAXELL);
		Rprintf("*j-1: %d MAXJ: %d\n", *ell-1, MAXJ);
		Rprintf("k: %d MAXK: %d\n", k, MAXK);
		Rprintf("d: %d d+DOFFSET %d; MAXD: %d\n", d, d+DOFFSET, MAXD);
		*/

		/*
		 * Storage here
		 */

		if (((*ell-1) < MAXELL) && ((*j-1) < MAXJ) && (k < MAXK) &&
			(d+DOFFSET >= 0) && (d+DOFFSET < MAXD))	{
			if (ValExists[*ell-1][*j-1][k][d+DOFFSET])	{
			    rh = ThmStore[*ell-1][*j-1][k][d+DOFFSET];
			    /*
			    Rprintf("Found: %d %d %d %d\n", *ell, *j, k, d);
			    */
			    ++nfound;

			    }

			else	{
				CPkPlj(PsiJ+*(linPsiJ+k), lvPsiJ+k,
					psil, lpsil, psij, lpsij,
					&d, &rh, verbose, error);
				ThmStore[*ell-1][*j-1][k][d+DOFFSET] = rh;
				ValExists[*ell-1][*j-1][k][d+DOFFSET] = TRUE;
				/*
			        Rprintf("Store: %d %d %d %d\n", *ell, *j, k, d);
				*/
				++nstored;
				}
			}
		else	/* Storage does not exist for this range of parms */
			{
			/*
			Rprintf("Out of range: calculating raw\n");
			*/
			CPkPlj(PsiJ+*(linPsiJ+k), lvPsiJ+k,
				psil, lpsil, psij, lpsij,
				&d, &rh, verbose, error);
			++noutside;
			}


		if (*error != 0)
			return;

		if (*verbose>=1)
			Rprintf("k: %d; S[k]: %lf; rh: %lf\n", k, *(S+k), rh);

		bigsum += *(S+k) * rh;
		}

	*ans = 2.0*bigsum*bigsum;
	}

}

/* CstarIcov
 *
 * Computes the covariance of the running mean smooth of the wavelet
 * periodogram. It calculates the covariance between \tilde{I}_{ell, nz}
 * and \tilde{I}_{j, nz} where ell and j are scales and nz is the location.
 * The quantity \tilde{I} is the smoothed periodogram (smoothed using a
 * simple running mean of window width of s (formula (13) in the tech
 * report version of the paper).
 *
 * Parameters are
 *
 * ell: one scale
 * j: other scale
 * nz: location
 * s: smoothing bandwidth
 * TT: length of series
 * IIvec: the entire smoothed wavelet periodogram (see call to CstarIcov
 * 	for example of how this can be passed to this routine from R in the
 * 	Rvarlacf function). Essentially, the whole smoothed wavelet periodgram
 * 	matrix is passed as a vector to this routine.
 * Svec: as for IIvec, except for the evolutionary wavelet estimator.
 * J: number of scales associated with IIvec and Svec
 *  PsiJ: the autocorrelation wavelet, in one big vector, finest scale first
 *  lPsiJ: length of PsiJ vector
 *  linPsiJ: array that contains entry points into PsiJ. First entry in
 *	linPsiJ corresponds to finest scale (actually, 0),
 *      Second entry in linPsiJ corresponds to second scale, etc.
 *  lvPsiJ: length of each scale of autocorrelation coefficients
 *  psil, lpsil: discrete wavelet at scale l, and its length
 *  psij, lpsij: discrete wavelet at scale j, and its length
 *  verbose: if 1, some messages, if 2, loads of messages
 *  ans: the answer to return
 *  error: the error code (which will be non zero if error occurs)
 */


void CstarIcov(int *ell, int *j, int *nz, int *s, int *TT,
		double *IIvec, double *Svec, int *J,
		double *PsiJ, int *lPsiJ, int *linPsiJ, int *lvPsiJ,
		double *psil, int *lpsil,
		double *psij, int *lpsij,
		int *verbose,
		double *ans, int *error)
{
int minlim;
int maxlim;
int u, v;
int avpos;
double LocalAns;
double denom;
void CcovIxscale();	/* Calculates the covariance in Thm 1 of Nason2013*/

*ans = 0.0;
*error = 0;


minlim = max(*nz - *s, 1);
maxlim = min(*nz + *s, *TT);

for(u=minlim; u <= maxlim; ++u)	{

	for(v=minlim; v <= maxlim; ++v)	{

		LocalAns = 0.0;

		avpos = (int)(((double)u+(double)v)/2.0);

		CcovIxscale(ell, j, &u, &v,
			IIvec + *J*(avpos-1),
			Svec + *J*(avpos-1),
			J,
			PsiJ, lPsiJ, linPsiJ, lvPsiJ,
			psil, lpsil,
			psij, lpsij,
			verbose,
			&LocalAns, error);

		if (*error != 0)
			return;

		*ans += LocalAns;
		}
	}

/*
Rprintf("nfound %lf, ncomputed: %lf; computed percent %lf\n", nfound, ncomputed, 100*(double)ncomputed/(double)(nfound+ncomputed));
*/

/*
 * Denominator should really be the number that we calculated
 */

denom = 2.0*(double)*s + 1.0;
denom = denom*denom;

*ans /= denom;
}
