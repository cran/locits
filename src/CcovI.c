#include <R.h>

#define ACCESS(array, nrows, col, row)	*(array+(col)*(nrows)+row)
#define max(a,b)        ((a) > (b) ? (a) : (b))
#define min(a,b)        ((a) > (b) ? (b) : (a))


void Cvarip2(int *i, int *ll,
		double *S, int *lS,
		double *Pmat, int *ncP, int *nrP,
		int *lP,
		double *ans)
{
int r;
int rul;
int TWOimONE;
int m,n;
double covAA,covAB,covBB;
double tmpans;

void CcovI();

covAA = covAB = covBB = 0.0;

/*
 * Compute first r upper limit
 */

TWOimONE = (1 << (*i-1));

rul = TWOimONE-1;

m = 1;

for(r=0; r<= rul; ++r)	{

	n = r+1;

	CcovI(S, lS, 
		&m, &n, ll,
		Pmat, ncP, nrP,
		lP,
		&tmpans);

	if (r == 0)
		covAA += (TWOimONE-r)*tmpans;
	else
		covAA += (TWOimONE-r)*tmpans*2;
	}

for(r=-rul; r <= rul; ++r)	{

	n = 1 + r + TWOimONE;

	CcovI(S, lS, 
		&m, &n, ll,
		Pmat, ncP, nrP,
		lP,
		&tmpans);

	covAB += (TWOimONE - abs(r))*tmpans;
	}

covBB = covAA;

/*
Rprintf("covAA=%lf, covAV=%lf\n", covAA, covAB);
*/

*ans = covAA - 2.0*covAB + covBB;

*ans = pow(2.0, (double)(-*i))*(*ans);
}


void CcovI(double *S, int *lS, 
		int *m, int *n, int *ll,
		double *Pmat, int *ncP, int *nrP,
		int *lP,
		double *bigans)
{
int i,k;
int Nll, Nk;
int mintau, maxtau;
double ans;
double tmpa, tmpb;
int iix;

/* Pmat is the C array version of Pmat in R. So, the C array indices
 * are always one less than those in R. */


iix = *ncP/2;

*bigans = 0.0;

Nll = (*(lP+*ll-1) -1)/2;

for(k=0; k < *lS; ++k)	{

	Nk = (*(lP+k) - 1)/2;

	mintau = max(-Nll+*n-*m, -Nk);
	maxtau = min(Nll+*n-*m, Nk);

	ans = 0.0;

	/*
	Rprintf("mintau:maxtau %d %d\n", mintau, maxtau);
	*/

	if (mintau <= maxtau)
		for(i=mintau; i<=maxtau; ++i)	{
			tmpa = ACCESS(Pmat, *nrP, i+iix, k);
			tmpb = ACCESS(Pmat, *nrP, *m-*n+i+iix, *ll-1);

			/*
			Rprintf("\t P[%d,%d]=%lf, P[%d, %d]=%lf\n", 
					i+iix, k, tmpa,
					*m-*n+i+iix, +*ll-1, tmpb);	
			*/

			ans += tmpa*tmpb;
			}

	/*
	Rprintf("%d %lf\n", k, ans);
	*/

	*bigans += *(S+k) * ans;

}

/*
for(i=0; i<*nrP; ++i)
Rprintf("Pmat[4096, %d] = %lf\n", i, ACCESS(Pmat, *nrP, 4096, i));
*/

*bigans = 2*(*bigans)*(*bigans);
}
