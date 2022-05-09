#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */

extern void CstarIcov(int *ell, int *j, int *nz, int *s, int *TT,
                double *IIvec, double *Svec, int *J,
                double *PsiJ, int *lPsiJ, int *linPsiJ, int *lvPsiJ,
                double *psil, int *lpsil,
                double *psij, int *lpsij,
                int *truedenom,
                int *verbose,
                double *ans, int *error);

extern void Cvarip2(int *i, int *ll,
                double *S, int *lS,
                double *Pmat, int *ncP, int *nrP,
                int *lP,
                double *ans);


extern void initThmStore(int *error);

extern void StoreStatistics(double *lfound, double *lstored, double *loutside);

static const R_CMethodDef CEntries[] = {
    {"CstarIcov",       (DL_FUNC) &CstarIcov,       20},
    {"Cvarip2",         (DL_FUNC) &Cvarip2,          9},
    {"initThmStore",    (DL_FUNC) &initThmStore,     1},
    {"StoreStatistics", (DL_FUNC) &StoreStatistics,  3},
    {NULL, NULL, 0}
};

void R_init_locits(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
