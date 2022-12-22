#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(divset_f)(int *alg, int *iv, int liv, int lv, double *v);

//' divset_c
//'
//' c interface to divset fortran function
//'
//' @param alg integer indicating type of algorithm (1 = optim, 2=est)
//' @return returns an integer array and a real array
extern SEXP divset_c(SEXP alg, SEXP iv, SEXP v){
  const int liv = LENGTH(iv);
  const int lv = LENGTH(v);

  /* SEXP ret_iv = PROTECT(allocVector(INTSXP, liv));
  SEXP ret_v  = PROTECT(allocVector(REALSXP, lv));*/
  SEXP out = PROTECT(allocVector(VECSXP,2));
  F77_CALL(divset_f)(INTEGER(alg), INTEGER(iv), liv, lv, REAL(v));
  /* ret_iv = iv;
  ret_v = v;  */
  SET_VECTOR_ELT(out,0,iv);
  SET_VECTOR_ELT(out,1,v);
  UNPROTECT(1);
  return(out);
}

void F77_NAME(drglg) (double *d, double *dr, int *iv, int *liv, \
              int *lv, int *n, int *nd, int *nn, int *p, int *ps, \
              int lrhoi, int lrhor, double *r, double *rd, \
              double *v, double *x, int *rhoi, double *rhor, \
                      int *i_itsum);

extern SEXP drglg_c (SEXP d, SEXP dr, SEXP iv, SEXP liv, SEXP lv, \
                     SEXP n, SEXP nd, \
                     SEXP nn, SEXP p, SEXP ps, SEXP r, SEXP rd, \
                     SEXP v, SEXP x, SEXP rhoi, SEXP rhor, \
                     SEXP i_itsum) {
  const int lrhoi = LENGTH(rhoi);
  const int lrhor  = LENGTH(rhor);
  /*  The big question is:  How many things need to be returned?*/
  /* SEXP ret_iv = PROTECT(allocVector(INTSXP, liv));
  SEXP ret_v  = PROTECT(allocVector(REALSXP, lv));*/
  SEXP out = PROTECT(allocVector(VECSXP,10));
  F77_CALL(drglg)(REAL(d), REAL(dr), INTEGER(iv), INTEGER(liv), \
           INTEGER(lv), \
           INTEGER (n), INTEGER(nd), INTEGER(nn), INTEGER(p), \
           INTEGER(ps), lrhoi, lrhor, REAL(r), REAL(rd), \
           REAL(v), REAL(x), INTEGER(rhoi), REAL(rhor), \
           INTEGER(i_itsum));
  /* ret_iv = iv;
  ret_v = v;  */
  SET_VECTOR_ELT(out,0,d);
  SET_VECTOR_ELT(out,1,dr);
  SET_VECTOR_ELT(out,2,iv);
  SET_VECTOR_ELT(out,3,r);
  SET_VECTOR_ELT(out,4,rd);
  SET_VECTOR_ELT(out,5,v);
  SET_VECTOR_ELT(out,6,x);
  SET_VECTOR_ELT(out,7,rhoi);
  SET_VECTOR_ELT(out,8,rhor);
  SET_VECTOR_ELT(out,9,i_itsum);
  UNPROTECT(1);
  return(out);
}

static const R_CallMethodDef CallEntries[] = {
  {"divset_c",   (DL_FUNC) &divset_c,  3},
  {"drglg_c",   (DL_FUNC) &drglg_c,  17},
  {NULL,         NULL,               0}
};

void R_init_bgw (DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
