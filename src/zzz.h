#include <R.h>
#include <Rinternals.h>

SEXP invmlogit(SEXP R_z, SEXP R_nrow, SEXP R_ncol);
SEXP lp_c_raw(SEXP R_lp_vec, SEXP R_nrow, SEXP R_ncol, SEXP R_Codon_id,
		SEXP R_naa);
