#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum,
		double *scale_exp, int *flag_out_range);

/* Inverse of multivariate logit.
 * R_z: double(nrow, ncol), exponent terms of -(eta phi t).
 */
SEXP invmlogit(SEXP R_z, SEXP R_nrow, SEXP R_ncol){
	// Inverse multinomial logit function, vectorized calculates
	// probabilities from logodds, using last value as reference level.
	SEXP R_p;
	double *C_z, *C_p, *C_tmp;
        double tmp_total_sum, tmp_scale_exp;
	int tmp_flag_out_range;
	int nrow, ncol, n_category, i, j;

	C_z = REAL(R_z);
	nrow = INTEGER(R_nrow)[0];
	ncol = INTEGER(R_ncol)[0];
	n_category = ncol + 1;

	// Allocate storage for returning.
	PROTECT(R_p = allocMatrix(REALSXP, nrow, n_category));
	C_p = REAL(R_p);

	// For adjusting overflow.
	C_tmp = (double *) malloc(n_category * sizeof(double));
	if(C_tmp == NULL){
		error("Memory allocation fails!\n");
	}

	for(i = 0; i < nrow; i++){
		// Make a small copy from C_z.
		for(j = 0; j < ncol; j++){
			*(C_tmp + j) = *(C_z + j * nrow);
		}
		// R_p has one extra column than R_z.
		*(C_tmp + ncol) = 0.0;

		// Call stable_exp to updae C_tmp.
		e_step_with_stable_exp(&n_category, C_tmp, &tmp_total_sum,
			&tmp_scale_exp, &tmp_flag_out_range);

		// Copy results back to C_p;
		for(j = 0; j < n_category; j++){
			*(C_p + j * nrow) = *(C_tmp + j);
		}

		// Move pointers.
		C_p++;
		C_z++;
	}

	free(C_tmp);
	UNPROTECT(1);
	return(R_p);
} /* End of invmlogit(). */
