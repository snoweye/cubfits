#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

/* log posterior vector for rocnsef model.
 * R_lp_vec: double[R_nrow, R_ncol] proportional probability observing codons
 * R_Codon_id: integer[sum(R_naa)] observed codon id
 * R_naa: interer[total_sequences] total codon counts for each sequence.
 */
SEXP lp_c_raw(SEXP R_lp_vec, SEXP R_nrow, SEXP R_ncol, SEXP R_Codon_id,
		SEXP R_naa){
	SEXP R_lp_c_raw;
	double *C_lp_c_raw, *C_lp_vec;
	int *C_ncol, *C_Codon_id, *C_naa;
	int i, j, total_sequences, counter;

	C_lp_vec = REAL(R_lp_vec);
	C_ncol = INTEGER(R_ncol);
	C_Codon_id = INTEGER(R_Codon_id);
	C_naa = INTEGER(R_naa);		// Sum of all R_naa should be R_nrow.

	total_sequences = LENGTH(R_naa);
	PROTECT(R_lp_c_raw = allocVector(REALSXP, total_sequences));

	// Set initials to 0.
	C_lp_c_raw = REAL(R_lp_c_raw);
	for(i = 0; i < total_sequences; i++){
		*C_lp_c_raw = 0;
		C_lp_c_raw++;
	}

	for(j = 0; j < *C_ncol; j++){
		// Rewind to the head of pointers.
		C_lp_c_raw = REAL(R_lp_c_raw);
		C_Codon_id = INTEGER(R_Codon_id);
		C_naa = INTEGER(R_naa);

		for(i = 0; i < total_sequences; i++){
			counter = *C_naa;
			// Skip sequence which has 0 codon for given AA.
			while(counter > 0){
				// Only add value for the matched codon.
				if(*C_Codon_id == j){
					*C_lp_c_raw = *C_lp_c_raw + *C_lp_vec;
				}
				C_lp_vec++;
				C_Codon_id++;
				counter--;
			}
			C_lp_c_raw++;
			C_naa++;
		}
	}

	UNPROTECT(1);
	return(R_lp_c_raw);
} /* End of lp_c_raw(). */

