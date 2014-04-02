// Copied from "phyclust/src/phyclust_em_step.c".

#include <math.h>

/* Gamma[n][k]
   = pi_k * L_k(X_n) / sum_i(pi_i * L_i(X_n))
   = 1 / sum_i(pi_i * L_i(X_n) / (pi_k * L_k))
   = 1 / sum_i(exp(log(pi_i) + log(L_i(X_n)) - log(pi_k) - log(L_k(X_n))))
   This function return stable exponential values with a scale exponeitial value and flag.
   If flag = 1, scale_exp will be used to adjuste results eveywhere. Otherwise scale_exp = 0.
   *total_sum is for logL_observed.
*/
void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum, double *scale_exp, int *flag_out_range){
	int k;
	double tmp_exp, max_exp;

	*total_sum = 0.0;
	*scale_exp = 0.0;
	*flag_out_range = 0;
	max_exp = a_Z_normalized[0];
	for(k = 1; k < *K; k++){
		if(a_Z_normalized[k] > max_exp){
			max_exp = a_Z_normalized[k];
		}
	}

	/* tmp_exp = HUGE_VAL for overflow and 0 for underflow.
	 *   e.g. max_exp is large when parameters are near the boundary and is tiny when too many products.
	 *        errno = ERANGE, only when tmp_exp is too huge (HUGE_VAL).
	 *        errno = 0, when tmp_exp is too small. "Avoid to use ERANGE to check overflow or underflow!"
	 * Scale max_exp by 2 such that close to +0 or -0 until errno is not ERANGE or tmp_exp is not HUGE_VAL. */
/* BUG!
	errno = 0;
	tmp_exp = exp(max_exp);
	if(tmp_exp == HUGE_VAL){
		*flag_out_range = 1;
		*scale_exp = (tmp_exp == HUGE_VAL) ? max_exp : -max_exp;
		do{
			errno = 0;
			*scale_exp *= 0.5;
			tmp_exp = exp(*scale_exp);
		} while(errno == ERANGE);
		*scale_exp = max_exp - *scale_exp;
	}
*/
	tmp_exp = exp(max_exp);
	if(tmp_exp == HUGE_VAL || tmp_exp == 0.0){
		*flag_out_range = 1;
		*scale_exp = (tmp_exp == HUGE_VAL) ? max_exp : -max_exp;
		do{
			*scale_exp *= 0.5;
			tmp_exp = exp(*scale_exp);
		} while(tmp_exp == HUGE_VAL);
		*scale_exp = max_exp - *scale_exp;
		/* The *scale_exp is always greater than 0.
		 * If max_exp > 0 and too large, then shift all to left and computable.
		 *   c = max_exp - *scale_exp > 0, and all should minus c.
		 * If max_exp < 0 and too small, then shift all to right and positive computable.
		 *   c = max_exp - *scale_exp < max_exp < 0, and all should minus c. */
	}
	if(*flag_out_range){
		for(k = 0; k < *K; k++){
			a_Z_normalized[k] -= *scale_exp;
		}
	}

	*total_sum = 0.0;
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = exp(a_Z_normalized[k]);
		*total_sum += a_Z_normalized[k];
	}
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = a_Z_normalized[k] / *total_sum;
	}
} /* End of e_step_with_stable_exp(); */
