// Copied from "phyclust/src/phyclust_em_step.c".

#include <math.h>
#include <float.h>

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


/*	Subtract the maximum from everything, guaranteeing values in (-inf,0]
		This may cause some values to go so low that when they are
		exponentiated, they go to 0. This isn't a huge problem,
		if the distance between the two phi values is 100 orders
		of magnitude, the precision lost is irrelevant.			*/

	*total_sum = 0.0;
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = exp(a_Z_normalized[k] - max_exp);
		*total_sum += a_Z_normalized[k];
	}
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = a_Z_normalized[k] / *total_sum;
	}
} /* End of e_step_with_stable_exp(); */
