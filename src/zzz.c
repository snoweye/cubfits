#include <R.h>
#include <R_ext/Rdynload.h>

#include "zzz.h"

static const R_CallMethodDef callMethods[] = {
	{"invmlogit", (DL_FUNC) &invmlogit, 3},
	{"lp_c_raw", (DL_FUNC) &lp_c_raw, 5},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_cubfits(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_cubfits(). */
