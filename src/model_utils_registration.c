
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "model_fitters.h"

static R_NativePrimitiveArgType fit_model_t[15]={INTSXP, REALSXP, REALSXP, INTSXP,INTSXP,INTSXP, 
						 REALSXP,REALSXP,REALSXP,REALSXP, INTSXP,REALSXP, 
						 INTSXP,INTSXP,REALSXP};
static R_NativePrimitiveArgType fit_model_with_allocation_t[13]={INTSXP, REALSXP, REALSXP, INTSXP,
								 INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,
								 REALSXP, INTSXP,REALSXP,INTSXP};
R_CMethodDef cMethods[] = {
  {"fit_model", (DL_FUNC) &fit_model, 15, fit_model_t},
  {"fit_model_with_allocation", (DL_FUNC) &fit_model_with_allocation, 13, fit_model_with_allocation_t},
  {NULL, NULL, 0}
};
          
void R_init_modelUtils(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);

  R_RegisterCCallable("modelUtils", "model_allocate_int", (DL_FUNC) &model_allocate_int);
  R_RegisterCCallable("modelUtils", "model_allocate_double", (DL_FUNC) &model_allocate_double);
  R_RegisterCCallable("modelUtils", "fit_model", (DL_FUNC) &fit_model);
  R_RegisterCCallable("modelUtils", "fit_model_with_allocation", (DL_FUNC) &fit_model_with_allocation);

  /** set the default value of the debug level. **/
  R_MODEL_DEBUG_LEVEL = 10;
}
          
void
R_unload_modelUtils(DllInfo *info)
{
  /* Release resources. */
}
