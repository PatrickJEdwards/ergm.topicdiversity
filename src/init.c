#include <stdlib.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
  {NULL, NULL, 0}
};

void R_init_ergm_topicdiversity(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
