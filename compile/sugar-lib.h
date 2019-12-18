#ifndef __SUGAR_LIB_H
#define __SUGAR_LIB_H

#ifdef __SUGAR_MEX
#include "mex.h"
#define printf mexPrintf
#define malloc mxMalloc
#define free mxFree
#endif

#endif /* __SUGAR_LIB_H */
