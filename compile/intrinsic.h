#ifndef __INTRINSIC_H
#define __INTRINSIC_H

/*
 * Evaluate a SUGAR intrinsic.  Used only by evalexpr
 */

#include "evalexpr.h"

int eval_intrinsic_call(code_t* node, char* name, val_t* val_out, 
		        val_t* argval, int nargs);

#endif /* __INTRINSIC_H */
