#ifndef __EVALEXPR_H
#define __EVALEXPR_H

/*
 * Evaluate a SUGAR expression.
 *
 * Functions:
 *  eval_name -- create a string corresponding to a fully indexed name
 *               (but don't include scoping info)
 *  eval_expr -- evaluate an expression
 *
 * The only reason these functions take a memory pool as an argument is
 * so they have a place to create new string objects if need be.
 */

#include "parse.h"

/* Store a variable value */
typedef struct {
    char type;    /* 'd' == double, 's' == string, or 'u' == undef */
    union {
        double d;
	char* s;
    } val;
} val_t;

void eval_name(mempool_t pool, char* buf, code_t* node);
void eval_expr(mempool_t pool, val_t* val_out, code_t* node);

#endif /* __EVALEXPR_H */
