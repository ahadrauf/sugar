#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "codestack.h"
#include "codegen.h"
#include "evalexpr.h"
#include "intrinsic.h"

/* This is really obnoxious, but MSC apparently thinks it doesn't need to
 * implement the standard C headers.  In particular, it doesn't define M_PI.
 * Meh.
 */
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

static void sugarf_negate(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = -argval[0].val.d;
}

static void sugarf_add(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = argval[0].val.d + argval[1].val.d;
}

static void sugarf_sub(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = argval[0].val.d - argval[1].val.d;
}

static void sugarf_mul(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = argval[0].val.d * argval[1].val.d;
}

static void sugarf_div(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = argval[0].val.d / argval[1].val.d;
}

static void sugarf_pow(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = pow(argval[0].val.d, argval[1].val.d);
}

static void sugarf_not(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d == 0);
}

static void sugarf_eqd(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d == argval[1].val.d);
}

static void sugarf_eqs(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (strcmp(argval[0].val.s, argval[1].val.s) == 0);
}

static void sugarf_neqd(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d != argval[1].val.d);
}

static void sugarf_neqs(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (strcmp(argval[0].val.s, argval[1].val.s) != 0);
}

static void sugarf_lt(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d < argval[1].val.d);
}

static void sugarf_gt(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d > argval[1].val.d);
}

static void sugarf_leq(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d <= argval[1].val.d);
}

static void sugarf_geq(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].val.d >= argval[1].val.d);
}

static void sugarf_cos(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = cos(argval[0].val.d);
}

static void sugarf_sin(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = sin(argval[0].val.d);
}

static void sugarf_tan(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = tan(argval[0].val.d);
}

static void sugarf_sqrt(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = sqrt(argval[0].val.d);
}

static void sugarf_deg(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = M_PI * (argval[0].val.d / 180);
}

static void sugarf_mod(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = fmod(argval[0].val.d, argval[1].val.d);
}

static void sugarf_assert(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = 0;
    if (argval[0].val.d == 0) {
        add_error_trace_atcode(node, argval[1].val.s);
    }
}

static void sugarf_isdef(val_t* val_out, val_t* argval, code_t* node)
{
    val_out->type = 'd';
    val_out->val.d = (argval[0].type != 'u');
}

struct sugar_function {
    char* name;
    char* signature;
    void (*eval)(val_t* val_out, val_t* argval, code_t* node);
};

static struct sugar_function func_list[] = {
    {"unary-", "d", sugarf_negate},
    {"+", "dd", sugarf_add},
    {"-", "dd", sugarf_sub},
    {"*", "dd", sugarf_mul},
    {"/", "dd", sugarf_div},
    {"^", "dd", sugarf_pow},
    {"!", "d",  sugarf_not},
    {"==", "dd", sugarf_eqd},
    {"==", "ss", sugarf_eqs},
    {"!=", "dd", sugarf_neqd},
    {"!=", "ss", sugarf_neqs},
    {"<",  "dd", sugarf_lt},
    {">",  "dd", sugarf_gt},
    {"<=", "dd", sugarf_leq},
    {">=", "dd", sugarf_geq},
    {"cos", "d", sugarf_cos},
    {"sin", "d", sugarf_sin},
    {"tan", "d", sugarf_tan},
    {"sqrt", "d", sugarf_sqrt},
    {"deg", "d", sugarf_deg},
    {"mod", "dd", sugarf_mod},
    {"assert", "ds", sugarf_assert},
    {"isdef", "?", sugarf_isdef},
    {NULL, NULL, NULL}
};

/* Check if a signature matches, doing any casts to make it match
 * if possible
 */
static int check_signature(val_t* argval, int nargs, char* sign)
{
    int i;

    /* Check to make sure args match */
    for (i = 0; i < nargs && sign[i] != '\0'; ++i) {
        if (argval[i].type != sign[i] && sign[i] != '?') {
            return 0;
        }
    }

    /* Make sure number of args matches */
    if (i < nargs || sign[i] != '\0')
        return 0;

    return 1;
}

/* Evaluate a call to some intrinsic
 */
int eval_intrinsic_call(code_t* node, char* name, val_t* val_out, 
		        val_t* argval, int nargs)
{
    struct sugar_function* func = func_list;
    int matched_name = 0;
    
    while (func->name != NULL) {
        if (strcmp(func->name, name) == 0) {
	    matched_name = 1;
            if (check_signature(argval, nargs, func->signature)) {
                (func->eval)(val_out, argval, node);
		return 1;
            } 
        }
	++func;
    }

    if (matched_name) {
	char buf[128];
        val_out->type = 'e';
        sprintf(buf, "Invalid operands in %s", name);
        add_error_trace_atcode(node, buf);
	return 1;
    }

    return 0;
}


