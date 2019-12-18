#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "codestack.h"
#include "codegen.h"
#include "evalexpr.h"
#include "intrinsic.h"

static void eval_name_idx(mempool_t pool, char* buf, code_t* node);
static void eval_call(mempool_t pool, val_t* val_out, code_t* node);
static int eval_arguments(mempool_t pool, code_t* args, val_t* argval, int n);

/* Evaluate a name (without scoping info)
 */
void eval_name(mempool_t pool, char* buf, code_t* node)
{
    code_t* idx;

    strcpy(buf, node->v.name.name);

    if (node->v.name.indices) {
        strcat(buf, "(");
        idx = node->v.name.indices;
        while (idx->next != NULL) {
          eval_name_idx(pool, buf, idx);
          strcat(buf, ",");
          idx = idx->next;
        }
        eval_name_idx(pool, buf, idx);
        strcat(buf, ")");
    }
}

/* Cat an index onto a name-in-progress.  We don't like non-integer indices
 */
static void eval_name_idx(mempool_t pool, char* buf, code_t* node)
{
    val_t idx_val;
    char numbuf[128];

    eval_expr(pool, &idx_val, node);
    if (idx_val.type != 'd' || idx_val.val.d != (int) idx_val.val.d)
        add_error_trace_atcode(node, "Non-integer index!");
    else
        sprintf(numbuf, "%d", (int) idx_val.val.d);
    strcat(buf, numbuf);
}

/* Evaluate an expression
 */
void eval_expr(mempool_t pool, val_t* val_out, code_t* node)
{
    switch (node->tag) {
    case DOUBLE_NODE:
        val_out->type = 'd';
        val_out->val.d = node->v.dval;
        break;
    case STRING_NODE:
	val_out->type = 's';
	val_out->val.s = node->v.sval;
	break;
    case UNDEF_NODE:
	val_out->type = 'u';
	break;
    case CALL_NODE:
        eval_call(pool, val_out, node);
        break;
    case VAR_NODE:
        *val_out = *frame_var_entry(node->v.var->v.def.id, 
			            node->v.var->v.def.scope);
        break;
    default:
        assert(0);
    }

}

/* Evaluate an argument list
 */
static int eval_arguments(mempool_t pool, code_t* args, val_t* argval, int n)
{
    int nargs;

    /* Evaluate the arguments, and exit if any of them are erroneous */
    code_t* arg = args;
    for (nargs = 0; arg != NULL && nargs < n; ++nargs) {
        eval_expr(pool, &(argval[nargs]), arg);
        if (argval[nargs].type == 'e')
            return -1;
        arg = arg->next;
    }

    if (arg != NULL) {
        add_error_trace_atcode(arg, "Too many arguments");
        return -1;
    }

    return nargs;
}

/* Evaluate a function call.  Currently I only have sin and cos implemented.
 * There should probably be a "punt" case which calls Matlab with appropriate
 * arguments if we don't recognize the function name, and only returns
 * function not found if Matlab can't evaluate it, either.
 */
static void eval_call(mempool_t pool, val_t* val_out, code_t* node)
{
    if (strcmp(node->v.call.name, "&") == 0) {
        code_t* args = node->v.call.args;

	val_out->type = 'd';
	val_out->val.d = 1;

	while (args != NULL) {
	    val_t argval;
            eval_expr(pool, &argval, args);
            if (argval.type == 'e') {
	        val_out->type = 'e';
		return;
	    } else if (argval.type != 'd') {
                val_out->type = 'e';
                add_error_trace_atcode(node, "Invalid operand in &");
		return;
	    } else if (argval.val.d == 0) {
                val_out->val.d = 0;
		break;
	    }
	    args = args->next;
	}
	
    } else if (strcmp(node->v.call.name, "|") == 0) {
        code_t* args = node->v.call.args;

	val_out->type = 'd';
	val_out->val.d = 0;

	while (args != NULL) {
	    val_t argval;
            eval_expr(pool, &argval, args);
            if (argval.type == 'e') {
	        val_out->type = 'e';
		return;
	    } else if (argval.type != 'd') {
                val_out->type = 'e';
                add_error_trace_atcode(node, "Invalid operand in |");
		return;
	    } else if (argval.val.d != 0) {
                val_out->val.d = 1;
		break;
	    }
	    args = args->next;
	}

    } else if (strcmp(node->v.call.name, "cond") == 0) {
        code_t* cond     = node->v.call.args;
	code_t* ifstmt   = (cond != NULL ? cond->next : NULL);
	code_t* elsestmt = (ifstmt != NULL ? ifstmt->next : NULL);

	val_t condval;

	if (elsestmt == NULL || elsestmt->next != NULL) {
	    val_out->type = 'e';
            add_error_trace_atcode(node, "Wrong number of operands in cond");
	    return;
	}

        eval_expr(pool, &condval, cond);
	if (condval.type == 'e') {
	    val_out->type = 'e';
	    return;
	} else if (condval.type != 'd') {
	    val_out->type = 'e';
            add_error_trace_atcode(node, "Wrong operand type in cond");
	    return;
	}

	if (condval.val.d != 0) 
	    eval_expr(pool, val_out, ifstmt);
	else
	    eval_expr(pool, val_out, elsestmt);

    } else if (strcmp(node->v.call.name, "print") == 0) {
        code_t* args = node->v.call.args;

	val_out->type = 'd';
	val_out->val.d = 0;

	while (args != NULL) {
	    val_t argval;
            eval_expr(pool, &argval, args);
            if (argval.type == 's') 
                printf("%s", argval.val.s);
            else if (argval.type == 'd') 
                printf("%g", argval.val.d);
            else if (argval.type == 'u') 
                printf("%s", "*undef*");
            else 
                printf("%s", "*error*");
	    args = args->next;
	}
	printf("\n");

    } else {
        val_t argval[15];  /* Somewhat arbitrary limitation... */
        int nargs;

        nargs = eval_arguments(pool, node->v.call.args, argval, 16);
        if (nargs < 0) {
            val_out->type = 'e';
            return;
        }

        /* Evaluate the actual functions */
        if (eval_intrinsic_call(node, node->v.call.name, 
				val_out, argval, nargs)) {
           /* Already done whatever we're going to do */ 

        } else {
#ifdef __SUGAR_MEX
            mxArray** prhs = (mxArray**) mxMalloc(nargs * sizeof(mxArray*));
            mxArray* result;
            int i;

            for (i = 0; i < nargs; ++i) {
                if (argval[i].type == 's') {
                    prhs[i] = mxCreateString(argval[i].val.s);
#ifdef MATLAB6
                } else if (argval[i].type == 'd') {
                    prhs[i] = mxCreateScalarDouble(argval[i].val.d);
                }
#else
                } else if (argval[i].type == 'd') {
                    prhs[i] = mxCreateDoubleMatrix(1,1,mxREAL);
                    *mxGetPr(prhs[i]) = argval[i].val.d;
                }
#endif
            }

            mexCallMATLAB(1, &result, nargs, prhs, node->v.call.name);

            for (i = 0; i < nargs; ++i) {
                mxDestroyArray(prhs[i]);
            }

            if (mxIsNumeric(result)) {
                val_out->type = 'd';
                val_out->val.d = mxGetScalar(result);
            } else if (mxIsChar(result)) {
                int buflen = mxGetM(result) + 1;
                val_out->type = 's';
                val_out->val.s = (char*) 
                    mempool_get(pool, buflen * sizeof(char));
                mxGetString(result, val_out->val.s, buflen);
            } else {
                char buf[256];
                val_out->type = 'e';
                sprintf(buf, "Invalid return type from %s", 
				node->v.call.name);
                add_error_trace_atcode(node, buf);
            }

            mxDestroyArray(result);
            mxFree(prhs);
#else            
            char buf[256];
            sprintf(buf, "Unrecognized function %s", node->v.call.name);
            add_error_trace_atcode(node, buf);
#endif
        }
    }
}

