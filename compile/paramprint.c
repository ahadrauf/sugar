#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "paramprint.h"

static void print_expression(FILE* fp, code_t* node);
static void print_call(FILE* fp, code_t* node);

void print_params(FILE* fp, code_t* code)
{
    while (code != NULL) {
	if (code->tag == PARAM_NODE) {
	    fprintf(fp, "%s", code->v.def.name);
	    if (code->v.def.value) {
	        fprintf(fp, " ");
		print_expression(fp, code->v.def.value);
	    }
	    fprintf(fp, "\n");
	}
        code = code->next;
    }
}

/* Print an expression
 */
static void print_expression(FILE* fp, code_t* node)
{
    switch (node->tag) {
    case INT_NODE:
	fprintf(fp, "%d", node->v.ival);
        break;
    case DOUBLE_NODE:
	fprintf(fp, "%.18g", node->v.dval);
        break;
    case STRING_NODE:
	fprintf(fp, "\"%s\"", node->v.sval);
	break;
    case CALL_NODE:
        print_call(fp, node);
        break;
    case VAR_NODE:
	fprintf(fp, "%s", node->v.var->v.def.name);
        break;
    default:
        assert(0);
    }

}

/* Print a function call.
 */
static void print_call(FILE* fp, code_t* node)
{
    code_t* args = node->v.call.args;

    if (strcmp(node->v.call.name, "unary-") == 0) {
        fprintf(fp, "(-");
	print_expression(fp, args);
	fprintf(fp, ")");
    } else if (strcmp(node->v.call.name, "!") == 0) {
        fprintf(fp, "(!");
	print_expression(fp, args);
	fprintf(fp, ")");
	
    } else if (strcmp(node->v.call.name, "+") == 0 ||
               strcmp(node->v.call.name, "-") == 0 ||
               strcmp(node->v.call.name, "*") == 0 ||
               strcmp(node->v.call.name, "/") == 0 ||
               strcmp(node->v.call.name, "^") == 0 ||
               strcmp(node->v.call.name, "&&") == 0 ||
               strcmp(node->v.call.name, "||") == 0 ||
               strcmp(node->v.call.name, "==") == 0 ||
               strcmp(node->v.call.name, "!=") == 0 ||
               strcmp(node->v.call.name, "<") == 0 ||
               strcmp(node->v.call.name, ">") == 0 ||
               strcmp(node->v.call.name, "<=") == 0 ||
               strcmp(node->v.call.name, ">=") == 0) {
        fprintf(fp, "(");
	print_expression(fp, args);
        fprintf(fp, "%s", node->v.call.name);
	print_expression(fp, args->next);
	fprintf(fp, ")");

    } else {
        fprintf(fp, "%s(", node->v.call.name);
        while (args) {
	    print_expression(fp, args);
	    if (args->next) 
	        fprintf(fp, ",");
	    args = args->next;
        }
        fprintf(fp, ")");
    }
}

