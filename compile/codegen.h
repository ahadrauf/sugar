#ifndef __CODEGEN_H
#define __CODEGEN_H

/*
 * Data structures and routines associated with the low-level representation
 * of the netlist.  This low-level representation does not retain any loops,
 * subnets, conditionals, etc.  Instead, it only describes the elements,
 * nodes, and variables that go into the device description.
 *
 * There are a few unobvious points about the low level representation.
 *
 * - Process parameters are appended to the tail of each element's parameter
 *   list.  Since the final list segments are shared, they should not be
 *   modified or overwritten.
 *   
 * - The variables are re-ordered during index assignment so that ungrounded
 *   degrees of freedom come first.  The first num_dof entries in the
 *   variable table correspond to ungrounded degrees of freedom, in the order
 *   they appear in the global system.  The remaining num_vars - num_dof
 *   entries correspond to grounded variables.
 *
 * Functions:
 *  gen_code -- Generate low-level representation from tree representation
 *  record_new_node -- Used by codestack to record a new node not seen before 
 */

#include <stdio.h>
#include "evalexpr.h"
#include "shash.h"

struct model_t;

/* Node (or branch) variable structure */
typedef struct var_ir {
    char* name;            /* Name of the variable */
    char type;             /* Type 'u' == ungrounded nodal variable, 
			           'g' == grounded nodal variable, 
				   'b' == branch variable */
    int owner_id;          /* Index of node (type == 'u' or 'g') or
			      element (type == 'b') with which variable
			      is associated. */
    int index;             /* Global index of the variable (-1 == grounded) */
    struct var_ir* next; 
} var_ir;

/* Element parameter structure */
typedef struct parameter_ir {
    char* name;                 /* Name */
    val_t value;                /* Assigned value */
    struct parameter_ir* next;
} parameter_ir;

/* Element structure */
typedef struct element_ir {
    int id;                     /* Element number */
    char* name;                 /* Element name */
    char* model;                /* Model name */
    char* process_name;         /* Process name */
    int* node_ids;              /* Indices of nodes */
    int num_nodes;              /* Number of node_ids */
    int* var_ids;               /* Indices of variables */
    int num_vars;               /* Number of var_ids */
    struct parameter_ir* params;/* Element parameters */
    double Q[9];                /* Rotation to local, column major order */
    void* data;                 /* Storage used by the model functions */
    struct model_t* modelfun;   /* Model function structure */
    struct element_ir* next;
} element_ir;

/* Node structure */
typedef struct node_ir {
    int id;                 /* Node identifier */
    char* name;             /* Name of the node */
    int* elt_ids;           /* Elements that use this node */
    int num_elts;           /* Number of elt_ids */
    double pos[3];          /* Undisplaced location of mechanical node */
    struct node_ir* next;
} node_ir;

/* The elements, nodes, and variables that go into a netlist.
 * The first num_dof variables are ungrounded; the remainder
 * are grounded */
typedef struct netlist_ir {
    element_ir** elt_table;
    int num_elts;
    node_ir** node_table;
    int num_nodes;
    var_ir** var_table;
    int num_vars;
    int num_dof;
} netlist_ir;

/* Generate low-level representation from tree representation */
netlist_ir* gen_code(mempool_t code_pool, code_t* code, shash_t params);

/* Used by codestack to record a new node not seen before */
void record_new_node(int id, char* name);

/* Report a generation error with stack trace */
void add_error_trace_atcode(code_t* code, char* errmsg);
void add_warning_trace_atcode(code_t* code, char* errmsg);

#endif /* __CODEGEN_H */
