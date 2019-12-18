#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "lexer.h"
#include "parse.h"
#include "evalexpr.h"
#include "shash.h"
#include "codestack.h"
#include "codegen.h"

#include "position.h"
#include "varassign.h"
#include "modelmgr.h"

#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

static mempool_t code_pool;
static shash_t param_hash;
static element_ir** elt_tail;
static node_ir** node_tail;

static int element_count;
static int node_count;

static void gen_model_init(element_ir* elts);
static netlist_ir* gen_tabulate(node_ir* node_list, element_ir* elt_list);
static void gen_statement(code_t* node);
static void gen_statements(code_t* node);
static void gen_def(code_t* node);
static void gen_param(code_t* node);
static parameter_ir* cat_process(parameter_ir* params, code_t* process);
static void gen_forloop(code_t* node);
static void gen_ifstmt(code_t* node);
static void gen_element(code_t* node);
static void gen_subnet_element(code_t* node);

static void gen_matlab_element(element_ir* elt, code_t* node);
static void gen_element_name(element_ir* elt, code_t* node);
static void gen_element_nodes(element_ir* elt, code_t* node_names);
static parameter_ir* gen_param_list(code_t* params);
static void gen_rotation(element_ir* elt, code_t* args);
static int is_rotation(char* s);

void add_error_trace_atcode(code_t* code, char* errmsg)
{
    add_error_atcode(code, errmsg);
    frame_trace();
}

void add_warning_trace_atcode(code_t* code, char* errmsg)
{
    add_warning_atcode(code, errmsg);
    frame_trace();
}

/* Set up code generation, someone else does the work
 */
netlist_ir* gen_code(mempool_t pool, code_t* code, shash_t params)
{
    netlist_ir* netlist = NULL;
    element_ir* elt_list = NULL;
    node_ir* node_list = NULL;
    
    code_pool = pool;
    param_hash = params;
    frame_stack_create(code_pool, get_num_ids());

    elt_tail = &elt_list;
    node_tail = &node_list;
    
    /* Initialize element and node counters */
    element_count = 0;
    node_count = 0;
    
    /* Generate the netlist */
    gen_statements(code);
    netlist = gen_tabulate(node_list, elt_list);
    
    frame_stack_destroy();

    if (no_errors()) {
        gen_model_init(elt_list);
        position_nodes(netlist); 
        assign_var_indices(pool, netlist);
    }
    return netlist;
}


/* Initialize all model functions */
static void gen_model_init(element_ir* elts)
{
    while (elts != NULL) {
        elts->modelfun = model_lookup(elts->model);
        if (elts->modelfun && elts->modelfun->init_data) {
            elts->data = (elts->modelfun->init_data)(code_pool, elts);
        } else {
            elts->data = elts;
        }
        elts = elts->next;
    }
}


/* Create indexed arrays for the elements and nodes; generate links
 * from nodes back to referring elements
 */
static netlist_ir* gen_tabulate(node_ir* node_list, element_ir* elt_list)
{
    netlist_ir* netlist;
    node_ir* node;
    element_ir* elt;
    int i, j;

    node_ir** node_table;
    element_ir** elt_table;

    netlist = (netlist_ir*)
            mempool_get(code_pool, sizeof(netlist_ir));

    /* Build node array */
    node = node_list;
    node_table = (node_ir**)
            mempool_get(code_pool, node_count * sizeof(node_ir*));
    for (i = 0; i < node_count; ++i) {
        node_table[i] = node;
        node = node->next;
    }
    
    /* Build element array and count links needed */
    elt = elt_list;
    elt_table = (element_ir**)
            mempool_get(code_pool, element_count * sizeof(element_ir*));
    for (i = 0; i < element_count; ++i) {
        elt_table[i] = elt;
        for (j = 0; j < elt->num_nodes; ++j)
            node_table[elt->node_ids[j]]->num_elts++;
        elt = elt->next;
    }

    /* Allocate elt_ids fields */
    for (i = 0; i < node_count; ++i) 
        node_table[i]->elt_ids = (int *)
                mempool_get(code_pool, node_table[i]->num_elts*sizeof(int));

    /* Bucket in the links */
    for (i = 0; i < element_count; ++i) {
        elt = elt_table[i];
        for (j = 0; j < elt->num_nodes; ++j)
            *node_table[elt->node_ids[j]]->elt_ids++ = i;
    }

    /* Reset elt_ids fields */
    for (i = 0; i < node_count; ++i) 
        node_table[i]->elt_ids -= node_table[i]->num_elts;

    /* Pack everything into the netlist struct for return */
    netlist->elt_table = elt_table;
    netlist->num_elts = element_count;
    netlist->node_table = node_table;
    netlist->num_nodes = node_count;

    return netlist;
}


/* Loop through the statements and generate code (duh) */
static void gen_statements(code_t* code)
{
    while (code != NULL) {
        gen_statement(code);
        code = code->next;
    }
}


/* Switch to the appropriate subroutine */
static void gen_statement(code_t* node)
{
    val_t temp_val;
    switch (node->tag) {
    case CALL_NODE:
        eval_expr(code_pool, &temp_val, node);
	break;
    case DEF_NODE:
        gen_def(node);
        break;
    case PARAM_NODE:
        gen_param(node);
        break;
    case BLOCK_NODE:
        gen_statements(node->v.block);
        break;
    case ELEMENT_NODE:
        gen_element(node);
        break;
    case FOR_NODE:
        gen_forloop(node);
        break;
    case IF_NODE:
        gen_ifstmt(node);
        break;
    }
}


/* Handle a normal variable def node by assigning the appropriate entry
 * in the variable table.
 */
static void gen_def(code_t* node)
{
    assert(node->v.def.value != NULL);
    eval_expr(code_pool, frame_var_entry(node->v.def.id, node->v.def.scope), 
              node->v.def.value);
}


/* Generate a parameter node... like a def node, but with
 * the possibility of having overridden parameters.
 * I haven't dealt with them here yet.
 */
static void gen_param(code_t* node)
{
    val_t* value = frame_var_entry(node->v.def.id, node->v.def.scope);

    if (param_hash != NULL && 
            shash_get(param_hash, node->v.def.name, value) != NULL) {
        /* Don't need to do anything! */
    } else if (node->v.def.value) {
        eval_expr(code_pool, value, node->v.def.value);
    } else {
        char buf[256];
        sprintf(buf, "Undefined parameter %s has no default",
                node->v.def.name);
        add_error_trace_atcode(node, buf);
    }
}


/* Link process parameters onto the end of a parameter list
 */
static parameter_ir* cat_process(parameter_ir* params, code_t* process)
{
    parameter_ir* tail = params;
    parameter_ir* param_list;
    
    param_list = gen_param_list(process->v.process.defs);
    
    if (tail == NULL)
        return param_list;
    
    while (tail->next != NULL)
        tail = tail->next;

    tail->next = param_list;
    
    return params;
}


/* Evaluate the lower and upper bounds, then do the loop,
 * re-binding the index variable at each iteration.
 * If the bounds are not integer, we'll get snippy.
 */
static void gen_forloop(code_t* node)
{
    val_t lower;
    val_t upper;
    val_t* idx;
    int i;
    int lbound;
    int ubound;

    eval_expr(code_pool, &lower, node->v.forloop.lower);
    eval_expr(code_pool, &upper, node->v.forloop.upper);

    if (lower.type != 'd' || lower.val.d != (int) lower.val.d || 
            upper.type != 'd' || upper.val.d != (int) upper.val.d) {
        add_error_trace_atcode(node, "Loop does not have integer bounds!");
        return;
    }

    lbound = (int) lower.val.d;
    ubound = (int) upper.val.d;
    idx = frame_var_entry(node->v.forloop.idx->v.def.id,
                          node->v.forloop.idx->v.def.scope);
    idx->type = 'd';

    for (i = lbound; i <= ubound; ++i) {
        idx->val.d = i;
        gen_statements(node->v.forloop.body->v.block);
    }
}


/* Evaluate the condition and do the appropriate clause
 */
static void gen_ifstmt(code_t* node)
{
    val_t cond;
    int which_case;

    eval_expr(code_pool, &cond, node->v.ifstmt.cond);

    if (cond.type == 'd')
        which_case = (cond.val.d != 0);
    else {
        add_error_trace_atcode(node, "Invalid expression for conditional");
        return;
    }

    if (which_case) {
        gen_statements(node->v.ifstmt.then_clause);
    } else if (node->v.ifstmt.else_clause) {
        gen_statements(node->v.ifstmt.else_clause);
    }
}

/* Pick out whether this is a "normal" element or a subnet,
 * and call the appropriate subroutine.
 */
static void gen_element(code_t* node)
{
    if (node->v.element.model_node)
        gen_subnet_element(node);
    else {
        element_ir* elt = (element_ir*) 
                mempool_get(code_pool, sizeof(element_ir));
        gen_matlab_element(elt, node);
        elt->next = NULL;
        *elt_tail = elt;
        elt_tail = &(elt->next);
    }
}


/* Generate a subnet instance.
 */
static void gen_subnet_element(code_t* node)
{
    code_t* subnet = node->v.element.model_node;
    char* parent_name = parent_process->v.process.name;
    code_t* parent_defs = parent_process->v.process.defs;

    /* "Push" inherited process on top of parent process */
    parent_process->v.process.name = node->v.element.process->v.process.name;
    parent_process->v.process.defs = node->v.element.process->v.process.defs;

    /* Generate subnet definitions, bind formal to actual node names,
     * then generate actual subnet code */
    frame_push(node, cat_process(NULL, 
                        node->v.element.process));
    gen_statements(subnet->v.subnet.body);
    frame_pop();

    /* "Pop" inherited process from parent process stack */
    parent_process->v.process.name = parent_name;
    parent_process->v.process.defs = parent_defs;
}


/* Generate a "normal" element internal representation.
 */
static void gen_matlab_element(element_ir* elt, code_t* node)
{
    parameter_ir* params;

    elt->id = ++element_count;
    elt->model = mempool_strdup(code_pool, node->v.element.model);
    gen_element_name(elt, node->v.element.name);
    gen_element_nodes(elt, node->v.element.nodes);
    elt->params = gen_param_list(node->v.element.params);
    
    params = elt->params;
    while (params != NULL) {
        if (params->value.type == 'u') {
	    char buf[128];
	    sprintf(buf, "Parameter %s undefined", params->name);
            add_error_trace_atcode(node, buf);
	}
	params = params->next;
    }
    
    if (node->v.element.process->v.process.defs != NULL) {
        elt->params = cat_process(elt->params, 
                        node->v.element.process);
    }
    gen_rotation(elt, node->v.element.params);
}


/* Write out a fully qualified and indexed name (for an element).
 */
static void gen_element_name(element_ir* elt, code_t* node)
{
    char buf[512];
    *buf = '\0';
    eval_name(code_pool, buf, node);
    elt->name = frame_scoped_name(code_pool, buf);
}


/* Generate global index assignments.
 */
static void gen_element_nodes(element_ir* elt, code_t* node_names)
{
    code_t* node;
    int node_index;

    elt->num_nodes = 0;
    for (node = node_names; node != NULL; node = node->next)
        elt->num_nodes++;

    elt->node_ids = mempool_get(code_pool, sizeof(int) * elt->num_nodes);
    node_index = 0;

    for (node = node_names; node != NULL; node = node->next) {
        
        int node_id;
        char buf[512];
        *buf = '\0';

        eval_name(code_pool, buf, node);
        node_id = frame_node_index(buf);
        elt->node_ids[node_index++] = node_id-1;
    }
}



/* TODO: This is still sort of kludgy
 * If this is a new node, ... */
void record_new_node(int id, char* name)
{
    node_ir* new_node = (node_ir*)
                mempool_get(code_pool, sizeof(node_ir));
    new_node->id = id;
    new_node->name = frame_scoped_name(code_pool, name);
    new_node->elt_ids = NULL;
    new_node->num_elts = 0;
    new_node->pos[0] = new_node->pos[1] = new_node->pos[2] = 0;
    new_node->next = NULL;
    *node_tail = new_node;
    node_tail = &(new_node->next);
    ++node_count;
}


/* Given an existing rotation and a list of arguments which might contain
 * o[xyz] bindings, add on a new rotation.
 */
static void gen_rotation(element_ir* elt, code_t* args)
{
    double ox = 0, oy = 0, oz = 0;
    double Q[9];

    /* Extract ox, oy, oz.  
     */
    while (args != NULL) {
        char* param_name = args->v.def.name;
        val_t val;
        if (is_rotation(param_name)) {
            if (param_name[1] == 'x') {
                eval_expr(code_pool, &val, args->v.def.value);
                if (val.type == 'd')
	            ox = val.val.d;
		else
		    add_error_trace_atcode(args, "ox must be a number");
#ifndef USE_DEGREES
		if ((ox > 2*M_PI || ox < -2*M_PI) && 
                      (360.0/ox) == (int) (360.0/ox))
		    add_warning_trace_atcode(args, "ox should be in radians");
#endif
            } else if (param_name[1] == 'y') {
                eval_expr(code_pool, &val, args->v.def.value);
                if (val.type == 'd')
	            oy = val.val.d;
		else
		    add_error_trace_atcode(args, "oy must be a number");
#ifndef USE_DEGREES
		if ((oy > 2*M_PI || oy < -2*M_PI) && 
                      (360.0/oy) == (int) (360.0/oy))
		    add_warning_trace_atcode(args, "oy should be in radians");
#endif
            } else if (param_name[1] == 'z') {
                eval_expr(code_pool, &val, args->v.def.value);
                if (val.type == 'd')
	            oz = val.val.d;
		else
		    add_error_trace_atcode(args, "oz must be a number");
#ifndef USE_DEGREES
		if ((oz > 2*M_PI || oz < -2*M_PI) && 
                      (360.0/oz) == (int) (360.0/oz))
		    add_warning_trace_atcode(args, "oz should be in radians");
#endif
            }
        }
        args = args->next;
    }

    frame_rot2local(Q, ox, oy, oz);
    memcpy(elt->Q, Q, 9*sizeof(double));
}


/* Generate parameter list representation for non-rotation params
 */
static parameter_ir* gen_param_list(code_t* params)
{
    parameter_ir* params_ir = NULL;
    parameter_ir* param_ptr;
    parameter_ir** param_pptr = &params_ir;

    while (params != NULL) {
        char* param_name = params->v.def.name;
        if (!is_rotation(param_name)) {
                
            param_ptr = (parameter_ir*)
                    mempool_get(code_pool, sizeof(parameter_ir));
            param_ptr->name = mempool_strdup(code_pool, param_name);
            eval_expr(code_pool, &(param_ptr->value), params->v.def.value);
            param_ptr->next = NULL;

	    if (param_ptr->value.type != 'u') {
                *param_pptr = param_ptr;
                param_pptr = &(param_ptr->next);
	    }
        }
        params = params->next;
    }

    return params_ir;
}


/* Helper routine: check if string is the name of a rotation param
 */
static int is_rotation(char* s)
{
    return (s[0] == 'o' && (s[1] == 'x' || s[1] == 'y' || s[1] == 'z') &&
            s[2] == '\0');
}

