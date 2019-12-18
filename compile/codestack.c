#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "shash.h"
#include "codestack.h"
#include "codegen.h"
#include "lexer.h"

/* This is really obnoxious, but MSC apparently thinks it doesn't need to
 * implement the standard C headers.  In particular, it doesn't define M_PI.
 * Meh.
 */
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif


#define NODE_TABLE_SIZE 255

typedef struct frame_t {
    struct frame_t* prev_frame;     /* Previous stack frame             */
    char*           name_prefix;    /* Prefix for creating scoped names */
    code_t*         active_element; /* Child subnet in progress         */
    val_t*          var_table;      /* Table of variable values         */
    shash_t         node_table;     /* Table of node numbers            */
    double          Q[9];           /* Rotation to subnet local coords  */
} frame_t;

static mempool_t framepool;
static mempool_t code_pool;
static frame_t* current_frame = NULL;
static frame_t* global_frame;
static int current_node_index;

static frame_t* frame_alloc(int var_table_size);
static void frame_init_names(code_t* element);
static void frame_bind_names(frame_t* new_frame, code_t* subnet, 
                             code_t* element);
static void frame_bind_args(frame_t* new_frame, code_t* formals, 
		            code_t* element, parameter_ir* process_params);
static void frame_subnet_rotation(frame_t* new_frame, code_t* args);


/* Helper routine: check if string is the name of a rotation param
 */
static int is_rotation(char* s)
{
    return (s[0] == 'o' && (s[1] == 'x' || s[1] == 'y' || s[1] == 'z') &&
            s[2] == '\0');
}

/* Allocate frame pool and global stack frame, and set the change of
 * coordinates for the top level to be the identity.
 */
void frame_stack_create(mempool_t pool, int num_global_vars)
{
    code_pool = pool;
    framepool = mempool_create(MEMPOOL_DEFAULT_SPAN);

    global_frame = frame_alloc(num_global_vars);
    current_frame = global_frame;

    current_frame->name_prefix = "";
    memset(current_frame->Q, 0, 9*sizeof(double));
    current_frame->Q[0] = current_frame->Q[4] = current_frame->Q[8] = 1;

    current_node_index = 1;
}

/* Destroy the stack.
 */
void frame_stack_destroy(void)
{
    mempool_destroy(framepool);
    current_frame = NULL;
    global_frame = NULL;
}

/* Allocate a new frame with variable table of the given size and
 * put the new frame at the top of the stack.
 */
static frame_t* frame_alloc(int var_table_size)
{
    frame_t* new_frame = (frame_t*) mempool_geth(framepool, sizeof(frame_t));
    new_frame->var_table = (val_t*) 
            mempool_get(framepool, sizeof(val_t) * var_table_size);
    new_frame->node_table =
            shash_pcreate(NODE_TABLE_SIZE, sizeof(int), framepool);
    new_frame->prev_frame = NULL;
    return new_frame;
}

/* Run through the actual names passed to a subnet instance in order
 * to ensure space is created for them and indices assigned *before*
 * a new stack frame is allocated.
 */
static void frame_init_names(code_t* element)
{
    code_t* actual_name = element->v.element.nodes;
    while (actual_name != NULL) {
        char buf[256];
        eval_name(code_pool, buf, actual_name);
        frame_node_index(buf);
        actual_name = actual_name->next;
    }
}

/* Bind node indices for actual node names in the parent environment
 * to formal node names in the new environment.  Also set the new name
 * prefix.
 */
static void frame_bind_names(frame_t* new_frame, code_t* subnet, 
                             code_t* element)
{
    code_t* actual_name = element->v.element.nodes;
    code_t* formal_name = subnet->v.subnet.nodes;
    char buf[256];

    while (actual_name != NULL && formal_name != NULL) {

        int node_index;
        
        eval_name(code_pool, buf, actual_name);
        node_index = frame_node_index(buf);

        eval_name(code_pool, buf, formal_name);
        shash_set(new_frame->node_table, buf, &node_index);
        
        actual_name = actual_name->next;
        formal_name = formal_name->next;
    }

    if (actual_name != NULL || formal_name != NULL) {
        char buf[256];
        sprintf(buf, "Mismatched number of nodes in call to %s", 
               subnet->v.subnet.model);
        add_error_trace_atcode(element, buf);
    }

    eval_name(code_pool, buf, element->v.element.name);
    new_frame->name_prefix = (char*) 
	    mempool_get(framepool, strlen(current_frame->name_prefix) +
			strlen(buf) + 2);
    strcpy(new_frame->name_prefix, current_frame->name_prefix);
    strcat(new_frame->name_prefix, buf);
    strcat(new_frame->name_prefix, ".");
}

/* Bind actual parameter values to the formal arguments in the
 * variable table for the new frame.
 */
static void frame_bind_args(frame_t* new_frame, code_t* formals, 
		            code_t* element, parameter_ir* process_params)
{
    code_t* actuals = element->v.element.params;
    code_t* formal = formals;
    code_t* actual;
    val_t* var_table = new_frame->var_table;
    char buf[256];
    
    while (formal != NULL) {

	int argument_bound = 0;   /* Is the argument bound yet?   */
	int argument_unk = 0;     /* Is the value an "undefined"? */

        /* --- Find matching actual */
        actual = actuals;
        while (actual != NULL && 
               strcmp(formal->v.def.name, actual->v.def.name) != 0)
            actual = actual->next;

        /* --- If there was a match, bind the value in */
        if (actual != NULL) {
            eval_expr(code_pool, &(var_table[formal->v.def.id]), 
                      actual->v.def.value);

	    argument_bound = (var_table[formal->v.def.id].type != 'u');
	} 

        /* --- If no match, use default (if it exists) 
	 *     NOTE: The default should be evaluated in the new stack
	 *       frame, since it may refer to previous arguments.
	 */
	if (!argument_bound && formal->v.def.value) {
	    frame_t* saved_current_frame = current_frame;

	    current_frame = new_frame;
            eval_expr(code_pool, &(var_table[formal->v.def.id]), 
                      formal->v.def.value);
	    current_frame = saved_current_frame;
	    
	    argument_bound = 1;
	    argument_unk = (var_table[formal->v.def.id].type == 'u');
        } 
	
        /* --- At last resort, try the process info */
	if (!argument_bound || argument_unk) {

            parameter_ir* actual_param = process_params;
            while (actual_param != NULL && 
                   strcmp(formal->v.def.name, actual_param->name) != 0) {
	        actual_param = actual_param->next;
	    }

            if (actual_param != NULL) {
		var_table[formal->v.def.id] = actual_param->value;
		argument_bound = 1;
		argument_unk = 0;
            } 
        }

	if (!argument_bound) {
            sprintf(buf, "Unbound argument %s in subnet instance",
                    formal->v.def.name);
            add_error_trace_atcode(element, buf);
        }

        formal = formal->next;
    }
}

/* Given an existing rotation and a list of arguments which might contain
 * o[xyz] bindings, add on a new rotation.
 */
static void frame_subnet_rotation(frame_t* new_frame, code_t* args)
{
    double ox = 0, oy = 0, oz = 0;

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
		    add_warning_atcode(args, "ox should be in radians");
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
		    add_warning_atcode(args, "oy should be in radians");
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
		    add_warning_atcode(args, "oz should be in radians");
#endif
            }
        }

        args = args->next;
    }

    frame_rot2local(new_frame->Q, ox, oy, oz);
}

void frame_push(code_t* element, parameter_ir* process_params)
{
    code_t* subnet = element->v.element.model_node;
    frame_t* new_frame;

    /* -- Instantiate info for all actual nodes *before* creating
     *    new frame (don't want them destroyed on free!) */
    frame_init_names(element);

    /* --- Allocate frame */
    new_frame = frame_alloc(subnet->v.subnet.def_count);
    new_frame->prev_frame = current_frame;

    /* --- Bind names */
    frame_bind_names(new_frame, subnet, element);


    /* --- Bind actuals to formals */
    frame_bind_args(new_frame, subnet->v.subnet.params, 
                    element, process_params);

    /* --- Generate rotation */
    frame_subnet_rotation(new_frame, element->v.element.params);

    /* --- Link onto top of stack */
    current_frame->active_element = element;
    current_frame = new_frame;
}

/* Pop a frame from the top of the stack.
 */
void frame_pop(void)
{
    frame_t* dead_frame = current_frame;
    current_frame = current_frame->prev_frame;
    mempool_freeh(dead_frame); 
}

/* Print a stack trace in case of error
 */
void frame_trace(void)
{
    if (current_frame != global_frame) {
        frame_t* trace_frame = current_frame->prev_frame;
        while (trace_frame != NULL) {
	    code_t* element = trace_frame->active_element;
	    add_error_stack(element->fileno, element->lineno,
	                    element->v.element.model);
	    trace_frame = trace_frame->prev_frame;
        }
    }
}

/* Allocate and create a fully scoped text version of a name.
 */
char* frame_scoped_name(mempool_t pool, char* name)
{
    char* scoped_name;

    scoped_name = (char*) 
	    mempool_get(pool, strlen(current_frame->name_prefix) +
			strlen(name) + 1);
    strcpy(scoped_name, current_frame->name_prefix);
    strcat(scoped_name, name);
    return scoped_name;
}

/* Look up the index of a node in the local table.  If no such node,
 * assign a new index.
 */
int frame_node_index(char* name)
{
    int* node_index = shash_index(current_frame->node_table, name);
    if (*node_index == 0) {
        *node_index = current_node_index++;
        record_new_node(*node_index, name);
    }

    return *node_index;
}

int get_node_count()
{
    return current_node_index-1;
}

/* Look up the entry for a variable in the given scope (local or global)
 */
val_t* frame_var_entry (int var_id, int scope)
{
    return ((scope == 0) ? &(global_frame->var_table[var_id]) :
                           &(current_frame->var_table[var_id]));
}

/* Generates rotation to local (rx*rz*ry) in column major arrangement.
 * Conveniently, this is rotation to global in row-major.
 */
void frame_rot2local(double* Q, double ox, double oy, double oz)
{
    double cx, sx, cy, sy, cz, sz;

#ifdef USE_DEGREES    
    ox = (ox/180)*M_PI;
    oy = (oy/180)*M_PI;
    oz = (oz/180)*M_PI;
#endif
    
    cx = cos(ox);  cy =  cos(oy);  cz = cos(oz);
    sx = sin(ox);  sy = -sin(oy);  sz = sin(oz);

    memcpy(Q, current_frame->Q, 9*sizeof(double));
    
    #define Qij(i,j) (Q[(i)-1+((j)-1)*3])
    
    #define rotv(c,s,i1,i2,j) \
        { double v1 = c*Qij(i1,j) + s*Qij(i2,j); \
          double v2 = -s*Qij(i1,j) + c*Qij(i2,j); \
          Qij(i1,j) = v1; Qij(i2,j) = v2; }
        
    #define rotm(c,s,i1,i2) \
        rotv(c,s,i1,i2,1) \
        rotv(c,s,i1,i2,2) \
        rotv(c,s,i1,i2,3) 

    rotm(cy,sy,1,3)
    rotm(cz,sz,1,2)
    rotm(cx,sx,2,3)

    #undef rotm
    #undef rotv
    #undef Qij
}
