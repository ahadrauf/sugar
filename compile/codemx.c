#include "mex.h"
#include "matrix.h"

#include <string.h>
#include <assert.h>
#include "codegen.h"
#include "codemx.h"

#define ELT_NAME_FIELD    0
#define ELT_MODEL_FIELD   1
#define ELT_NODES_FIELD   2
#define ELT_PARAMS_FIELD  3
#define ELT_ROT_FIELD     4
#define ELT_VARS_FIELD    5
#define ELT_NUM_FIELDS    6

#define NODE_NAME_FIELD   0
#define NODE_ELTS_FIELD   1
#define NODE_POS_FIELD    2
#define NODE_NUM_FIELDS   3

#define VAR_NAME_FIELD    0
#define VAR_TYPE_FIELD    1
#define VAR_OWNER_FIELD   2
#define VAR_NUM_FIELDS    3

#define NETLIST_ELTS_FIELD  0
#define NETLIST_NODES_FIELD 1
#define NETLIST_VARS_FIELD  2
#define NETLIST_DOF_FIELD   3
#define NETLIST_NUM_FIELDS  4

static void fill_param_mx(mxArray* result, int field_num, parameter_ir* param);
static void fill_element_mx(mxArray* result, int entry_num, element_ir* elt, 
                            int dof);
static void fill_node_mx(mxArray* result, int entry_num, node_ir* node);
	
static void fill_var_mx(mxArray* result, int entry_num, var_ir* var);
static mxArray* var_to_mx(var_ir** var_table, int num_vars);

static void fill_param_mx(mxArray* result, int field_num, parameter_ir* param)
{
    mxArray* param_result;
    if (param->value.type == 'd') {
       param_result = mxCreateDoubleMatrix(1,1,mxREAL);
       *mxGetPr(param_result) = param->value.val.d;
    } else if (param->value.type == 's') {
       param_result = mxCreateString(param->value.val.s); 
    } else {
        /*assert(0);*/
        mxAssert(0, "Unknown param type");
    }
    mxSetFieldByNumber(result,0,field_num,param_result);
}

/* Convert parameter list to Matlab.
 * Note: This is a little trickier than some of the other conversions,
 *   simply because the low-level representation of the netlist can
 *   have repeated entries for a particular parameter, while the
 *   Matlab structure cannot.  Later entries (typically process entries)
 *   may be ignored.
 */
mxArray* params_to_mx(parameter_ir* params)
{
    int dims[2] = {1,1};
    mxArray* result;
    const char** field_names;
    char* is_unique;
    parameter_ir* param;
    int param_count, unique_count;
    int i, j;

    /* Count the number of parameters and use that info for space
     * allocation.
     */
    param_count = 0; 
    param = params;
    while (param != NULL) {
        ++param_count;
	param = param->next;
    }

    field_names = (const char**) mxMalloc(param_count * sizeof(char*));
    is_unique = (char*) mxMalloc(param_count * sizeof(char));

    /* Record the names of previously unseen ("unique") parameters.
     * Indices of repeat parameters are noted, but a corresponding
     * entry is not added to the field_names list
     */
    param = params;
    unique_count = 0;
    for (i = 0; i < param_count; ++i) {
	is_unique[i] = 1;
	for (j = 0; j < unique_count; ++j) {
	    if (strcmp(field_names[j], param->name) == 0) {
		is_unique[i] = 0;
		break;
	    }
	}
	if (is_unique[i]) 
	    field_names[unique_count++] = param->name;
	param = param->next;
    }

    /* Fill in the parameter fields, skipping over any repeated
     * definitions.
     */
    result = mxCreateStructArray(2, dims, unique_count, field_names);
    param = params;
    j = 0;
    for (i = 0; i < param_count; ++i) {
	if (is_unique[i]) 
            fill_param_mx(result, j++, param);
	param = param->next;
    }

    mxFree(is_unique);
    mxFree(field_names);
    return result;
}

/* Fill in a particular element's data in the Matlab struct array
 */
static void fill_element_mx(mxArray* result, int entry_num, element_ir* elt,
                            int dof)
{
    mxArray* field_value;
    double* field_ptr;
    int i, j;
    char buf[256];

    field_value = mxCreateString(elt->name);
    mxSetFieldByNumber(result,entry_num, ELT_NAME_FIELD, field_value);
    
    strcpy(buf, "MF_");
    strcat(buf, elt->model);
    field_value = mxCreateString(buf);
    mxSetFieldByNumber(result,entry_num, ELT_MODEL_FIELD, field_value);
   
    field_value = mxCreateDoubleMatrix(1, elt->num_nodes, mxREAL);
    field_ptr = mxGetPr(field_value);
    for (j = 0; j < elt->num_nodes; ++j)
        field_ptr[j] = elt->node_ids[j]+1;
    mxSetFieldByNumber(result,entry_num, ELT_NODES_FIELD, field_value);

    field_value = params_to_mx(elt->params);
    mxSetFieldByNumber(result,entry_num, ELT_PARAMS_FIELD, field_value);

    field_value = mxCreateDoubleMatrix(3, 3, mxREAL);
    field_ptr = mxGetPr(field_value);
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
	    field_ptr[3*i+j] = elt->Q[3*j+i];
    mxSetFieldByNumber(result,entry_num, ELT_ROT_FIELD, field_value);

    field_value = mxCreateDoubleMatrix(1, elt->num_vars, mxREAL);
    field_ptr = mxGetPr(field_value);
    for (i = 0; i < elt->num_vars; ++i)
        field_ptr[i] = ((elt->var_ids[i] < dof) ? elt->var_ids[i]+1 : 0);
    mxSetFieldByNumber(result,entry_num, ELT_VARS_FIELD, field_value);
}

/* Transfer element table to Matlab structure array
 */
mxArray* element_to_mx(element_ir** elt_table, int num_elts, int dof)
{
    const char* field_names[] = 
        {"name", "model", "node_ids", "parameter", "R", "var_ids"};
    int dims[2];
    int i;
    mxArray* result;

    /* Create a cell structure array to put stuff into */
    dims[0] = num_elts;
    dims[1] = 1;
    result = mxCreateStructArray(2, dims, ELT_NUM_FIELDS, field_names);

    /* Fill in the cell structure */
    for (i = 0; i < num_elts; ++i)
        fill_element_mx(result, i, elt_table[i], dof);
    
    return result;
}

/* Fill in the data for a node structure
 */
static void fill_node_mx(mxArray* result, int entry_num, node_ir* node)
{
    mxArray* field_value;
    double* field_ptr;
    int j;

    field_value = mxCreateString(node->name);
    mxSetFieldByNumber(result,entry_num, NODE_NAME_FIELD, field_value);
    
    field_value = mxCreateDoubleMatrix(1, node->num_elts, mxREAL);
    field_ptr = mxGetPr(field_value);
    for (j = 0; j < node->num_elts; ++j)
        field_ptr[j] = node->elt_ids[j]+1;
    mxSetFieldByNumber(result,entry_num, NODE_ELTS_FIELD, field_value);

    field_value = mxCreateDoubleMatrix(3, 1, mxREAL);
    memcpy(mxGetPr(field_value), node->pos, 3 * sizeof(double));
    mxSetFieldByNumber(result,entry_num, NODE_POS_FIELD, field_value);
}

/* Convert a table of nodes to a Matlab structure
 */
mxArray* node_to_mx(node_ir** node_table, int num_nodes)
{
    const char* field_names[] = {"name", "elt_ids", "pos"};
    int dims[2];
    int i;
    mxArray* result;

    /* Create a cell structure array to put stuff into */
    dims[0] = num_nodes;
    dims[1] = 1;
    result = mxCreateStructArray(2, dims, NODE_NUM_FIELDS, field_names);

    /* Fill in the cell structure */
    for (i = 0; i < num_nodes; ++i)
        fill_node_mx(result, i, node_table[i]);
    
    return result;
}

/* Fill in the data for a var structure
 */
static void fill_var_mx(mxArray* result, int entry_num, var_ir* var)
{
    mxArray* field_value;
    double* field_ptr;
    int j;
    char type[2] = {0, 0};

    field_value = mxCreateString(var->name);
    mxSetFieldByNumber(result,entry_num, VAR_NAME_FIELD, field_value);
  
    type[0] = var->type; 
    field_value = mxCreateString(type);
    mxSetFieldByNumber(result,entry_num, VAR_TYPE_FIELD, field_value);

#ifdef MATLAB6
    field_value = mxCreateScalarDouble(var->owner_id + 1);
#else
    field_value = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field_value) = var->owner_id + 1;
#endif
    mxSetFieldByNumber(result,entry_num, VAR_OWNER_FIELD, field_value);
}


/* Transfer variable table to Matlab structure array
 */
static mxArray* var_to_mx(var_ir** var_table, int num_vars)
{
    const char* field_names[] = 
        {"name", "type", "owner"};
    int dims[2];
    int i;
    mxArray* result;

    /* Create a cell structure array to put stuff into */
    dims[0] = num_vars;
    dims[1] = 1;
    result = mxCreateStructArray(2, dims, VAR_NUM_FIELDS, field_names);

    /* Fill in the cell structure */
    for (i = 0; i < num_vars; ++i)
        fill_var_mx(result, i, var_table[i]);
    
    return result;
}

/* Put together a Matlab netlist structure
 */
mxArray* netlist_to_mx(netlist_ir* netlist)
{
    const char* field_names[] = {"elements", "nodes", "vars", "dof"};
    int dims[2] = {1,1};
    mxArray* result;
    mxArray* field_value;

    result = mxCreateStructArray(2, dims, NETLIST_NUM_FIELDS, field_names);

    field_value = element_to_mx(netlist->elt_table, netlist->num_elts,
                                netlist->num_dof);
    mxSetFieldByNumber(result,0, NETLIST_ELTS_FIELD, field_value);

    field_value = node_to_mx(netlist->node_table, netlist->num_nodes);
    mxSetFieldByNumber(result,0, NETLIST_NODES_FIELD, field_value);

    field_value = var_to_mx(netlist->var_table, netlist->num_vars);
    mxSetFieldByNumber(result,0, NETLIST_VARS_FIELD, field_value);

#ifdef MATLAB6
    field_value = mxCreateScalarDouble(netlist->num_dof);
#else
    field_value = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field_value) = netlist->num_dof;
#endif

    mxSetFieldByNumber(result,0, NETLIST_DOF_FIELD, field_value);
    return result;
}
