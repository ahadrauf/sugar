#include "shash.h"
#include "modelmgr.h"
#include "models.h"

#ifdef __SUGAR_MEX
  #include "mex.h"
  #include "matrix.h"
  #include "codemx.h"
#endif /* __SUGAR_MEX */

#define MODEL_TABLE_SIZE 127
static shash_t model_table;

#ifdef __SUGAR_MEX
static void* matlab_init_data(mempool_t data_pool, element_ir* elt);
static void matlab_find_position(void* data, mempool_t pool,
                                 element_pos_t* elt_pos);
static void matlab_post_position(void* data, node_ir** node_table);
static void matlab_add_vars(void* data, var_ir** tmp_var_table,
                            mempool_t code_pool);

static model_t matlab_model = {
    matlab_init_data, 
    matlab_find_position, 
    matlab_post_position,
    NULL,
    NULL,
    matlab_add_vars,
    NULL,
    NULL
};
#endif /* __SUGAR_MEX */


void model_manager_init()
{
    model_table = shash_create(MODEL_TABLE_SIZE, sizeof(model_t));
#ifndef __SUGAR_MEX
/* TODO */
    anchor_register();
    beams_register();
    force_register();
    pos_register();
    gap_register();
#endif
}


void model_register(char* name, model_t* model)
{
#ifdef __SUGAR_MEX
    if (model->init_data == NULL)
        model->init_data = matlab_init_data;
    if (model->find_position == NULL)
        model->find_position = matlab_find_position;
    if (model->display == NULL)
        model->find_position = matlab_find_position;
#endif /* __SUGAR_MEX */
    shash_set(model_table, name, model);
}


model_t* model_lookup(char* name)
{
    if (shash_contains(model_table, name))
        return (model_t*) shash_index(model_table, name);
    else
#ifdef __SUGAR_MEX
        return &matlab_model;
#else
        return NULL;
#endif
}


void model_manager_shutdown()
{
    shash_destroy(model_table);
}


#ifdef __SUGAR_MEX

static void* matlab_init_data(mempool_t data_pool, element_ir* elt)
{
    return elt;
}


static void matlab_find_position(void* data, mempool_t pool, 
                                 element_pos_t* elt_pos)
{
    char command_name[128];
    mxArray* result;
    mxArray* prhs[3];
    int i, j, nelements;
    double* field_value;
    element_ir* elt = (element_ir*) data;

    strcpy(command_name, "MF_");
    strcat(command_name, elt->model);

    prhs[0] = mxCreateString("pos");
    
    prhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
    field_value = mxGetPr(prhs[1]);
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            field_value[i*3+j] = elt->Q[j*3+i];

    prhs[2] = params_to_mx(elt->params);
    
    mexCallMATLAB(1, &result, 3, prhs, command_name);

    nelements = mxGetNumberOfElements(result);
    elt_pos->num_mech_nodes = nelements / 3;
    elt_pos->relpos = (double*) 
            mempool_get(pool, nelements * sizeof(double));
    memcpy(elt_pos->relpos, mxGetPr(result), nelements * sizeof(double));

    mxDestroyArray(result);
    mxDestroyArray(prhs[2]); 
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[0]);
}

/* TODO: Finish up */
static void matlab_post_position(void* data, node_ir** node_table)
{
    char command_name[128];
    mxArray* result;
    mxArray* prhs[3];
    int i, j, nelements;
    double* field_value;
    element_ir* elt = (element_ir*) data;
    node_ir** local_node_table;

    strcpy(command_name, "MF_");
    strcat(command_name, elt->model);

    prhs[0] = mxCreateString("postpos");
    
    prhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
    field_value = mxGetPr(prhs[1]);
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            field_value[i*3+j] = elt->Q[j*3+i];

    prhs[2] = params_to_mx(elt->params);
    prhs[3] = mxCreateDoubleMatrix(0, 0, mxREAL);
    prhs[4] = mxCreateDoubleMatrix(0, 0, mxREAL);

    local_node_table = mxMalloc(elt->num_nodes * sizeof(node_ir*));
    for (i = 0; i < elt->num_nodes; ++i) {
        local_node_table[i] = node_table[elt->node_ids[i]];
    }
    prhs[5] = node_to_mx(local_node_table, elt->num_nodes);

    mexCallMATLAB(1, &result, 5, prhs, command_name);

/*
    nelements = mxGetNumberOfElements(result);
    elt_pos->num_mech_nodes = nelements / 3;
    elt_pos->relpos = (double*) 
            mempool_get(pool, nelements * sizeof(double));
    memcpy(elt_pos->relpos, mxGetPr(result), nelements * sizeof(double));
*/

    mxDestroyArray(result);
    mxFree(local_node_table);
    mxDestroyArray(prhs[4]); 
    mxDestroyArray(prhs[3]); 
    mxDestroyArray(prhs[2]); 
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[0]);
}

static void matlab_add_vars(void* data, var_ir** tmp_var_table, 
                            mempool_t code_pool)
{
    char command_name[128];
    mxArray* result;
    mxArray* prhs[1];
    mxArray* node;
    mxArray* ground;
    mxArray* branch;
    element_ir* elt = (element_ir*) data;

    int num_vars = 0;
    int var_ids[64]; /* TODO: Hope you don't have many nodal vars */

    strcpy(command_name, "MF_");
    strcat(command_name, elt->model);

    prhs[0] = mxCreateString("vars");
    
    mexCallMATLAB(1, &result, 1, prhs, command_name);

    node = mxGetField(result, 0, "dynamic");
    ground = mxGetField(result, 0, "ground");
    branch = mxGetField(result, 0, "branch");

    /* I *think* this is an M-by-2 cell array, but I'm not sure. */

    if (ground) {
        int M, i;
        M = mxGetM(ground);

        for (i = 0; i < M; ++i) {
            mxArray* which_node;
            mxArray* vars;
            int node_num, node_id, N, j;

            mxAssert(mxIsCell(ground), "Not a cell!\n");
            mxAssert(mxGetN(ground) == 2, "Not two columns!\n");

            which_node = mxGetCell(ground, i);
            mxAssert(mxIsNumeric(which_node), "Node not a number!\n");

            vars = mxGetCell(ground, M+i);
            mxAssert(mxIsCell(vars), "Vars not a cell!\n");

            node_num = (int) mxGetScalar(which_node);
            mxAssert(node_num > 0 && node_num <= elt->num_nodes,
                     "Node number out of range");

            node_id = elt->node_ids[node_num-1];
            N = mxGetN(vars);
 
            for (j = 0; j < N; ++j) {
                char name[128];
                mxGetString(mxGetCell(vars, j), name, sizeof(name));
                merge_var(tmp_var_table, node_id, name, 'g');
            }
        }
    }

    if (node) {
        int M = mxGetM(node);
        int i;

        for (i = 0; i < M; ++i) {
            mxArray* which_node;
            mxArray* vars;
            int node_num, node_id, N, j;

            mxAssert(mxIsCell(node), "Not a cell!\n");
            mxAssert(mxGetN(node) == 2, "Not two columns!\n");

            which_node = mxGetCell(node, i);
            mxAssert(mxIsNumeric(which_node), "Node not a number!\n");

            vars = mxGetCell(node, M+i);
            mxAssert(mxIsCell(vars), "Vars not a cell!\n");

            node_num = (int) mxGetScalar(which_node);
            mxAssert(node_num > 0 && node_num <= elt->num_nodes,
                     "Node number out of range");

            node_id = elt->node_ids[node_num-1];
            N = mxGetN(vars);

            for (j = 0; j < N; ++j) {
                char name[128];
                mxGetString(mxGetCell(vars, j), name, sizeof(name));
                var_ids[num_vars++] = 
                        merge_var(tmp_var_table, node_id, name, 'u');
            }
        }
    }

    if (branch) {
        int N = mxGetN(branch);
        int j;

        for (j = 0; j < N; ++j) {
            char name[128];
            mxGetString(mxGetCell(branch, j), name, sizeof(name));
            var_ids[num_vars++] = 
                    merge_var(tmp_var_table, elt->id, name, 'b');
        }
    }

    elt->num_vars = num_vars;
    elt->var_ids = (int*) mempool_cget(code_pool, num_vars * sizeof(int));
    memcpy(elt->var_ids, var_ids, num_vars * sizeof(int));

    mxDestroyArray(result);
    mxDestroyArray(prhs[0]);
}


#endif /* __SUGAR_MEX */

