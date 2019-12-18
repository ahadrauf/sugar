#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "evalexpr.h"
#include "shash.h"
#include "codestack.h"
#include "codegen.h"
#include "varassign.h"
#include "modelmgr.h"

#ifdef __SUGAR_MEX
  #include "codemx.h"
#endif

static mempool_t code_pool;
static mempool_t temp_pool;
static int var_count;
static int ground_count;
static int num_nodes;

/* Assign variable indices for the entries in a netlist
 */
void assign_var_indices(mempool_t pool, netlist_ir* netlist)
{
    int i;
    var_ir** tmp_var_table;
    var_ir** var_table;
    int* mapping;
    int current_unground, current_ground;

    code_pool = pool;
    temp_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);

    tmp_var_table = (var_ir**) mempool_cget(temp_pool, sizeof(var_ir*) *
                            (netlist->num_nodes + netlist->num_elts));

    /* -- Get the contributions from each element to the variable lists */
    
    var_count = 0;
    ground_count = 0;
    num_nodes = netlist->num_nodes;

    for (i = 0; i < netlist->num_elts; ++i) {
        element_ir* elt = netlist->elt_table[i];
        if (elt->modelfun && elt->modelfun->add_vars)
            (elt->modelfun->add_vars)(elt->data, tmp_var_table, code_pool);
    }

    /* -- Build the variable table, reordering with grounds at end */
   
    mapping = (int*) mempool_cget(temp_pool, var_count * sizeof(int));
    var_table = (var_ir**) mempool_cget(code_pool, var_count * 
                                        sizeof(var_ir*));

    current_unground = 0;
    current_ground = var_count - ground_count;

    for (i = 0; i < netlist->num_nodes + netlist->num_elts; ++i) {
        var_ir* var = tmp_var_table[i];
        while (var != NULL) {
            if (var->type == 'g') {
                mapping[var->index] = current_ground;
                var_table[current_ground] = var;
                var->index = current_ground++;
            } else {
                mapping[var->index] = current_unground;
                var_table[current_unground] = var;
                var->index = current_unground++;
            }
            var = var->next;
        }
    }
    
    /* -- Map from old variable ordering to new one */

    for (i = 0; i < netlist->num_elts; ++i) {
        element_ir* elt = netlist->elt_table[i];
        int j;
        for (j = 0; j < elt->num_vars; ++j) {
            elt->var_ids[j] = mapping[elt->var_ids[j]];
        }
    }

    /* -- Save and clean up */

    netlist->var_table = var_table;
    netlist->num_vars = var_count;
    netlist->num_dof = var_count - ground_count;

    mempool_destroy(temp_pool);
}

/* Merge a variable into the currently accumulating list
 */
int merge_var(var_ir** var_table, int owner_id, char* name, char type)
{
    var_ir* new_var;
    var_ir** var_list;

    if (type == 'b')
        var_list = &(var_table[num_nodes + owner_id]);
    else
        var_list = &(var_table[owner_id]);
    
    /* If the variable exists, update the status (if necessary) and
     * return a pointer to the structure.
     */
    while (*var_list != NULL) {
        if (strcmp(name, (*var_list)->name) == 0) {
            if (type == 'g' && (*var_list)->type != 'g') {
                (*var_list)->type = 'g';
                ++ground_count;
            }
            return (*var_list)->index;
        }
        var_list = &((*var_list)->next);
    }

    /* Otherwise, create a new variable entry and add it in 
     */
    new_var = (var_ir*)
        mempool_cget(code_pool, sizeof(var_ir));
    new_var->name = mempool_strdup(code_pool, name);
    new_var->type = type;
    new_var->owner_id = owner_id;
    new_var->index = var_count++;
    new_var->next = NULL;
    *var_list = new_var;
    if (type == 'g')
        ++ground_count;

    return new_var->index;
}

