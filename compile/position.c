#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "evalexpr.h"
#include "shash.h"
#include "codestack.h"
#include "codegen.h"

#ifdef __SUGAR_MEX
  #include "codemx.h"
#endif

#include "modelmgr.h"

static mempool_t local_pool;

static int find_unpositioned(char* colors, netlist_ir* netlist);
static int position_using_element(element_ir* elt, element_pos_t* elt_pos,
                                  char* colors, node_ir** node_table);
static void initialize_positioning(char* colors, element_pos_t* elt_pos,
                                   netlist_ir* netlist);

static void get_elt_relpos(element_ir* elt, element_pos_t* elt_pos);


void position_nodes(netlist_ir* netlist)
{
    node_ir** node_table = netlist->node_table;
    element_ir** elt_table = netlist->elt_table;
    int num_nodes = netlist->num_nodes;
    int num_elts = netlist->num_elts;
   
    element_pos_t* elt_pos; 
    char* colors;
    int* stack;
    int stack_top;
    int current_node;
    int i;

    local_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);

    elt_pos = (element_pos_t*) 
            mempool_cget(local_pool, num_elts * sizeof(element_pos_t));
    colors = (char*) mempool_cget(local_pool, num_nodes * sizeof(char));
    initialize_positioning(colors, elt_pos, netlist);
   
    stack = (int*) mempool_get(local_pool, num_nodes * sizeof(int));
    stack_top = 0;

    current_node = find_unpositioned(colors, netlist);
    while (current_node != -1) {

        /* Position an otherwise unpositioned node at the origin
         * and push it onto the stack
         */
        node_table[current_node]->pos[0] = 0;
        node_table[current_node]->pos[1] = 0;
        node_table[current_node]->pos[2] = 0;
        colors[current_node] = 'g';
        stack[stack_top++] = current_node;
       
        /* Do a tree search, positioning nodes along the way
         */
        while (stack_top != 0) {
            int j, k;
                
            current_node = stack[--stack_top];
            colors[current_node] = 'b';

            for (j = 0; j < node_table[current_node]->num_elts; ++j) {
                int current_element = node_table[current_node]->elt_ids[j];
                element_ir* elt = elt_table[current_element];
                if (position_using_element(elt, &elt_pos[current_element],
                                           colors, node_table)) {
                    for (k = 0; k < elt->num_nodes; ++k) {
                        if (colors[elt->node_ids[k]] == 'w') {
                            stack[stack_top++] = elt->node_ids[k];
                            colors[elt->node_ids[k]] = 'g';
                        }
                    }
                }
            }
        }
       
        /* If there are left-over nodes, warn the user that there is
         * probably a mistake and give them positions.
         */
        current_node = find_unpositioned(colors, netlist);
        if (current_node >= 0)
            printf("Warning: Not all mechanical nodes positioned!\n");
    }

#ifndef __SUGAR_MEX
    /* TODO: Incorporate even when in Matlab */

    /* Now do a post-positioning sweep */
    for (i = 0; i < num_elts; ++i) {
        element_ir* elt = elt_table[i];
	if (elt->modelfun && elt->modelfun->post_position)
	    (elt->modelfun->post_position)(elt->data, netlist->node_table);
    }
#endif

    mempool_destroy(local_pool);
}


static int find_unpositioned(char* colors, netlist_ir* netlist)
{
    int i;
    for (i = 0; i < netlist->num_nodes; ++i)
        if (colors[i] == 'w')
            return i;

    return -1;
}


static int position_using_element(element_ir* elt, element_pos_t* elt_pos,
                                  char* colors, node_ir** node_table)
{
    int marked_index = -1;
    int all_marked = 1;
    int i;

    for (i = 0; i < elt_pos->num_mech_nodes; ++i) {
        all_marked = (all_marked && colors[elt->node_ids[i]] != 'w');
        if (colors[elt->node_ids[i]] == 'b')
            marked_index = i;
    }

    if (!all_marked) {

        double* rp = elt_pos->relpos;
        double shift[3];

        shift[0] = node_table[elt->node_ids[marked_index]]->pos[0];
        shift[1] = node_table[elt->node_ids[marked_index]]->pos[1];
        shift[2] = node_table[elt->node_ids[marked_index]]->pos[2];
            
        shift[0] -= rp[marked_index*3 + 0];
        shift[1] -= rp[marked_index*3 + 1];
        shift[2] -= rp[marked_index*3 + 2];

        for (i = 0; i < elt_pos->num_mech_nodes; ++i) {
            if (colors[elt->node_ids[i]] == 'w') {
                node_table[elt->node_ids[i]]->pos[0] = rp[i*3+0] + shift[0];
                node_table[elt->node_ids[i]]->pos[1] = rp[i*3+1] + shift[1];
                node_table[elt->node_ids[i]]->pos[2] = rp[i*3+2] + shift[2];
            }
        }
        return 1;
    }
    return 0;
}


static void initialize_positioning(char* colors, element_pos_t* elt_pos,
                                   netlist_ir* netlist)
{
    int i, j;
    for (i = 0; i < netlist->num_elts; ++i) {
        get_elt_relpos(netlist->elt_table[i], &(elt_pos[i]));
        for (j = 0; j < elt_pos[i].num_mech_nodes; ++j) {
            int node_id = netlist->elt_table[i]->node_ids[j];
            colors[node_id] = 'w';
        }
    }
}

static void get_elt_relpos(element_ir* elt, element_pos_t* elt_pos)
{
    if (elt->modelfun && elt->modelfun->find_position) {
        (elt->modelfun->find_position)(elt->data, local_pool, elt_pos);
    } else {
        elt_pos->num_mech_nodes = 0;
    }
}

void post_positioning(netlist_ir* netlist)
{
    int i;
    for (i = 0; i < netlist->num_elts; ++i) {
        element_ir* elt = netlist->elt_table[i];
	if (elt->modelfun && elt->modelfun->post_position) {
	    (elt->modelfun->post_position)(elt->data, netlist->node_table);
	}
    }
}
