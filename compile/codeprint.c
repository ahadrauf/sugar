#include <assert.h>
#include <string.h>
#include <math.h>

#include "sugar-lib.h"
#include "codegen.h"
#include "codeprint.h"

void print_nodes(FILE* fp, node_ir** node_table, int num_nodes);
void print_elements(FILE* fp, element_ir** elt_table, int num_elts);
void print_matlab_element(FILE* fp, element_ir* elt);
void print_element_nodes(FILE* fp, element_ir* elt);
void print_element_params(FILE* fp, element_ir* elt);
void print_rotation(FILE* fp, element_ir* elt);


/* Print out the netlist
 */
void print_code(FILE* fp, netlist_ir* netlist)
{
    /* Print header */
    fprintf(fp, "function [net] = netlist(param);\n\n");
    fprintf(fp, "elt_id = 1;\n");
    fprintf(fp, "if nargin == 0\n  param = [];\nend\n\n");

    print_elements(fp, netlist->elt_table, netlist->num_elts);
    
    /* Generate the footer */
    if (netlist->num_elts == 0) {
        fprintf(fp, "  error('Netlist contains no elements!');\n");
    } else {
        print_nodes(fp, netlist->node_table, netlist->num_nodes);
        fprintf(fp, "node_count = %d;\n", netlist->num_nodes);
        fprintf(fp, "net = parse_enrich3(net.elements, nodes, node_count);\n");
    }
}


/* Dump out the nodes and their information
 */
void print_nodes(FILE* fp, node_ir** node_table, int num_nodes)
{
    int i, j;
    for (i = 0; i < num_nodes; ++i) {
        fprintf(fp, "nodes(%d).name = '%s';\n", 
                 node_table[i]->id, node_table[i]->name);
        fprintf(fp, "nodes(%d).elt_ids = [", node_table[i]->id);
        if (node_table[i]->num_elts > 0) {
            for (j = 0; j < node_table[i]->num_elts-1; ++j)
                fprintf(fp, "%d,", node_table[i]->elt_ids[j]+1);
            fprintf(fp, "%d", node_table[i]->elt_ids[j]+1);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "nodes(%d).pos = [%.18e;\n\t %.18e;\n\t %.18e];\n",
                        node_table[i]->id,
                        node_table[i]->pos[0],
                        node_table[i]->pos[1],
                        node_table[i]->pos[2]);
    }
}


/* Print elements
 */
void print_elements(FILE* fp, element_ir** elt_table, int num_elts)
{
    int i;
    for (i = 0; i < num_elts; ++i)
        print_matlab_element(fp, elt_table[i]);
}


/* Print global node index assignments.
 */
void print_element_nodes(FILE* fp, element_ir* elt)
{
    int node_index = 0;
    fprintf(fp, "elt.node_ids = [");
    for (node_index = 0; node_index < elt->num_nodes; ++node_index) {
        fprintf(fp, "%d ", elt->node_ids[node_index]+1);
    }
    fprintf(fp, "];\n");
}


/* Print variable index assignments
 */
void print_element_vars(FILE* fp, element_ir* elt)
{
    int var_index = 0;
    fprintf(fp, "elt.var_ids = [");
    for (var_index = 0; var_index < elt->num_vars; ++var_index) {
        fprintf(fp, "%d ", elt->var_ids[var_index]+1);
    }
    fprintf(fp, "];\n");
}


/* Print a normal element
 */
void print_matlab_element(FILE* fp, element_ir* elt)
{
    fprintf(fp, "elt = [];\n");
    fprintf(fp, "elt.name = '%s';\n", elt->name);
    fprintf(fp, "elt.model = 'MF_%s';\n", elt->model);
    print_element_nodes(fp, elt);
    print_element_vars(fp, elt);
    print_element_params(fp, elt);

    fprintf(fp, "elt.R = ");
    print_rotation(fp, elt);

    fprintf(fp, "net.elements(%d) = elt;\n\n", elt->id);
}


/* Print out parameter list 
 */
void print_element_params(FILE* fp, element_ir* elt)
{
    parameter_ir* params = elt->params;
    fprintf(fp, "elt.parameter = [];\n");
    while (params != NULL) {
        val_t val = params->value;
        fprintf(fp, "elt.parameter.%s = ", params->name);
        if (val.type == 'd')
            fprintf(fp, "%.18e", val.val.d);
        else
            fprintf(fp, "'%s'", val.val.s);
        fprintf(fp, ";\n");
        params = params->next;
    }
}


/* Print the (transposed) rotation to get local->global 
 */
void print_rotation(FILE* fp, element_ir* elt)
{
    double* Q = elt->Q;

    fprintf(fp, "...\n"
                "    [%.18e,%.18e,%.18e;\n"
                "     %.18e,%.18e,%.18e;\n"
                "     %.18e,%.18e,%.18e];\n",
            Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7], Q[8]);
}

