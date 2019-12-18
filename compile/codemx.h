#ifndef __CODEMX_H
#define __CODEMX_H

#include "mex.h"
#include "matrix.h"
#include "codegen.h"

mxArray* params_to_mx(parameter_ir* params);
mxArray* element_to_mx(element_ir** elt_table, int num_elts, int dof);
mxArray* node_to_mx(node_ir** node_table, int num_nodes);
mxArray* netlist_to_mx(netlist_ir* netlist);

#endif /* __CODE_MX_H */
