#ifndef __VARASSIGN_H
#define __VARASSIGN_H

struct var_ir;

void assign_var_indices(mempool_t pool, netlist_ir* netlist);
int merge_var(var_ir** var_table, int owner_id, char* name, char type);

#endif /* __VARASSIGN_H */
