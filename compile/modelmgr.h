#ifndef __MODELMGR_H
#define __MODELMGR_H

#include "mempool.h"
#include "codegen.h"

typedef struct element_pos {
    int num_mech_nodes;
    double* relpos;
} element_pos_t;

typedef struct model_t {
    void* (*init_data)(mempool_t data_pool, element_ir* elt);
    void (*find_position)(void* data, mempool_t pool, element_pos_t* elt_pos);
    void (*post_position)(void* data, node_ir** node_table);
    void (*display)(void* data, node_ir** node_table, double* dqlocal);
    void (*writegeom)(void* xdrsp, void* data, node_ir** node_table);
    void (*add_vars)(void* data, var_ir** tmp_var_table, mempool_t code_pool);
    void (*add_K)(void* data, double* K);
    void (*add_F)(void* data, double* F);
} model_t;

void model_manager_init();
void model_register(char* name, model_t* model);
model_t* model_lookup(char* name);
void model_manager_shutdown();

#endif /* __MODELMGR_H */
