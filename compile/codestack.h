#ifndef __CODESTACK_H
#define __CODESTACK_H

/*
 * Manage run-time (?) stack for data structure generation
 *
 * frame_push -- Create a stack frame for a new subnet instance
 * frame_pop  -- Pop the current stack frame
 *
 * frame_node_index -- Get the index of a node (or create a new one)
 *                     Node names are entirely local, except so far as they
 *                     can alias external nodes.
 * get_node_count   -- Get current number of nodes
 * frame_var_entry  -- Get a pointer to a variable in the value table
 * frame_rot2local  -- Generate a rotation to local coordinates for an
 *                     element, starting with the local coordinates for the
 *                     subnet.
 */

#include "parse.h"
#include "evalexpr.h"
#include "codegen.h"

void    frame_stack_create (mempool_t code_pool, int num_global_vars);
void    frame_stack_destroy(void);

void    frame_push      (code_t* element, parameter_ir* process_params);
void    frame_pop       (void);
void    frame_trace     ();

char*   frame_scoped_name(mempool_t pool, char* name);
int     frame_node_index(char* name);
int     get_node_count();
val_t*  frame_var_entry (int var_id, int scope);
void    frame_rot2local (double* Qbuf, double ox, double oy, double oz);

#endif /* __CODESTACK_H */
