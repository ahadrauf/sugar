#ifndef __CODEPRINT_H
#define __CODEPRINT_H

/* Create a Matlab file which will generate the data structures for
 * a netlist instance.  This module should be phased out soon, since
 * the codemx module creates the data structures directly.
 *
 * print_code -- prints the routine for a netlist to file
 */

#include <stdio.h>
#include "codegen.h"

void print_code(FILE* fp, netlist_ir* netlist);

#endif /* __CODEPRINT_H */
