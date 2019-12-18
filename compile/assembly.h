#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "codegen.h"

void assemble_globalm(element_ir* elt, double* Mlocal, double* M, int ldM);
void assemble_globalv(element_ir* elt, double* vlocal, double* v);
void assemble_localv(element_ir* elt, double* vlocal, double* v);

void rotv2local(double* Q, int n, double* v);
void rotv2global(double* Q, int n, double* v);
void rotm2global(double* Q, int n, double* M, int ldM);

#endif /* __ASSEMBLY_H */

