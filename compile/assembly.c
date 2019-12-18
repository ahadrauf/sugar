#include "assembly.h"

void assemble_globalm(element_ir* elt, double* Mlocal, double* M, int ldM)
{
    int i, j;
    int nvars = elt->num_vars;
    int* vars = elt->var_ids;
    
    for (j = 0; j < nvars; ++j)
        for (i = 0; i < nvars; ++i)
	    M[vars[j]*ldM + vars[i]] += Mlocal[j*nvars + i];
}

void assemble_globalv(element_ir* elt, double* vlocal, double* v)
{
    int i;
    int nvars = elt->num_vars;
    int* vars = elt->var_ids;
    
    for (i = 0; i < nvars; ++i)
        v[vars[i]] += vlocal[i];
}

void assemble_localv(element_ir* elt, double* vlocal, double* v)
{
    int i;
    int nvars = elt->num_vars;
    int* vars = elt->var_ids;
    
    for (i = 0; i < nvars; ++i)
        vlocal[i] = (v ? v[vars[i]] : 0);
}

#define rotv2locals(Q, n, v, stride) \
{ \
    double vloc[3]; \
    vloc[0] = v[0]; \
    vloc[1] = v[stride]; \
    vloc[2] = ((n == 2) ? 0 : v[2*stride]); \
\
    v[0*stride] = Q[0]*vloc[0] + Q[3]*vloc[1] + Q[6]*vloc[2]; \
    v[1*stride] = Q[1]*vloc[0] + Q[4]*vloc[1] + Q[7]*vloc[2]; \
    if (n != 2) \
        v[2*stride] = Q[2]*vloc[0] + Q[5]*vloc[1] + Q[8]*vloc[2]; \
}

#define rotv2globals(Q, n, v, stride) \
{ \
    double vloc[3]; \
    vloc[0] = v[0]; \
    vloc[1] = v[stride]; \
    vloc[2] = ((n == 2) ? 0 : v[2*stride]); \
\
    v[0*stride] = Q[0]*vloc[0] + Q[1]*vloc[1] + Q[2]*vloc[2]; \
    v[1*stride] = Q[3]*vloc[0] + Q[4]*vloc[1] + Q[5]*vloc[2]; \
    if (n != 2) \
        v[2*stride] = Q[6]*vloc[0] + Q[7]*vloc[1] + Q[8]*vloc[2]; \
}

void rotv2local(double* Q, int n, double* v)
{
    rotv2locals(Q, n, v, 1)
}

void rotv2global(double* Q, int n, double* v)
{
    rotv2globals(Q, n, v, 1)
}

void rotm2global(double* Q, int n, double* M, int ldM)
{
    /* Apply Q' */
    rotv2globals(Q, n, (&(M[0*ldM])), 1)
    rotv2globals(Q, n, (&(M[1*ldM])), 1)
    rotv2globals(Q, n, (&(M[2*ldM])), 1)

    /* Apply Q */
    rotv2globals(Q, n, (&(M[0])), ldM)
    rotv2globals(Q, n, (&(M[1])), ldM)
    rotv2globals(Q, n, (&(M[2])), ldM)
}

