#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "assembly.h"
#include "codegen.h"
#include "modelmgr.h"

int lusolve(int n, int lda, double* A, double* b)
{
    extern void dgesv_(int* n, int* nrhs, double* A, int* lda,
		       int* ipiv, double* B, int* ldb, int* info);

    int nrhs = 1;
    int* ipiv;
    int info;

    ipiv = (int*) malloc(n * sizeof(int));
    dgesv_(&n, &nrhs, A, &lda, ipiv, b, &n, &info);
    free(ipiv);
    
    return info;
}

double* dcsolve(mempool_t pool, netlist_ir* net)
{
    double* x = (double*)
	    mempool_cget(pool, net->num_vars * sizeof(double));

    double* K = (double*)
	    mempool_cgeth(pool, net->num_vars * net->num_vars * sizeof(double));
    double* F = (double*)
	    mempool_cget(pool, net->num_vars * sizeof(double));
    double* local = (double*)
	    mempool_cgeth(pool, net->num_vars * net->num_vars * sizeof(double));

    int i,j;

    for (i = 0; i < net->num_elts; ++i) {
        element_ir* elt = net->elt_table[i];
	if (elt->modelfun) {
	    if (elt->modelfun->add_K) {
	        (elt->modelfun->add_K)(elt->data, local);
		assemble_globalm(elt, local, K, net->num_vars);
	    }
	    if (elt->modelfun->add_F) {
		(elt->modelfun->add_F)(elt->data, local);
		assemble_globalv(elt, local, F);
	    }
	}
    }

    memcpy(x, F, net->num_vars * sizeof(double));
    i = lusolve(net->num_dof, net->num_vars, K, x);
    assert(i == 0);

    for (i = 0; i < net->num_dof; ++i) {
        double residual = -F[i];
	for (j = 0; j < net->num_dof; ++j)
	    residual += K[j*net->num_vars + i] * x[j];
	printf("%g\t%g\t%g\n", x[i], F[i], residual);
    }
    
    mempool_freeh(K);

    return x;
}
