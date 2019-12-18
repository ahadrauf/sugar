#include <stdio.h>
#include <rpc/xdr.h>

#include "codegen.h"
#include "modelmgr.h"
#include "writegeom.h"

void writebeam(void* xdrsp, double* Q, double L, double W, double H,
               double px, double py, double pz, int* vindex)
{
    XDR* xdrs = (XDR*) xdrsp;
    float ftmp;
    int i, type;

    type = 1;
    xdr_int(xdrs, &type);

    for (i = 0; i < 9; ++i) {
	ftmp = (float) Q[i];
        xdr_float(xdrs, &ftmp);
    }
    
    ftmp = (float) L;
    xdr_float(xdrs, &ftmp);
    ftmp = (float) W;
    xdr_float(xdrs, &ftmp);
    ftmp = (float) H;
    xdr_float(xdrs, &ftmp);
    
    ftmp = (float) px;
    xdr_float(xdrs, &ftmp);
    ftmp = (float) py;
    xdr_float(xdrs, &ftmp);
    ftmp = (float) pz;
    xdr_float(xdrs, &ftmp);

    for (i = 0; i < 12; ++i) 
        xdr_int(xdrs, &vindex[i]);
}

void writegeom(FILE* fp, netlist_ir* netlist, double* q)
{
    XDR xdrs;
    int i, nitems;
    element_ir** elt_table = netlist->elt_table;
    node_ir** node_table = netlist->node_table;
    
    xdrstdio_create(&xdrs, fp, XDR_ENCODE);

    nitems = 0;
    for (i = 0; i < netlist->num_elts; ++i) {
        if (elt_table[i]->modelfun && elt_table[i]->modelfun->writegeom)
	    ++nitems;
    }

    xdr_int(&xdrs, &nitems);

    for (i = 0; i < netlist->num_elts; ++i) {
        if (elt_table[i]->modelfun && elt_table[i]->modelfun->writegeom)
            (elt_table[i]->modelfun->writegeom)(&xdrs, elt_table[i]->data, 
						node_table);
    }

    xdr_int(&xdrs, &netlist->num_dof);

    if (q != NULL) {
        float qtemp;
	for (i = 0; i < netlist->num_dof; ++i) {
	    qtemp = (float) q[i];
	    xdr_float(&xdrs, &qtemp);
	}
    } else {
        float qtemp = 0;
	for (i = 0; i < netlist->num_dof; ++i) 
	    xdr_float(&xdrs, &qtemp);
    }
    
    xdr_destroy(&xdrs);
}

