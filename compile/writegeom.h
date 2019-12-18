#ifndef __WRITEGEOM_H
#define __WRITEGEOM_H

void writebeam(void* xdrsp, double* Q, double L, double W, double H,
               double px, double py, double pz, int* vindex);
void writegeom(FILE* fp, netlist_ir* net, double* q);

#endif /* __WRITEGEOM_H */
