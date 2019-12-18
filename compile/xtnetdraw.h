#ifndef __XFNETDRAW_H
#define __XFNETDRAW_H

#include "codegen.h"

void display_device(int argc, char** argv, netlist_ir* net, double* dq);
void draw_beam_old(double x1, double y1, double x2, double y2, double w);
void draw_beam(double* R, double x, double y, double z,
	       double l, double w, double h,
	       double dx1, double dy1, double dz1,
	       double drx1, double dry1, double drz1,
	       double dx2, double dy2, double dz2,
	       double drx2, double dry2, double drz2);

#endif /* __XFNETDRAW_H */
