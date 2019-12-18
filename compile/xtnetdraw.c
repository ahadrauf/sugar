#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "codegen.h"
#include "modelmgr.h"
#include "assembly.h"

#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <X11/Core.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/Command.h>

#define FRAME_SPACE 10
#define FRAME_IN_SPACE 5
#define FRAMEW 510
#define FRAMEH 510
#define BUTTON_SPACE 20

#define FRAME_INW (FRAMEW - 2*FRAME_IN_SPACE)
#define FRAME_INH (FRAMEH - 2*FRAME_IN_SPACE)

#define min(a,b)  ((a < b) ? a : b)
#define max(a,b)  ((a > b) ? a : b)

static double Q[9];

static double xmin, ymin, scl;
static netlist_ir* netlist;
static double* dq;

static Display* disp;
static Drawable win;
static GC gc;

static int spin_on;
static Widget drawing;

#define SPIN_DELAY 150L

static void hermite_cubic(double* fout, int len, double L, 
		          double f0, double f00, double fL, double fLL)
{
    double f0L = (fL - f0) / L;
    double f00L = (f0L - f00) / L;
    double f0LL = (fLL - f0L) / L;
    double f00LL = (f0LL - f00L) / L;
    int i;
    for (i = 0; i < len; ++i) {
        double s = i*L / (len-1);
	fout[i] = f0 + s*(f00 + s*(f00L + (s-L)*f00LL));
    }
}

void draw_beam(double* R, double x, double y, double z,
	       double l, double w, double h,
	       double dx1, double dy1, double dz1,
	       double drx1, double dry1, double drz1,
	       double dx2, double dy2, double dz2,
	       double drx2, double dry2, double drz2)
{
    #define NPTS 10
	
    XPoint xp[4*NPTS];
    double pts[12*NPTS];
    double* p;
    XPoint* q;
    int i;

    double yout[NPTS];
    double zout[NPTS];

    hermite_cubic(yout, NPTS, l, dy1, drz1, dy2, drz2);
    hermite_cubic(zout, NPTS, l, dz1,-dry1, dz2,-dry2);

    /* Form points in reference frame */
    p = pts;
    for (i = 0; i < NPTS; ++i) {
	double s = ((double) i)/(NPTS-1);
        p[0] = p[3] = p[6] = p[9] = s*l + s*dx1 + (1-s)*dx2;
	
	p[1 ] =  w/2 + yout[i];  p[2 ] =  h/2 + zout[i];
	p[4 ] = -w/2 + yout[i];  p[5 ] =  h/2 + zout[i];
	p[7 ] =  w/2 + yout[i];  p[8 ] = -h/2 + zout[i];
	p[10] = -w/2 + yout[i];  p[11] = -h/2 + zout[i];
	p += 12;
    }
    
    /* Rotate to global, then to view; then scale and translate */
    p = pts;
    q = xp;
    for (i = 0; i < 4*NPTS; ++i) {
	double v[3];
	
	v[0] = x + R[0]*p[0] + R[1]*p[1] + R[2]*p[2];
	v[1] = y + R[3]*p[0] + R[4]*p[1] + R[5]*p[2];
	v[2] = z + R[6]*p[0] + R[7]*p[1] + R[8]*p[2];
	
	p[0] = Q[0]*v[0] + Q[3]*v[1] + Q[6]*v[2];
	p[1] = Q[1]*v[0] + Q[4]*v[1] + Q[7]*v[2];
	p[2] = Q[2]*v[0] + Q[5]*v[1] + Q[8]*v[2];
	
        p[0] = scl*(p[0]-xmin); p[1] = scl*(p[1]-ymin);
        q->x = (int)(FRAME_INW * p[0]) + FRAME_IN_SPACE;
        q->y = (int)(FRAME_INH * (1-p[1])) + FRAME_IN_SPACE;

	p += 3;
	q++;
    }
    
    /* Lucy, you got some 'splayin' to do! */
    q = xp;
    for (i = 0; i < NPTS-1; ++i) {
	    
        XDrawLine(disp, win, gc, q[0].x, q[0].y, q[1].x, q[1].y);
        XDrawLine(disp, win, gc, q[1].x, q[1].y, q[2].x, q[2].y);
        XDrawLine(disp, win, gc, q[2].x, q[2].y, q[3].x, q[3].y);
        XDrawLine(disp, win, gc, q[3].x, q[3].y, q[0].x, q[0].y);
	
        XDrawLine(disp, win, gc, q[0].x, q[0].y, q[4].x, q[4].y);
        XDrawLine(disp, win, gc, q[1].x, q[1].y, q[5].x, q[5].y);
        XDrawLine(disp, win, gc, q[2].x, q[2].y, q[6].x, q[6].y);
        XDrawLine(disp, win, gc, q[3].x, q[3].y, q[7].x, q[7].y);

	q += 4;
    }

    XDrawLine(disp, win, gc, q[0].x, q[0].y, q[1].x, q[1].y);
    XDrawLine(disp, win, gc, q[1].x, q[1].y, q[2].x, q[2].y);
    XDrawLine(disp, win, gc, q[2].x, q[2].y, q[3].x, q[3].y);
    XDrawLine(disp, win, gc, q[3].x, q[3].y, q[0].x, q[0].y);
}


static void rotX(double radians)
{
    double Qnew[9];
    double c = cos(radians);
    double s = sin(radians);

    memcpy(Qnew, Q, 9*sizeof(double));
    Qnew[1] =  c*Q[1] + s*Q[2];
    Qnew[2] = -s*Q[1] + c*Q[2];
    Qnew[4] =  c*Q[4] + s*Q[5];
    Qnew[5] = -s*Q[4] + c*Q[5];
    Qnew[7] =  c*Q[7] + s*Q[8];
    Qnew[8] = -s*Q[7] + c*Q[8];
    memcpy(Q, Qnew, 9*sizeof(double));
}


static void rotY(double radians)
{
    double Qnew[9];
    double c = cos(radians);
    double s = sin(radians);

    memcpy(Qnew, Q, 9*sizeof(double));
    Qnew[0] =  c*Q[0] + s*Q[2];
    Qnew[2] = -s*Q[0] + c*Q[2];
    Qnew[3] =  c*Q[3] + s*Q[5];
    Qnew[5] = -s*Q[3] + c*Q[5];
    Qnew[6] =  c*Q[6] + s*Q[8];
    Qnew[8] = -s*Q[6] + c*Q[8];
    memcpy(Q, Qnew, 9*sizeof(double));
}


static void rotZ(double radians)
{
    double Qnew[9];
    double c = cos(radians);
    double s = sin(radians);

    memcpy(Qnew, Q, 9*sizeof(double));
    Qnew[0] =  c*Q[0] + s*Q[1];
    Qnew[1] = -s*Q[0] + c*Q[1];
    Qnew[3] =  c*Q[3] + s*Q[4];
    Qnew[4] = -s*Q[3] + c*Q[4];
    Qnew[6] =  c*Q[6] + s*Q[7];
    Qnew[7] = -s*Q[6] + c*Q[7];
    memcpy(Q, Qnew, 9*sizeof(double));
}


static void draw_netlist(Widget w)
{
    int i;
    element_ir** elt_table = netlist->elt_table;
    node_ir** node_table = netlist->node_table;

    disp = XtDisplay(w);
    win = XtWindow(w);

    gc = XCreateGC(disp, win, 0, NULL);
    XSetForeground(disp, gc, 1);
    XSetBackground(disp, gc, 0);

    for (i = 0; i < netlist->num_elts; ++i) {
        if (elt_table[i]->modelfun && elt_table[i]->modelfun->display) {
	    static double local[256]; /* TODO -- allocate dynamically */
	    assemble_localv(elt_table[i], local, dq);
            (elt_table[i]->modelfun->display)(elt_table[i]->data, 
					      node_table, local);
	}
    }

    XFreeGC(disp, gc);
}

static void get_scaling()
{
    int i;
    double xmax, ymax;
    double xrng, yrng;

    xmin = ymin = 0;
    xmax = ymax = 0;
    
    for (i = 0; i < netlist->num_nodes; ++i) {
        xmin = min(xmin, netlist->node_table[i]->pos[0]);
        ymin = min(ymin, netlist->node_table[i]->pos[1]);
        xmax = max(xmax, netlist->node_table[i]->pos[0]);
        ymax = max(ymax, netlist->node_table[i]->pos[1]);
    }

    xrng = xmax-xmin;
    yrng = ymax-ymin;
    
    if (xrng < yrng) {
        scl = 1/yrng;
        xmin -= (yrng-xrng)/2;
    } else {
        scl = 1/xrng;
        ymin -= (xrng-yrng)/2;
    }
}

static void redisplay_event(Widget w, XtPointer client, XExposeEvent* ev,
                            Boolean* continue_dispatch)
{
    if (ev->count != 0)
        return;

    XClearWindow(XtDisplay(w), XtWindow(w));
    draw_netlist(w);
}

static void keypress_event(Widget w, XtPointer client, XKeyEvent* ev,
                           Boolean* continue_dispatch)
{
    double inc = (M_PI / 18);
    switch (ev->keycode) {
	case 98:    /* Up    */
	    rotX(-inc);
	    break;

	case 104:   /* Down  */
	    rotX(inc);
	    break;
	    
	case 100:   /* Left  */
	    rotY(inc);
	    break;
	    
	case 102:   /* Right */
	    rotY(-inc);
	    break;

	case 99:    /* PgUp  */
	    rotZ(-inc);
	    break;

	case 105:   /* PgDn  */
	    rotZ(inc);
	    break;
    }

    XClearWindow(XtDisplay(w), XtWindow(w));
    draw_netlist(w);
}

static void timer_func(XtPointer closure, XtIntervalId* id)
{
    if (spin_on) {
        double inc = (M_PI / 18);
        rotY(inc);
        XClearWindow(XtDisplay(drawing), XtWindow(drawing));
        draw_netlist(drawing);
    }

    XtAddTimeOut(SPIN_DELAY, timer_func, closure);
}

static void spin_func(Widget w, XtPointer client, XtPointer call)
{
    spin_on = !spin_on;
}

static void quit_func(Widget w, XtPointer client, XtPointer call)
{
    exit(0);
}

void display_device(int argc, char** argv, netlist_ir* net, double* dqarg)
{ 
    Widget toplevel;
    Widget box;
    Widget quit;
    Widget spin;

    int n;
    Arg wargs[10];

    netlist = net;
    dq = dqarg;

    get_scaling();
    memset(Q, 0, 9*sizeof(double));
    Q[0] = Q[4] = Q[8] = 1;

    toplevel = XtInitialize(argv[0],"drawing", NULL, 0,
                            &argc, argv);

    box = XtCreateManagedWidget("box", boxWidgetClass,
                                toplevel, NULL, 0);

    drawing = XtCreateManagedWidget("drawing",coreWidgetClass,
                                    box, NULL, 0);

    quit = XtCreateManagedWidget("quit", commandWidgetClass, box, NULL, 0);
    XtAddCallback(quit, XtNcallback, quit_func, NULL);

    spin_on = 0;
    spin = XtCreateManagedWidget("spin", commandWidgetClass, box, NULL, 0);
    XtAddCallback(spin, XtNcallback, spin_func, NULL);

    n = 0;
    XtSetArg(wargs[n], XtNheight, FRAMEH); n++;
    XtSetArg(wargs[n], XtNwidth, FRAMEW); n++;
    XtSetValues(drawing, wargs, n);

    XtAddEventHandler(drawing, ExposureMask, FALSE,
                      (XtEventHandler) redisplay_event, NULL);
    XtAddEventHandler(drawing, KeyPressMask, FALSE,
		      (XtEventHandler) keypress_event, NULL);

    XtAddTimeOut(SPIN_DELAY, timer_func, NULL);

    XtRealizeWidget(toplevel);

    XtMainLoop();

}

