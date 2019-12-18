
#include <string.h>
#include <assert.h>
#include <math.h>
#include "modelmgr.h"
#include "models.h"
#include "varassign.h"

#include "assembly.h"
#include "writegeom.h"
#ifdef USE_XT
#include "xtnetdraw.h"
#endif /* USE_XT */

/* ---------- */

static void test_symm(double* K, int N)
{
    int i, j;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < i; ++j) {
            if (fabs(K[i*N+j] - K[j*N+i]) >= 2e-16)
                printf("Asymmetry at (%d, %d) [%g != %g]\n", i, j,
                       K[i*N+j], K[j*N+i]);
        }
    }
}


static int get_double_param(double* val, char* name, parameter_ir* param)
{
    while (param != NULL) {
        if (strcmp(param->name, name) == 0) {
	    if (param->value.type == 'd') {
	        *val = param->value.val.d;
	        return 0;
	    } else
		return -1;
	}
	param = param->next;
    }
    return -2;
}


/* ----- Anchor model ----- */

static void* anchor_init_data(mempool_t data_pool, element_ir* elt);
static void anchor_add_vars(void* data, var_ir** tmp_var_table, mempool_t pool);
static void anchor_display(void* data, node_ir** node_table, double* dq);

struct anchor_params {
    double l, w, h;
    int node_ids[1];
    int* var_ids; 
    double* Q;
};

void anchor_register()
{
    model_t anchor = {
	anchor_init_data,
	NULL,
	NULL,
	anchor_display,
	NULL,
	anchor_add_vars,
	NULL,
	NULL
    };

    model_register("anchor", &anchor);
}

static void* anchor_init_data(mempool_t data_pool, element_ir* elt)
{
    struct anchor_params* aparams = (struct anchor_params*)
            mempool_cget(data_pool, sizeof(struct anchor_params));

    get_double_param(&(aparams->l), "l", elt->params);
    get_double_param(&(aparams->w), "w", elt->params);
    get_double_param(&(aparams->h), "h", elt->params);

    aparams->Q = elt->Q;
    aparams->node_ids[0] = elt->node_ids[0];

    elt->var_ids = mempool_cget(data_pool, 6*sizeof(int));
    elt->num_vars = 6;
    aparams->var_ids = elt->var_ids;

    return aparams;
}

static void anchor_add_vars(void* data, var_ir** var_table, mempool_t pool)
{
    struct anchor_params* aparams = (struct anchor_params*) data;
    int n1 = aparams->node_ids[0];

    aparams->var_ids[0] = merge_var(var_table, n1, "x", 'g');
    aparams->var_ids[1] = merge_var(var_table, n1, "y", 'g');
    aparams->var_ids[2] = merge_var(var_table, n1, "z", 'g');
    aparams->var_ids[3] = merge_var(var_table, n1, "rx", 'g');
    aparams->var_ids[4] = merge_var(var_table, n1, "ry", 'g');
    aparams->var_ids[5] = merge_var(var_table, n1, "rz", 'g');
}

static void anchor_display(void* data, node_ir** node_table, double* dq)
{
#ifdef USE_XT
    struct anchor_params* aparams = (struct anchor_params*) data;
    node_ir* n1 = node_table[aparams->node_ids[0]];

    draw_beam(aparams->Q, n1->pos[0], n1->pos[1], n1->pos[2],
              aparams->l, aparams->w, aparams->h,
	      0, 0, 0,  0, 0, 0,
	      0, 0, 0,  0, 0, 0);
#endif /* USE_XT*/
}


/* ----- Force model ----- */

static void* force_init_data(mempool_t data_pool, element_ir* elt);
static void force_add_vars(void* data, var_ir** tmp_var_table, mempool_t pool);
static void force_F(void* data, double* F);

struct force_params {
    double F;
    int node_ids[1];
    int* var_ids;
    double* Q;
};

void force_register()
{
    model_t force = {
        force_init_data,
        NULL,
        NULL,
        NULL,
	NULL,
        force_add_vars,
        NULL,
        force_F
    };

    model_register("f2d", &force);
    model_register("f3d", &force);
}

static void* force_init_data(mempool_t data_pool, element_ir* elt)
{
    struct force_params* fparams = (struct force_params*)
            mempool_cget(data_pool, sizeof(struct force_params));

    get_double_param(&(fparams->F), "F", elt->params);

    fparams->Q = elt->Q;
    fparams->node_ids[0] = elt->node_ids[0];

    elt->var_ids = mempool_cget(data_pool, 3*sizeof(int));
    elt->num_vars = 3;
    fparams->var_ids = elt->var_ids;

    return fparams;
}


static void force_add_vars(void* data, var_ir** var_table, mempool_t pool)
{
    struct force_params* fparams = (struct force_params*) data;
    
    int n1 = fparams->node_ids[0];

    fparams->var_ids[0] = merge_var(var_table, n1, "x", 'u');
    fparams->var_ids[1] = merge_var(var_table, n1, "y", 'u');
    fparams->var_ids[2] = merge_var(var_table, n1, "z", 'u');
}


static void force_F(void* data, double* F)
{
    struct force_params* fparams = (struct force_params*) data;

    F[0] = fparams->F * fparams->Q[0];
    F[1] = fparams->F * fparams->Q[1];
    F[2] = fparams->F * fparams->Q[2];
}


/* ----- Pos model ----- */

static void* pos_init_data(mempool_t data_pool, element_ir* elt);
static void pos_find_position(void* data, mempool_t pool, 
                                element_pos_t* elt_pos);

struct pos_params {
    double x, y, z;
    int node_ids[2];
    double* Q;
};

void pos_register()
{
    model_t pos = {
        pos_init_data,
        pos_find_position,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
        NULL
    };
 
    model_register("pos", &pos);
}

static void* pos_init_data(mempool_t data_pool, element_ir* elt)
{
    struct pos_params* pparams = (struct pos_params*)
            mempool_cget(data_pool, sizeof(struct pos_params));

    pparams->x = 0;
    pparams->y = 0;
    pparams->z = 0;

    get_double_param(&(pparams->x), "x", elt->params);
    get_double_param(&(pparams->y), "y", elt->params);
    get_double_param(&(pparams->z), "z", elt->params);

    pparams->Q = elt->Q;
    pparams->node_ids[0] = elt->node_ids[0];
    pparams->node_ids[1] = elt->node_ids[1];

    elt->var_ids = NULL;
    elt->num_vars = 0;

    return pparams;
}

static void pos_find_position(void* data, mempool_t pool, 
                                element_pos_t* elt_pos)
{
    struct pos_params* pparams = (struct pos_params*) data;

    elt_pos->num_mech_nodes = 2;
    elt_pos->relpos = (double*) 
            mempool_cget(pool, 6 * sizeof(double));

    elt_pos->relpos[0] = 0;
    elt_pos->relpos[1] = 0;
    elt_pos->relpos[2] = 0;
    
    elt_pos->relpos[3] = pparams->x;
    elt_pos->relpos[4] = pparams->y;
    elt_pos->relpos[5] = pparams->z;

    rotv2global(pparams->Q, 3, &(elt_pos->relpos[3]));
}




/* ----- Beam models ----- */

static void* beams_init_data(mempool_t data_pool, element_ir* elt);
static void beams_find_position(void* data, mempool_t pool, 
                                element_pos_t* elt_pos);
static void beams_post_position(void* data, node_ir** node_table);

static void beam2d_display(void* data, node_ir** node_table, double* dq);
static void beam2d_writegeom(void* xdrsp, void* data, node_ir** node_table);
static void beam2d_add_vars(void* data, var_ir** tmp_var_table, mempool_t pool);
static void beam2d_K(void* data, double* K);

static void beam3d_display(void* data, node_ir** node_table, double* dq);
static void beam3d_writegeom(void* xdrsp, void* data, node_ir** node_table);
static void beam3d_add_vars(void* data, var_ir** tmp_var_table, mempool_t pool);
static void beam3d_K(void* data, double* K);


struct beam_params {
    double l, w, h, E, Nu;
    int node_ids[2];
    int* var_ids;
    double* Q;
};


void beams_register()
{
    model_t beam2d = {
        beams_init_data,
        beams_find_position,
        beams_post_position,
        beam2d_display,
        beam2d_writegeom,
        beam2d_add_vars,
        beam2d_K,
        NULL
    };
 
    model_t beam3d = {
        beams_init_data,
        beams_find_position,
        beams_post_position,
        beam3d_display,
        beam3d_writegeom,
        beam3d_add_vars,
        beam3d_K,
        NULL
    };
    
    model_register("beam2d", &beam2d);
    model_register("beam3d", &beam3d);
}


static void* beams_init_data(mempool_t data_pool, element_ir* elt)
{
    struct beam_params* bparams = (struct beam_params*)
            mempool_cget(data_pool, sizeof(struct beam_params));

    bparams->l = 0;
    bparams->w = 0;
    bparams->h = 0;

    get_double_param(&(bparams->l), "l", elt->params);
    get_double_param(&(bparams->w), "w", elt->params);
    get_double_param(&(bparams->h), "h", elt->params);
    get_double_param(&(bparams->E), "Youngsmodulus", elt->params);
    get_double_param(&(bparams->Nu), "Poisson", elt->params);

    bparams->Q = elt->Q;
    bparams->node_ids[0] = elt->node_ids[0];
    bparams->node_ids[1] = elt->node_ids[1];

    elt->var_ids = mempool_cget(data_pool, 12*sizeof(int));
    elt->num_vars = 12;
    bparams->var_ids = elt->var_ids;

    return bparams;
}


static void beams_find_position(void* data, mempool_t pool, 
                                element_pos_t* elt_pos)
{
    struct beam_params* bparams = (struct beam_params*) data;

    if (bparams->l == 0) {
        elt_pos->num_mech_nodes = 0;
	return;
    }

    elt_pos->num_mech_nodes = 2;
    elt_pos->relpos = (double*) 
            mempool_cget(pool, 6 * sizeof(double));

    elt_pos->relpos[0] = 0;
    elt_pos->relpos[1] = 0;
    elt_pos->relpos[2] = 0;
    
    elt_pos->relpos[3] = bparams->l * bparams->Q[0];
    elt_pos->relpos[4] = bparams->l * bparams->Q[3];
    elt_pos->relpos[5] = bparams->l * bparams->Q[6];
}


static void beams_post_position(void* data, node_ir** node_table)
{
    #define norm(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    #define normalize(a) \
    do { \
	double l = norm(a); \
	a[0] /= l; a[1] /= l; a[2] /= l; \
    } while (0)
	
    struct beam_params* bparams = (struct beam_params*) data;
    node_ir* n1 = node_table[bparams->node_ids[0]];
    node_ir* n2 = node_table[bparams->node_ids[1]];

    double cutoff = 1e-8;
    double Rnew[9];
    double* rx = &(Rnew[0]);
    double* ry = &(Rnew[3]);
    double* rz = &(Rnew[6]);
    int i;

    if (bparams->l != 0)
        return;

    /* Reference x direction is along the beam 
     * (in closest enclosing coordinate system)
     */ 
    for (i = 0; i < 3; ++i)
        rx[i] = n2->pos[i] - n1->pos[i];
    rotv2local(bparams->Q, 3, rx);
    bparams->l = norm(rx);
    for (i = 0; i < 3; ++i)
        rx[i] /= bparams->l;

    /* Reference y and z are a little trickier */

    if (rx[0]*rx[0] + rx[1]*rx[1] > cutoff*cutoff) {

        /* If possible, the reference y should lie in the subnet x-y plane,
         * orthogonal to the projection of the reference x.
         */ 

	ry[0] = -rx[1];
	ry[1] =  rx[0];
	ry[2] =  0;
	normalize(ry);

        rz[0] = 0 - rx[2]*rx[0];
	rz[1] = 0 - rx[2]*rx[1];
	rz[2] = 1 - rx[2]*rx[2];
	normalize(rz);

    } else {

        /* If the beam sticks nearly exactly out-of-plane, make the
         * reference y direction line in the subnet x-z plane
         */

	ry[0] =  rx[2];
	ry[1] =  0;
	ry[2] = -rx[0];
	normalize(ry);

	rz[0] = 0 - rx[1]*rx[0];
	rz[1] = 1 - rx[1]*rx[1];
	rz[2] = 0 - rx[1]*rx[2];
	normalize(rz);
	
    }

    /* Now construct the rotation from local to global coordinates */
    rotv2global(bparams->Q, 3, rx);
    rotv2global(bparams->Q, 3, ry);
    rotv2global(bparams->Q, 3, rz);

    /* And finally copy it over */
    bparams->Q[0] = rx[0];
    bparams->Q[3] = rx[1];
    bparams->Q[6] = rx[2];

    bparams->Q[1] = ry[0];
    bparams->Q[4] = ry[1];
    bparams->Q[7] = ry[2];

    bparams->Q[2] = rz[0];
    bparams->Q[5] = rz[1];
    bparams->Q[8] = rz[2];

    #undef normalize
    #undef norm
}

static void beam2d_display(void* data, node_ir** node_table, double* dq)
{
#ifdef USE_XT
    struct beam_params* bparams = (struct beam_params*) data;
    node_ir* n1 = node_table[bparams->node_ids[0]];

    rotv2local(bparams->Q, 2, dq+0);
    rotv2local(bparams->Q, 2, dq+3);
    draw_beam(bparams->Q, n1->pos[0], n1->pos[1], n1->pos[2],
              bparams->l, bparams->w, bparams->h,
	      dq[0], dq[1], 0,  0, 0, dq[2],
	      dq[3], dq[4], 0,  0, 0, dq[5]);
#endif /* USE_XT*/
}


static void beam2d_writegeom(void* xdrsp, void* data, node_ir** node_table)
{
/*
    struct beam_params* bparams = (struct beam_params*) data;
    node_ir* n1 = node_table[bparams->node_ids[0]];
    int vindex[12] = {-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1};

    writebeam(xdrsp, bparams->Q, 
              bparams->l, bparams->w, bparams->h,
	      n1->pos[0], n1->pos[1], n1->pos[2], 
	      vindex);
*/
}


static void beam2d_add_vars(void* data, var_ir** var_table, mempool_t pool)
{
    struct beam_params* bparams = (struct beam_params*) data;
    int n1 = bparams->node_ids[0];
    int n2 = bparams->node_ids[1];

    bparams->var_ids[0] = merge_var(var_table, n1, "x", 'u');
    bparams->var_ids[1] = merge_var(var_table, n1, "y", 'u');
    bparams->var_ids[2] = merge_var(var_table, n1, "rz", 'u');

    bparams->var_ids[3] = merge_var(var_table, n2, "x", 'u');
    bparams->var_ids[4] = merge_var(var_table, n2, "y", 'u');
    bparams->var_ids[5] = merge_var(var_table, n2, "rz", 'u');

    bparams->var_ids[6] = merge_var(var_table, n1, "rx", 'g');
    bparams->var_ids[7] = merge_var(var_table, n1, "ry", 'g');
    bparams->var_ids[8] = merge_var(var_table, n1, "z", 'g');

    bparams->var_ids[9] = merge_var(var_table, n2, "rx", 'g');
    bparams->var_ids[10] = merge_var(var_table, n2, "ry", 'g');
    bparams->var_ids[11] = merge_var(var_table, n2, "z", 'g');
}


static void beam2d_K(void* data, double* K)
{
    struct beam_params* bparams = (struct beam_params*) data;
    
    double E = bparams->E;
    double l = bparams->l;
    double w = bparams->w;
    double h = bparams->h;
    
    double A = w*h;
    double I = (w*w*w)*h/12;
    double c = E*I/(l*l*l);

    #define Kij(i,j) K[(j-1)*12+(i-1)]

    memset(K, 0, 144*sizeof(double));

    /* Fill in the 11 block, first */
    Kij(1,1) = E*A/l;
    Kij(2,2) = 12*c;
    Kij(3,3) = 4*c*l*l;
    Kij(2,3) = Kij(3,2) = 6*c*l;

    /* Now fill in the 22 block */
    Kij(4,4) = E*A/l;
    Kij(5,5) = 12*c;
    Kij(6,6) = 4*c*l*l;
    Kij(5,6) = Kij(6,5) = -6*c*l;

    /* Now the off-diagonal blocks */
    Kij(4,1) = Kij(1,4) = -E*A/l;
    Kij(5,2) = Kij(2,5) = -12*c;
    Kij(6,3) = Kij(3,6) = 2*c*l*l;
    Kij(6,2) = Kij(2,6) = 6*c*l;
    Kij(5,3) = Kij(3,5) = -6*c*l;

    rotm2global(bparams->Q, 2, &Kij(1,1), 12);
    rotm2global(bparams->Q, 2, &Kij(4,4), 12);
    rotm2global(bparams->Q, 2, &Kij(1,4), 12);
    rotm2global(bparams->Q, 2, &Kij(4,1), 12);

    #undef Kij

    test_symm(K,12);
}


static void beam3d_display(void* data, node_ir** node_table, double* dq)
{
#ifdef USE_XT
    struct beam_params* bparams = (struct beam_params*) data;
    node_ir* n1 = node_table[bparams->node_ids[0]];

    rotv2local(bparams->Q, 3, dq+0);
    rotv2local(bparams->Q, 3, dq+3);
    rotv2local(bparams->Q, 3, dq+6);
    rotv2local(bparams->Q, 3, dq+9);
    draw_beam(bparams->Q, n1->pos[0], n1->pos[1], n1->pos[2],
              bparams->l, bparams->w, bparams->h,
	      dq[0], dq[1], dq[2],  dq[3], dq[4], dq[5],
	      dq[6], dq[7], dq[8],  dq[9], dq[10], dq[11]);
#endif /* USE_XT*/
}


static void beam3d_writegeom(void* xdrsp, void* data, node_ir** node_table)
{
	/*
    struct beam_params* bparams = (struct beam_params*) data;
    node_ir* n1 = node_table[bparams->node_ids[0]];
    int vindex[12] = {-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1};

    writebeam(xdrsp, bparams->Q, 
              bparams->l, bparams->w, bparams->h,
	      n1->pos[0], n1->pos[1], n1->pos[2], 
	      vindex);
	      */
}


static void beam3d_add_vars(void* data, var_ir** var_table, mempool_t pool)
{
    struct beam_params* bparams = (struct beam_params*) data;
    int n1 = bparams->node_ids[0];
    int n2 = bparams->node_ids[1];

    bparams->var_ids[0]  = merge_var(var_table, n1, "x", 'u');
    bparams->var_ids[1]  = merge_var(var_table, n1, "y", 'u');
    bparams->var_ids[2]  = merge_var(var_table, n1, "z", 'u');
    bparams->var_ids[3]  = merge_var(var_table, n1, "rx", 'u');
    bparams->var_ids[4]  = merge_var(var_table, n1, "ry", 'u');
    bparams->var_ids[5]  = merge_var(var_table, n1, "rz", 'u');

    bparams->var_ids[6]  = merge_var(var_table, n2, "x", 'u');
    bparams->var_ids[7]  = merge_var(var_table, n2, "y", 'u');
    bparams->var_ids[8]  = merge_var(var_table, n2, "z", 'u');
    bparams->var_ids[9]  = merge_var(var_table, n2, "rx", 'u');
    bparams->var_ids[10] = merge_var(var_table, n2, "ry", 'u');
    bparams->var_ids[11] = merge_var(var_table, n2, "rz", 'u');
}


static void beam3d_K(void* data, double* K)
{
    struct beam_params* bparams = (struct beam_params*) data;
    
    double E = bparams->E;
    double Nu = bparams->Nu;
    double L = bparams->l;
    double W = bparams->w;
    double H = bparams->h;
    
    double A = W*H;
    double G = E / (2*(1+Nu));

    double asp = (W>H ? H/W : W/H);
    double J = (W>H) ? 
        W*H*H*H*(16/3-3.36*asp*(1-asp*asp*asp*asp/12))/16 :
        H*W*W*W*(16/3-3.36*asp*(1-asp*asp*asp*asp/12))/16;

    double Iy = W*H*H*H/12;
    double Iz = H*W*W*W/12;

    double a = E * A / L;
    double b = 12 * E * Iz / (L*L*L);
    double c = 12 * E * Iy / (L*L*L);
    double d = G * J / L;
    double e = 4 * E * Iy / L;
    double f = 2 * E * Iy / L;
    double g = 4 * E * Iz / L;
    double h = 2 * E * Iz / L;
    double i = 6 * E * Iy / (L*L);
    double j = 6 * E * Iz / (L*L);

    int n, m;
   
    #define Kij(i,j) K[(j-1)*12+(i-1)]

    memset(K, 0, 144*sizeof(double));

    /* Fill in the 11 block, first */
    Kij(1,1) = a;
    Kij(2,2) = b;
    Kij(3,3) = c;
    Kij(4,4) = d;
    Kij(5,5) = e;
    Kij(6,6) = g;
    Kij(6,2) = Kij(2,6) = j;
    Kij(5,3) = Kij(3,5) = -i;

    /* Now fill in the 22 block */
    Kij(7,7) = a;
    Kij(8,8) = b;
    Kij(9,9) = c;
    Kij(10,10) = d;
    Kij(11,11) = e;
    Kij(12,12) = g;
    Kij(12,8) = Kij(8,12) = -j;
    Kij(11,9) = Kij(9,11) = i;

    /* Now the off-diagonal blocks */
    Kij(7,1) = Kij(1,7) = -a;
    Kij(8,2) = Kij(2,8) = -b;
    Kij(9,3) = Kij(3,9) = -c;
    Kij(10,4) = Kij(4,10) = -d;
    Kij(11,5) = Kij(5,11) = f;
    Kij(12,6) = Kij(6,12) = h;
    Kij(12,2) = Kij(2,12) = j;
    Kij(11,3) = Kij(3,11) = -i;
    Kij(9,5) = Kij(5,9) = i;
    Kij(8,6) = Kij(6,8) = -j;

    for (n = 0; n < 4; ++n)
        for (m = 0; m < 4; ++m)
            rotm2global(bparams->Q, 3, &Kij(3*n+1,3*m+1), 12);

    #undef Kij

    test_symm(K,12);
}


/* ----- Gap model ----- */

static void* gap_init_data(mempool_t data_pool, element_ir* elt);
static void gap_find_position(void* data, mempool_t pool, 
                              element_pos_t* elt_pos);

struct gap_params {
    double l, w1, w2, gap, permittivity;
    int node_ids[4];
    double* Q;
};

void gap_register()
{
    model_t gap = {
        gap_init_data,
        gap_find_position,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
        NULL
    };
 
    model_register("gap2dforce", &gap);
    model_register("gap3dforce", &gap);
}

static void* gap_init_data(mempool_t data_pool, element_ir* elt)
{
    struct gap_params* gparams = (struct gap_params*)
            mempool_cget(data_pool, sizeof(struct gap_params));

    gparams->l = 0;
    gparams->w1 = 0;
    gparams->w2 = 0;
    gparams->gap = 0;
    gparams->permittivity = 0;

    get_double_param(&(gparams->l), "l", elt->params);
    get_double_param(&(gparams->w1), "w1", elt->params);
    get_double_param(&(gparams->w2), "w2", elt->params);
    get_double_param(&(gparams->gap), "gap", elt->params);
    get_double_param(&(gparams->permittivity), "permittivity", elt->params);

    gparams->Q = elt->Q;
    gparams->node_ids[0] = elt->node_ids[0];
    gparams->node_ids[1] = elt->node_ids[1];
    gparams->node_ids[2] = elt->node_ids[2];
    gparams->node_ids[3] = elt->node_ids[3];

    elt->var_ids = NULL;
    elt->num_vars = 0;

    return gparams;
}

static void gap_find_position(void* data, mempool_t pool, 
                              element_pos_t* elt_pos)
{
    struct gap_params* gparams = (struct gap_params*) data;
    double l = gparams->l;
    double g = gparams->gap + (gparams->w1 + gparams->w2)/2;

    elt_pos->num_mech_nodes = 4;
    elt_pos->relpos = (double*) 
            mempool_cget(pool, 12 * sizeof(double));

    elt_pos->relpos[0]  =  0;
    elt_pos->relpos[1]  =  0;
    elt_pos->relpos[2]  =  0;
    
    elt_pos->relpos[3]  =  l;
    elt_pos->relpos[4]  =  0;
    elt_pos->relpos[5]  =  0;

    elt_pos->relpos[6]  =  0;
    elt_pos->relpos[7]  = -g;
    elt_pos->relpos[8]  =  0;

    elt_pos->relpos[9]  =  l;
    elt_pos->relpos[10] = -g;
    elt_pos->relpos[11] =  0;
    
    rotv2global(gparams->Q, 3, &(elt_pos->relpos[3]));
    rotv2global(gparams->Q, 3, &(elt_pos->relpos[6]));
    rotv2global(gparams->Q, 3, &(elt_pos->relpos[9]));
}



