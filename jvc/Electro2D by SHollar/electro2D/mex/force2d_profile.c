static char mc_version[] = "MATLAB Compiler 1.2.1 Jan 15 1999 infun";
/*
 *  MATLAB Compiler: 1.2.1
 *  Date: Jan 15 1999
 *  Arguments: -Z -i -r force2d_profile 
 */
#ifndef ARRAY_ACCESS_INLINING
#error You must use the -inline option when compiling MATLAB compiler generated code with MEX or MBUILD
#endif
#ifndef MATLAB_COMPILER_GENERATED_CODE
#define MATLAB_COMPILER_GENERATED_CODE
#endif

#include <math.h>
#include "mex.h"
#include "mcc.h"

/* static array S0_ (1 x 4) int, line 75 */
static int S0r_[] =
{
	0, 0, 0, 0, 
};
static mxArray S0_ = mccCINIT( mccINT, 1, 4, S0r_, 0 );
/* static array S1_ (1 x 4) int, line 81 */
static int S1r_[] =
{
	0, 0, 0, 0, 
	
};
static mxArray S1_ = mccCINIT( mccINT, 1, 4, S1r_, 0 );
/* static array S2_ (1 x 19) text, line 107: 'matrix calculations' */
static unsigned short S2__r_[] =
{
        109,   97,  116,  114,  105,  120,   32,   99,
         97,  108,   99,  117,  108,   97,  116,  105,
        111,  110,  115,
};
static mxArray S2_ = mccCINIT( mccTEXT,  1, 19, S2__r_, 0);
/* static array S3_ (1 x 24) text, line 149: 'calculate force matrices' */
static unsigned short S3__r_[] =
{
         99,   97,  108,   99,  117,  108,   97,  116,
        101,   32,  102,  111,  114,   99,  101,   32,
        109,   97,  116,  114,  105,   99,  101,  115,
};
static mxArray S3_ = mccCINIT( mccTEXT,  1, 24, S3__r_, 0);
/* static array S4_ (1 x 5) text, line 159: 'sub_x' */
static unsigned short S4__r_[] =
{
        115,  117,   98,   95,  120,
};
static mxArray S4_ = mccCINIT( mccTEXT,  1, 5, S4__r_, 0);
/* static array S5_ (1 x 5) text, line 160: 'sub_y' */
static unsigned short S5__r_[] =
{
        115,  117,   98,   95,  121,
};
static mxArray S5_ = mccCINIT( mccTEXT,  1, 5, S5__r_, 0);
/* static array S6_ (1 x 9) text, line 161: 'y_discont' */
static unsigned short S6__r_[] =
{
        121,   95,  100,  105,  115,   99,  111,  110,
        116,
};
static mxArray S6_ = mccCINIT( mccTEXT,  1, 9, S6__r_, 0);
/* static array S7_ (1 x 26) text, line 207: 'calculate potential matri...' */
static unsigned short S7__r_[] =
{
         99,   97,  108,   99,  117,  108,   97,  116,
        101,   32,  112,  111,  116,  101,  110,  116,
        105,   97,  108,   32,  109,   97,  116,  114,
        105,  120,
};
static mxArray S7_ = mccCINIT( mccTEXT,  1, 26, S7__r_, 0);
/* static array S8_ (1 x 23) text, line 255: 'invert potential matrix' */
static unsigned short S8__r_[] =
{
        105,  110,  118,  101,  114,  116,   32,  112,
        111,  116,  101,  110,  116,  105,   97,  108,
         32,  109,   97,  116,  114,  105,  120,
};
static mxArray S8_ = mccCINIT( mccTEXT,  1, 23, S8__r_, 0);
/***************** Compiler Assumptions ****************
 * M-File: C:/WINDOWS/Desktop/U C BERKELEY/CAD MEMS/Electro2D by SHollar/electro2D/mex/force2d_profile.m
 *
 *       F_xx        	real vector/matrix
 *       F_yy        	real vector/matrix
 *       I0_         	integer scalar temporary
 *       I1_         	integer scalar temporary
 *       IM0_        	integer vector/matrix temporary
 *       IM1_        	integer vector/matrix temporary
 *       L           	real vector/matrix
 *       R0_         	real scalar temporary
 *       R1_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       RM2_        	real vector/matrix temporary
 *       RM3_        	real vector/matrix temporary
 *       RM4_        	real vector/matrix temporary
 *       RM5_        	real vector/matrix temporary
 *       RM6_        	real vector/matrix temporary
 *       RM7_        	real vector/matrix temporary
 *       RM8_        	real vector/matrix temporary
 *       R_par       	real vector/matrix
 *       R_perp      	real vector/matrix
 *       S0_         	<constant>
 *       S1_         	<constant>
 *       X_trans     	real vector/matrix
 *       Y_cos       	real vector/matrix
 *       Y_trans     	real vector/matrix
 *       abs         	<function>
 *       adding_matrix	real vector/matrix
 *       alpha       	real vector/matrix
 *       atan        	<function>
 *       cos         	<function>
 *       cos_angle   	real vector/matrix
 *       cos_rot_theta	real vector/matrix
 *       costh       	real vector/matrix
 *       del_theta   	real vector/matrix
 *       del_x       	real vector/matrix
 *       del_y       	real vector/matrix
 *       delta       	real scalar
 *       delta_approx	real vector/matrix
 *       denom_minusL	real vector/matrix
 *       denom_plusL 	real vector/matrix
 *       dissect     	<function>
 *       div_sinth   	real vector/matrix
 *       epsilon     	real vector/matrix
 *       force2d_profile	<function being defined>
 *       inv         	<function>
 *       less_zero_minusL	real vector/matrix
 *       less_zero_plusL	real vector/matrix
 *       mseg1       	real vector/matrix
 *       mseg2       	real vector/matrix
 *       n           	integer vector/matrix
 *       n           	integer scalar  => n_1
 *       ones        	<function>
 *       permitivity 	real vector/matrix
 *       pi          	<function>
 *       pot_angle   	real vector/matrix
 *       reallog     	<function>
 *       rot_theta   	real vector/matrix
 *       s_matrix    	real vector/matrix
 *       s_plus_matrix	real vector/matrix
 *       seg1        	real vector/matrix
 *       seg2        	real vector/matrix
 *       sigma_charge	real vector/matrix
 *       sin         	<function>
 *       sin_angle   	real vector/matrix
 *       sin_rot_theta	real vector/matrix
 *       sinth       	real vector/matrix
 *       size        	<function>
 *       size_seg1   	integer scalar
 *       size_seg2   	integer scalar
 *       size_volt   	integer scalar
 *       size_volt1  	integer scalar
 *       size_volt2  	integer scalar
 *       sub_force_matrices	<function>
 *       sub_x       	real vector/matrix
 *       sub_x_sum   	real vector/matrix
 *       sub_y       	real vector/matrix
 *       sub_y_sum   	real vector/matrix
 *       tot_charge  	real vector/matrix
 *       total_segment_number	integer scalar
 *       volt_m      	real vector/matrix
 *       volt_m1     	real vector/matrix
 *       volt_m2     	real vector/matrix
 *       volt_mseg1  	real vector/matrix
 *       volt_mseg2  	real vector/matrix
 *       voltage     	real vector/matrix
 *       x_m         	real vector/matrix
 *       x_p         	real vector/matrix
 *       x_pot       	real vector/matrix
 *       x_value     	real vector/matrix
 *       y_discont   	real vector/matrix
 *       y_force_discon	real vector/matrix
 *       y_pot       	real vector/matrix
 *       y_value     	real vector/matrix
 *       z_axis      	real vector/matrix
 *       z_depth     	real vector/matrix
 *******************************************************/

void
mexFunction(
    int nlhs_,
    mxArray *plhs_[],
    int nrhs_,
    const mxArray *prhs_[]
)
{
   mxArray *Mplhs_[1];
   mxArray *Mprhs_[2];
   

   if (nrhs_ > 5 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }

   if (nlhs_ > 5 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   mcmSetLineNumber(0);
   {
      mxArray F_xx;
      mxArray F_yy;
      mxArray tot_charge;
      mxArray mseg1;
      mxArray mseg2;
      mxArray seg1;
      mxArray seg2;
      mxArray delta_approx;
      mxArray permitivity;
      mxArray z_depth;
      mxArray epsilon;
      double delta = 0.0;
      mxArray z_axis;
      mxArray adding_matrix;
      mxArray n;
      mxArray volt_mseg1;
      int n_1 = 0;
      mxArray volt_mseg2;
      mxArray volt_m1;
      mxArray volt_m2;
      int size_seg1 = 0;
      int size_seg2 = 0;
      int size_volt2 = 0;
      int size_volt1 = 0;
      int total_segment_number = 0;
      mxArray del_theta;
      mxArray costh;
      mxArray sinth;
      mxArray rot_theta;
      mxArray cos_rot_theta;
      mxArray sin_rot_theta;
      mxArray del_x;
      mxArray del_y;
      mxArray X_trans;
      mxArray Y_trans;
      mxArray Y_cos;
      mxArray R_par;
      mxArray R_perp;
      mxArray div_sinth;
      mxArray L;
      mxArray denom_plusL;
      mxArray denom_minusL;
      mxArray less_zero_plusL;
      mxArray less_zero_minusL;
      mxArray alpha;
      mxArray sub_x_sum;
      mxArray sub_x;
      mxArray sub_y_sum;
      mxArray sub_y;
      mxArray y_force_discon;
      mxArray y_discont;
      mxArray volt_m;
      int size_volt = 0;
      mxArray x_value;
      mxArray y_value;
      mxArray pot_angle;
      mxArray sin_angle;
      mxArray cos_angle;
      mxArray x_pot;
      mxArray y_pot;
      mxArray x_p;
      mxArray x_m;
      mxArray s_matrix;
      mxArray s_plus_matrix;
      mxArray voltage;
      mxArray sigma_charge;
      int I0_ = 0;
      int I1_ = 0;
      mxArray IM0_;
      mxArray IM1_;
      mxArray RM0_;
      mxArray RM1_;
      mxArray RM2_;
      mxArray RM3_;
      mxArray RM4_;
      double R0_ = 0.0;
      double R1_ = 0.0;
      mxArray RM5_;
      mxArray RM6_;
      mxArray RM7_;
      mxArray RM8_;
      
      mccRealInit(seg1);
      mccImport(&seg1, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(seg2);
      mccImport(&seg2, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(delta_approx);
      mccImport(&delta_approx, ((nrhs_>2) ? prhs_[2] : 0), 0, 0);
      mccRealInit(permitivity);
      mccImport(&permitivity, ((nrhs_>3) ? prhs_[3] : 0), 0, 0);
      mccRealInit(z_depth);
      mccImport(&z_depth, ((nrhs_>4) ? prhs_[4] : 0), 0, 0);
      mccRealInit(F_xx);
      mccRealInit(F_yy);
      mccRealInit(tot_charge);
      mccRealInit(mseg1);
      mccRealInit(mseg2);
      mccRealInit(epsilon);
      mccRealInit(z_axis);
      mccRealInit(adding_matrix);
      mccIntInit(n);
      mccRealInit(volt_mseg1);
      mccRealInit(volt_mseg2);
      mccRealInit(volt_m1);
      mccRealInit(volt_m2);
      mccRealInit(del_theta);
      mccRealInit(costh);
      mccRealInit(sinth);
      mccRealInit(rot_theta);
      mccRealInit(cos_rot_theta);
      mccRealInit(sin_rot_theta);
      mccRealInit(del_x);
      mccRealInit(del_y);
      mccRealInit(X_trans);
      mccRealInit(Y_trans);
      mccRealInit(Y_cos);
      mccRealInit(R_par);
      mccRealInit(R_perp);
      mccRealInit(div_sinth);
      mccRealInit(L);
      mccRealInit(denom_plusL);
      mccRealInit(denom_minusL);
      mccRealInit(less_zero_plusL);
      mccRealInit(less_zero_minusL);
      mccRealInit(alpha);
      mccRealInit(sub_x_sum);
      mccRealInit(sub_x);
      mccRealInit(sub_y_sum);
      mccRealInit(sub_y);
      mccRealInit(y_force_discon);
      mccRealInit(y_discont);
      mccRealInit(volt_m);
      mccRealInit(x_value);
      mccRealInit(y_value);
      mccRealInit(pot_angle);
      mccRealInit(sin_angle);
      mccRealInit(cos_angle);
      mccRealInit(x_pot);
      mccRealInit(y_pot);
      mccRealInit(x_p);
      mccRealInit(x_m);
      mccRealInit(s_matrix);
      mccRealInit(s_plus_matrix);
      mccRealInit(voltage);
      mccRealInit(sigma_charge);
      mccIntInit(IM0_);
      mccIntInit(IM1_);
      mccRealInit(RM0_);
      mccRealInit(RM1_);
      mccRealInit(RM2_);
      mccRealInit(RM3_);
      mccRealInit(RM4_);
      mccRealInit(RM5_);
      mccRealInit(RM6_);
      mccRealInit(RM7_);
      mccRealInit(RM8_);
      
      
      /* % function [F_xx, F_yy, tot_charge, mseg1, mseg2] = force2d_profile(seg1, seg2, ... */
      /* % delta_approx, permitivity, z_depth); */
      
      /* % force2d_profile calculates the force and charge distributions  */
      /* % on multiple conducting bodies.  Assuming the surfaces of the conducting  */
      /* % bodies are specified in seg1 and seg2, force2d_profile is a function that */
      /* % returns matrices of force and an array of charge.  The model is based  */
      /* % on the theory of a quasi-static electric field derived from Maxwell's */
      /* % equations.   */
      
      /* % F_xx (microNewtons/Coulomb^2) = force matrix in the X direction. */
      /* % To find the force acting on the ith element in mseg1, mseg1(i), by  */
      /* % the jth element in mseg2, mseg2(j), do this: */
      /* % jth_seg2 = j + size(mseg1,1) */
      /* % ith_seg1 = i */
      /* % Force	= tot_charge(ith_seg1)*F_xx(ith_seg1,jth_seg2)*tot_charge(jth_seg2) */
      
      /* % F_yy (microNewtons/Coulomb^2) = force matrix in the Y direction.  See F_xx */
      
      /* % tot_charge (Coulomb) = an array that contains the charge on each line */
      /* % as described by [mseg1 ; mseg2].  For example the charge on line segment */
      /* % mseg1(i) is tot_charge(i).  The charge on line segment mseg2(j) is  */
      /* % tot_charge((size(mseg1,1))+j). */
      
      /* % mseg1 (micron, micron, radians, voltage) = for finite element analysis */
      /* % seg1 is broken up into smaller segments of length delta_approx.  The  */
      /* % array of line segments comprising seg1 are mseg1.   */
      
      /* % The format of mseg1 is as follows: */
      /* % mseg1(i,:) = [ center_x center_y angle voltage ] where center_x, center_y */
      /* % represent the center of the line segment, angle is the angle the line segment */
      /* % makes with respect to the x axis, and voltage is the voltage of the  */
      /* % line segment. */
      
      /* % mseg2 (micron, micron, radians, voltage) = for finite element analysis, */
      /* % seg2 is broken up into smaller segments of length delta_approx.  The */
      /* % format is the same as mseg1. */
      
      /* % seg1 (microns, Volts) = array of line segments defining  */
      /* % a conducting body. */
      /* % Each line is defined as [x1 y1 x2 y2 Voltage] where  (x1,y1)  */
      /* % and (x2,y2) define the endpoints of the line and Voltage defines */
      /* % the voltage on the line. */
      
      /* % seg2 = same as seg1 */
      
      /* % delta_approx (microns) = seg1 and seg2 are broken up  */
      /* % into smaller line segments to approximate the ideal  */
      /* % behavior.  The length of these smaller line segments  */
      /* % is delta approx.  The smaller the delta_approx the */
      /* % more accurate the results.  Nevertheless, simulation  */
      /* % time is on order O(N^2) where N=1/delta_approx (ie  */
      /* % the smaller the delta approx, the longer it takes  */
      /* % to simulate). */
      
      /* % permitivity of free space is 8.854e-6 (mks units) * 1e-6.   */
      /* % Permitivity value is scaled by 1e6 to avoid truncation error */
      /* % in matlab. */
      
      /* % z_depth = since this is only a 2D model, z_depth defines the third */
      /* % dimension (ie the depth of conducting bodies out of the plane)  */
      
      
      /* epsilon = permitivity; % program calculates based on 1e6 scale factor */
      if(mccNOTSET(&permitivity))
      {
         mexErrMsgTxt( "variable permitivity undefined, line 67" );
      }
      mccCopy(&epsilon, &permitivity);
      /* delta=1e-15; % truncation error factor */
      delta = 1e-15;
      /* z_axis=z_depth;%microns */
      if(mccNOTSET(&z_depth))
      {
         mexErrMsgTxt( "variable z_depth undefined, line 69" );
      }
      mccCopy(&z_axis, &z_depth);
      
      
      /* % this section takes in seg1 and seg2 and breaks the lines up  */
      /* % into smaller line segments of length delta approx */
      
      /* adding_matrix=[0 0 0 0]; */
      mccCopy(&adding_matrix, &S0_);
      mccSTRING(&adding_matrix) = 0;
      /* for	n=1:size(seg1,1) */
      if(mccNOTSET(&seg1))
      {
         mexErrMsgTxt( "variable seg1 undefined, line 76" );
      }
      I1_ = mccGetDimensionSize(&seg1, 1);
      mccIntColon2(&IM0_, 1, I1_);
      mccCreateEmpty( &n );
      for (I0_ = 0; I0_ < mccN(&IM0_); ++I0_)
      {
         mccForCol(&n, &IM0_, I0_);
         /* adding_matrix=[adding_matrix;dissect(seg1(n,:),delta_approx)]; */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_seg1;
            int I_seg1=1, J_seg1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM1_) * mccN(&IM1_)), mccN(&seg1));
            mccAllocateMatrix(&RM0_, m_, n_);
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            J_seg1 = ((mccM(&seg1) != 1 || mccN(&seg1) != 1) ? mccM(&seg1) : 0);
            p_seg1 = mccPR(&seg1) + mccM(&seg1) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_seg1 += J_seg1)
               {
                  p_IM1_ = mccIPR(&IM1_);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_IM1_+=I_IM1_)
                  {
                     *p_RM0_ = p_seg1[((int)(*p_IM1_ - .5))];
                  }
               }
            }
         }
         if(mccNOTSET(&delta_approx))
         {
            mexErrMsgTxt( "variable delta_approx undefined, line 77" );
         }
         Mprhs_[0] = &RM0_;
         Mprhs_[1] = &delta_approx;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "dissect", 77);
         mccCatenateRows(&adding_matrix, &adding_matrix, &RM1_);
         /* end */
      }
      /* volt_mseg1=adding_matrix(2:size(adding_matrix,1),:); */
      I1_ = mccGetDimensionSize(&adding_matrix, 1);
      mccIntColon2(&IM0_, 2, I1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_volt_mseg1;
         int I_volt_mseg1=1;
         double *p_adding_matrix;
         int I_adding_matrix=1, J_adding_matrix;
         int *p_IM0_;
         int I_IM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM0_) * mccN(&IM0_)), mccN(&adding_matrix));
         mccAllocateMatrix(&volt_mseg1, m_, n_);
         I_volt_mseg1 = (mccM(&volt_mseg1) != 1 || mccN(&volt_mseg1) != 1);
         p_volt_mseg1 = mccPR(&volt_mseg1);
         J_adding_matrix = ((mccM(&adding_matrix) != 1 || mccN(&adding_matrix) != 1) ? mccM(&adding_matrix) : 0);
         p_adding_matrix = mccPR(&adding_matrix) + mccM(&adding_matrix) * 0;
         I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_adding_matrix += J_adding_matrix)
            {
               p_IM0_ = mccIPR(&IM0_);
               for (i_=0; i_<m_; ++i_, p_volt_mseg1+=I_volt_mseg1, p_IM0_+=I_IM0_)
               {
                  *p_volt_mseg1 = p_adding_matrix[((int)(*p_IM0_ - .5))];
               }
            }
         }
      }
      mccSTRING(&volt_mseg1) = mccSTRING(&adding_matrix);
      
      /* adding_matrix=[0 0 0 0]; */
      mccCopy(&adding_matrix, &S1_);
      mccSTRING(&adding_matrix) = 0;
      /* for	n=1:size(seg2,1) */
      if(mccNOTSET(&seg2))
      {
         mexErrMsgTxt( "variable seg2 undefined, line 82" );
      }
      I0_ = mccGetDimensionSize(&seg2, 1);
      for (I1_ = 1; I1_ <= I0_; I1_ = I1_ + 1)
      {
         n_1 = I1_;
         /* adding_matrix=[adding_matrix;dissect(seg2(n,:),delta_approx)]; */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM1_;
            int I_RM1_=1;
            double *p_seg2;
            int I_seg2=1, J_seg2;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&seg2));
            mccAllocateMatrix(&RM1_, m_, n_);
            I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
            p_RM1_ = mccPR(&RM1_);
            if (mccN(&seg2) == 1) { I_seg2 = J_seg2 = 0;}
            else { I_seg2 = 1; J_seg2=mccM(&seg2)-m_; }
            p_seg2 = mccPR(&seg2) + (n_1-1) + mccM(&seg2) * 0;
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_seg2 += J_seg2)
               {
                  for (i_=0; i_<m_; ++i_, p_RM1_+=I_RM1_, p_seg2+=I_seg2)
                  {
                     *p_RM1_ = *p_seg2;
                  }
               }
            }
         }
         if(mccNOTSET(&delta_approx))
         {
            mexErrMsgTxt( "variable delta_approx undefined, line 83" );
         }
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &delta_approx;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "dissect", 83);
         mccCatenateRows(&adding_matrix, &adding_matrix, &RM0_);
         /* end */
      }
      /* volt_mseg2=adding_matrix(2:size(adding_matrix,1),:); */
      I0_ = mccGetDimensionSize(&adding_matrix, 1);
      mccIntColon2(&IM0_, 2, I0_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_volt_mseg2;
         int I_volt_mseg2=1;
         double *p_adding_matrix;
         int I_adding_matrix=1, J_adding_matrix;
         int *p_IM0_;
         int I_IM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM0_) * mccN(&IM0_)), mccN(&adding_matrix));
         mccAllocateMatrix(&volt_mseg2, m_, n_);
         I_volt_mseg2 = (mccM(&volt_mseg2) != 1 || mccN(&volt_mseg2) != 1);
         p_volt_mseg2 = mccPR(&volt_mseg2);
         J_adding_matrix = ((mccM(&adding_matrix) != 1 || mccN(&adding_matrix) != 1) ? mccM(&adding_matrix) : 0);
         p_adding_matrix = mccPR(&adding_matrix) + mccM(&adding_matrix) * 0;
         I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_adding_matrix += J_adding_matrix)
            {
               p_IM0_ = mccIPR(&IM0_);
               for (i_=0; i_<m_; ++i_, p_volt_mseg2+=I_volt_mseg2, p_IM0_+=I_IM0_)
               {
                  *p_volt_mseg2 = p_adding_matrix[((int)(*p_IM0_ - .5))];
               }
            }
         }
      }
      mccSTRING(&volt_mseg2) = mccSTRING(&adding_matrix);
      
      /* % set volt_m1 and volt_m2 to refer to all segments */
      /* % need to get force profile matrix */
      
      /* mseg1 = volt_mseg1; */
      mccCopy(&mseg1, &volt_mseg1);
      mccSTRING(&mseg1) = mccSTRING(&volt_mseg1);
      /* mseg2 = volt_mseg2; */
      mccCopy(&mseg2, &volt_mseg2);
      mccSTRING(&mseg2) = mccSTRING(&volt_mseg2);
      
      /* volt_m1=[volt_mseg1;volt_mseg2]; */
      mccCatenateRows(&volt_m1, &volt_mseg1, &volt_mseg2);
      /* volt_m2=volt_m1; */
      mccCopy(&volt_m2, &volt_m1);
      
      /* size_seg1 = size(volt_mseg1,1); */
      size_seg1 = mccGetDimensionSize(&volt_mseg1, 1);
      /* size_seg2 = size(volt_mseg2,1); */
      size_seg2 = mccGetDimensionSize(&volt_mseg2, 1);
      
      /* size_volt2 = size(volt_m2,1); */
      size_volt2 = mccGetDimensionSize(&volt_m2, 1);
      /* size_volt1 = size(volt_m1,1); */
      size_volt1 = mccGetDimensionSize(&volt_m1, 1);
      
      /* % total segment number represents the total number of small segments */
      /* % size_volt1 is total segment numbers */
      /* total_segment_number=size_volt1 */
      total_segment_number = size_volt1;
      mccPrint (mccTempMatrix(total_segment_number, 0., mccINT, 0 ), "total_segment_number");
      
      
      /* 'matrix calculations' */
      mccPrint (&S2_, 0);
      
      /* % all calculations are done in matrix form.  This consumes quite a bit */
      /* % of memory but since matlab does fast matrix calculations, it is  */
      /* % very fast */
      
      /* % rotate line segments so that they lie on the x axis.  This is needed */
      /* % to do the force calculations */
      
      /* del_theta=ones(size_volt1,1)*volt_m2(:,3)'... */
      mccOnesMN(&IM0_, size_volt1, 1);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_volt_m2;
         int I_volt_m2=1, J_volt_m2;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m2), 1);
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         if (mccM(&volt_m2) == 1) { I_volt_m2 = J_volt_m2 = 0;}
         else { I_volt_m2 = 1; J_volt_m2=mccM(&volt_m2)-m_; }
         p_volt_m2 = mccPR(&volt_m2) + 0 + mccM(&volt_m2) * (3-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m2 += J_volt_m2)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_volt_m2+=I_volt_m2)
               {
                  *p_RM0_ = *p_volt_m2;
               }
            }
         }
      }
      mccConjTrans(&RM1_, &RM0_);
      mccRealMatrixMultiply(&RM2_, &IM0_, &RM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM3_;
         int I_RM3_=1;
         double *p_volt_m1;
         int I_volt_m1=1, J_volt_m1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m1), 1);
         mccAllocateMatrix(&RM3_, m_, n_);
         I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
         p_RM3_ = mccPR(&RM3_);
         if (mccM(&volt_m1) == 1) { I_volt_m1 = J_volt_m1 = 0;}
         else { I_volt_m1 = 1; J_volt_m1=mccM(&volt_m1)-m_; }
         p_volt_m1 = mccPR(&volt_m1) + 0 + mccM(&volt_m1) * (3-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m1 += J_volt_m1)
            {
               for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_volt_m1+=I_volt_m1)
               {
                  *p_RM3_ = *p_volt_m1;
               }
            }
         }
      }
      mccOnesMN(&IM1_, 1, size_volt2);
      mccRealMatrixMultiply(&RM4_, &RM3_, &IM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_del_theta;
         int I_del_theta=1;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM2_), mccN(&RM2_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&del_theta, m_, n_);
         I_del_theta = (mccM(&del_theta) != 1 || mccN(&del_theta) != 1);
         p_del_theta = mccPR(&del_theta);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_del_theta+=I_del_theta, p_RM2_+=I_RM2_, p_RM4_+=I_RM4_)
               {
                  *p_del_theta = (*p_RM2_ - *p_RM4_);
               }
            }
         }
      }
      /* costh=cos(del_theta); */
      mccCos(&costh, &del_theta);
      /* sinth=sin(del_theta); */
      mccSin(&sinth, &del_theta);
      /* clear del_theta */
      /* This comment replaces a call on clear */
      
      
      /* rot_theta=volt_m1(:,3)*ones(1,size_volt2); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_volt_m1;
         int I_volt_m1=1, J_volt_m1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m1), 1);
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (mccM(&volt_m1) == 1) { I_volt_m1 = J_volt_m1 = 0;}
         else { I_volt_m1 = 1; J_volt_m1=mccM(&volt_m1)-m_; }
         p_volt_m1 = mccPR(&volt_m1) + 0 + mccM(&volt_m1) * (3-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m1 += J_volt_m1)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_volt_m1+=I_volt_m1)
               {
                  *p_RM4_ = *p_volt_m1;
               }
            }
         }
      }
      mccOnesMN(&IM1_, 1, size_volt2);
      mccRealMatrixMultiply(&rot_theta, &RM4_, &IM1_);
      /* cos_rot_theta=cos(rot_theta); */
      mccCos(&cos_rot_theta, &rot_theta);
      /* sin_rot_theta=sin(rot_theta); */
      mccSin(&sin_rot_theta, &rot_theta);
      /* clear rot_theta */
      /* This comment replaces a call on clear */
      
      /* del_x=ones(size_volt1,1)*volt_m2(:,1)'... */
      mccOnesMN(&IM1_, size_volt1, 1);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_volt_m2;
         int I_volt_m2=1, J_volt_m2;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m2), 1);
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (mccM(&volt_m2) == 1) { I_volt_m2 = J_volt_m2 = 0;}
         else { I_volt_m2 = 1; J_volt_m2=mccM(&volt_m2)-m_; }
         p_volt_m2 = mccPR(&volt_m2) + 0 + mccM(&volt_m2) * (1-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m2 += J_volt_m2)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_volt_m2+=I_volt_m2)
               {
                  *p_RM4_ = *p_volt_m2;
               }
            }
         }
      }
      mccConjTrans(&RM3_, &RM4_);
      mccRealMatrixMultiply(&RM2_, &IM1_, &RM3_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM1_;
         int I_RM1_=1;
         double *p_volt_m1;
         int I_volt_m1=1, J_volt_m1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m1), 1);
         mccAllocateMatrix(&RM1_, m_, n_);
         I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
         p_RM1_ = mccPR(&RM1_);
         if (mccM(&volt_m1) == 1) { I_volt_m1 = J_volt_m1 = 0;}
         else { I_volt_m1 = 1; J_volt_m1=mccM(&volt_m1)-m_; }
         p_volt_m1 = mccPR(&volt_m1) + 0 + mccM(&volt_m1) * (1-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m1 += J_volt_m1)
            {
               for (i_=0; i_<m_; ++i_, p_RM1_+=I_RM1_, p_volt_m1+=I_volt_m1)
               {
                  *p_RM1_ = *p_volt_m1;
               }
            }
         }
      }
      mccOnesMN(&IM0_, 1, size_volt2);
      mccRealMatrixMultiply(&RM0_, &RM1_, &IM0_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_del_x;
         int I_del_x=1;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_RM0_;
         int I_RM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM2_), mccN(&RM2_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
         mccAllocateMatrix(&del_x, m_, n_);
         I_del_x = (mccM(&del_x) != 1 || mccN(&del_x) != 1);
         p_del_x = mccPR(&del_x);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_del_x+=I_del_x, p_RM2_+=I_RM2_, p_RM0_+=I_RM0_)
               {
                  *p_del_x = (*p_RM2_ - *p_RM0_);
               }
            }
         }
      }
      /* del_y=ones(size_volt1,1)*volt_m2(:,2)'... */
      mccOnesMN(&IM0_, size_volt1, 1);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_volt_m2;
         int I_volt_m2=1, J_volt_m2;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m2), 1);
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         if (mccM(&volt_m2) == 1) { I_volt_m2 = J_volt_m2 = 0;}
         else { I_volt_m2 = 1; J_volt_m2=mccM(&volt_m2)-m_; }
         p_volt_m2 = mccPR(&volt_m2) + 0 + mccM(&volt_m2) * (2-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m2 += J_volt_m2)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_volt_m2+=I_volt_m2)
               {
                  *p_RM0_ = *p_volt_m2;
               }
            }
         }
      }
      mccConjTrans(&RM1_, &RM0_);
      mccRealMatrixMultiply(&RM2_, &IM0_, &RM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM3_;
         int I_RM3_=1;
         double *p_volt_m1;
         int I_volt_m1=1, J_volt_m1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m1), 1);
         mccAllocateMatrix(&RM3_, m_, n_);
         I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
         p_RM3_ = mccPR(&RM3_);
         if (mccM(&volt_m1) == 1) { I_volt_m1 = J_volt_m1 = 0;}
         else { I_volt_m1 = 1; J_volt_m1=mccM(&volt_m1)-m_; }
         p_volt_m1 = mccPR(&volt_m1) + 0 + mccM(&volt_m1) * (2-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m1 += J_volt_m1)
            {
               for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_volt_m1+=I_volt_m1)
               {
                  *p_RM3_ = *p_volt_m1;
               }
            }
         }
      }
      mccOnesMN(&IM1_, 1, size_volt2);
      mccRealMatrixMultiply(&RM4_, &RM3_, &IM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_del_y;
         int I_del_y=1;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM2_), mccN(&RM2_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&del_y, m_, n_);
         I_del_y = (mccM(&del_y) != 1 || mccN(&del_y) != 1);
         p_del_y = mccPR(&del_y);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_del_y+=I_del_y, p_RM2_+=I_RM2_, p_RM4_+=I_RM4_)
               {
                  *p_del_y = (*p_RM2_ - *p_RM4_);
               }
            }
         }
      }
      
      /* X_trans=cos_rot_theta.*del_x... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_X_trans;
         int I_X_trans=1;
         double *p_cos_rot_theta;
         int I_cos_rot_theta=1;
         double *p_del_x;
         int I_del_x=1;
         double *p_sin_rot_theta;
         int I_sin_rot_theta=1;
         double *p_del_y;
         int I_del_y=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_rot_theta), mccN(&cos_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&del_x), mccN(&del_x));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_rot_theta), mccN(&sin_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&del_y), mccN(&del_y));
         mccAllocateMatrix(&X_trans, m_, n_);
         I_X_trans = (mccM(&X_trans) != 1 || mccN(&X_trans) != 1);
         p_X_trans = mccPR(&X_trans);
         I_cos_rot_theta = (mccM(&cos_rot_theta) != 1 || mccN(&cos_rot_theta) != 1);
         p_cos_rot_theta = mccPR(&cos_rot_theta);
         I_del_x = (mccM(&del_x) != 1 || mccN(&del_x) != 1);
         p_del_x = mccPR(&del_x);
         I_sin_rot_theta = (mccM(&sin_rot_theta) != 1 || mccN(&sin_rot_theta) != 1);
         p_sin_rot_theta = mccPR(&sin_rot_theta);
         I_del_y = (mccM(&del_y) != 1 || mccN(&del_y) != 1);
         p_del_y = mccPR(&del_y);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_X_trans+=I_X_trans, p_cos_rot_theta+=I_cos_rot_theta, p_del_x+=I_del_x, p_sin_rot_theta+=I_sin_rot_theta, p_del_y+=I_del_y)
               {
                  *p_X_trans = ((*p_cos_rot_theta * (double) *p_del_x) + (*p_sin_rot_theta * (double) *p_del_y));
               }
            }
         }
      }
      /* Y_trans=cos_rot_theta.*del_y... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Y_trans;
         int I_Y_trans=1;
         double *p_cos_rot_theta;
         int I_cos_rot_theta=1;
         double *p_del_y;
         int I_del_y=1;
         double *p_sin_rot_theta;
         int I_sin_rot_theta=1;
         double *p_del_x;
         int I_del_x=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_rot_theta), mccN(&cos_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&del_y), mccN(&del_y));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_rot_theta), mccN(&sin_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&del_x), mccN(&del_x));
         mccAllocateMatrix(&Y_trans, m_, n_);
         I_Y_trans = (mccM(&Y_trans) != 1 || mccN(&Y_trans) != 1);
         p_Y_trans = mccPR(&Y_trans);
         I_cos_rot_theta = (mccM(&cos_rot_theta) != 1 || mccN(&cos_rot_theta) != 1);
         p_cos_rot_theta = mccPR(&cos_rot_theta);
         I_del_y = (mccM(&del_y) != 1 || mccN(&del_y) != 1);
         p_del_y = mccPR(&del_y);
         I_sin_rot_theta = (mccM(&sin_rot_theta) != 1 || mccN(&sin_rot_theta) != 1);
         p_sin_rot_theta = mccPR(&sin_rot_theta);
         I_del_x = (mccM(&del_x) != 1 || mccN(&del_x) != 1);
         p_del_x = mccPR(&del_x);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Y_trans+=I_Y_trans, p_cos_rot_theta+=I_cos_rot_theta, p_del_y+=I_del_y, p_sin_rot_theta+=I_sin_rot_theta, p_del_x+=I_del_x)
               {
                  *p_Y_trans = ((*p_cos_rot_theta * (double) *p_del_y) - (*p_sin_rot_theta * (double) *p_del_x));
               }
            }
         }
      }
      /* clear del_x del_y */
      /* This comment replaces a call on clear */
      
      /* Y_cos=Y_trans.*costh; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Y_cos;
         int I_Y_cos=1;
         double *p_Y_trans;
         int I_Y_trans=1;
         double *p_costh;
         int I_costh=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&Y_trans), mccN(&Y_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&costh), mccN(&costh));
         mccAllocateMatrix(&Y_cos, m_, n_);
         I_Y_cos = (mccM(&Y_cos) != 1 || mccN(&Y_cos) != 1);
         p_Y_cos = mccPR(&Y_cos);
         I_Y_trans = (mccM(&Y_trans) != 1 || mccN(&Y_trans) != 1);
         p_Y_trans = mccPR(&Y_trans);
         I_costh = (mccM(&costh) != 1 || mccN(&costh) != 1);
         p_costh = mccPR(&costh);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Y_cos+=I_Y_cos, p_Y_trans+=I_Y_trans, p_costh+=I_costh)
               {
                  *p_Y_cos = (*p_Y_trans * (double) *p_costh);
               }
            }
         }
      }
      /* R_par=X_trans.*costh+Y_trans.*sinth; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_R_par;
         int I_R_par=1;
         double *p_X_trans;
         int I_X_trans=1;
         double *p_costh;
         int I_costh=1;
         double *p_Y_trans;
         int I_Y_trans=1;
         double *p_sinth;
         int I_sinth=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&X_trans), mccN(&X_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&costh), mccN(&costh));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&Y_trans), mccN(&Y_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sinth), mccN(&sinth));
         mccAllocateMatrix(&R_par, m_, n_);
         I_R_par = (mccM(&R_par) != 1 || mccN(&R_par) != 1);
         p_R_par = mccPR(&R_par);
         I_X_trans = (mccM(&X_trans) != 1 || mccN(&X_trans) != 1);
         p_X_trans = mccPR(&X_trans);
         I_costh = (mccM(&costh) != 1 || mccN(&costh) != 1);
         p_costh = mccPR(&costh);
         I_Y_trans = (mccM(&Y_trans) != 1 || mccN(&Y_trans) != 1);
         p_Y_trans = mccPR(&Y_trans);
         I_sinth = (mccM(&sinth) != 1 || mccN(&sinth) != 1);
         p_sinth = mccPR(&sinth);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_R_par+=I_R_par, p_X_trans+=I_X_trans, p_costh+=I_costh, p_Y_trans+=I_Y_trans, p_sinth+=I_sinth)
               {
                  *p_R_par = ((*p_X_trans * (double) *p_costh) + (*p_Y_trans * (double) *p_sinth));
               }
            }
         }
      }
      /* R_perp=-X_trans.*sinth+Y_cos; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_R_perp;
         int I_R_perp=1;
         double *p_X_trans;
         int I_X_trans=1;
         double *p_sinth;
         int I_sinth=1;
         double *p_Y_cos;
         int I_Y_cos=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&X_trans), mccN(&X_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sinth), mccN(&sinth));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&Y_cos), mccN(&Y_cos));
         mccAllocateMatrix(&R_perp, m_, n_);
         I_R_perp = (mccM(&R_perp) != 1 || mccN(&R_perp) != 1);
         p_R_perp = mccPR(&R_perp);
         I_X_trans = (mccM(&X_trans) != 1 || mccN(&X_trans) != 1);
         p_X_trans = mccPR(&X_trans);
         I_sinth = (mccM(&sinth) != 1 || mccN(&sinth) != 1);
         p_sinth = mccPR(&sinth);
         I_Y_cos = (mccM(&Y_cos) != 1 || mccN(&Y_cos) != 1);
         p_Y_cos = mccPR(&Y_cos);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_R_perp+=I_R_perp, p_X_trans+=I_X_trans, p_sinth+=I_sinth, p_Y_cos+=I_Y_cos)
               {
                  *p_R_perp = (((-*p_X_trans) * (double) *p_sinth) + *p_Y_cos);
               }
            }
         }
      }
      /* div_sinth=(2*(sinth>=0)-1)./(abs(sinth)+delta); */
      mccAbs(&RM4_, &sinth);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_div_sinth;
         int I_div_sinth=1;
         double *p_sinth;
         int I_sinth=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sinth), mccN(&sinth));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&div_sinth, m_, n_);
         I_div_sinth = (mccM(&div_sinth) != 1 || mccN(&div_sinth) != 1);
         p_div_sinth = mccPR(&div_sinth);
         I_sinth = (mccM(&sinth) != 1 || mccN(&sinth) != 1);
         p_sinth = mccPR(&sinth);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_div_sinth+=I_div_sinth, p_sinth+=I_sinth, p_RM4_+=I_RM4_)
               {
                  *p_div_sinth = (((2 * (double) ( (*p_sinth >= 0) && !mccREL_NAN(*p_sinth) )) - 1) / (double) (*p_RM4_ + delta));
               }
            }
         }
      }
      
      
      /* % once rotated we can now calculate the force contribution from each segment */
      /* % in seg1 to each segment in seg2.  Matrices are of the size  */
      /* % size_volt_m1 X size_volt_m2 */
      
      /* 'calculate force matrices' */
      mccPrint (&S3_, 0);
      
      /* L=delta_approx/2; */
      if(mccNOTSET(&delta_approx))
      {
         mexErrMsgTxt( "variable delta_approx undefined, line 151" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_L;
         int I_L=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&L, m_, n_);
         I_L = (mccM(&L) != 1 || mccN(&L) != 1);
         p_L = mccPR(&L);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_L+=I_L, p_delta_approx+=I_delta_approx)
               {
                  *p_L = (*p_delta_approx / (double) 2);
               }
            }
         }
      }
      /* denom_plusL=Y_trans+L*sinth; */
      mccRealMatrixMultiply(&RM4_, &L, &sinth);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_denom_plusL;
         int I_denom_plusL=1;
         double *p_Y_trans;
         int I_Y_trans=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&Y_trans), mccN(&Y_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&denom_plusL, m_, n_);
         I_denom_plusL = (mccM(&denom_plusL) != 1 || mccN(&denom_plusL) != 1);
         p_denom_plusL = mccPR(&denom_plusL);
         I_Y_trans = (mccM(&Y_trans) != 1 || mccN(&Y_trans) != 1);
         p_Y_trans = mccPR(&Y_trans);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_denom_plusL+=I_denom_plusL, p_Y_trans+=I_Y_trans, p_RM4_+=I_RM4_)
               {
                  *p_denom_plusL = (*p_Y_trans + *p_RM4_);
               }
            }
         }
      }
      /* denom_minusL=Y_trans-L*sinth; */
      mccRealMatrixMultiply(&RM4_, &L, &sinth);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_denom_minusL;
         int I_denom_minusL=1;
         double *p_Y_trans;
         int I_Y_trans=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&Y_trans), mccN(&Y_trans));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&denom_minusL, m_, n_);
         I_denom_minusL = (mccM(&denom_minusL) != 1 || mccN(&denom_minusL) != 1);
         p_denom_minusL = mccPR(&denom_minusL);
         I_Y_trans = (mccM(&Y_trans) != 1 || mccN(&Y_trans) != 1);
         p_Y_trans = mccPR(&Y_trans);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_denom_minusL+=I_denom_minusL, p_Y_trans+=I_Y_trans, p_RM4_+=I_RM4_)
               {
                  *p_denom_minusL = (*p_Y_trans - *p_RM4_);
               }
            }
         }
      }
      /* less_zero_plusL=(2*(denom_plusL>=0)-1); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_less_zero_plusL;
         int I_less_zero_plusL=1;
         double *p_denom_plusL;
         int I_denom_plusL=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&denom_plusL), mccN(&denom_plusL));
         mccAllocateMatrix(&less_zero_plusL, m_, n_);
         I_less_zero_plusL = (mccM(&less_zero_plusL) != 1 || mccN(&less_zero_plusL) != 1);
         p_less_zero_plusL = mccPR(&less_zero_plusL);
         I_denom_plusL = (mccM(&denom_plusL) != 1 || mccN(&denom_plusL) != 1);
         p_denom_plusL = mccPR(&denom_plusL);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_less_zero_plusL+=I_less_zero_plusL, p_denom_plusL+=I_denom_plusL)
               {
                  *p_less_zero_plusL = ((2 * (double) ( (*p_denom_plusL >= 0) && !mccREL_NAN(*p_denom_plusL) )) - 1);
               }
            }
         }
      }
      /* less_zero_minusL=(2*(denom_minusL>=0)-1); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_less_zero_minusL;
         int I_less_zero_minusL=1;
         double *p_denom_minusL;
         int I_denom_minusL=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&denom_minusL), mccN(&denom_minusL));
         mccAllocateMatrix(&less_zero_minusL, m_, n_);
         I_less_zero_minusL = (mccM(&less_zero_minusL) != 1 || mccN(&less_zero_minusL) != 1);
         p_less_zero_minusL = mccPR(&less_zero_minusL);
         I_denom_minusL = (mccM(&denom_minusL) != 1 || mccN(&denom_minusL) != 1);
         p_denom_minusL = mccPR(&denom_minusL);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_less_zero_minusL+=I_less_zero_minusL, p_denom_minusL+=I_denom_minusL)
               {
                  *p_less_zero_minusL = ((2 * (double) ( (*p_denom_minusL >= 0) && !mccREL_NAN(*p_denom_minusL) )) - 1);
               }
            }
         }
      }
      
      /* alpha=-delta_approx/2; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_alpha;
         int I_alpha=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&alpha, m_, n_);
         I_alpha = (mccM(&alpha) != 1 || mccN(&alpha) != 1);
         p_alpha = mccPR(&alpha);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_alpha+=I_alpha, p_delta_approx+=I_delta_approx)
               {
                  *p_alpha = ((-*p_delta_approx) / (double) 2);
               }
            }
         }
      }
      /* sub_force_matrices; */
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 0, Mprhs_, "sub_force_matrices", 158);
      /* sub_x_sum=sub_x; */
      mccUndefVariable(&sub_x_sum, &S4_);
      /* sub_y_sum=sub_y; */
      mccUndefVariable(&sub_y_sum, &S5_);
      /* y_force_discon=-y_discont; */
      mccUndefVariable(&y_discont, &S6_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_y_force_discon;
         int I_y_force_discon=1;
         double *p_y_discont;
         int I_y_discont=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_discont), mccN(&y_discont));
         mccAllocateMatrix(&y_force_discon, m_, n_);
         I_y_force_discon = (mccM(&y_force_discon) != 1 || mccN(&y_force_discon) != 1);
         p_y_force_discon = mccPR(&y_force_discon);
         I_y_discont = (mccM(&y_discont) != 1 || mccN(&y_discont) != 1);
         p_y_discont = mccPR(&y_discont);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_y_force_discon+=I_y_force_discon, p_y_discont+=I_y_discont)
               {
                  *p_y_force_discon = (-*p_y_discont);
               }
            }
         }
      }
      
      /* alpha= delta_approx/2; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_alpha;
         int I_alpha=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&alpha, m_, n_);
         I_alpha = (mccM(&alpha) != 1 || mccN(&alpha) != 1);
         p_alpha = mccPR(&alpha);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_alpha+=I_alpha, p_delta_approx+=I_delta_approx)
               {
                  *p_alpha = (*p_delta_approx / (double) 2);
               }
            }
         }
      }
      
      /* % sub_force_matrices in a .m file that actually calculates the force */
      /* % on each segment */
      
      /* sub_force_matrices; */
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 0, Mprhs_, "sub_force_matrices", 168);
      /* sub_x_sum=sub_x_sum-sub_x; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_sub_x_sum;
         int I_sub_x_sum=1;
         double *p_1sub_x_sum;
         int I_1sub_x_sum=1;
         double *p_sub_x;
         int I_sub_x=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_x_sum), mccN(&sub_x_sum));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_x), mccN(&sub_x));
         mccAllocateMatrix(&sub_x_sum, m_, n_);
         I_sub_x_sum = (mccM(&sub_x_sum) != 1 || mccN(&sub_x_sum) != 1);
         p_sub_x_sum = mccPR(&sub_x_sum);
         I_1sub_x_sum = (mccM(&sub_x_sum) != 1 || mccN(&sub_x_sum) != 1);
         p_1sub_x_sum = mccPR(&sub_x_sum);
         I_sub_x = (mccM(&sub_x) != 1 || mccN(&sub_x) != 1);
         p_sub_x = mccPR(&sub_x);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_sub_x_sum+=I_sub_x_sum, p_1sub_x_sum+=I_1sub_x_sum, p_sub_x+=I_sub_x)
               {
                  *p_sub_x_sum = (*p_1sub_x_sum - *p_sub_x);
               }
            }
         }
      }
      /* sub_y_sum=sub_y_sum-sub_y; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_sub_y_sum;
         int I_sub_y_sum=1;
         double *p_1sub_y_sum;
         int I_1sub_y_sum=1;
         double *p_sub_y;
         int I_sub_y=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_y_sum), mccN(&sub_y_sum));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_y), mccN(&sub_y));
         mccAllocateMatrix(&sub_y_sum, m_, n_);
         I_sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_sub_y_sum = mccPR(&sub_y_sum);
         I_1sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_1sub_y_sum = mccPR(&sub_y_sum);
         I_sub_y = (mccM(&sub_y) != 1 || mccN(&sub_y) != 1);
         p_sub_y = mccPR(&sub_y);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_sub_y_sum+=I_sub_y_sum, p_1sub_y_sum+=I_1sub_y_sum, p_sub_y+=I_sub_y)
               {
                  *p_sub_y_sum = (*p_1sub_y_sum - *p_sub_y);
               }
            }
         }
      }
      /* y_force_discon=y_force_discon+y_discont; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_y_force_discon;
         int I_y_force_discon=1;
         double *p_1y_force_discon;
         int I_1y_force_discon=1;
         double *p_y_discont;
         int I_y_discont=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_force_discon), mccN(&y_force_discon));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_discont), mccN(&y_discont));
         mccAllocateMatrix(&y_force_discon, m_, n_);
         I_y_force_discon = (mccM(&y_force_discon) != 1 || mccN(&y_force_discon) != 1);
         p_y_force_discon = mccPR(&y_force_discon);
         I_1y_force_discon = (mccM(&y_force_discon) != 1 || mccN(&y_force_discon) != 1);
         p_1y_force_discon = mccPR(&y_force_discon);
         I_y_discont = (mccM(&y_discont) != 1 || mccN(&y_discont) != 1);
         p_y_discont = mccPR(&y_discont);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_y_force_discon+=I_y_force_discon, p_1y_force_discon+=I_1y_force_discon, p_y_discont+=I_y_discont)
               {
                  *p_y_force_discon = (*p_1y_force_discon + *p_y_discont);
               }
            }
         }
      }
      
      /* sub_y_sum=sub_y_sum+y_force_discon... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_sub_y_sum;
         int I_sub_y_sum=1;
         double *p_1sub_y_sum;
         int I_1sub_y_sum=1;
         double *p_y_force_discon;
         int I_y_force_discon=1;
         double *p_less_zero_plusL;
         int I_less_zero_plusL=1;
         double *p_less_zero_minusL;
         int I_less_zero_minusL=1;
         double *p_1less_zero_plusL;
         int I_1less_zero_plusL=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_y_sum), mccN(&sub_y_sum));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_force_discon), mccN(&y_force_discon));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&less_zero_plusL), mccN(&less_zero_plusL));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&less_zero_minusL), mccN(&less_zero_minusL));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&less_zero_plusL), mccN(&less_zero_plusL));
         mccAllocateMatrix(&sub_y_sum, m_, n_);
         I_sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_sub_y_sum = mccPR(&sub_y_sum);
         I_1sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_1sub_y_sum = mccPR(&sub_y_sum);
         I_y_force_discon = (mccM(&y_force_discon) != 1 || mccN(&y_force_discon) != 1);
         p_y_force_discon = mccPR(&y_force_discon);
         I_less_zero_plusL = (mccM(&less_zero_plusL) != 1 || mccN(&less_zero_plusL) != 1);
         p_less_zero_plusL = mccPR(&less_zero_plusL);
         I_less_zero_minusL = (mccM(&less_zero_minusL) != 1 || mccN(&less_zero_minusL) != 1);
         p_less_zero_minusL = mccPR(&less_zero_minusL);
         I_1less_zero_plusL = (mccM(&less_zero_plusL) != 1 || mccN(&less_zero_plusL) != 1);
         p_1less_zero_plusL = mccPR(&less_zero_plusL);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_sub_y_sum+=I_sub_y_sum, p_1sub_y_sum+=I_1sub_y_sum, p_y_force_discon+=I_y_force_discon, p_less_zero_plusL+=I_less_zero_plusL, p_less_zero_minusL+=I_less_zero_minusL, p_1less_zero_plusL+=I_1less_zero_plusL)
               {
                  *p_sub_y_sum = (*p_1sub_y_sum + ((*p_y_force_discon * (double) ( ((*p_less_zero_plusL * (double) *p_less_zero_minusL) < 0) && !mccREL_NAN((*p_less_zero_plusL * (double) *p_less_zero_minusL)) )) * (double) *p_1less_zero_plusL));
               }
            }
         }
      }
      
      /* % include epsilon plus factors here */
      
      /* sub_x_sum=-1/2/pi/epsilon*sub_x_sum; */
      R0_ = mcmPi();
      R1_ = (( -1 / (double) 2) / (double) R0_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_epsilon;
         int I_epsilon=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&epsilon), mccN(&epsilon));
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         I_epsilon = (mccM(&epsilon) != 1 || mccN(&epsilon) != 1);
         p_epsilon = mccPR(&epsilon);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_epsilon+=I_epsilon)
               {
                  *p_RM4_ = (R1_ / (double) *p_epsilon);
               }
            }
         }
      }
      mccRealMatrixMultiply(&sub_x_sum, &RM4_, &sub_x_sum);
      /* sub_y_sum=-1/2/pi/epsilon*sub_y_sum; */
      R1_ = mcmPi();
      R0_ = (( -1 / (double) 2) / (double) R1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_epsilon;
         int I_epsilon=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&epsilon), mccN(&epsilon));
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         I_epsilon = (mccM(&epsilon) != 1 || mccN(&epsilon) != 1);
         p_epsilon = mccPR(&epsilon);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_epsilon+=I_epsilon)
               {
                  *p_RM4_ = (R0_ / (double) *p_epsilon);
               }
            }
         }
      }
      mccRealMatrixMultiply(&sub_y_sum, &RM4_, &sub_y_sum);
      
      /* % now we have to rotate back to original position */
      
      /* F_xx=cos_rot_theta.*sub_x_sum-sin_rot_theta... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_F_xx;
         int I_F_xx=1;
         double *p_cos_rot_theta;
         int I_cos_rot_theta=1;
         double *p_sub_x_sum;
         int I_sub_x_sum=1;
         double *p_sin_rot_theta;
         int I_sin_rot_theta=1;
         double *p_sub_y_sum;
         int I_sub_y_sum=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_rot_theta), mccN(&cos_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_x_sum), mccN(&sub_x_sum));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_rot_theta), mccN(&sin_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_y_sum), mccN(&sub_y_sum));
         mccAllocateMatrix(&F_xx, m_, n_);
         I_F_xx = (mccM(&F_xx) != 1 || mccN(&F_xx) != 1);
         p_F_xx = mccPR(&F_xx);
         I_cos_rot_theta = (mccM(&cos_rot_theta) != 1 || mccN(&cos_rot_theta) != 1);
         p_cos_rot_theta = mccPR(&cos_rot_theta);
         I_sub_x_sum = (mccM(&sub_x_sum) != 1 || mccN(&sub_x_sum) != 1);
         p_sub_x_sum = mccPR(&sub_x_sum);
         I_sin_rot_theta = (mccM(&sin_rot_theta) != 1 || mccN(&sin_rot_theta) != 1);
         p_sin_rot_theta = mccPR(&sin_rot_theta);
         I_sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_sub_y_sum = mccPR(&sub_y_sum);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_F_xx+=I_F_xx, p_cos_rot_theta+=I_cos_rot_theta, p_sub_x_sum+=I_sub_x_sum, p_sin_rot_theta+=I_sin_rot_theta, p_sub_y_sum+=I_sub_y_sum)
               {
                  *p_F_xx = ((*p_cos_rot_theta * (double) *p_sub_x_sum) - (*p_sin_rot_theta * (double) *p_sub_y_sum));
               }
            }
         }
      }
      /* F_yy=sin_rot_theta.*sub_x_sum+cos_rot_theta... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_F_yy;
         int I_F_yy=1;
         double *p_sin_rot_theta;
         int I_sin_rot_theta=1;
         double *p_sub_x_sum;
         int I_sub_x_sum=1;
         double *p_cos_rot_theta;
         int I_cos_rot_theta=1;
         double *p_sub_y_sum;
         int I_sub_y_sum=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_rot_theta), mccN(&sin_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_x_sum), mccN(&sub_x_sum));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_rot_theta), mccN(&cos_rot_theta));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sub_y_sum), mccN(&sub_y_sum));
         mccAllocateMatrix(&F_yy, m_, n_);
         I_F_yy = (mccM(&F_yy) != 1 || mccN(&F_yy) != 1);
         p_F_yy = mccPR(&F_yy);
         I_sin_rot_theta = (mccM(&sin_rot_theta) != 1 || mccN(&sin_rot_theta) != 1);
         p_sin_rot_theta = mccPR(&sin_rot_theta);
         I_sub_x_sum = (mccM(&sub_x_sum) != 1 || mccN(&sub_x_sum) != 1);
         p_sub_x_sum = mccPR(&sub_x_sum);
         I_cos_rot_theta = (mccM(&cos_rot_theta) != 1 || mccN(&cos_rot_theta) != 1);
         p_cos_rot_theta = mccPR(&cos_rot_theta);
         I_sub_y_sum = (mccM(&sub_y_sum) != 1 || mccN(&sub_y_sum) != 1);
         p_sub_y_sum = mccPR(&sub_y_sum);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_F_yy+=I_F_yy, p_sin_rot_theta+=I_sin_rot_theta, p_sub_x_sum+=I_sub_x_sum, p_cos_rot_theta+=I_cos_rot_theta, p_sub_y_sum+=I_sub_y_sum)
               {
                  *p_F_yy = ((*p_sin_rot_theta * (double) *p_sub_x_sum) + (*p_cos_rot_theta * (double) *p_sub_y_sum));
               }
            }
         }
      }
      
      /* % multiple by depth to scale in 3rd dimension */
      /* % F_xx is a force matrix from every line segment in seg1 */
      /* % to every line segment in seg2.  When F_xx is multiplied by */
      /* % charge densities from seg1 and seg2, the total force is obtained. */
      /* % ie. total force in x = (charge density for seg1)*F_xx*(charge density */
      /* % for seg2) */
      
      /* F_xx=F_xx*z_axis; */
      mccRealMatrixMultiply(&F_xx, &F_xx, &z_axis);
      /* F_yy=F_yy*z_axis; */
      mccRealMatrixMultiply(&F_yy, &F_yy, &z_axis);
      
      /* % clear up some memory */
      
      /* clear cos_rot_theta sin_rot_theta sub_x_sum ... */
      /* This comment replaces a call on clear */
      
      /* 'calculate potential matrix' */
      mccPrint (&S7_, 0);
      
      /* % now we have to determine the charge density on each line segment. */
      /* % we now the potential on each line segment (initial condition). */
      /* % What we do now is similar to the above force calculation.  We find */
      /* % the potential contribution from each line segment to every other line */
      /* % segment.  The matrices here are even larger (total_line_segments X  */
      /* % total_line_segments). */
      
      /* volt_m = volt_m1;		%some vector as volt_m1 */
      mccCopy(&volt_m, &volt_m1);
      /* size_volt=size(volt_m,1); */
      size_volt = mccGetDimensionSize(&volt_m, 1);
      
      /* % again we have to rotate the frame such that the line segment is  */
      /* % parallel to the X axis and centered at the origin */
      
      /* x_value=volt_m(:,1)*ones(1,size_volt)... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (1-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_volt_m+=I_volt_m)
               {
                  *p_RM4_ = *p_volt_m;
               }
            }
         }
      }
      mccOnesMN(&IM1_, 1, size_volt);
      mccRealMatrixMultiply(&RM3_, &RM4_, &IM1_);
      mccOnesMN(&IM0_, size_volt, 1);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&RM2_, m_, n_);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (1-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_volt_m+=I_volt_m)
               {
                  *p_RM2_ = *p_volt_m;
               }
            }
         }
      }
      mccConjTrans(&RM1_, &RM2_);
      mccRealMatrixMultiply(&RM0_, &IM0_, &RM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_x_value;
         int I_x_value=1;
         double *p_RM3_;
         int I_RM3_=1;
         double *p_RM0_;
         int I_RM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
         mccAllocateMatrix(&x_value, m_, n_);
         I_x_value = (mccM(&x_value) != 1 || mccN(&x_value) != 1);
         p_x_value = mccPR(&x_value);
         I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
         p_RM3_ = mccPR(&RM3_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_x_value+=I_x_value, p_RM3_+=I_RM3_, p_RM0_+=I_RM0_)
               {
                  *p_x_value = (*p_RM3_ - *p_RM0_);
               }
            }
         }
      }
      /* y_value=volt_m(:,2)*ones(1,size_volt)... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (2-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_volt_m+=I_volt_m)
               {
                  *p_RM0_ = *p_volt_m;
               }
            }
         }
      }
      mccOnesMN(&IM0_, 1, size_volt);
      mccRealMatrixMultiply(&RM1_, &RM0_, &IM0_);
      mccOnesMN(&IM1_, size_volt, 1);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&RM2_, m_, n_);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (2-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_volt_m+=I_volt_m)
               {
                  *p_RM2_ = *p_volt_m;
               }
            }
         }
      }
      mccConjTrans(&RM3_, &RM2_);
      mccRealMatrixMultiply(&RM4_, &IM1_, &RM3_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_y_value;
         int I_y_value=1;
         double *p_RM1_;
         int I_RM1_=1;
         double *p_RM4_;
         int I_RM4_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM1_), mccN(&RM1_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM4_), mccN(&RM4_));
         mccAllocateMatrix(&y_value, m_, n_);
         I_y_value = (mccM(&y_value) != 1 || mccN(&y_value) != 1);
         p_y_value = mccPR(&y_value);
         I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
         p_RM1_ = mccPR(&RM1_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_y_value+=I_y_value, p_RM1_+=I_RM1_, p_RM4_+=I_RM4_)
               {
                  *p_y_value = (*p_RM1_ - *p_RM4_);
               }
            }
         }
      }
      /* pot_angle=volt_m(:,3)*ones(1,size_volt); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (3-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_volt_m+=I_volt_m)
               {
                  *p_RM4_ = *p_volt_m;
               }
            }
         }
      }
      mccOnesMN(&IM1_, 1, size_volt);
      mccRealMatrixMultiply(&pot_angle, &RM4_, &IM1_);
      /* sin_angle=sin(pot_angle); */
      mccSin(&sin_angle, &pot_angle);
      /* cos_angle=cos(pot_angle); */
      mccCos(&cos_angle, &pot_angle);
      /* x_pot=x_value.*cos_angle+y_value.*sin_angle; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_x_pot;
         int I_x_pot=1;
         double *p_x_value;
         int I_x_value=1;
         double *p_cos_angle;
         int I_cos_angle=1;
         double *p_y_value;
         int I_y_value=1;
         double *p_sin_angle;
         int I_sin_angle=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_value), mccN(&x_value));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_angle), mccN(&cos_angle));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_value), mccN(&y_value));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_angle), mccN(&sin_angle));
         mccAllocateMatrix(&x_pot, m_, n_);
         I_x_pot = (mccM(&x_pot) != 1 || mccN(&x_pot) != 1);
         p_x_pot = mccPR(&x_pot);
         I_x_value = (mccM(&x_value) != 1 || mccN(&x_value) != 1);
         p_x_value = mccPR(&x_value);
         I_cos_angle = (mccM(&cos_angle) != 1 || mccN(&cos_angle) != 1);
         p_cos_angle = mccPR(&cos_angle);
         I_y_value = (mccM(&y_value) != 1 || mccN(&y_value) != 1);
         p_y_value = mccPR(&y_value);
         I_sin_angle = (mccM(&sin_angle) != 1 || mccN(&sin_angle) != 1);
         p_sin_angle = mccPR(&sin_angle);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_x_pot+=I_x_pot, p_x_value+=I_x_value, p_cos_angle+=I_cos_angle, p_y_value+=I_y_value, p_sin_angle+=I_sin_angle)
               {
                  *p_x_pot = ((*p_x_value * (double) *p_cos_angle) + (*p_y_value * (double) *p_sin_angle));
               }
            }
         }
      }
      /* y_pot=y_value.*cos_angle-x_value.*sin_angle; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_y_pot;
         int I_y_pot=1;
         double *p_y_value;
         int I_y_value=1;
         double *p_cos_angle;
         int I_cos_angle=1;
         double *p_x_value;
         int I_x_value=1;
         double *p_sin_angle;
         int I_sin_angle=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_value), mccN(&y_value));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&cos_angle), mccN(&cos_angle));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_value), mccN(&x_value));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&sin_angle), mccN(&sin_angle));
         mccAllocateMatrix(&y_pot, m_, n_);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         I_y_value = (mccM(&y_value) != 1 || mccN(&y_value) != 1);
         p_y_value = mccPR(&y_value);
         I_cos_angle = (mccM(&cos_angle) != 1 || mccN(&cos_angle) != 1);
         p_cos_angle = mccPR(&cos_angle);
         I_x_value = (mccM(&x_value) != 1 || mccN(&x_value) != 1);
         p_x_value = mccPR(&x_value);
         I_sin_angle = (mccM(&sin_angle) != 1 || mccN(&sin_angle) != 1);
         p_sin_angle = mccPR(&sin_angle);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_y_pot+=I_y_pot, p_y_value+=I_y_value, p_cos_angle+=I_cos_angle, p_x_value+=I_x_value, p_sin_angle+=I_sin_angle)
               {
                  *p_y_pot = ((*p_y_value * (double) *p_cos_angle) - (*p_x_value * (double) *p_sin_angle));
               }
            }
         }
      }
      
      /* x_p=x_pot-delta_approx/2; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_x_p;
         int I_x_p=1;
         double *p_x_pot;
         int I_x_pot=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_pot), mccN(&x_pot));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&x_p, m_, n_);
         I_x_p = (mccM(&x_p) != 1 || mccN(&x_p) != 1);
         p_x_p = mccPR(&x_p);
         I_x_pot = (mccM(&x_pot) != 1 || mccN(&x_pot) != 1);
         p_x_pot = mccPR(&x_pot);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_x_p+=I_x_p, p_x_pot+=I_x_pot, p_delta_approx+=I_delta_approx)
               {
                  *p_x_p = (*p_x_pot - (*p_delta_approx / (double) 2));
               }
            }
         }
      }
      /* x_m=x_pot+delta_approx/2; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_x_m;
         int I_x_m=1;
         double *p_x_pot;
         int I_x_pot=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_pot), mccN(&x_pot));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&x_m, m_, n_);
         I_x_m = (mccM(&x_m) != 1 || mccN(&x_m) != 1);
         p_x_m = mccPR(&x_m);
         I_x_pot = (mccM(&x_pot) != 1 || mccN(&x_pot) != 1);
         p_x_pot = mccPR(&x_pot);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_x_m+=I_x_m, p_x_pot+=I_x_pot, p_delta_approx+=I_delta_approx)
               {
                  *p_x_m = (*p_x_pot + (*p_delta_approx / (double) 2));
               }
            }
         }
      }
      
      /* % below is the actual calculation of the potnential contribution */
      /* % from one segment to another. */
      
      /* s_matrix=(x_m/2.*log(x_m.^2+y_pot.^2+delta)... */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM4_;
         int I_RM4_=1;
         double *p_x_m;
         int I_x_m=1;
         double *p_y_pot;
         int I_y_pot=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_m), mccN(&x_m));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         mccAllocateMatrix(&RM4_, m_, n_);
         I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
         p_RM4_ = mccPR(&RM4_);
         I_x_m = (mccM(&x_m) != 1 || mccN(&x_m) != 1);
         p_x_m = mccPR(&x_m);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_x_m+=I_x_m, p_y_pot+=I_y_pot)
               {
                  *p_RM4_ = ((mcmRealPowerInt(*p_x_m, 2) + mcmRealPowerInt(*p_y_pot, 2)) + delta);
               }
            }
         }
      }
      mccLog(&RM3_, &RM4_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM2_;
         int I_RM2_=1;
         double *p_x_p;
         int I_x_p=1;
         double *p_y_pot;
         int I_y_pot=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_p), mccN(&x_p));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         mccAllocateMatrix(&RM2_, m_, n_);
         I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
         p_RM2_ = mccPR(&RM2_);
         I_x_p = (mccM(&x_p) != 1 || mccN(&x_p) != 1);
         p_x_p = mccPR(&x_p);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_x_p+=I_x_p, p_y_pot+=I_y_pot)
               {
                  *p_RM2_ = ((mcmRealPowerInt(*p_x_p, 2) + mcmRealPowerInt(*p_y_pot, 2)) + delta);
               }
            }
         }
      }
      mccLog(&RM1_, &RM2_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_x_m;
         int I_x_m=1;
         double *p_y_pot;
         int I_y_pot=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_m), mccN(&x_m));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         I_x_m = (mccM(&x_m) != 1 || mccN(&x_m) != 1);
         p_x_m = mccPR(&x_m);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_x_m+=I_x_m, p_y_pot+=I_y_pot)
               {
                  *p_RM0_ = (*p_x_m / (double) (*p_y_pot + delta));
               }
            }
         }
      }
      mccAtan(&RM5_, &RM0_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM6_;
         int I_RM6_=1;
         double *p_x_p;
         int I_x_p=1;
         double *p_y_pot;
         int I_y_pot=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_p), mccN(&x_p));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         mccAllocateMatrix(&RM6_, m_, n_);
         I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
         p_RM6_ = mccPR(&RM6_);
         I_x_p = (mccM(&x_p) != 1 || mccN(&x_p) != 1);
         p_x_p = mccPR(&x_p);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_x_p+=I_x_p, p_y_pot+=I_y_pot)
               {
                  *p_RM6_ = (*p_x_p / (double) (*p_y_pot + delta));
               }
            }
         }
      }
      mccAtan(&RM7_, &RM6_);
      R0_ = mcmPi();
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM8_;
         int I_RM8_=1;
         double *p_x_m;
         int I_x_m=1;
         double *p_RM3_;
         int I_RM3_=1;
         double *p_x_p;
         int I_x_p=1;
         double *p_RM1_;
         int I_RM1_=1;
         double *p_y_pot;
         int I_y_pot=1;
         double *p_RM5_;
         int I_RM5_=1;
         double *p_1y_pot;
         int I_1y_pot=1;
         double *p_RM7_;
         int I_RM7_=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_m), mccN(&x_m));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&x_p), mccN(&x_p));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM1_), mccN(&RM1_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&y_pot), mccN(&y_pot));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&RM7_), mccN(&RM7_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&RM8_, m_, n_);
         I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
         p_RM8_ = mccPR(&RM8_);
         I_x_m = (mccM(&x_m) != 1 || mccN(&x_m) != 1);
         p_x_m = mccPR(&x_m);
         I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
         p_RM3_ = mccPR(&RM3_);
         I_x_p = (mccM(&x_p) != 1 || mccN(&x_p) != 1);
         p_x_p = mccPR(&x_p);
         I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
         p_RM1_ = mccPR(&RM1_);
         I_y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_y_pot = mccPR(&y_pot);
         I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
         p_RM5_ = mccPR(&RM5_);
         I_1y_pot = (mccM(&y_pot) != 1 || mccN(&y_pot) != 1);
         p_1y_pot = mccPR(&y_pot);
         I_RM7_ = (mccM(&RM7_) != 1 || mccN(&RM7_) != 1);
         p_RM7_ = mccPR(&RM7_);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_x_m+=I_x_m, p_RM3_+=I_RM3_, p_x_p+=I_x_p, p_RM1_+=I_RM1_, p_y_pot+=I_y_pot, p_RM5_+=I_RM5_, p_1y_pot+=I_1y_pot, p_RM7_+=I_RM7_, p_delta_approx+=I_delta_approx)
               {
                  *p_RM8_ = ((((((((*p_x_m / (double) 2) * (double) *p_RM3_) - ((*p_x_p / (double) 2) * (double) *p_RM1_)) + (*p_y_pot * (double) *p_RM5_)) - (*p_1y_pot * (double) *p_RM7_)) + *p_delta_approx) / (double) 2) / (double) R0_);
               }
            }
         }
      }
      mccRealRightDivide(&s_matrix, &RM8_, &epsilon);
      
      /* s_matrix=s_matrix'; */
      mccConjTrans(&s_matrix, &s_matrix);
      
      /* % s_plus_matrix adds an extra line.  This essentiall forces the sum */
      /* % of the charge of the entire system to be zero. */
      
      /* s_plus_matrix=[s_matrix ones(size_volt,1)]; */
      mccOnesMN(&IM1_, size_volt, 1);
      mccCatenateColumns(&s_plus_matrix, &s_matrix, &IM1_);
      /* s_plus_matrix=[s_plus_matrix; ones(1,size_volt+1)]; */
      mccOnesMN(&IM1_, 1, (size_volt+1));
      mccCatenateRows(&s_plus_matrix, &s_plus_matrix, &IM1_);
      /* s_plus_matrix(size_volt+1,size_volt+1)=0; */
      mccPR(&s_plus_matrix)[((size_volt+1)-1) + mccM(&s_plus_matrix)*((size_volt+1)-1)] = 0;
      /* voltage=volt_m(:,4); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_voltage;
         int I_voltage=1;
         double *p_volt_m;
         int I_volt_m=1, J_volt_m;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&volt_m), 1);
         mccAllocateMatrix(&voltage, m_, n_);
         I_voltage = (mccM(&voltage) != 1 || mccN(&voltage) != 1);
         p_voltage = mccPR(&voltage);
         if (mccM(&volt_m) == 1) { I_volt_m = J_volt_m = 0;}
         else { I_volt_m = 1; J_volt_m=mccM(&volt_m)-m_; }
         p_volt_m = mccPR(&volt_m) + 0 + mccM(&volt_m) * (4-1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_volt_m += J_volt_m)
            {
               for (i_=0; i_<m_; ++i_, p_voltage+=I_voltage, p_volt_m+=I_volt_m)
               {
                  *p_voltage = *p_volt_m;
               }
            }
         }
      }
      /* voltage(size_volt+1)=0; */
      mccPR(&voltage)[((size_volt+1)-1)] = 0;
      
      /* 'invert potential matrix' */
      mccPrint (&S8_, 0);
      /* % since voltage = s_plus_matrix*sigma_charge */
      /* % we can invert the matrix to find the charge density */
      /* % (sigma_charge) */
      
      /* sigma_charge=inv(s_plus_matrix)*voltage; */
      Mprhs_[0] = &s_plus_matrix;
      Mplhs_[0] = &RM8_;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "inv", 260);
      mccRealMatrixMultiply(&sigma_charge, &RM8_, &voltage);
      
      /* % volt_m=[volt_m sigma_charge(1:size_volt)];   */
      
      /* clear s_plus_matrix s_matrix cos_angle ... */
      /* This comment replaces a call on clear */
      
      /* % remove last line of sigma_charge (doesn't represent charge) */
      /* tot_charge=sigma_charge(1:size(sigma_charge,1)-1); */
      I0_ = mccGetDimensionSize(&sigma_charge, 1);
      mccIntColon2(&IM1_, 1, (I0_-1));
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_tot_charge;
         int I_tot_charge=1;
         double *p_sigma_charge;
         int *p_IM1_;
         int I_IM1_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM1_), mccN(&IM1_), &sigma_charge);
         mccAllocateMatrix(&tot_charge, m_, n_);
         I_tot_charge = (mccM(&tot_charge) != 1 || mccN(&tot_charge) != 1);
         p_tot_charge = mccPR(&tot_charge);
         I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
         p_IM1_ = mccIPR(&IM1_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_tot_charge+=I_tot_charge, p_IM1_+=I_IM1_)
               {
                  *p_tot_charge = mccPR(&sigma_charge)[((int)(*p_IM1_ - .5))];
               }
            }
         }
      }
      
      
      mccReturnFirstValue(&plhs_[0], &F_xx);
      mccReturnValue(&plhs_[1], &F_yy);
      mccReturnValue(&plhs_[2], &tot_charge);
      mccReturnValue(&plhs_[3], &mseg1);
      mccReturnValue(&plhs_[4], &mseg2);
   }
   return;
}
