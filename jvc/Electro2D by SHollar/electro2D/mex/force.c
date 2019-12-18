static char mc_version[] = "MATLAB Compiler 1.2.1 Jan 15 1999 infun";
/*
 *  MATLAB Compiler: 1.2.1
 *  Date: Jan 15 1999
 *  Arguments: -Z -i -r force 
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

/***************** Compiler Assumptions ****************
 * M-File: C:/WINDOWS/Desktop/U C BERKELEY/CAD MEMS/Electro2D by SHollar/electro2D/mex/force.m
 *
 *       FF_X        	real vector/matrix
 *       FF_Y        	real vector/matrix
 *       F_xx        	real vector/matrix
 *       F_yy        	real vector/matrix
 *       IM0_        	integer vector/matrix temporary
 *       IM1_        	integer vector/matrix temporary
 *       IM2_        	integer vector/matrix temporary
 *       RM0_        	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       RM2_        	real vector/matrix temporary
 *       charge_1    	real vector/matrix
 *       charge_2    	real vector/matrix
 *       delta_approx	real scalar
 *       force       	<function being defined>
 *       force2d_profile	<function>
 *       load_and_parse	<function>
 *       mseg1       	real vector/matrix
 *       mseg2       	real vector/matrix
 *       ones        	<function>
 *       permitivity 	real scalar
 *       plot_electro2d	<function>
 *       s1          	real vector/matrix
 *       s2          	real vector/matrix
 *       seg1        	real vector/matrix
 *       seg2        	real vector/matrix
 *       size        	<function>
 *       size_seg1   	integer scalar
 *       size_seg2   	integer scalar
 *       tot_charge  	real vector/matrix
 *       z_depth     	integer scalar
 *       zeros       	<function>
 *******************************************************/

void
mexFunction(
    int nlhs_,
    mxArray *plhs_[],
    int nrhs_,
    const mxArray *prhs_[]
)
{
   mxArray *Mplhs_[5];
   mxArray *Mprhs_[5];
   

   if (nrhs_ > 2 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }

   if (nlhs_ > 0 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   mcmSetLineNumber(0);
   {
      mxArray s1;
      mxArray s2;
      mxArray seg1;
      mxArray seg2;
      double delta_approx = 0.0;
      double permitivity = 0.0;
      int z_depth = 0;
      mxArray F_xx;
      mxArray F_yy;
      mxArray tot_charge;
      mxArray mseg1;
      mxArray mseg2;
      int size_seg1 = 0;
      int size_seg2 = 0;
      mxArray charge_1;
      mxArray charge_2;
      mxArray FF_X;
      mxArray FF_Y;
      mxArray IM0_;
      mxArray IM1_;
      mxArray IM2_;
      mxArray RM0_;
      mxArray RM1_;
      mxArray RM2_;
      
      mccRealInit(s1);
      mccImport(&s1, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(s2);
      mccImport(&s2, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(seg1);
      mccRealInit(seg2);
      mccRealInit(F_xx);
      mccRealInit(F_yy);
      mccRealInit(tot_charge);
      mccRealInit(mseg1);
      mccRealInit(mseg2);
      mccRealInit(charge_1);
      mccRealInit(charge_2);
      mccRealInit(FF_X);
      mccRealInit(FF_Y);
      mccIntInit(IM0_);
      mccIntInit(IM1_);
      mccIntInit(IM2_);
      mccRealInit(RM0_);
      mccRealInit(RM1_);
      mccRealInit(RM2_);
      
      
      /* seg1 = load_and_parse(s1); */
      if(mccNOTSET(&s1))
      {
         mexErrMsgTxt( "variable s1 undefined, line 6" );
      }
      Mprhs_[0] = &s1;
      Mplhs_[0] = &seg1;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "load_and_parse", 6);
      /* seg2 = load_and_parse(s2); */
      if(mccNOTSET(&s2))
      {
         mexErrMsgTxt( "variable s2 undefined, line 7" );
      }
      Mprhs_[0] = &s2;
      Mplhs_[0] = &seg2;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "load_and_parse", 7);
      
      /* delta_approx = 0.125; 		%number of microns that the segments get divided. */
      delta_approx = 0.125;
      /* permitivity = 8.854e-6;	%mks units * 1e-6 (avoids truncation error). */
      permitivity = 8.854e-6;
      /* z_depth = 200; 				%microns, thickness. */
      z_depth = 200;
      
      /* [F_xx, F_yy, tot_charge, mseg1, mseg2] = force2d_profile(seg1,seg2,... */
      Mprhs_[0] = &seg1;
      Mprhs_[1] = &seg2;
      Mprhs_[2] = mccTempMatrix(delta_approx, 0., mccREAL, 0 );
      Mprhs_[3] = mccTempMatrix(permitivity, 0., mccREAL, 0 );
      Mprhs_[4] = mccTempMatrix(z_depth, 0., mccINT, 0 );
      Mplhs_[0] = &F_xx;
      Mplhs_[1] = &F_yy;
      Mplhs_[2] = &tot_charge;
      Mplhs_[3] = &mseg1;
      Mplhs_[4] = &mseg2;
      mccCallMATLAB(5, Mplhs_, 5, Mprhs_, "force2d_profile", 14);
      
      /* size_seg1 = size(mseg1,1); */
      size_seg1 = mccGetDimensionSize(&mseg1, 1);
      /* size_seg2 = size(mseg2,1); */
      size_seg2 = mccGetDimensionSize(&mseg2, 1);
      
      /* % charge_1 is the charge vector for seg1 */
      /* charge_1=[ones(size_seg1,1) ; zeros(size_seg2,1)].*tot_charge; */
      mccOnesMN(&IM0_, size_seg1, 1);
      mccZerosMN(&IM1_, size_seg2, 1);
      mccCatenateRows(&IM2_, &IM0_, &IM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_charge_1;
         int I_charge_1=1;
         int *p_IM2_;
         int I_IM2_=1;
         double *p_tot_charge;
         int I_tot_charge=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&IM2_), mccN(&IM2_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&tot_charge), mccN(&tot_charge));
         mccAllocateMatrix(&charge_1, m_, n_);
         I_charge_1 = (mccM(&charge_1) != 1 || mccN(&charge_1) != 1);
         p_charge_1 = mccPR(&charge_1);
         I_IM2_ = (mccM(&IM2_) != 1 || mccN(&IM2_) != 1);
         p_IM2_ = mccIPR(&IM2_);
         I_tot_charge = (mccM(&tot_charge) != 1 || mccN(&tot_charge) != 1);
         p_tot_charge = mccPR(&tot_charge);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_charge_1+=I_charge_1, p_IM2_+=I_IM2_, p_tot_charge+=I_tot_charge)
               {
                  *p_charge_1 = (((int)*p_IM2_) * (double) *p_tot_charge);
               }
            }
         }
      }
      /* % charge_2 is the charge vector for seg2 */
      /* charge_2=[zeros(size_seg1,1) ; ones(size_seg2,1)].*tot_charge; */
      mccZerosMN(&IM2_, size_seg1, 1);
      mccOnesMN(&IM1_, size_seg2, 1);
      mccCatenateRows(&IM0_, &IM2_, &IM1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_charge_2;
         int I_charge_2=1;
         int *p_IM0_;
         int I_IM0_=1;
         double *p_tot_charge;
         int I_tot_charge=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&tot_charge), mccN(&tot_charge));
         mccAllocateMatrix(&charge_2, m_, n_);
         I_charge_2 = (mccM(&charge_2) != 1 || mccN(&charge_2) != 1);
         p_charge_2 = mccPR(&charge_2);
         I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
         p_IM0_ = mccIPR(&IM0_);
         I_tot_charge = (mccM(&tot_charge) != 1 || mccN(&tot_charge) != 1);
         p_tot_charge = mccPR(&tot_charge);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_charge_2+=I_charge_2, p_IM0_+=I_IM0_, p_tot_charge+=I_tot_charge)
               {
                  *p_charge_2 = (((int)*p_IM0_) * (double) *p_tot_charge);
               }
            }
         }
      }
      
      /* % calculate total force (sum of force contributions from segment1 */
      /* % onto segment2 */
      /* % value is in micro newtons */
      /* FF_X=charge_2'*F_xx'*charge_1 */
      mccConjTrans(&RM0_, &charge_2);
      mccConjTrans(&RM1_, &F_xx);
      mccRealMatrixMultiply(&RM2_, &RM0_, &RM1_);
      mccRealMatrixMultiply(&FF_X, &RM2_, &charge_1);
      mccPrint (&FF_X, "FF_X");
      /* FF_Y=charge_2'*F_yy'*charge_1	 */
      mccConjTrans(&RM2_, &charge_2);
      mccConjTrans(&RM1_, &F_yy);
      mccRealMatrixMultiply(&RM0_, &RM2_, &RM1_);
      mccRealMatrixMultiply(&FF_Y, &RM0_, &charge_1);
      mccPrint (&FF_Y, "FF_Y");
      
      /* plot_electro2d(seg1,seg2, delta_approx); */
      Mprhs_[0] = &seg1;
      Mprhs_[1] = &seg2;
      Mprhs_[2] = mccTempMatrix(delta_approx, 0., mccREAL, 0 );
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 3, Mprhs_, "plot_electro2d", 30);
   }
   return;
}
