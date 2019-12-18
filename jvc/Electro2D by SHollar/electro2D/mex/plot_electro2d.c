static char mc_version[] = "MATLAB Compiler 1.2.1 Jan 15 1999 infun";
/*
 *  MATLAB Compiler: 1.2.1
 *  Date: Jan 15 1999
 *  Arguments: -Z -i -r plot_electro2d 
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

/* static array S0_ (1 x 4) int, line 17 */
static int S0r_[] =
{
	0, 0, 0, 0, 
};
static mxArray S0_ = mccCINIT( mccINT, 1, 4, S0r_, 0 );
/* static array S1_ (1 x 4) int, line 23 */
static int S1r_[] =
{
	0, 0, 0, 0, 
	
};
static mxArray S1_ = mccCINIT( mccINT, 1, 4, S1r_, 0 );
/* static array S2_ (1 x 2) int, line 45 */
static int S2r_[] =
{
	1, -1, 
};
static mxArray S2_ = mccCINIT( mccINT, 1, 2, S2r_, 0 );
/* static array S3_ (1 x 2) int, line 45 */
static int S3r_[] =
{
	1, -1, 
};
static mxArray S3_ = mccCINIT( mccINT, 1, 2, S3r_, 0 );
/* static array S4_ (1 x 2) int, line 55 */
static int S4r_[] =
{
	1, -1, 
};
static mxArray S4_ = mccCINIT( mccINT, 1, 2, S4r_, 0 );
/* static array S5_ (1 x 2) int, line 55 */
static int S5r_[] =
{
	1, -1, 
	
};
static mxArray S5_ = mccCINIT( mccINT, 1, 2, S5r_, 0 );
/***************** Compiler Assumptions ****************
 * M-File: C:/WINDOWS/Desktop/U C BERKELEY/CAD MEMS/Electro2D by SHollar/electro2D/mex/plot_electro2d.m
 *
 *       I0_         	integer scalar temporary
 *       I1_         	integer scalar temporary
 *       IM0_        	integer vector/matrix temporary
 *       IM1_        	integer vector/matrix temporary
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
 *       S0_         	<constant>
 *       S1_         	<constant>
 *       S2_         	<constant>
 *       S3_         	<constant>
 *       S4_         	<constant>
 *       S5_         	<constant>
 *       adding_matrix	real vector/matrix
 *       clf         	<function>
 *       cos         	<function>
 *       delta_approx	real vector/matrix
 *       dissect     	<function>
 *       length      	real vector/matrix
 *       line        	<function>
 *       n           	integer vector/matrix
 *       n           	integer scalar  => n_1
 *       plot_electro2d	<function being defined>
 *       plot_matrix 	real vector/matrix
 *       seg1        	real vector/matrix
 *       seg2        	real vector/matrix
 *       sin         	<function>
 *       size        	<function>
 *       size_volt1  	integer scalar
 *       size_volt2  	integer scalar
 *       theta       	real vector/matrix
 *       theta       	real scalar  => theta_1
 *       total_segment_number	integer scalar
 *       volt_m1     	real vector/matrix
 *       volt_m2     	real vector/matrix
 *       x_c         	real vector/matrix
 *       x_c         	real scalar  => x_c_1
 *       y_c         	real vector/matrix
 *       y_c         	real scalar  => y_c_1
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
   

   if (nrhs_ > 3 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }

   if (nlhs_ > 0 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   mcmSetLineNumber(0);
   {
      mxArray seg1;
      mxArray seg2;
      mxArray delta_approx;
      mxArray adding_matrix;
      mxArray n;
      mxArray volt_m1;
      mxArray volt_m2;
      int size_volt1 = 0;
      int size_volt2 = 0;
      int total_segment_number = 0;
      mxArray length;
      mxArray plot_matrix;
      mxArray theta;
      mxArray x_c;
      mxArray y_c;
      int n_1 = 0;
      double theta_1 = 0.0;
      double x_c_1 = 0.0;
      double y_c_1 = 0.0;
      int I0_ = 0;
      int I1_ = 0;
      mxArray IM0_;
      mxArray IM1_;
      mxArray RM0_;
      mxArray RM1_;
      mxArray RM2_;
      mxArray RM3_;
      mxArray RM4_;
      mxArray RM5_;
      mxArray RM6_;
      mxArray RM7_;
      double R0_ = 0.0;
      double R1_ = 0.0;
      
      mccRealInit(seg1);
      mccImport(&seg1, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(seg2);
      mccImport(&seg2, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(delta_approx);
      mccImport(&delta_approx, ((nrhs_>2) ? prhs_[2] : 0), 0, 0);
      mccRealInit(adding_matrix);
      mccIntInit(n);
      mccRealInit(volt_m1);
      mccRealInit(volt_m2);
      mccRealInit(length);
      mccRealInit(plot_matrix);
      mccRealInit(theta);
      mccRealInit(x_c);
      mccRealInit(y_c);
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
      
      
      /* % function plot_electro2d(seg1,seg2,delta_approx); */
      
      /* % plot_electro2d plots the approximation of seg1  */
      /* % and seg2.  Each line segment in seg1 and seg2 is */
      /* % approximated as an array of line segments of  */
      /* % length delta_approx.  */
      
      /* % see also electro2d */
      
      
      /* % clear original plot figure */
      /* clf; */
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 0, Mprhs_, "clf", 14);
      
      /* % break up segments into smaller segments of length delta_approx */
      /* adding_matrix=[0 0 0 0]; */
      mccCopy(&adding_matrix, &S0_);
      /* for	n=1:size(seg1,1) */
      if(mccNOTSET(&seg1))
      {
         mexErrMsgTxt( "variable seg1 undefined, line 18" );
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
            mexErrMsgTxt( "variable delta_approx undefined, line 19" );
         }
         Mprhs_[0] = &RM0_;
         Mprhs_[1] = &delta_approx;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "dissect", 19);
         mccCatenateRows(&adding_matrix, &adding_matrix, &RM1_);
         /* end */
      }
      /* volt_m1=adding_matrix(2:size(adding_matrix,1),:); */
      I1_ = mccGetDimensionSize(&adding_matrix, 1);
      mccIntColon2(&IM0_, 2, I1_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_volt_m1;
         int I_volt_m1=1;
         double *p_adding_matrix;
         int I_adding_matrix=1, J_adding_matrix;
         int *p_IM0_;
         int I_IM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM0_) * mccN(&IM0_)), mccN(&adding_matrix));
         mccAllocateMatrix(&volt_m1, m_, n_);
         I_volt_m1 = (mccM(&volt_m1) != 1 || mccN(&volt_m1) != 1);
         p_volt_m1 = mccPR(&volt_m1);
         J_adding_matrix = ((mccM(&adding_matrix) != 1 || mccN(&adding_matrix) != 1) ? mccM(&adding_matrix) : 0);
         p_adding_matrix = mccPR(&adding_matrix) + mccM(&adding_matrix) * 0;
         I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_adding_matrix += J_adding_matrix)
            {
               p_IM0_ = mccIPR(&IM0_);
               for (i_=0; i_<m_; ++i_, p_volt_m1+=I_volt_m1, p_IM0_+=I_IM0_)
               {
                  *p_volt_m1 = p_adding_matrix[((int)(*p_IM0_ - .5))];
               }
            }
         }
      }
      
      /* adding_matrix=[0 0 0 0]; */
      mccCopy(&adding_matrix, &S1_);
      /* for	n=1:size(seg2,1) */
      if(mccNOTSET(&seg2))
      {
         mexErrMsgTxt( "variable seg2 undefined, line 24" );
      }
      I0_ = mccGetDimensionSize(&seg2, 1);
      mccIntColon2(&IM0_, 1, I0_);
      mccCreateEmpty( &n );
      for (I1_ = 0; I1_ < mccN(&IM0_); ++I1_)
      {
         mccForCol(&n, &IM0_, I1_);
         /* adding_matrix=[adding_matrix;dissect(seg2(n,:),delta_approx)]; */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM1_;
            int I_RM1_=1;
            double *p_seg2;
            int I_seg2=1, J_seg2;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM1_) * mccN(&IM1_)), mccN(&seg2));
            mccAllocateMatrix(&RM1_, m_, n_);
            I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
            p_RM1_ = mccPR(&RM1_);
            J_seg2 = ((mccM(&seg2) != 1 || mccN(&seg2) != 1) ? mccM(&seg2) : 0);
            p_seg2 = mccPR(&seg2) + mccM(&seg2) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_seg2 += J_seg2)
               {
                  p_IM1_ = mccIPR(&IM1_);
                  for (i_=0; i_<m_; ++i_, p_RM1_+=I_RM1_, p_IM1_+=I_IM1_)
                  {
                     *p_RM1_ = p_seg2[((int)(*p_IM1_ - .5))];
                  }
               }
            }
         }
         if(mccNOTSET(&delta_approx))
         {
            mexErrMsgTxt( "variable delta_approx undefined, line 25" );
         }
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &delta_approx;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "dissect", 25);
         mccCatenateRows(&adding_matrix, &adding_matrix, &RM0_);
         /* end */
      }
      /* volt_m2=adding_matrix(2:size(adding_matrix,1),:); */
      I0_ = mccGetDimensionSize(&adding_matrix, 1);
      mccIntColon2(&IM0_, 2, I0_);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_volt_m2;
         int I_volt_m2=1;
         double *p_adding_matrix;
         int I_adding_matrix=1, J_adding_matrix;
         int *p_IM0_;
         int I_IM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM0_) * mccN(&IM0_)), mccN(&adding_matrix));
         mccAllocateMatrix(&volt_m2, m_, n_);
         I_volt_m2 = (mccM(&volt_m2) != 1 || mccN(&volt_m2) != 1);
         p_volt_m2 = mccPR(&volt_m2);
         J_adding_matrix = ((mccM(&adding_matrix) != 1 || mccN(&adding_matrix) != 1) ? mccM(&adding_matrix) : 0);
         p_adding_matrix = mccPR(&adding_matrix) + mccM(&adding_matrix) * 0;
         I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_adding_matrix += J_adding_matrix)
            {
               p_IM0_ = mccIPR(&IM0_);
               for (i_=0; i_<m_; ++i_, p_volt_m2+=I_volt_m2, p_IM0_+=I_IM0_)
               {
                  *p_volt_m2 = p_adding_matrix[((int)(*p_IM0_ - .5))];
               }
            }
         }
      }
      
      /* size_volt1 = size(volt_m1,1); */
      size_volt1 = mccGetDimensionSize(&volt_m1, 1);
      /* size_volt2 = size(volt_m2,1); */
      size_volt2 = mccGetDimensionSize(&volt_m2, 1);
      
      /* % total number of segments of length delta_approx */
      /* total_segment_number=size_volt1+size_volt2 */
      total_segment_number = (size_volt1 + size_volt2);
      mccPrint (mccTempMatrix(total_segment_number, 0., mccINT, 0 ), "total_segment_number");
      
      
      /* % plot all lines */
      /* length=delta_approx/2; */
      if(mccNOTSET(&delta_approx))
      {
         mexErrMsgTxt( "variable delta_approx undefined, line 37" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_length;
         int I_length=1;
         double *p_delta_approx;
         int I_delta_approx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta_approx), mccN(&delta_approx));
         mccAllocateMatrix(&length, m_, n_);
         I_length = (mccM(&length) != 1 || mccN(&length) != 1);
         p_length = mccPR(&length);
         I_delta_approx = (mccM(&delta_approx) != 1 || mccN(&delta_approx) != 1);
         p_delta_approx = mccPR(&delta_approx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_length+=I_length, p_delta_approx+=I_delta_approx)
               {
                  *p_length = (*p_delta_approx / (double) 2);
               }
            }
         }
      }
      /* plot_matrix=volt_m1; */
      mccCopy(&plot_matrix, &volt_m1);
      
      /* for n=1:size(plot_matrix,1); */
      I1_ = mccGetDimensionSize(&plot_matrix, 1);
      mccIntColon2(&IM0_, 1, I1_);
      mccCreateEmpty( &n );
      for (I0_ = 0; I0_ < mccN(&IM0_); ++I0_)
      {
         mccForCol(&n, &IM0_, I0_);
         /* theta=plot_matrix(n,3); */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_theta;
            int I_theta=1;
            double *p_plot_matrix;
            int I_plot_matrix=1, J_plot_matrix;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM1_) * mccN(&IM1_)), 1);
            mccAllocateMatrix(&theta, m_, n_);
            I_theta = (mccM(&theta) != 1 || mccN(&theta) != 1);
            p_theta = mccPR(&theta);
            J_plot_matrix = ((mccM(&plot_matrix) != 1 || mccN(&plot_matrix) != 1) ? mccM(&plot_matrix) : 0);
            p_plot_matrix = mccPR(&plot_matrix) + mccM(&plot_matrix) * (3-1);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_plot_matrix += J_plot_matrix)
               {
                  p_IM1_ = mccIPR(&IM1_);
                  for (i_=0; i_<m_; ++i_, p_theta+=I_theta, p_IM1_+=I_IM1_)
                  {
                     *p_theta = p_plot_matrix[((int)(*p_IM1_ - .5))];
                  }
               }
            }
         }
         /* x_c=plot_matrix(n,1); */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_x_c;
            int I_x_c=1;
            double *p_plot_matrix;
            int I_plot_matrix=1, J_plot_matrix;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM1_) * mccN(&IM1_)), 1);
            mccAllocateMatrix(&x_c, m_, n_);
            I_x_c = (mccM(&x_c) != 1 || mccN(&x_c) != 1);
            p_x_c = mccPR(&x_c);
            J_plot_matrix = ((mccM(&plot_matrix) != 1 || mccN(&plot_matrix) != 1) ? mccM(&plot_matrix) : 0);
            p_plot_matrix = mccPR(&plot_matrix) + mccM(&plot_matrix) * (1-1);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_plot_matrix += J_plot_matrix)
               {
                  p_IM1_ = mccIPR(&IM1_);
                  for (i_=0; i_<m_; ++i_, p_x_c+=I_x_c, p_IM1_+=I_IM1_)
                  {
                     *p_x_c = p_plot_matrix[((int)(*p_IM1_ - .5))];
                  }
               }
            }
         }
         /* y_c= plot_matrix(n,2); */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_y_c;
            int I_y_c=1;
            double *p_plot_matrix;
            int I_plot_matrix=1, J_plot_matrix;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM1_) * mccN(&IM1_)), 1);
            mccAllocateMatrix(&y_c, m_, n_);
            I_y_c = (mccM(&y_c) != 1 || mccN(&y_c) != 1);
            p_y_c = mccPR(&y_c);
            J_plot_matrix = ((mccM(&plot_matrix) != 1 || mccN(&plot_matrix) != 1) ? mccM(&plot_matrix) : 0);
            p_plot_matrix = mccPR(&plot_matrix) + mccM(&plot_matrix) * (2-1);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_plot_matrix += J_plot_matrix)
               {
                  p_IM1_ = mccIPR(&IM1_);
                  for (i_=0; i_<m_; ++i_, p_y_c+=I_y_c, p_IM1_+=I_IM1_)
                  {
                     *p_y_c = p_plot_matrix[((int)(*p_IM1_ - .5))];
                  }
               }
            }
         }
         /* line(length*cos(theta)*[1 -1] + x_c, ... */
         mccCos(&RM0_, &theta);
         mccRealMatrixMultiply(&RM1_, &length, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &S2_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_x_c;
            int I_x_c=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM2_), mccN(&RM2_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&x_c), mccN(&x_c));
            mccAllocateMatrix(&RM3_, m_, n_);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_x_c = (mccM(&x_c) != 1 || mccN(&x_c) != 1);
            p_x_c = mccPR(&x_c);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_RM2_+=I_RM2_, p_x_c+=I_x_c)
                  {
                     *p_RM3_ = (*p_RM2_ + *p_x_c);
                  }
               }
            }
         }
         mccSTRING(&RM3_) = 0;
         mccSin(&RM4_, &theta);
         mccRealMatrixMultiply(&RM5_, &length, &RM4_);
         mccRealMatrixMultiply(&RM6_, &RM5_, &S3_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM7_;
            int I_RM7_=1;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_y_c;
            int I_y_c=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM6_), mccN(&RM6_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&y_c), mccN(&y_c));
            mccAllocateMatrix(&RM7_, m_, n_);
            I_RM7_ = (mccM(&RM7_) != 1 || mccN(&RM7_) != 1);
            p_RM7_ = mccPR(&RM7_);
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_y_c = (mccM(&y_c) != 1 || mccN(&y_c) != 1);
            p_y_c = mccPR(&y_c);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM7_+=I_RM7_, p_RM6_+=I_RM6_, p_y_c+=I_y_c)
                  {
                     *p_RM7_ = (*p_RM6_ + *p_y_c);
                  }
               }
            }
         }
         mccLOG(&RM7_) = 0;
         mccSTRING(&RM7_) = 0;
         Mprhs_[0] = &RM3_;
         Mprhs_[1] = &RM7_;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 2, Mprhs_, "line", 45);
         /* end */
      }
      
      /* plot_matrix=volt_m2; */
      mccCopy(&plot_matrix, &volt_m2);
      
      /* for n=1:size(plot_matrix,1); */
      I0_ = mccGetDimensionSize(&plot_matrix, 1);
      for (I1_ = 1; I1_ <= I0_; I1_ = I1_ + 1)
      {
         n_1 = I1_;
         /* theta=plot_matrix(n,3); */
         theta_1 = (mccPR(&plot_matrix)[(3-1)*mccM(&plot_matrix) + (n_1-1)]);
         /* x_c=plot_matrix(n,1); */
         x_c_1 = (mccPR(&plot_matrix)[(1-1)*mccM(&plot_matrix) + (n_1-1)]);
         /* y_c= plot_matrix(n,2); */
         y_c_1 = (mccPR(&plot_matrix)[(2-1)*mccM(&plot_matrix) + (n_1-1)]);
         /* line(length*cos(theta)*[1 -1] + x_c, ... */
         R0_ = cos(theta_1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM7_;
            int I_RM7_=1;
            double *p_length;
            int I_length=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&length), mccN(&length));
            mccAllocateMatrix(&RM7_, m_, n_);
            I_RM7_ = (mccM(&RM7_) != 1 || mccN(&RM7_) != 1);
            p_RM7_ = mccPR(&RM7_);
            I_length = (mccM(&length) != 1 || mccN(&length) != 1);
            p_length = mccPR(&length);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM7_+=I_RM7_, p_length+=I_length)
                  {
                     *p_RM7_ = (*p_length * (double) R0_);
                  }
               }
            }
         }
         mccRealMatrixMultiply(&RM6_, &RM7_, &S4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_RM6_;
            int I_RM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM6_), mccN(&RM6_));
            mccAllocateMatrix(&RM5_, m_, n_);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM5_+=I_RM5_, p_RM6_+=I_RM6_)
                  {
                     *p_RM5_ = (*p_RM6_ + x_c_1);
                  }
               }
            }
         }
         mccSTRING(&RM5_) = 0;
         R1_ = sin(theta_1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_length;
            int I_length=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&length), mccN(&length));
            mccAllocateMatrix(&RM4_, m_, n_);
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            I_length = (mccM(&length) != 1 || mccN(&length) != 1);
            p_length = mccPR(&length);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_length+=I_length)
                  {
                     *p_RM4_ = (*p_length * (double) R1_);
                  }
               }
            }
         }
         mccRealMatrixMultiply(&RM3_, &RM4_, &S5_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = (*p_RM3_ + y_c_1);
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM5_;
         Mprhs_[1] = &RM2_;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 2, Mprhs_, "line", 55);
         /* end */
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
   }
   return;
}
