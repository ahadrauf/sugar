static char mc_version[] = "MATLAB Compiler 1.2.1 Jan 15 1999 infun";
/*
 *  MATLAB Compiler: 1.2.1
 *  Date: Jan 15 1999
 *  Arguments: -Z -i -r dissect 
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
 * M-File: C:/WINDOWS/Desktop/U C BERKELEY/CAD MEMS/Electro2D by SHollar/electro2D/mex/dissect.m
 *
 *       I0_         	integer scalar temporary
 *       R0_         	real scalar temporary
 *       R1_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       RM2_        	real vector/matrix temporary
 *       angle       	real scalar
 *       array_of_segments	real vector/matrix
 *       atan2       	<function>
 *       ceil        	<function>
 *       delta       	real vector/matrix
 *       dissect     	<function being defined>
 *       distance    	real scalar
 *       i           	integer scalar
 *       line        	real vector/matrix
 *       num_seg     	real vector/matrix
 *       realsqrt    	<function>
 *       x_mean      	real vector/matrix
 *       y_mean      	real vector/matrix
 *******************************************************/

void
mexFunction(
    int nlhs_,
    mxArray *plhs_[],
    int nrhs_,
    const mxArray *prhs_[]
)
{
   

   if (nrhs_ > 2 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }

   if (nlhs_ > 1 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   mcmSetLineNumber(0);
   {
      mxArray array_of_segments;
      mxArray line;
      mxArray delta;
      double distance = 0.0;
      mxArray num_seg;
      double angle = 0.0;
      int i = 0;
      mxArray x_mean;
      mxArray y_mean;
      mxArray RM0_;
      int I0_ = 0;
      double R0_ = 0.0;
      double R1_ = 0.0;
      mxArray RM1_;
      mxArray RM2_;
      
      mccRealInit(line);
      mccImport(&line, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(delta);
      mccImport(&delta, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(array_of_segments);
      mccRealInit(num_seg);
      mccRealInit(x_mean);
      mccRealInit(y_mean);
      mccRealInit(RM0_);
      mccRealInit(RM1_);
      mccRealInit(RM2_);
      
      
      /* % function [array_of_segments] = dissect(line,delta) */
      
      /* % dissect takes in a line of the form specified in electro2d and  */
      /* % dissects it into an array of segments of length delta.   */
      /* % This function is used by electro2d. */
      
      /* % See also electro2d, plot_electro2d */
      
      /* distance = sqrt((line(2)-line(4))^2+(line(1)-line(3))^2); */
      if(mccNOTSET(&line))
      {
         mexErrMsgTxt( "variable line undefined, line 11" );
      }
      if(mccNOTSET(&line))
      {
         mexErrMsgTxt( "variable line undefined, line 11" );
      }
      if(mccNOTSET(&line))
      {
         mexErrMsgTxt( "variable line undefined, line 11" );
      }
      if(mccNOTSET(&line))
      {
         mexErrMsgTxt( "variable line undefined, line 11" );
      }
      distance = sqrt((mcmRealPowerInt(((mccPR(&line)[(2-1)]) - (mccPR(&line)[(4-1)])), 2) + mcmRealPowerInt(((mccPR(&line)[(1-1)]) - (mccPR(&line)[(3-1)])), 2)));
      /* num_seg = ceil(distance/delta); */
      if(mccNOTSET(&delta))
      {
         mexErrMsgTxt( "variable delta undefined, line 12" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_delta;
         int I_delta=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&delta), mccN(&delta));
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
         p_RM0_ = mccPR(&RM0_);
         I_delta = (mccM(&delta) != 1 || mccN(&delta) != 1);
         p_delta = mccPR(&delta);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_delta+=I_delta)
               {
                  *p_RM0_ = (distance / (double) *p_delta);
               }
            }
         }
      }
      mccCeil(&num_seg, &RM0_);
      /* angle = atan2(line(4)-line(2),line(3)-line(1)); */
      angle = atan2(((mccPR(&line)[(4-1)]) - (mccPR(&line)[(2-1)])), ((mccPR(&line)[(3-1)]) - (mccPR(&line)[(1-1)])));
      
      /* for i=1:num_seg */
      for (I0_ = 1; I0_ <= *mccPR(&num_seg); I0_ = I0_ + 1)
      {
         i = I0_;
         /* x_mean=line(1)+(line(3)-line(1))*(i-0.5)*delta/distance; */
         R0_ = (mccPR(&line)[(1-1)]);
         R1_ = (((mccPR(&line)[(3-1)]) - (mccPR(&line)[(1-1)])) * (double) (i - 0.5));
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_x_mean;
            int I_x_mean=1;
            double *p_delta;
            int I_delta=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&delta), mccN(&delta));
            mccAllocateMatrix(&x_mean, m_, n_);
            I_x_mean = (mccM(&x_mean) != 1 || mccN(&x_mean) != 1);
            p_x_mean = mccPR(&x_mean);
            I_delta = (mccM(&delta) != 1 || mccN(&delta) != 1);
            p_delta = mccPR(&delta);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_x_mean+=I_x_mean, p_delta+=I_delta)
                  {
                     *p_x_mean = (R0_ + ((R1_ * (double) *p_delta) / (double) distance));
                  }
               }
            }
         }
         /* y_mean=line(2)+(line(4)-line(2))*(i-0.5)*delta/distance; */
         R1_ = (mccPR(&line)[(2-1)]);
         R0_ = (((mccPR(&line)[(4-1)]) - (mccPR(&line)[(2-1)])) * (double) (i - 0.5));
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_y_mean;
            int I_y_mean=1;
            double *p_delta;
            int I_delta=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&delta), mccN(&delta));
            mccAllocateMatrix(&y_mean, m_, n_);
            I_y_mean = (mccM(&y_mean) != 1 || mccN(&y_mean) != 1);
            p_y_mean = mccPR(&y_mean);
            I_delta = (mccM(&delta) != 1 || mccN(&delta) != 1);
            p_delta = mccPR(&delta);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_y_mean+=I_y_mean, p_delta+=I_delta)
                  {
                     *p_y_mean = (R1_ + ((R0_ * (double) *p_delta) / (double) distance));
                  }
               }
            }
         }
         /* array_of_segments(i,:)=[x_mean y_mean angle line(5)]; */
         mccCatenateColumns(&RM0_, &x_mean, &y_mean);
         mccCatenateColumns(&RM1_, &RM0_, mccTempMatrix(angle, 0., mccREAL, 0 ));
         mccCatenateColumns(&RM2_, &RM1_, mccTempVectorElement(&line, 5));
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_array_of_segments;
            int I_array_of_segments=1, J_array_of_segments;
            double *p_RM2_;
            int I_RM2_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM2_), mccN(&RM2_));
            if (!mccNOTSET(&array_of_segments)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&array_of_segments));
            if (mccN(&array_of_segments) == 1) { I_array_of_segments = J_array_of_segments = 0;}
            else { I_array_of_segments = 1; J_array_of_segments=mccM(&array_of_segments)-m_; }
            p_array_of_segments = mccPR(&array_of_segments) + (i-1) + mccM(&array_of_segments) * 0;
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_array_of_segments += J_array_of_segments)
               {
                  for (i_=0; i_<m_; ++i_, p_array_of_segments+=I_array_of_segments, p_RM2_+=I_RM2_)
                  {
                     *p_array_of_segments = *p_RM2_;
                  }
               }
            }
         }
         mccSTRING(&array_of_segments) = mccSTRING(&RM2_);
         /* end */
      }
      mccReturnFirstValue(&plhs_[0], &array_of_segments);
   }
   return;
}
