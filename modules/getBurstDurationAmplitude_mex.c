/* 21-08-19     increasing output array sizes to avoid crash in rare instances */
/* 03-07-19     creation */

#include "mex.h"
#include "math.h"

void burstAnalysis(double *envelope, int envSize, double SR, double threshold, 
        double minBurstDuration, double *dur, double *amp)   
{
  int count = 0;
  int nBurst = 0;
  double burstAmp = 0;
  mwSize k;

  
  for (k = 0; k < envSize; k++) 
  {
    if (envelope[k] >= threshold) 
    {
      count++;
      if (count == 1) 
      {
          nBurst++;
      }
      if (envelope[k] >= burstAmp) 
      {
         burstAmp = envelope[k];
      }
    } 
    else 
    {
      if (count > 0) 
      {
        /* If burst just finished */
        if (count / SR >= minBurstDuration) 
        {
            dur[(mwSize) nBurst-1] = count/SR;
            amp[(mwSize) nBurst-1] = burstAmp;
        }
        else
        {
            nBurst--;
        }
        count = 0;
        burstAmp = 0;
        
      }
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
                double *envelope = mxGetPr(prhs[0]);
                int envSize = mxGetScalar(prhs[1]);
                double SR = mxGetScalar(prhs[2]);
                double threshold = mxGetScalar(prhs[3]);   
                double minBurstDuration = mxGetScalar(prhs[4]);  
                
                double *dur;
                double *amp; 
                
				plhs[0] = mxCreateDoubleMatrix(1,(mwSize)round(envSize/2),mxREAL);/*initializes each element in the pr array to 0*/
				dur = mxGetPr(plhs[0]);
				plhs[1] = mxCreateDoubleMatrix(1,(mwSize)round(envSize/2),mxREAL);/*initializes each element in the pr array to 0*/
				amp = mxGetPr(plhs[1]);
				
				burstAnalysis(envelope, envSize, SR, threshold, minBurstDuration, dur, amp);
}

