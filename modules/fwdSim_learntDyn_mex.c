/* 01-09-19     first revision */

#include "mex.h"
#include "math.h"

void fwdSim(double minr, double maxr, double dr, int nr, double *rdot, double sigma, double theta_c, int nMax,
        double dt, double *randr, double *r, double *theta)
{
    /*  Euler method */
        
    mwSize i;
    double maxRdot = rdot[(mwSize) nr-1];
    double minRdot = rdot[0];
    for (i = 0; i < nMax-1; i++)
    {
        double vField;
        if (r[i]>=maxr)
        {
            vField = maxRdot;
        }
        else if (r[i]<=minr)
        {
            vField = minRdot;
        }
        else
        {
            double dToRi = r[i] - minr;
            mwSize I = floor((r[i] - minr)/dr);
            double dri = r[i] - (minr + I*dr);
            vField = rdot[I] + dri * (rdot[I+1]-rdot[I])/dr;
        }
                
        r[i+1] = r[i] + vField*dt + sigma*sqrt(dt)*randr[i];          
        r[i+1] = fabs(r[i+1]);
        theta[i+1] = theta[i] + 2*M_PI*theta_c * dt;
        
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double minr = mxGetScalar(prhs[0]);
    double maxr = mxGetScalar(prhs[1]);
    double dr = mxGetScalar(prhs[2]);
    int nr = mxGetScalar(prhs[3]);
    double *rdot = mxGetPr(prhs[4]);
    double sigma = mxGetScalar(prhs[5]);
    double theta_c = mxGetScalar(prhs[6]);
    int nMax =  mxGetScalar(prhs[7]);
    double dt = mxGetScalar(prhs[8]);
    double *randr = mxGetPr(prhs[9]);
    
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)nMax,mxREAL);/*initializes each element in the pr array to 0*/
    double *r = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)nMax,mxREAL);/*initializes each element in the pr array to 0*/
    double *theta = mxGetPr(plhs[1]);
    
    fwdSim(minr, maxr, dr, nr, rdot, sigma, theta_c, nMax, dt, randr, r, theta);
}

