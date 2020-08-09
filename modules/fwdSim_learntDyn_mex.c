/* 01-09-19     first revision */

#include "mex.h"
#include "math.h"

void fwdSim(double minr, double maxr, double dr, int nr, double *rdot, double sigma, double theta_c, int nMax,
        double dt, double *randr, double *r, double *theta)
{
    /*  Euler method */
    
//     mxArray *Arg[1];
//     Arg[0] = mxCreateString("starting in fwdSim c function");
//     mexCallMATLAB(0, NULL, 1, Arg, "myLogging_c");
//     _sleep(2000);
    
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
        
//         mxArray *Arg[1];
//         Arg[0] = mxCreateString("vField");
//         mexCallMATLAB(0, NULL, 1, Arg, "myLogging_c");
//         
//         char c[50]; //size of the number
//         sprintf(c, "%g", vField);
//         mxArray *Arg2[1];
//         Arg2[0] = mxCreateString(c);
//         mexCallMATLAB(0, NULL, 1, Arg2, "myLogging_c");
//         
//         mxArray *Arg3[1];
//         Arg3[0] = mxCreateString("randr");
//         mexCallMATLAB(0, NULL, 1, Arg3, "myLogging_c");
//         
//         
//         char c2[50]; //size of the number
//         sprintf(c2, "%g", randr[i+1]);
//         mxArray *Arg4[1];
//         Arg4[0] = mxCreateString(c2);
//         mexCallMATLAB(0, NULL, 1, Arg4, "myLogging_c");
//         
//         mxArray *Arg7[1];
//         Arg7[0] = mxCreateString("sigma");
//         mexCallMATLAB(0, NULL, 1, Arg7, "myLogging_c");
//         char c4[50]; //size of the number
//         sprintf(c4, "%g", sigma);
//         mxArray *Arg5[1];
//         Arg5[0] = mxCreateString(c4);
//         mexCallMATLAB(0, NULL, 1, Arg5, "myLogging_c");
//         
//         mxArray *Arg8[1];
//         Arg8[0] = mxCreateString("dt");
//         mexCallMATLAB(0, NULL, 1, Arg8, "myLogging_c");
//         char c5[50]; //size of the number
//         sprintf(c5, "%g", dt);
//         mxArray *Arg9[1];
//         Arg9[0] = mxCreateString(c5);
//         mexCallMATLAB(0, NULL, 1, Arg9, "myLogging_c");
        
        r[i+1] = r[i] + vField*dt + sigma*sqrt(dt)*randr[i];
        
//         mxArray *Arg10[1];
//         Arg10[0] = mxCreateString("r[i+1] before abs");
//         mexCallMATLAB(0, NULL, 1, Arg10, "myLogging_c");
//         char c6[50]; //size of the number
//         sprintf(c6, "%g", r[i+1]);
//         mxArray *Arg11[1];
//         Arg11[0] = mxCreateString(c6);
//         mexCallMATLAB(0, NULL, 1, Arg11, "myLogging_c");
        
                
        r[i+1] = fabs(r[i+1]);
//         float rfabsf = fabs(r[i+1]);
        theta[i+1] = theta[i] + 2*M_PI*theta_c * dt;
        
//         char c3[50]; //size of the number
//         sprintf(c3, "%g", r[i+1]);
//         mxArray *Arg12[1];
//         Arg12[0] = mxCreateString(c3);
//         mexCallMATLAB(0, NULL, 1, Arg12, "myLogging_c");
//         
//         mxArray *Arg13[1];
//         Arg13[0] = mxCreateString("fabsf");
//         mexCallMATLAB(0, NULL, 1, Arg13, "myLogging_c");
//         char c7[50]; //size of the number
//         sprintf(c7, "%g", rfabsf);
//         mxArray *Arg14[1];
//         Arg14[0] = mxCreateString(c7);
//         mexCallMATLAB(0, NULL, 1, Arg14, "myLogging_c");
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
//     mxArray *Arg[1];
//     Arg[0] = mxCreateString("done with fwdSim c function");
//     mexCallMATLAB(0, NULL, 1, Arg, "myLogging_c");
//     _sleep(2000);
}

