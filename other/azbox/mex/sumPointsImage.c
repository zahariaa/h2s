/* 
IM = sumPointsImage(allx,ally,x,y,r)

2013-03-21 AZ Created

MATLAB code:

[XI,YI] = meshgrid(allx,ally);
for i = 1:nx
   % Elliptical bound based on axis limits, meant for plotting circles
   IM = IM + single(((XI - x(i)).^2)/r(1)^2 + ((YI - y(i)).^2)/r(2)^2 < 1);
end
*/

#include <matrix.h>  /* Matlab matrices */
#include <mex.h>
#include <math.h>

void mexFunction(int            nlhs,	     /* Num return vals on lhs */
		           mxArray       *plhs[],     /* Matrices on lhs        */
		           int            nrhs,	     /* Num args on rhs        */
		           const mxArray *prhs[]      /* Matrices on rhs        */
		 )
  {
  const mxArray *arg0, *arg1, *arg2, *arg3, *arg4; /*    input pointers  */
  const mxArray *out0;                             /*   output pointers  */
  double *allx, *ally, *x, *y, *r;                 /*    input variables */
  double *IM;                                      /*   output variables */
  int nx, xdim, ydim, i, ix, iy, localval;         /* internal variables */
  
  /* Associate  input pointers */
  arg0 = mxDuplicateArray(prhs[0]);
  arg1 = mxDuplicateArray(prhs[1]);
  arg2 = mxDuplicateArray(prhs[2]);
  arg3 = mxDuplicateArray(prhs[3]);
  arg4 = mxDuplicateArray(prhs[4]);
  
  /* Get relevant dimensions of inputs */
  xdim = mxGetNumberOfElements(arg0);
  ydim = mxGetNumberOfElements(arg1);
  nx   = mxGetNumberOfElements(arg2);
  
  /* Associate output pointers */
  out0 = plhs[0] = mxCreateDoubleMatrix(ydim,xdim,mxREAL);  /* transposes matrix! */
  
  /* Associate pointers with internal variables */
  allx = mxGetPr(arg0);
  ally = mxGetPr(arg1);
  x    = mxGetPr(arg2);
  y    = mxGetPr(arg3);
  r    = mxGetPr(arg4);
  IM   = mxGetPr(out0);
  
  localval = 0;
  for (i=0; i<nx; i++)
  {
     for (iy=0; iy<ydim; iy++)
     {
        for (ix=0; ix<xdim; ix++)
        {
           localval  = pow((allx[ix] - x[i]),2)/pow(r[0],2) + pow((ally[iy] - y[i]),2)/pow(r[1],2);
           
           if (localval < 1)
           {
              IM[ix*ydim+iy]++;
           }
        }
     }
  }

  return;
  }