/*
Copyright 2016 Sotiris Papatheodorou

Licensed under the Apache License, Version 2.0 (the \"License\");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an \"AS IS\" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdio.h>
#include <math.h>
#include <mex.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
/* input and output arguments */
/* 	E: matrix double 4xp
	A: column vector doulbe 2x1
	B: column vector doulbe 2x1
	uncertA: double
	uncertB: double
	pf: double -> p: int
	max_line: double
*/
#define E_OUT			plhs[0]
#define a_OUT			plhs[1]
#define c_OUT			plhs[2]

#define A_IN			prhs[0]
#define B_IN			prhs[1]
#define uncertA_IN		prhs[2]
#define uncertB_IN		prhs[3]
#define p_IN			prhs[4]
#define max_line_IN		prhs[5]

double *A, *B, *E, uncertA, uncertB, pf, max_line ;

/* get the inputs */
A = mxGetPr(A_IN);
B = mxGetPr(B_IN);
uncertA = mxGetScalar(uncertA_IN);
uncertB = mxGetScalar(uncertB_IN);
pf = mxGetScalar(p_IN);
max_line = mxGetScalar(max_line_IN);






/* FUNCTION BEGINNING ---------------------------------------- */

int i, p;
p = (int) pf;
double *a, *c, x[p], y[p], dx, theta, tmp1, tmp2, fx, fy;

/* hyperbola parameters */
a_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
a = mxGetPr(a_OUT);
c_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
c = mxGetPr(c_OUT);

*a = (uncertA + uncertB) / 2;
*c = (sqrt( (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) )) / 2; /* euclidean distance of A and B */


/* if c <= a the disks overlap and there is no edge-cell */
if ( *c <= *a )
{
	/* create an empty matrix */
	E_OUT = mxCreateDoubleMatrix(0, 0, mxREAL);
	E = mxGetPr(E_OUT);
	return;
}


/* if a = 0 then the edge is the classic Voronoi edge, the perpendicular
bisector of A and B */
if ( *a == 0 )
{
	p = 4;
	x[0] = 0;
	y[0] = max_line;
	x[1] = 0;
	y[1] = -max_line;
	x[2] = max_line;
	y[2] = -max_line;
	x[3] = max_line;
	y[3] = max_line;
}

/* if a > 0 then the edge is a hyperbola branch */
else
{
	dx = (max_line - *a) / ((p-1) / 2);

	for (i=0; i < (p-1)/2; i++)
	{
		/* create x vector for hyperbola */
		x[i] = max_line - i*dx;
		x[p-i-1] = max_line - i*dx;
				
		/* calculate the y of the hyperbola branch */
		y[i] = -sqrt( (*c**c - *a**a) * (x[i]*x[i] / (*a**a) - 1) );
		y[p-i-1] = sqrt( (*c**c - *a**a) * (x[p-i-1]*x[p-i-1] / (*a**a) - 1) );
	}
	x[(p-1)/2] = *a;
	y[(p-1)/2] = 0;
}


E_OUT = mxCreateDoubleMatrix(4, p, mxREAL);
E = mxGetPr(E_OUT);
for (i=0; i < p; i++)
{
	E[0 + 4*i] = -x[i];
	E[1 + 4*i] = y[i];
	E[2 + 4*i] = x[i];
	E[3 + 4*i] = y[i];
}

/* Rotate */
theta = atan2( B[1]-A[1] , B[0]-A[0] );
for (i=0; i < p; i++)
{
	tmp1 = cos(theta)*E[0 + 4*i] - sin(theta)*E[1 + 4*i];
	tmp2 = sin(theta)*E[0 + 4*i] + cos(theta)*E[1 + 4*i];
	E[0 + 4*i] = tmp1;
	E[1 + 4*i] = tmp2;
	
	tmp1 = cos(theta)*E[2 + 4*i] - sin(theta)*E[3 + 4*i];
	tmp2 = sin(theta)*E[2 + 4*i] + cos(theta)*E[3 + 4*i];
	E[2 + 4*i] = tmp1;
	E[3 + 4*i] = tmp2;
}

/* Translate */
fx = (A[0] + B[0]) / 2;
fy = (A[1] + B[1]) / 2;
for (i=0; i < p; i++)
{
	
	E[0 + 4*i] = E[0 + 4*i] + fx;
	E[1 + 4*i] = E[1 + 4*i] + fy;
	E[2 + 4*i] = E[2 + 4*i] + fx;
	E[3 + 4*i] = E[3 + 4*i] + fy;
}

return;

}
