/**************************************************************************
 *
 * $Id: CN2.c,v 1.3 2011/03/26 09:20:41 patrick Exp $
 *
 * Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <stdio.h>
#include <math.h>
#include "mex.h"

void tridag(double *a, double *b, double *c, double *r, double *u, int n)
{
  // Solver for a vector u, the tridiagonal linear set given by
  // A u = r where A is tridiagonal a2...an, b1...bn, c1...cn-1
  double bet;
  double *gam = (double *)mxCalloc(n, sizeof(double) ); // Workspace vector
  char msg1[100];
  char msg2[100];
  int j;

  sprintf(msg1, " Error 1 in tridag");
  sprintf(msg2, " Error 2 in tridag");

  if (b[0] == 0.0)
    mexErrMsgTxt(msg1);

  u[0] = r[0]/(bet=b[0]);
  for (j=1; j<=n-1; j++) {
    // Decomposition and forward substitution
    gam[j] = c[j-1]/bet;
    bet = b[j]-a[j]*gam[j];
    if (bet == 0.0)
      mexErrMsgTxt(msg2);
    u[j] = (r[j]-a[j]*u[j-1])/bet;
  }
  for (j=n-2; j>=0; j--) {
    // Backsubstitution
    u[j] -= gam[j+1]*u[j+1];
  }
  mxFree(gam);
}

void cn2(double dx, double dt, double *u, double *D, double *dk_dv, int n)
{
  // Crank-Nicholson scheme (semi-implicit method)
  // with Dirichlet boundary condition
  int i;
  double *Di = (double *)mxCalloc(n-1, sizeof(double) );
  double *dk_dvi = (double *)mxCalloc(n-1, sizeof(double) );
  double *a = (double *)mxCalloc(n, sizeof(double) );
  double *b = (double *)mxCalloc(n, sizeof(double) );
  double *c = (double *)mxCalloc(n, sizeof(double) );
  double *u1 = (double *)mxCalloc(n, sizeof(double) );
  double alpha = dt/(2.0*dx*dx);

  for (i=0; i<n-1; i++) {
    Di[i] = 0.5*(D[i]+D[i+1]);
    dk_dvi[i] = 0.5*(dk_dv[i]+dk_dv[i+1]);
    Di[i] *= dk_dvi[i];
  }

  a[0]=a[n-1]=0.0;
  b[0]=b[n-1]=1.0;
  c[0]=c[n-1]=0.0;

  u1[0] = u[0];
  u1[n-1] = u[n-1];

  for (i=0; i<n-2; i++) {
    a[i+1] = -dk_dv[i+1]*alpha*Di[i];
    b[i+1] = 1.0+dk_dv[i+1]*alpha*(Di[i]+Di[i+1]);
    c[i+1] = -dk_dv[i+1]*alpha*Di[i+1];
  }

  for (i=0; i<n-2; i++) {
    u1[i+1] = (dk_dv[i+1]*alpha*Di[i]*u[i]+
               (1-dk_dv[i+1]*alpha*(Di[i]+Di[i+1]))*u[i+1]+
               dk_dv[i+1]*alpha*Di[i+1]*u[i+2]);
  }
  tridag(a, b, c, u1, u, n);

  mxFree(Di);
  mxFree(dk_dvi);
  mxFree(a);
  mxFree(b);
  mxFree(c);
  mxFree(u1);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i;
  /* INPUT */
  double dx=*mxGetPr(prhs[0]);
  double dt=*mxGetPr(prhs[1]);

  int n=mxGetM(prhs[2])*mxGetN(prhs[2]);
  double *u=mxGetPr(prhs[2]);
  double *D=mxGetPr(prhs[3]);
  double *dk_dv=mxGetPr(prhs[4]);

  /* OUTPUT */
  double *uout;
  plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
  uout = mxGetPr(plhs[0]);
  for (i=0; i<n; i++) uout[i]=u[i];
  cn2(dx, dt, uout, D, dk_dv, n);
}


