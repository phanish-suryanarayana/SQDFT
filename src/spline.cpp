/** \file spline.cpp
  \brief This file contains the functions that do 1d cubic spline interpolation.


*/

#include "headers.h"
#include "ds_sq.h"
#include "func_sq.h"

void ispline_gen(double *X1,double *Y1,int len1,double *X2,double *Y2,int len2,double *YD)
{
	int i,j;
	double A0=0.0,A1=0.0,A2=0.0,A3=0.0,x,dx,dy,p1=0.0,p2=0.0,p3=0.0;

	//check error
	if(X2[0]<X1[0] || X2[len2-1]>X1[len1-1])
	{
		printf("%lf, %lf, %lf, %lf\n",X1[0],X2[0],X1[len1-1],X2[len2-1]);
		printf("out of range in spline interpolation\n");
		exit(1);
	}

	// p1 is left endpoint of the interval
	// p2 is resampling position
	// p3 is right endpoint of interval
	// j is input index of current interval

	p3=X2[0]-1;    // force coefficient initialization
	for(i=j=0;i<len2;i++)
	{
		// check if in new interval
		p2=X2[i];
		if(p2>p3){
			// find interval which contains p2 
			for(;j<len1 && p2>X1[j];j++);
			if(p2<X1[j]) j--;
			p1=X1[j];  //update left endpoint 
			p3=X1[j+1]; // update right endpoint

			// compute spline coefficients
			dx = 1.0/(X1[j+1]-X1[j]);
			dy = (Y1[j+1]-Y1[j])*dx;
			A0 = Y1[j];
			A1 = YD[j];
			A2 = dx*(3.0*dy - 2.0*YD[j]-YD[j+1]);
			A3 = dx*dx*(-2.0*dy+YD[j] + YD[j+1]);  
		}
		// use Horner's rule to calculate cubic polynomial
		x =  p2-p1;
		Y2[i] = ((A3*x +A2)*x + A1)*x + A0;
	}

	return;  
}

void getYD_gen(double *X, double *Y, double *YD,int len)
{
	int i;
	double h0,h1,r0,r1,*A,*B,*C;
	A = new double [len](); // need to de-allocate later
	B = new double [len](); // need to de-allocate later
	C = new double [len](); // need to de-allocate later

	// init first row data
	h0 =  X[1]-X[0]; h1 = X[2]-X[1];
	r0 = (Y[1]-Y[0])/h0; r1=(Y[2]-Y[1])/h1;
	B[0] = h1*(h0+h1);
	C[0] = (h0+h1)*(h0+h1);
	YD[0] = r0*(3*h0*h1 + 2*h1*h1) + r1*h0*h0;

	//init tridiagonal bands A,B,C and column vector YD
	// YD will be used to return the derivatives
	for(i=1;i<len-1;i++) {
		h0 = X[i]-X[i-1]; h1=X[i+1]-X[i];
		r0 = (Y[i]-Y[i-1])/h0;  r1=(Y[i+1]-Y[i])/h1;
		A[i] = h1;
		B[i] = 2*(h0+h1);
		C[i] = h0;
		YD[i] = 3*(r0*h1 + r1*h0);
	}

	// last row 
	A[i] = (h0+h1)*(h0+h1);
	B[i] = h0*(h0+h1);
	YD[i] = r0*h1*h1 + r1*(3*h0*h1 + 2*h0*h0);

	// solve for tridiagonal matrix: YD = YD*inv(tridiag matrix)
	tridiag_gen(A,B,C,YD,len);

	delete [] A;
	delete [] B;
	delete [] C;

	return;
}


void tridiag_gen(double *A,double *B,double *C,double *D,int len)
{
	int i;
	double b, *F;

	F = new double [len](); // need to de-allocate later

	// Gauss elimination; forward substitution
	b = B[0];
	D[0] = D[0]/b;
	for(i=1;i<len;i++){
		F[i] = C[i-1]/b;
		b= B[i] - A[i]*F[i];
		if(b==0) {
			printf("Divide by zero in tridiag_gen\n"); 
			exit(1);
		}

		D[i] =(D[i] - D[i-1]*A[i])/b;
	}

	// backsubstitution 
	for(i=len-2;i >= 0;i--)
		D[i] -= (D[i+1]*F[i+1]);

	delete [] F;
	return;
}

