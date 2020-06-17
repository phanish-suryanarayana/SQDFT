/** \file anderson.cpp
  \brief This file contains the functions needed to find out the Anderson mixing update vector.


*/
#include "headers.h"
#include "ds_sq.h"
#include "func_sq.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef TEMP_TOL
#define TEMP_TOL 1e-12
#endif

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// function to compute Anderson update vector
void AndersonExtrapolation(double **DX, double **DF, double *phi_res_vec, double beta_mix, int anderson_history, int N, double *am_vec, double *mpi_time)
{
	// DX and DF are iterate and residual history matrices of size Nxm where N=no. of nodes in proc domain and m=anderson_history
	// am_vec is the final update vector of size Nx1
	// phi_res_vec is the residual vector at current iteration of fixed point method

	int i,j,k,ctr;
	int m=anderson_history;
	double tcomm1,tcomm2;

	// ------------------- First find DF'*DF m x m matrix (local and then AllReduce) ------------------- //
	double temp_sum=0;
	double *FtF; // DF'*DF matrix in column major format
	FtF = new double [m*m](); // need to de-allocate later
	ctr=0;
	for(j=0;j<m;j++)
	{
		for(i=0;i<m;i++)
		{
			temp_sum=0;
			for(k=0;k<N;k++)
			{
				temp_sum = temp_sum + DF[i][k]*DF[j][k];
			}
			FtF[ctr] = temp_sum; // (FtF)_ij element
			ctr = ctr + 1;
		}
	}

	tcomm1 = MPI_Wtime();   
	MPI_Allreduce(MPI_IN_PLACE, FtF, m*m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	*mpi_time+=tcomm2-tcomm1;

	// ------------------- Compute DF'*phi_res_vec ------------------- //
	double *Ftf; // DF'*phi_res_vec vector of size m x 1
	Ftf = new double [m](); // need to de-allocate later
	ctr=0;
	for(j=0;j<m;j++)
	{
		temp_sum=0;
		for(k=0;k<N;k++)
		{
			temp_sum = temp_sum + DF[j][k]*phi_res_vec[k];
		}
		Ftf[ctr] = temp_sum; // (Ftf)_j element
		ctr = ctr + 1;
	}

	tcomm1 = MPI_Wtime(); 
	MPI_Allreduce(MPI_IN_PLACE, Ftf, m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	*mpi_time+=tcomm2-tcomm1;

	// find pseudoinverse
	double *svec; // vector to store singular values
	svec = new double [m](); // need to de-allocate later
	PseudoInverseTimesVec(FtF,Ftf,svec,m); // svec=Y    
	for(j=0;j<m;j++)
	{
		Ftf[j]=svec[j];
	}

	// ------------------- Compute Anderson update vector am_vec ------------------- //
	for(k=0;k<N;k++)
	{
		temp_sum=0;
		for(j=0;j<m;j++)
		{
			temp_sum = temp_sum + (DX[j][k] + beta_mix*DF[j][k])*Ftf[j];
		}
		am_vec[k] = temp_sum; // (am_vec)_k element
	}

	// De-allocate memory
	delete [] FtF;
	delete [] Ftf;
	delete [] svec;
}

void PseudoInverseTimesVec(double *Ac,double *b,double *x,int m) // returns x=pinv(A)*b, matrix A is m x m  but given in column major format Ac (so m^2 x 1) and vector b is m x 1
{
	int i,j,k,ctr,jj;
	double **A,**U,**V,*w; // A matrix in column major format as Ac. w is the array of singular values
	A = new double* [m](); // need to de-allocate later
	U = new double* [m](); // need to de-allocate later
	V = new double* [m](); // need to de-allocate later
	for(k=0;k<m;k++)
	{
		A[k] = new double [m](); // need to de-allocate later
		U[k] = new double [m](); // need to de-allocate later
		V[k] = new double [m](); // need to de-allocate later
	}

	w = new double [m](); // need to de-allocate later

	// convert column major to matrix form
	ctr=0;
	for(j=0;j<m;j++) // column index
	{
		for(i=0;i<m;i++) // row index
		{
			A[i][j]=Ac[ctr]; // (A)_ij element
			U[i][j]=A[i][j];
			ctr = ctr + 1;
		}
	}

	// Perform SVD on matrix A=UWV'.
	SingularValueDecomp(U,m,m,w,V); // While input, U=A, and while output U=U. need to give output singular values which have been zeroed out if they are small

	// Find Pseudoinverse times vector (pinv(A)*b=(V*diag(1/wj)*U')*b)
	double s,*tmp;
	tmp = new double [m](); // need to de-allocate later
	for(j=0;j<m;j++) // calculate U'*b
	{
		s=0.0;
		if(w[j]) // nonzero result only if wj is nonzero
		{
			for(i=0;i<m;i++) s += U[i][j]*b[i];
			s /= w[j]; // This is the divide by wj
		}
		tmp[j]=s;
	}
	for(j=0;j<m;j++) // Matrix multiply by V to get answer
	{
		s=0.0;
		for(jj=0;jj<m;jj++) s += V[j][jj]*tmp[jj];
		x[j]=s;
	}

	for(k=0;k<m;k++)
	{
		delete [] A[k];
		delete [] U[k];
		delete [] V[k];
	}
	delete [] A;
	delete [] U;
	delete [] V;
	delete [] w;
	delete [] tmp;

}

void SingularValueDecomp(double **a,int m,int n, double *w, double **v) // a is matrix (array) size m x n, A=UWV'. U replaces "a" on output. w is an array of singular values, size 1 x n. V is output as matrix v of size n x n. 
{
	int flag,i,its,j,jj,k,l,nm,Max_its=250;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1 = new double [n](); // need to de-allocate later
	g=scale=anorm=0.0;
	// Householder reduction to bidiagonal form
	for(i=0;i<n;i++)
	{
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if(i<m)
		{
			for(k=i;k<m;k++) scale += fabs(a[k][i]);
			if(scale)
			{
				for(k=i;k<m;k++)
				{
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for(j=1;j<n;j++)
				{
					for(s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for(k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for(k=i;k<m;k++) a[k][i] *= scale;

			}
		}

		w[i]=scale *g;
		g=s=scale=0.0;
		if(i<=m-1 && i!=n-1)
		{
			for(k=l;k<n;k++) scale += fabs(a[i][k]);
			if(scale)
			{
				for(k=l;k<n;k++)
				{
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for(k=l;k<n;k++) rv1[k]=a[i][k]/h;
				for(j=l;j<m;j++)
				{
					for(s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
					for(k=l;k<n;k++) a[j][k] += s*rv1[k];
				}
				for(k=l;k<n;k++) a[i][k] *= scale;
			}
		}

		anorm = max(anorm,(fabs(w[i])+fabs(rv1[i])));

	} // end for loop over i

	// Accumulation of right-hand transformations
	for(i=n-1;i>=0;i--)
	{
		if(i<n-1)
		{
			if(g)
			{
				for(j=l;j<n;j++) // Double division to avoid possible underflow
					v[j][i]=(a[i][j]/a[i][l])/g;
				for(j=l;j<n;j++)
				{
					for(s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for(k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for(j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	} // end for loop over i

	// Accumulation of left-hand transformations
	for(i=min(m,n)-1;i>=0;i--)
	{
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if(g)
		{
			g=1.0/g;
			for(j=l;j<n;j++)
			{
				for(s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for(k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for(j=i;j<m;j++) a[j][i] *= g;
		}else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	} // end for over i

	// Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations
	for(k=n-1;k>=0;k--)
	{
		for(its=0;its<=Max_its;its++)
		{
			flag=1;
			for(l=k;l>=0;l--) // Test for splitting
			{
				nm=l-1; // Note that rv1[0] is always zero
				if((double)(fabs(rv1[l])+anorm) == anorm)
				{
					flag=0;
					break;
				}
				if((double)(fabs(w[nm])+anorm) == anorm) break;
			} // end for over l
			if(flag)
			{
				c=0.0; // Cancellation of rv1[1], if l>1
				s=1.0;
				for(i=l;i<=k;i++)
				{
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if((double)(fabs(f)+anorm)==anorm) break;
					g=w[i];
					h=pythag(f,g);//sqrt(f*f+g*g); //pythag
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=-f*h;
					for(j=0;j<m;j++)
					{
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if(l==k) // Convergence
			{
				if(z<0.0) // Singular value is made nonnegative
				{
					w[k] = -z;
					for(j=0;j<n;j++) v[j][k]= -v[j][k];
				}
				break;
			}
			if(its==Max_its){ printf("no convergence in %d svd iterations \n",Max_its);exit(1);}
			x=w[l]; // Shift from bottom 2-by-2 minor
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0); //sqrt(f*f+1.0); // pythag
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; // Next QR transformation
			for(j=l;j<=nm;j++)
			{
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);//sqrt(f*f+h*h); // pythag
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for(jj=0;jj<n;jj++)
				{
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z = pythag(f,h);//sqrt(f*f+h*h); // pythag
				//printf("f=%.14f,h=%.14f,z=%.14f \n",f,h,z);
				w[j]=z; // Rotation can be arbitrary if z=0
				if (z)
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(jj=0;jj<m;jj++)
				{
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;

		} // end for over its
	} //end for over k

	delete [] rv1;

	// on output a should be u. But for square matrix u and v are the same. so re-assign a as v.
	for(j=0;j<m;j++) 
	{
		for(i=0;i<m;i++)
			a[i][j]=v[i][j];
	}

	// zero out small singular values
	double wmin,wmax=0.0; // will be the maximum singular value obtained
	for(j=0;j<n;j++) if(w[j] > wmax) wmax=w[j];
	wmin = n*wmax*(2.22044605e-16); // using eps=1e-15
	//wmin = n*(5e-16);
	for(j=0;j<n;j++) if(w[j] < wmin) w[j]=0.0;

}

double pythag(double a, double b) // computes (a^2 + b^2)^0.5
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if(absa>absb) return absa*sqrt(1.0+(double)(absb*absb/(absa*absa)));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(double)(absa*absa/(absb*absb))));
}
