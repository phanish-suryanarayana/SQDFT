/** \file forces.cpp
  \brief This file contains functions that compute the forces. 

  Functions to compute both local and non-local component of the forces.

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

// function to compute local component of the forces
void ForcesLocal(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double ttime00,ttime0,ttime1,ttime2,ttime3,ttime4,timei=0.0,timeb=0.0,timev=0.0,timeg=0.0;
	double ***VlocJ,***VlocJ_tilda,***bJ; // can modify the code to do it without storing VlocJ 3d array
	int xstart=-1,xend=-1,ystart=-1,yend=-1,zstart=-1,zend=-1;// limits of intersection domain, if there is intersection then all the indices should be greater than or equal to zero, since all of them are in main domain only
	int nxVloc,nyVloc,nzVloc; // no. of nodes in each direction of intersection+FDn region
	int i,j,k,i0,j0,k0,a;
	double r; // radial distance
	double *DVloc; // derivative of pseudopotential
	double b_temp,Db_temp,f_temp,Db_tilda_temp,DV_temp,DV_tilda_temp,***bJ_tilda;//,*Forces_x,*Forces_y,*Forces_z;
	double rc_tilda = pSQ->rc_tilda; 

	for(k=0;k<pSQ->n_atm;k++)
	{
		pSQ->flocx[k]=0.0;
		pSQ->flocy[k]=0.0;
		pSQ->flocz[k]=0.0;

		pSQ->fcorrx[k]=0.0;
		pSQ->fcorry[k]=0.0;
		pSQ->fcorrz[k]=0.0;
	}

	// Loop over all the atoms in main+rb and compute VlocJ for those atoms whose rb-domain intersects
	// Need to compute VlocJ on intersected+FDn domain to apply FD stencil
	int atom_count=0;
	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{

		ttime00 = MPI_Wtime();

		DVloc = new double [pSQ->Psd[JJ_typ].size](); // need to de-allocate later
		getYD_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc,DVloc,pSQ->Psd[JJ_typ].size); //derivatives for spline

		ttime0 = MPI_Wtime();
		timei += ttime0-ttime00;

		for(int JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) // loop over main atoms of current type
		{	  
			for(int JJr=0;JJr<pSQ->Atm[JJ_typ].Rx[JJ].n_replica;JJr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rb-domain intersects the processor domain
			{

				xstart = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].start; xend = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].end;
				ystart = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].start; yend = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].end;
				zstart = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].start; zend = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].end;

				nxVloc = xend-xstart+1+2*pSQ->FDn+2*pSQ->FDn; //------// no. of nodes in x-direction of intersection+FDn+FDn_grad domain
				nyVloc = yend-ystart+1+2*pSQ->FDn+2*pSQ->FDn; //------// no. of nodes in y-direction of intersection+FDn+FDn_grad domain
				nzVloc = zend-zstart+1+2*pSQ->FDn+2*pSQ->FDn; //------// no. of nodes in z-direction of intersection+FDn+FDn_grad domain

				// allocate memory to store VlocJ on intersection+FDn domain
				VlocJ = new double** [nzVloc](); // need to de-allocate later
				bJ = new double** [nzVloc](); // need to de-allocate later
				VlocJ_tilda = new double** [nzVloc](); // need to de-allocate later
				bJ_tilda = new double** [nzVloc](); // need to de-allocate later
				if(VlocJ == NULL)
				{
					if(rank==0)
						cout << "Memory allocation failed in VlocJ"<< endl;
					exit(1);
				}
				for(k=0;k<nzVloc;k++)
				{
					VlocJ[k] = new double* [nyVloc](); // need to de-allocate later
					bJ[k] = new double* [nyVloc](); // need to de-allocate later
					VlocJ_tilda[k] = new double* [nyVloc](); // need to de-allocate later
					bJ_tilda[k] = new double* [nyVloc](); // need to de-allocate later
					if(VlocJ[k] == NULL)
					{
						if(rank==0)
							cout << "Memory allocation failed in VlocJ[k]"<< endl;
						exit(1);
					}

					for(j=0;j<nyVloc;j++)
					{
						VlocJ[k][j] = new double [nxVloc](); // need to de-allocate later
						bJ[k][j] = new double [nxVloc](); // need to de-allocate later
						VlocJ_tilda[k][j] = new double [nxVloc](); // need to de-allocate later
						bJ_tilda[k][j] = new double [nxVloc](); // need to de-allocate later
						if(VlocJ == NULL)
						{
							if(rank==0)
								cout << "Memory allocation failed in VlocJ[k][j]"<< endl;
							exit(1);
						}
					}
				}

				ttime1 = MPI_Wtime();

				// Now compute VlocJ in the intersection+FDn+FDn_grad domain
				i0 = xstart-pSQ->FDn-pSQ->FDn; //-------------
				j0 = ystart-pSQ->FDn-pSQ->FDn; //-------------
				k0 = zstart-pSQ->FDn-pSQ->FDn; //-------------


				for(k=0;k<nzVloc;k++) 
				{
					for(j=0;j<nyVloc;j++) 
					{
						for(i=0;i<nxVloc;i++) 
						{
							r = sqrt(    ((i0+i)*pSQ->delta  -  pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord)*((i0+i)*pSQ->delta  -  pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord)             +               ((j0+j)*pSQ->delta  -  pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord)*((j0+j)*pSQ->delta  -  pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord)                 +                ((k0+k)*pSQ->delta  -  pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord)*((k0+k)*pSQ->delta  -  pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord)   ); // radial distance to the current atom

							ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc, pSQ->Psd[JJ_typ].size, &r,&VlocJ[k][j][i],1,DVloc);

							if(fabs(r-0.0)<TEMP_TOL)
							{
								VlocJ[k][j][i] = pSQ->Psd[JJ_typ].rVloc[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
							}else
							{
								VlocJ[k][j][i] = VlocJ[k][j][i]/r;
							}

							if(r < rc_tilda)
							{
								VlocJ_tilda[k][j][i] = -pSQ->Atm[JJ_typ].Z*(9*pow(r,7) - 30*rc_tilda*pow(r,6) + 28*rc_tilda*rc_tilda*pow(r,5) - 14*pow(rc_tilda,5)*pow(r,2) + 12*pow(rc_tilda,7) )/(5*pow(rc_tilda,8));
							}else
							{
								VlocJ_tilda[k][j][i] = -pSQ->Atm[JJ_typ].Z/r;
							}

						}
					}
				}

				ttime2 = MPI_Wtime();
				timev += ttime2-ttime1;

				// Now compute bJ = -(1/4pi)*Laplacian(VlocJ) on intersect+FDn_grad domain by using FD stencil on points in intersection+FDn+FDn_grad region, bJ[k][j][i]
				for(k=pSQ->FDn;k<nzVloc-pSQ->FDn;k++)
				{
					for(j=pSQ->FDn;j<nyVloc-pSQ->FDn;j++)
					{
						for(i=pSQ->FDn;i<nxVloc-pSQ->FDn;i++)
						{//i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 
							// cout << pSQ->b[k][j][i]<<endl;
							// printf("%u,%u,%u \n",k,j,i);
							//cout<<pSQ->coeff_lap[0]<<endl;
							b_temp = VlocJ[k][j][i]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));
							for(a=1;a<=pSQ->FDn;a++)
							{
								b_temp += (VlocJ[k][j][i-a] + VlocJ[k][j][i+a] + VlocJ[k][j-a][i] + VlocJ[k][j+a][i] + VlocJ[k-a][j][i] + VlocJ[k+a][j][i])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 				 
							}

							bJ[k][j][i] = b_temp; // bJ(x), i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 
							b_temp = VlocJ_tilda[k][j][i]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));
							for(a=1;a<=pSQ->FDn;a++)
							{
								b_temp += (VlocJ_tilda[k][j][i-a] + VlocJ_tilda[k][j][i+a] + VlocJ_tilda[k][j-a][i] + VlocJ_tilda[k][j+a][i] + VlocJ_tilda[k-a][j][i] + VlocJ_tilda[k+a][j][i])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 				 
							}
							bJ_tilda[k][j][i] = b_temp; // bJ(x), i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 
						}
					}
				}

				ttime3 = MPI_Wtime();
				timeb += ttime3-ttime2;

				// Now compute grad(bJ) in the intersection domain by using FD stencil on points in intersection+FDn_grad region
				for(k=2*pSQ->FDn;k<nzVloc-2*pSQ->FDn;k++)
				{
					for(j=2*pSQ->FDn;j<nyVloc-2*pSQ->FDn;j++)
					{
						for(i=2*pSQ->FDn;i<nxVloc-2*pSQ->FDn;i++)
						{//i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 

							// x-direction
							Db_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_temp += (-bJ[k][j][i-a] + bJ[k][j][i+a])*pSQ->coeff_grad[a]; 				 
							}

							f_temp = Db_temp*(pSQ->phi[k+k0-pSQ->pnode_s[2]+pSQ->FDn][j+j0-pSQ->pnode_s[1]+pSQ->FDn][i+i0-pSQ->pnode_s[0]+pSQ->FDn]-VlocJ[k][j][i])*pow(pSQ->delta,3); 
							pSQ->flocx[atom_count]+=f_temp; // x-force on current atom

							Db_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_tilda_temp += (-bJ_tilda[k][j][i-a] + bJ_tilda[k][j][i+a])*pSQ->coeff_grad[a]; 				 
							}
							DV_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_temp += (-VlocJ[k][j][i-a] + VlocJ[k][j][i+a])*pSQ->coeff_grad[a]; 				 
							}
							DV_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_tilda_temp += (-VlocJ_tilda[k][j][i-a] + VlocJ_tilda[k][j][i+a])*pSQ->coeff_grad[a]; 				 
							}

							f_temp = 0.5*(Db_tilda_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]-VlocJ_tilda[k][j][i]) + Db_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]+VlocJ[k][j][i])  + (DV_tilda_temp - DV_temp)*(pSQ->b_tilda[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]] + pSQ->b[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]])  +  bJ[k][j][i]*DV_temp  -  bJ_tilda[k][j][i]*DV_tilda_temp)*pow(pSQ->delta,3); 
							pSQ->fcorrx[atom_count]+=f_temp; // x-force on current atom



							// y-direction
							Db_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_temp += (-bJ[k][j-a][i] + bJ[k][j+a][i])*pSQ->coeff_grad[a]; 
							}
							f_temp = Db_temp*(pSQ->phi[k+k0-pSQ->pnode_s[2]+pSQ->FDn][j+j0-pSQ->pnode_s[1]+pSQ->FDn][i+i0-pSQ->pnode_s[0]+pSQ->FDn]-VlocJ[k][j][i])*pow(pSQ->delta,3); 
							pSQ->flocy[atom_count]+=f_temp; // y-force on current atom

							Db_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_tilda_temp += (-bJ_tilda[k][j-a][i] + bJ_tilda[k][j+a][i])*pSQ->coeff_grad[a]; 				 
							}
							DV_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_temp += (-VlocJ[k][j-a][i] + VlocJ[k][j+a][i])*pSQ->coeff_grad[a]; 
							}
							DV_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_tilda_temp += (-VlocJ_tilda[k][j-a][i] + VlocJ_tilda[k][j+a][i])*pSQ->coeff_grad[a]; 				 
							}

							f_temp = 0.5*(Db_tilda_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]-VlocJ_tilda[k][j][i]) + Db_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]+VlocJ[k][j][i])  + (DV_tilda_temp - DV_temp)*(pSQ->b_tilda[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]] + pSQ->b[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]])  +  bJ[k][j][i]*DV_temp  -  bJ_tilda[k][j][i]*DV_tilda_temp)*pow(pSQ->delta,3); 
							pSQ->fcorry[atom_count]+=f_temp; // y-force on current atom

							// z-direction
							Db_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_temp += (-bJ[k-a][j][i] + bJ[k+a][j][i])*pSQ->coeff_grad[a]; 		 
							}
							f_temp = Db_temp*(pSQ->phi[k+k0-pSQ->pnode_s[2]+pSQ->FDn][j+j0-pSQ->pnode_s[1]+pSQ->FDn][i+i0-pSQ->pnode_s[0]+pSQ->FDn]-VlocJ[k][j][i])*pow(pSQ->delta,3); 
							pSQ->flocz[atom_count]+=f_temp; // z-force on current atom

							Db_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_tilda_temp += (-bJ_tilda[k-a][j][i] + bJ_tilda[k+a][j][i])*pSQ->coeff_grad[a]; 				 
							}
							DV_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_temp += (-VlocJ[k-a][j][i] + VlocJ[k+a][j][i])*pSQ->coeff_grad[a]; 		 
							}
							DV_tilda_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								DV_tilda_temp += (-VlocJ_tilda[k-a][j][i] + VlocJ_tilda[k+a][j][i])*pSQ->coeff_grad[a]; 				 
							}

							f_temp = 0.5*(Db_tilda_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]-VlocJ_tilda[k][j][i]) + Db_temp*(pSQ->Vc[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]+VlocJ[k][j][i])  + (DV_tilda_temp - DV_temp)*(pSQ->b_tilda[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]] + pSQ->b[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]])  +  bJ[k][j][i]*DV_temp  -  bJ_tilda[k][j][i]*DV_tilda_temp)*pow(pSQ->delta,3); 
							pSQ->fcorrz[atom_count]+=f_temp; // z-force on current atom

						}
					}
				}

				ttime4 = MPI_Wtime();
				timeg += ttime4-ttime3;

				// de-allocate memory
				for(k=0;k<nzVloc;k++)
				{
					for(j=0;j<nyVloc;j++)
					{
						delete [] VlocJ[k][j];
						delete [] bJ[k][j];
						delete [] VlocJ_tilda[k][j];
						delete [] bJ_tilda[k][j];
					}
					delete [] VlocJ[k];
					delete [] bJ[k];
					delete [] VlocJ_tilda[k];
					delete [] bJ_tilda[k];
				}
				delete [] VlocJ;
				delete [] bJ;
				delete [] VlocJ_tilda;
				delete [] bJ_tilda;

			} // end for loop over replica atoms of current main atom (includes the main domain atom too)
			atom_count += 1;
		} // end for loop over main atoms of current type
		delete [] DVloc;
	} // end for loop over atom types

	if(rank==0)
	{
		printf("Time stats: Interp=%.4f, VJ calc = %.4f, bJ calc = %.4f, gradbJ calc = %.4f \n",timei,timev,timeb,timeg);
	}

}

// function to compute non-local component of the forces
void ForcesNonLocal(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k,ii,jj,kk,a,atom_count,JJ,JJ_typ;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//if(rank==0)
	//printf("Computing non-local component of the forces... ");

	for(k=0;k<pSQ->n_atm;k++)
	{
		pSQ->fnlocx[k]=0.0;
		pSQ->fnlocy[k]=0.0;
		pSQ->fnlocz[k]=0.0;
	}

	double ***vec,***Hv,***t0,***t1,***t2,***DMcol,***gradx_DM,***grady_DM,***gradz_DM;

	gradx_DM = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	grady_DM = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	gradz_DM = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	vec = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	Hv = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t0 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t1 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t2 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	for(k=0;k<2*pSQ->nloc+1;k++)
	{
		gradx_DM[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		grady_DM[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		gradz_DM[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		vec[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		Hv[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t0[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t1[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t2[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		for(j=0;j<2*pSQ->nloc+1;j++)
		{
			gradx_DM[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			grady_DM[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			gradz_DM[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			vec[k][j]=new double [2*pSQ->nloc+1](); // need to de-allocate later
			Hv[k][j]=new double [2*pSQ->nloc+1](); // need to de-allocate later
			t0[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			t1[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			t2[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
		}
	}
	DMcol = new double** [2*pSQ->nloc+1+2*pSQ->FDn](); // need to de-allocate later
	for(k=0;k<2*pSQ->nloc+1+2*pSQ->FDn;k++)
	{
		DMcol[k] = new double* [2*pSQ->nloc+1+2*pSQ->FDn](); // need to de-allocate later
		for(j=0;j<2*pSQ->nloc+1+2*pSQ->FDn;j++)
		{
			DMcol[k][j]=new double [2*pSQ->nloc+1+2*pSQ->FDn](); // need to de-allocate later
		}
	}

	// loop over finite difference nodes in the processor domain
	int nd=0,nq; 
	for(k=pSQ->nloc+1-1;k<=pSQ->np_z+pSQ->nloc-1;k++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(j=pSQ->nloc+1-1;j<=pSQ->np_y+pSQ->nloc-1;j++) // no. of proc nodes in each direction
		{
			for(i=pSQ->nloc+1-1;i<=pSQ->np_x+pSQ->nloc-1;i++) // no. of proc nodes in each direction
			{
				/// For each FD node, loop over quadrature order to find Chebyshev expansion components, Use the HsubTimesVec function this.
				for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
				{
					for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
					{
						for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
						{
							t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=0.0;
							if(ii==0 && jj==0 && kk==0)
								t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=1.0;
						}
					}
				}
				HsubTimesVec(pSQ,t0,i,j,k,Hv); // Hv=Hsub*t0. Here t0 is the vector and i,j,k are the indices of node in proc domain to which the Hsub corresponds to and the indices are w.r.t proc+Rcut domain
				for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
				{
					for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
					{
						for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
						{
							t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=(Hv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]-pSQ->chi[nd]*t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc])/pSQ->zee[nd];

							DMcol[kk+pSQ->nloc+pSQ->FDn][jj+pSQ->nloc+pSQ->FDn][ii+pSQ->nloc+pSQ->FDn] = 0.0;
						}
					}
				}

				// loop over quadrature order
				for(nq=0;nq<=pSQ->npl;nq++)
				{
					if(nq==0)
					{
						for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
						{
							for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
							{
								for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
								{    
									DMcol[kk+pSQ->nloc+pSQ->FDn][jj+pSQ->nloc+pSQ->FDn][ii+pSQ->nloc+pSQ->FDn] += (double)(t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]*pSQ->Ci[nd][nq]/(pow(pSQ->delta,3)));

								}
							}
						}

					}else if(nq==1)
					{
						for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
						{
							for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
							{
								for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
								{    
									DMcol[kk+pSQ->nloc+pSQ->FDn][jj+pSQ->nloc+pSQ->FDn][ii+pSQ->nloc+pSQ->FDn] += (double)(t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]*pSQ->Ci[nd][nq]/(pow(pSQ->delta,3)));
								}
							}
						}

					}else
					{
						HsubTimesVec(pSQ,t1,i,j,k,Hv); // Hv=Hsub*t1
						for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
						{
							for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
							{
								for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
								{
									t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=(2*(Hv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]-pSQ->chi[nd]*t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc])/pSQ->zee[nd]) - t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];		      

									DMcol[kk+pSQ->nloc+pSQ->FDn][jj+pSQ->nloc+pSQ->FDn][ii+pSQ->nloc+pSQ->FDn] += (double)(t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]*pSQ->Ci[nd][nq]/(pow(pSQ->delta,3)));

									t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];
									t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];
								}
							}
						}

					} // end if conditions

				} // end loop over nq


				// find gradient of DM cols
				for(kk=pSQ->FDn;kk<2*pSQ->nloc+1+pSQ->FDn;kk++)
				{
					for(jj=pSQ->FDn;jj<2*pSQ->nloc+1+pSQ->FDn;jj++)
					{
						for(ii=pSQ->FDn;ii<2*pSQ->nloc+1+pSQ->FDn;ii++)
						{
							gradx_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] = 0.0; 		
							grady_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] = 0.0; 
							gradz_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] = 0.0; 
						}
					}
				}

				for(kk=pSQ->FDn;kk<2*pSQ->nloc+1+pSQ->FDn;kk++)
				{
					for(jj=pSQ->FDn;jj<2*pSQ->nloc+1+pSQ->FDn;jj++)
					{
						for(ii=pSQ->FDn;ii<2*pSQ->nloc+1+pSQ->FDn;ii++)
						{
							// x-direction
							for(a=1;a<=pSQ->FDn;a++)
							{
								gradx_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] += (-DMcol[kk][jj][ii-a] + DMcol[kk][jj][ii+a])*pSQ->coeff_grad[a]; 				 
							}
							// y-direction
							for(a=1;a<=pSQ->FDn;a++)
							{
								grady_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] += (-DMcol[kk][jj-a][ii] + DMcol[kk][jj+a][ii])*pSQ->coeff_grad[a]; 				 
							}
							// z-direction
							for(a=1;a<=pSQ->FDn;a++)
							{
								gradz_DM[kk-pSQ->FDn][jj-pSQ->FDn][ii-pSQ->FDn] += (-DMcol[kk-a][jj][ii] + DMcol[kk+a][jj][ii])*pSQ->coeff_grad[a]; 				 
							}

						}
					}
				}

				// Loop over all main atoms to compute contribution of current FD node to the non-local force component for that atom
				atom_count = 0;
				for(JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
				{

					for(JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) // loop over main atoms of current type
					{

						if(pSQ->Atm[JJ_typ].Rx[JJ].n_replica_Rcut > 0)
						{

							VnlocsubTimesVec_J(pSQ,gradx_DM,i,j,k,Hv,JJ,JJ_typ); // Hv = VnlocsubJ x gradx_DM
							pSQ->fnlocx[atom_count] += -4*Hv[pSQ->nloc][pSQ->nloc][pSQ->nloc]*pow(pSQ->delta,3);

							VnlocsubTimesVec_J(pSQ,grady_DM,i,j,k,Hv,JJ,JJ_typ); // Hv = VnlocsubJ x grady_DM
							pSQ->fnlocy[atom_count] += -4*Hv[pSQ->nloc][pSQ->nloc][pSQ->nloc]*pow(pSQ->delta,3);

							VnlocsubTimesVec_J(pSQ,gradz_DM,i,j,k,Hv,JJ,JJ_typ); // Hv = VnlocsubJ x gradz_DM
							pSQ->fnlocz[atom_count] += -4*Hv[pSQ->nloc][pSQ->nloc][pSQ->nloc]*pow(pSQ->delta,3);

						}


						atom_count += 1;
					} // end for loop over main atoms of current type

				} // end for loop over atom types

				nd = nd+1;

			}
		}
	} // end loop over FD nodes in proc domain

	for(k=0;k<2*pSQ->nloc+1;k++)
	{
		for(j=0;j<2*pSQ->nloc+1;j++)
		{
			delete [] vec[k][j];
			delete [] gradx_DM[k][j];
			delete [] grady_DM[k][j];
			delete [] gradz_DM[k][j];
			delete [] Hv[k][j];
			delete [] t0[k][j];
			delete [] t1[k][j];
			delete [] t2[k][j];
		}
		delete [] gradx_DM[k];
		delete [] grady_DM[k];
		delete [] gradz_DM[k];
		delete [] vec[k];
		delete [] Hv[k];
		delete [] t0[k];
		delete [] t1[k];
		delete [] t2[k];
	}
	delete [] gradx_DM;
	delete [] grady_DM;
	delete [] gradz_DM;
	delete [] vec;
	delete [] Hv;
	delete [] t0;
	delete [] t1;
	delete [] t2;

	for(k=0;k<2*pSQ->nloc+1+2*pSQ->FDn;k++)
	{
		for(j=0;j<2*pSQ->nloc+1+2*pSQ->FDn;j++)
		{
			delete [] DMcol[k][j];
		}
		delete [] DMcol[k];
	}
	delete [] DMcol;

}

void ForcesTotal(DS_SQ* pSQ)
{
	int rank,k,j,i;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double tt0,tt1,tt2,tt2m,tt1m;
	tt0 = MPI_Wtime();

	double *fmag;
	fmag = new double [pSQ->n_atm]();
	double *Max_forces_ind,*Max_forces,*Max_forces_atyp,*Max_forces_aind; int nmaxf=5; // to print max "nmaxf" forces
	nmaxf = min(nmaxf,pSQ->n_atm);
	Max_forces_ind = new double [nmaxf]();
	Max_forces = new double [nmaxf]();
	Max_forces_atyp = new double [nmaxf]();
	Max_forces_aind = new double [nmaxf]();

	tt1 = MPI_Wtime();
	ForcesLocal(pSQ);
	tt2 = MPI_Wtime();

	if(rank==0)
	{printf("Time taken for local forces = %.2f sec \n",tt2-tt1);cout << " "<<endl;}

	tt1 = MPI_Wtime();
	ForcesNonLocal(pSQ);
	tt2 = MPI_Wtime();

	if(rank==0)
	{printf("Time taken for non-local forces = %.2f sec \n",tt2-tt1);cout << " "<<endl;}
	if(pSQ->Correction == 2)
		OverlapCorrection_forces(pSQ);

	for(k=0;k<pSQ->n_atm;k++)
	{
		pSQ->fx[k] = pSQ->flocx[k] + pSQ->fnlocx[k] + pSQ->fcorrx[k];
		pSQ->fy[k] = pSQ->flocy[k] + pSQ->fnlocy[k] + pSQ->fcorry[k];
		pSQ->fz[k] = pSQ->flocz[k] + pSQ->fnlocz[k] + pSQ->fcorrz[k];
		pSQ->forces[k] = pSQ->fx[k]; 
		pSQ->forces[k+pSQ->n_atm] = pSQ->fy[k];
		pSQ->forces[k+2*pSQ->n_atm] = pSQ->fz[k];
	}
	tt1 = MPI_Wtime();
	MPI_Allreduce(MPI_IN_PLACE, pSQ->forces, 3*pSQ->n_atm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tt2 = MPI_Wtime();

	pSQ->forcs_time=tt2-tt0;
	pSQ->tforcs_mpi=tt2-tt1;

	if(rank==0)
	{printf("Time taken for MPI_Allreduce of forces = %.2f sec \n",tt2-tt1);cout << " "<<endl;}

	if(rank==0)
	{  

		printf("Atomic forces: \n");
		for(k=0;k<pSQ->n_atm;k++)
		{	  
			if(pSQ->n_atm<=8)
			{
				printf("  %.12f    %.12f    %.12f \n",pSQ->forces[k],pSQ->forces[k+pSQ->n_atm],pSQ->forces[k+2*pSQ->n_atm]);
			}else
			{
				if(k==0)
					printf("Total number of atoms is greater than 8. Not printing atom forces.\n");
			}

		}

		printf("Total forces time      : %.2f (%.2f) sec \n",pSQ->forcs_time,pSQ->tforcs_mpi);
	}

	/////////////////////////////////////////////////////////////
	tt1m = MPI_Wtime();
	// Find out maximum forces
	double fmax=-1.0;int max_ind=0,chk,atyp=0,aind=0; //,latconst=1.0
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)	    
		{
			fmag[k]=sqrt(pSQ->forces[k]*pSQ->forces[k] + pSQ->forces[k+pSQ->n_atm]*pSQ->forces[k+pSQ->n_atm] + pSQ->forces[k+2*pSQ->n_atm]*pSQ->forces[k+2*pSQ->n_atm]); // magnitude of the force
			if(fmag[k]>fmax)
			{
				max_ind=k; // store index of max forces
				fmax = fmag[k];
				atyp=ii;
				aind=jj;
			}
			k+=1;
		}
	}
	Max_forces[0]=fmax;
	Max_forces_ind[0]=max_ind;
	Max_forces_atyp[0]=atyp;
	Max_forces_aind[0]=aind;

	for(k=1;k<nmaxf;k++)
	{
		fmax=-1.0;
		j=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
			{
				if(fmag[j]>=fmax && fmag[j]<=Max_forces[k-1])
				{
					chk=0;
					for(i=0;i<k;i++)
					{
						if(j==Max_forces_ind[i])
							chk=1;
					}
					if(chk==0) // not stored before
					{
						max_ind=j;
						fmax=fmag[j];
						atyp=ii;
						aind=jj;		  
					}
				} //end if
				j+=1;
			}
		}
		Max_forces[k]=fmax;
		Max_forces_ind[k]=max_ind;
		Max_forces_atyp[k]=atyp;
		Max_forces_aind[k]=aind;
	}
	tt2m = MPI_Wtime();
	if(rank==0)
	{printf("Time taken to find max %d forces = %.2f sec \n",nmaxf,tt2m-tt1m);cout << " "<<endl;}

	//************************* Symmetrize forces ******************************//
	// Make net force zero
	if(rank==0)
		printf("Symmetrizing forces. \n");

	double fsum_x=0.0,fsum_y=0.0,fsum_z=0.0;
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)	    
		{
			fsum_x += pSQ->forces[k];
			fsum_y += pSQ->forces[k+pSQ->n_atm];
			fsum_z += pSQ->forces[k+2*pSQ->n_atm];
			k+=1;
		}
	}

	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)	    
		{
			pSQ->forces[k] = pSQ->forces[k] - (fsum_x/pSQ->n_atm);
			pSQ->forces[k+pSQ->n_atm] = pSQ->forces[k+pSQ->n_atm] - (fsum_y/pSQ->n_atm);
			pSQ->forces[k+2*pSQ->n_atm] = pSQ->forces[k+2*pSQ->n_atm] - (fsum_z/pSQ->n_atm);
			k+=1;
		}
	}
	//*************************************************************************//
	delete [] fmag;
	delete [] Max_forces_ind;
	delete [] Max_forces;
	delete [] Max_forces_atyp;
	delete [] Max_forces_aind;

}



