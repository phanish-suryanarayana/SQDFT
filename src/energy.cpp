/** \file energy.cpp
  \brief This file contains functions to compute energy.


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

// function to total energy
void EvaluateTotalEnergy(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double Etot=0.0,Eele=0.0,Ex=0.0,Ec=0.0,Exc=0.0, Exc_del=0.0;
	pSQ->mpi_time=0.0; 
	double tcomm1,tcomm2;

	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				Eele = Eele + 0.5*(pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]*(pSQ->b[k][j][i]-pSQ->rho[k][j][i]))*pSQ->delta*pSQ->delta*pSQ->delta; 
			}
		}
	} 

	// Exchange correlation energy, Exc = Ex + Ec
	double A,alpha1,beta1,beta2,beta3,beta4,C2,rhoi,p;
	p = 1.0;
	A = 0.031091 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;
	C2 = 0.73855876638202 ;
	Exc=0.0;
	Exc_del=0.0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				rhoi = pSQ->rho[k][j][i];
				Ex = -C2*pow(rhoi,4.0/3.0);

				if(rhoi==0)
				{
					Ec=0.0;
				}
				else
				{
					Ec = pow(0.75/(M_PI*rhoi),(1.0/3.0));
					Ec = -2.0*A*(1.0+alpha1*Ec)*log(1.0+1.0/(2.0*A*(beta1*pow(Ec,0.5)+beta2*Ec+beta3*pow(Ec,1.5)+beta4*pow(Ec,(p+1.0))))); 
				}
				Exc = Exc + (Ec*rhoi+Ex)*pSQ->delta*pSQ->delta*pSQ->delta;

				Exc_del = Exc_del + -1.0*pSQ->Vxc[k][j][i]*rhoi*pSQ->delta*pSQ->delta*pSQ->delta;

			}
		}
	} 

	Etot = Eele + Exc + Exc_del;

	tcomm1 = MPI_Wtime();
	double Etot_global;
	MPI_Allreduce(&Etot, &Etot_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	pSQ->Etot = Etot_global - pSQ->Eself;
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;

	pSQ->Etot = pSQ->Etot + pSQ->Ebs + pSQ->Eent + pSQ->Ecorr; 

}

// function to compute replusive energy correction for atoms having rc-overlap
void OverlapCorrection(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double ***VlocJ,***bJ,***VlocI; // can modify the code to do it without storing VlocJ 3d array
	int xstart=-1,xend=-1,ystart=-1,yend=-1,zstart=-1,zend=-1;// limits of intersection domain, if there is intersection then all the indices should be greater than or equal to zero, since all of them are in main domain only
	int Xstart=-1,Xend=-1,Ystart=-1,Yend=-1,Zstart=-1,Zend=-1;
	int xl,yl,zl,xr,yr,zr;
	int nxVloc,nyVloc,nzVloc; // no. of nodes in each direction of intersection+FDn region
	int nXVloc,nYVloc,nZVloc;
	int i,j,k,i0,j0,k0,a,indx,indy,indz;
	double r; // radial distance
	double *DVloc,*DVlocI; // derivative of pseudopotential
	double b_temp,distance,Ecorr=0.0,E_zz_exact=0.0,E_zz_reform=0.0,factor; 
	int atomI_count,atomJ_count;

	atomJ_count=0;
	// Loop over all the atoms in main+rb and compute VlocJ for those atoms whose rb-domain intersects
	// Need to compute VlocJ on intersected+FDn domain to apply FD stencil
	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{
		DVloc = new double [pSQ->Psd[JJ_typ].size](); // need to de-allocate later
		getYD_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc,DVloc,pSQ->Psd[JJ_typ].size); //derivatives for spline

		for(int JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) // loop over main atoms of current type
		{	  
			for(int JJr=0;JJr<pSQ->Atm[JJ_typ].Rx[JJ].n_replica;JJr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rb-domain intersects the processor domain
			{
				xstart = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].start; xend = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].end;
				ystart = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].start; yend = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].end;
				zstart = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].start; zend = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].end;

				nxVloc = xend-xstart+1+2*pSQ->FDn; // no. of nodes in x-direction of intersection+FDn domain
				nyVloc = yend-ystart+1+2*pSQ->FDn; // no. of nodes in y-direction of intersection+FDn domain
				nzVloc = zend-zstart+1+2*pSQ->FDn; // no. of nodes in z-direction of intersection+FDn domain

				// allocate memory to store VlocJ on intersection+FDn domain
				VlocJ = new double** [nzVloc](); // need to de-allocate later
				if(VlocJ == NULL)
				{
					if(rank==0)
						cout << "Memory allocation failed in VlocJ"<< endl;
					exit(1);
				}
				for(k=0;k<nzVloc;k++)
				{
					VlocJ[k] = new double* [nyVloc](); // need to de-allocate later
					if(VlocJ[k] == NULL)
					{
						if(rank==0)
							cout << "Memory allocation failed in VlocJ[k]"<< endl;
						exit(1);
					}

					for(j=0;j<nyVloc;j++)
					{
						VlocJ[k][j] = new double [nxVloc](); // need to de-allocate later
						if(VlocJ == NULL)
						{
							if(rank==0)
								cout << "Memory allocation failed in VlocJ[k][j]"<< endl;
							exit(1);
						}
					}
				}
				bJ = new double** [nzVloc-2*pSQ->FDn](); // need to de-allocate later
				for(k=0;k<nzVloc-2*pSQ->FDn;k++)
				{
					bJ[k] = new double* [nyVloc-2*pSQ->FDn](); // need to de-allocate later
					for(j=0;j<nyVloc-2*pSQ->FDn;j++)
					{
						bJ[k][j] = new double [nxVloc-2*pSQ->FDn](); // need to de-allocate later
					}
				}
				// Now compute VlocJ in the intersection+FDn domain
				i0 = xstart-pSQ->FDn;
				j0 = ystart-pSQ->FDn;
				k0 = zstart-pSQ->FDn;

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

						}
					}
				}


				// Now compute bJ = -(1/4pi)*Laplacian(VlocJ) by using FD stencil on points in intersection region, b[k][j][i]
				for(k=0;k<nzVloc-2*pSQ->FDn;k++) 
				{
					for(j=0;j<nyVloc-2*pSQ->FDn;j++) 
					{
						for(i=0;i<nxVloc-2*pSQ->FDn;i++) 
						{

							b_temp = VlocJ[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));

							for(a=1;a<=pSQ->FDn;a++)
							{
								b_temp += (VlocJ[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn-a] + VlocJ[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn+a] + VlocJ[k+pSQ->FDn][j+pSQ->FDn-a][i+pSQ->FDn] + VlocJ[k+pSQ->FDn][j+pSQ->FDn+a][i+pSQ->FDn] + VlocJ[k+pSQ->FDn-a][j+pSQ->FDn][i+pSQ->FDn] + VlocJ[k+pSQ->FDn+a][j+pSQ->FDn][i+pSQ->FDn])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 			 
							}

							bJ[k][j][i] = b_temp; // bJ(x), i,j,k indices are w.r.t main domain hence subtract processor domains starting node to get local indexing
						}
					}
				}

				// ------------------------------------------------------------------------ //
				// 2nd LOOP OVER ATOMS IN PROC+rb domain
				atomI_count=0;
				for(int II_typ=0;II_typ<pSQ->n_typ;II_typ++) // loop over atom types
				{
					DVlocI = new double [pSQ->Psd[II_typ].size](); // need to de-allocate later
					getYD_gen(pSQ->Psd[II_typ].RadialGrid, pSQ->Psd[II_typ].rVloc,DVlocI,pSQ->Psd[II_typ].size); //derivatives for spline
					for(int II=0;II<pSQ->Atm[II_typ].natm;II++) // loop over main atoms of current type
					{	  
						for(int IIr=0;IIr<pSQ->Atm[II_typ].Rx[II].n_replica;IIr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rb-domain intersects the processor domain
						{
							distance = sqrt(pow((pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord),2) +
									pow((pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord),2) +
									pow((pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord),2));

							// Check if this I atom overlaps with outerloop's J atom (i.e. check if distance between I and J atoms is less than rz-->dist for Z/r)  
							if(distance < (pSQ->Psd[II_typ].rz+pSQ->Psd[JJ_typ].rz) && (atomI_count!=atomJ_count)) // true => overlap exits && exclude J atom itself
							{				 
								// need to find VlocI in the intersection of intersection regions of J and I atom rb-domains with proc domain
								// first find Xstart and Xend nodes of the interesection of intersection regions
								xl = pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].start; xr = pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].end;
								yl = pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].start; yr = pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].end;
								zl = pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].start; zr = pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].end;

								Xstart=-1;Xend=-1;Ystart=-1;Yend=-1;Zstart=-1;Zend=-1;

								if(xl >= xstart && xl <= xend)
									Xstart=xl;
								else if(xstart >= xl && xstart <= xr)
									Xstart=xstart;

								if(xr >= xstart && xr <= xend)
									Xend=xr;
								else if(xend >= xl && xend <= xr)
									Xend=xend;

								if(yl >= ystart && yl <= yend)
									Ystart=yl;
								else if(ystart >= yl && ystart <= yr)
									Ystart=ystart;

								if(yr >= ystart && yr <= yend)
									Yend=yr;
								else if(yend >= yl && yend <= yr)
									Yend=yend;

								if(zl >= zstart && zl <= zend)
									Zstart=zl;
								else if(zstart >= zl && zstart <= zr)
									Zstart=zstart;

								if(zr >= zstart && zr <= zend)
									Zend=zr;
								else if(zend >= zl && zend <= zr)
									Zend=zend;

								if((Xstart!=-1)&&(Xend!=-1)&&(Ystart!=-1)&&(Yend!=-1)&&(Zstart!=-1)&&(Zend!=-1)) // overlap exits (should be always true, since rz-already overlap)
								{
									nXVloc = Xend-Xstart+1; // no. of nodes in x-direction of intersection domain
									nYVloc = Yend-Ystart+1; // no. of nodes in y-direction of intersection domain
									nZVloc = Zend-Zstart+1; // no. of nodes in z-direction of intersection domain

									// allocate memory to store VlocI on intersection+FDn domain
									VlocI = new double** [nZVloc](); // need to de-allocate later
									for(k=0;k<nZVloc;k++)
									{
										VlocI[k] = new double* [nYVloc](); // need to de-allocate later
										for(j=0;j<nYVloc;j++)
										{
											VlocI[k][j] = new double [nXVloc](); // need to de-allocate later
										}
									}

									// Now compute VlocI in the intersection domain
									i0 = Xstart;
									j0 = Ystart;
									k0 = Zstart;

									factor=0.5;  // replica outside main but inside main+rb

									for(k=0;k<nZVloc;k++) 
									{
										for(j=0;j<nYVloc;j++) 
										{
											for(i=0;i<nXVloc;i++) 
											{
												r = sqrt(    ((i0+i)*pSQ->delta  -  pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord)*((i0+i)*pSQ->delta  -  pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord)             +               ((j0+j)*pSQ->delta  -  pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord)*((j0+j)*pSQ->delta  -  pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord)                 +                ((k0+k)*pSQ->delta  -  pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord)*((k0+k)*pSQ->delta  -  pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord)   ); // radial distance to the current atom

												ispline_gen(pSQ->Psd[II_typ].RadialGrid, pSQ->Psd[II_typ].rVloc, pSQ->Psd[II_typ].size, &r,&VlocI[k][j][i],1,DVlocI);

												if(fabs(r-0.0)<TEMP_TOL)
												{
													VlocI[k][j][i] = pSQ->Psd[II_typ].rVloc[1]/pSQ->Psd[II_typ].RadialGrid[1];
												}else
												{
													VlocI[k][j][i] = VlocI[k][j][i]/r;
												}	

												// compute energy correction from I-J interaction,  Ecorr = (E_ZZ)exact - (E_ZZ)reformulated, for only neighbors having rz-overlap
												indx = Xstart-xstart+i;
												indy = Ystart-ystart+j;
												indz = Zstart-zstart+k;
												E_zz_reform += factor*bJ[indz][indy][indx]*VlocI[k][j][i]*pSQ->delta*pSQ->delta*pSQ->delta; // (E_ZZ)reformulated

											}
										}
									}


									// de-allocate memory
									for(k=0;k<nZVloc;k++)
									{
										for(j=0;j<nYVloc;j++)
										{
											delete [] VlocI[k][j];
										}
										delete [] VlocI[k];
									}
									delete [] VlocI;


								} // end if on overlap of intersection regions

								// ***************************
								// Exact Repulsive energy -- Need to compute only for J atoms that are strictly within processor domain


								if(pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[0]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[0]+1)*pSQ->delta-TEMP_TOL))
									if(pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[1]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[1]+1)*pSQ->delta-TEMP_TOL))
										if(pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[2]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[2]+1)*pSQ->delta-TEMP_TOL))
										{
											E_zz_exact += (double)pSQ->Atm[II_typ].Z*pSQ->Atm[JJ_typ].Z*factor/distance; // exact repulsive energy between I and J nuclei, need factor 0.5 to account for double counting, (E_ZZ)exact = 0.5*ZI*ZJ/R_IJ (J is over main atoms and I is all replicas within main+rb)
										}
								// If an atom lies in between the border of proc domain and (right) border of main domain (due to mesh error), we still need to include that for E_zz_exact
								// TODO -- > done after including pnode_e+1 in above if condition

								// ***************************

							} // end if on rz-overlap check

							atomI_count+=1;
						}
					}
					delete [] DVlocI;
				} // end 2nd loop over atoms


				// ------------------------------------------------------------------------ //

				atomJ_count+=1;

				// de-allocate memory
				for(k=0;k<nzVloc;k++)
				{
					for(j=0;j<nyVloc;j++)
					{
						delete [] VlocJ[k][j];
					}
					delete [] VlocJ[k];
				}
				delete [] VlocJ;

				for(k=0;k<nzVloc-2*pSQ->FDn;k++)
				{
					for(j=0;j<nyVloc-2*pSQ->FDn;j++)
					{
						delete [] bJ[k][j];
					}
					delete [] bJ[k];
				}
				delete [] bJ;


			} // end for loop over replica atoms of current main atom (includes the main domain atom too)	  
		} // end for loop over main atoms of current type
		delete [] DVloc;
	} // end for loop over atom types

	Ecorr = E_zz_exact - E_zz_reform; // Ecorr = (E_ZZ)exact - (E_ZZ)reformulated, for only neighbors having rz-overlap

	double Ecorr_global;
	MPI_Allreduce(&Ecorr, &Ecorr_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	Ecorr = Ecorr_global;
	pSQ->Ecorr=Ecorr;
	if(rank==0)
	{
		printf("Energy Correction (Ha/atom) : %.12f\n",Ecorr/pSQ->n_atm);
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// function to compute replusive energy correction for atoms having rc-overlap
void OverlapCorrection_forces(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double ***VlocJ,***bJ,***gradbx,***gradby,***gradbz,***VlocI; // can modify the code to do it without storing VlocJ 3d array

	int xstart=-1,xend=-1,ystart=-1,yend=-1,zstart=-1,zend=-1;// limits of intersection domain, if there is intersection then all the indices should be greater than or equal to zero, since all of them are in main domain only
	int Xstart=-1,Xend=-1,Ystart=-1,Yend=-1,Zstart=-1,Zend=-1;
	int xl,yl,zl,xr,yr,zr;

	int nxVloc,nyVloc,nzVloc; // no. of nodes in each direction of intersection+FDn region
	int nXVloc,nYVloc,nZVloc;
	int i,j,k,i0,j0,k0,a,indx,indy,indz;
	double r; // radial distance
	double *DVloc,*DVlocI; // derivative of pseudopotential
	double b_temp,Db_temp,distance,factor=0.5;
	factor=factor;
	int atomI_count,atomJ_count,atom_count;
	double *fx_zz_reform,*fy_zz_reform,*fz_zz_reform,*fx_zz_exact,*fy_zz_exact,*fz_zz_exact;
	fx_zz_reform = new double [pSQ->n_atm]();
	fy_zz_reform = new double [pSQ->n_atm]();
	fz_zz_reform = new double [pSQ->n_atm]();
	fx_zz_exact = new double [pSQ->n_atm]();
	fy_zz_exact = new double [pSQ->n_atm]();
	fz_zz_exact = new double [pSQ->n_atm]();

	for(k=0;k<pSQ->n_atm;k++)
	{
		pSQ->fcorrx[k]=0.0;
		pSQ->fcorry[k]=0.0;
		pSQ->fcorrz[k]=0.0;
	}

	// Loop over all the atoms in main+rb and compute VlocJ for those atoms whose rb-domain intersects
	// Need to compute VlocJ on intersected+FDn domain to apply FD stencil
	atomJ_count=0;
	atom_count=0;
	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{
		DVloc = new double [pSQ->Psd[JJ_typ].size](); // need to de-allocate later
		getYD_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc,DVloc,pSQ->Psd[JJ_typ].size); //derivatives for spline

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
				for(k=0;k<nzVloc;k++)
				{
					VlocJ[k] = new double* [nyVloc](); // need to de-allocate later
					bJ[k] = new double* [nyVloc](); // need to de-allocate later
					for(j=0;j<nyVloc;j++)
					{
						VlocJ[k][j] = new double [nxVloc](); // need to de-allocate later
						bJ[k][j] = new double [nxVloc](); // need to de-allocate later
					}
				}

				gradbx = new double** [nzVloc-4*pSQ->FDn](); // need to de-allocate later
				gradby = new double** [nzVloc-4*pSQ->FDn](); // need to de-allocate later
				gradbz = new double** [nzVloc-4*pSQ->FDn](); // need to de-allocate later
				for(k=0;k<nzVloc-4*pSQ->FDn;k++)
				{
					gradbx[k] = new double* [nyVloc-4*pSQ->FDn](); // need to de-allocate later
					gradby[k] = new double* [nyVloc-4*pSQ->FDn](); // need to de-allocate later
					gradbz[k] = new double* [nyVloc-4*pSQ->FDn](); // need to de-allocate later

					for(j=0;j<nyVloc-4*pSQ->FDn;j++)
					{
						gradbx[k][j] = new double [nxVloc-4*pSQ->FDn](); // need to de-allocate later
						gradby[k][j] = new double [nxVloc-4*pSQ->FDn](); // need to de-allocate later
						gradbz[k][j] = new double [nxVloc-4*pSQ->FDn](); // need to de-allocate later
					}
				}

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

						}
					}
				}

				// Now compute bJ = -(1/4pi)*Laplacian(VlocJ) on intersect+FDn_grad domain by using FD stencil on points in intersection+FDn+FDn_grad region, bJ[k][j][i]
				for(k=pSQ->FDn;k<nzVloc-pSQ->FDn;k++)
				{
					for(j=pSQ->FDn;j<nyVloc-pSQ->FDn;j++)
					{
						for(i=pSQ->FDn;i<nxVloc-pSQ->FDn;i++)
						{//i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 
							b_temp = VlocJ[k][j][i]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));
							for(a=1;a<=pSQ->FDn;a++)
							{
								b_temp += (VlocJ[k][j][i-a] + VlocJ[k][j][i+a] + VlocJ[k][j-a][i] + VlocJ[k][j+a][i] + VlocJ[k-a][j][i] + VlocJ[k+a][j][i])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 				 
							}
							bJ[k][j][i] = b_temp; // bJ(x), i,j,k indices are w.r.t intersection+FDn+FDN_grad domain 

						}
					}
				}

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
							gradbx[k-2*pSQ->FDn][j-2*pSQ->FDn][i-2*pSQ->FDn]=Db_temp;

							// y-direction
							Db_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_temp += (-bJ[k][j-a][i] + bJ[k][j+a][i])*pSQ->coeff_grad[a]; 
							}
							gradby[k-2*pSQ->FDn][j-2*pSQ->FDn][i-2*pSQ->FDn]=Db_temp;

							// z-direction
							Db_temp = 0;
							for(a=1;a<=pSQ->FDn;a++)
							{
								Db_temp += (-bJ[k-a][j][i] + bJ[k+a][j][i])*pSQ->coeff_grad[a]; 		 
							}
							gradbz[k-2*pSQ->FDn][j-2*pSQ->FDn][i-2*pSQ->FDn]=Db_temp;

						}
					}
				}




				// ------------------------------------------------------------------------ //
				// 2nd LOOP OVER ATOMS IN PROC+rb domain
				atomI_count=0;
				for(int II_typ=0;II_typ<pSQ->n_typ;II_typ++) // loop over atom types
				{
					DVlocI = new double [pSQ->Psd[II_typ].size](); // need to de-allocate later
					getYD_gen(pSQ->Psd[II_typ].RadialGrid, pSQ->Psd[II_typ].rVloc,DVlocI,pSQ->Psd[II_typ].size); //derivatives for spline
					for(int II=0;II<pSQ->Atm[II_typ].natm;II++) // loop over main atoms of current type
					{	  
						for(int IIr=0;IIr<pSQ->Atm[II_typ].Rx[II].n_replica;IIr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rb-domain intersects the processor domain
						{
							distance = sqrt(pow((pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord),2) +
									pow((pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord),2) +
									pow((pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord - pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord),2));

							// Check if this I atom overlaps with outerloop's J atom (i.e. check if distance between I and J atoms is less than rz-->dist for Z/r)  
							if(distance < (pSQ->Psd[II_typ].rz+pSQ->Psd[JJ_typ].rz) && (atomI_count!=atomJ_count)) // true => overlap exits && exclude J atom itself (i.e distance=0)
							{				 
								// need to find VlocI in the intersection of intersection regions of J and I atom rb-domains with proc domain
								// first find Xstart and Xend nodes of the interesection of intersection regions
								xl = pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].start; xr = pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].end;
								yl = pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].start; yr = pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].end;
								zl = pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].start; zr = pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].end;

								Xstart=-1;Xend=-1;Ystart=-1;Yend=-1;Zstart=-1;Zend=-1;

								if(xl >= xstart && xl <= xend)
									Xstart=xl;
								else if(xstart >= xl && xstart <= xr)
									Xstart=xstart;

								if(xr >= xstart && xr <= xend)
									Xend=xr;
								else if(xend >= xl && xend <= xr)
									Xend=xend;

								if(yl >= ystart && yl <= yend)
									Ystart=yl;
								else if(ystart >= yl && ystart <= yr)
									Ystart=ystart;

								if(yr >= ystart && yr <= yend)
									Yend=yr;
								else if(yend >= yl && yend <= yr)
									Yend=yend;

								if(zl >= zstart && zl <= zend)
									Zstart=zl;
								else if(zstart >= zl && zstart <= zr)
									Zstart=zstart;

								if(zr >= zstart && zr <= zend)
									Zend=zr;
								else if(zend >= zl && zend <= zr)
									Zend=zend;

								if((Xstart!=-1)&&(Xend!=-1)&&(Ystart!=-1)&&(Yend!=-1)&&(Zstart!=-1)&&(Zend!=-1)) // overlap exits (should be always true, since rz-already overlap)
								{
									nXVloc = Xend-Xstart+1; // no. of nodes in x-direction of intersection domain
									nYVloc = Yend-Ystart+1; // no. of nodes in y-direction of intersection domain
									nZVloc = Zend-Zstart+1; // no. of nodes in z-direction of intersection domain

									// allocate memory to store VlocI on intersection+FDn domain
									VlocI = new double** [nZVloc](); // need to de-allocate later
									for(k=0;k<nZVloc;k++)
									{
										VlocI[k] = new double* [nYVloc](); // need to de-allocate later
										for(j=0;j<nYVloc;j++)
										{
											VlocI[k][j] = new double [nXVloc](); // need to de-allocate later
										}
									}


									// Now compute VlocI in the intersection domain
									i0 = Xstart;
									j0 = Ystart;
									k0 = Zstart;

									factor=0.5; // replica outside main but inside main+rb

									for(k=0;k<nZVloc;k++) 
									{
										for(j=0;j<nYVloc;j++) 
										{
											for(i=0;i<nXVloc;i++) 
											{
												r = sqrt(    ((i0+i)*pSQ->delta  -  pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord)*((i0+i)*pSQ->delta  -  pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord)             +               ((j0+j)*pSQ->delta  -  pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord)*((j0+j)*pSQ->delta  -  pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord)                 +                ((k0+k)*pSQ->delta  -  pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord)*((k0+k)*pSQ->delta  -  pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord)   ); // radial distance to the current atom

												ispline_gen(pSQ->Psd[II_typ].RadialGrid, pSQ->Psd[II_typ].rVloc, pSQ->Psd[II_typ].size, &r,&VlocI[k][j][i],1,DVlocI);

												if(fabs(r-0.0)<TEMP_TOL)
												{
													VlocI[k][j][i] = pSQ->Psd[II_typ].rVloc[1]/pSQ->Psd[II_typ].RadialGrid[1];
												}else
												{
													VlocI[k][j][i] = VlocI[k][j][i]/r;
												}	

												// compute energy correction from I-J interaction,  Ecorr = (E_ZZ)exact - (E_ZZ)reformulated, for only neighbors having rz-overlap
												indx = Xstart-xstart+i;
												indy = Ystart-ystart+j;
												indz = Zstart-zstart+k;

												fx_zz_reform[atom_count]+=gradbx[indz][indy][indx]*VlocI[k][j][i]*pow(pSQ->delta,3); // x-force on current atom
												fy_zz_reform[atom_count]+=gradby[indz][indy][indx]*VlocI[k][j][i]*pow(pSQ->delta,3); // y-force on current atom
												fz_zz_reform[atom_count]+=gradbz[indz][indy][indx]*VlocI[k][j][i]*pow(pSQ->delta,3); // z-force on current atom


											}
										}
									}


									// de-allocate memory
									for(k=0;k<nZVloc;k++)
									{
										for(j=0;j<nYVloc;j++)
										{
											delete [] VlocI[k][j];
										}
										delete [] VlocI[k];
									}
									delete [] VlocI;


								} // end if on overlap of intersection regions

								// ***************************
								// Exact Repulsive energy -- Need to compute only for J atoms that are strictly within processor domain


								if(pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[0]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[0]+1)*pSQ->delta-TEMP_TOL))
									if(pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[1]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[1]+1)*pSQ->delta-TEMP_TOL))
										if(pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord >= (pSQ->pnode_s[2]*pSQ->delta-TEMP_TOL) && pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord < ((pSQ->pnode_e[2]+1)*pSQ->delta-TEMP_TOL))
										{

											fx_zz_exact[atom_count] += (double)pSQ->Atm[II_typ].Z*pSQ->Atm[JJ_typ].Z*(-pSQ->Atm[II_typ].Rx[II].ProcAtm[IIr].coord + pSQ->Atm[JJ_typ].Rx[JJ].ProcAtm[JJr].coord)/pow(distance,3);
											fy_zz_exact[atom_count] += (double)pSQ->Atm[II_typ].Z*pSQ->Atm[JJ_typ].Z*(-pSQ->Atm[II_typ].Ry[II].ProcAtm[IIr].coord + pSQ->Atm[JJ_typ].Ry[JJ].ProcAtm[JJr].coord)/pow(distance,3);
											fz_zz_exact[atom_count] += (double)pSQ->Atm[II_typ].Z*pSQ->Atm[JJ_typ].Z*(-pSQ->Atm[II_typ].Rz[II].ProcAtm[IIr].coord + pSQ->Atm[JJ_typ].Rz[JJ].ProcAtm[JJr].coord)/pow(distance,3);
										}

								// ***************************

							} // end if on rz-overlap check

							atomI_count+=1;
						}
					}
					delete [] DVlocI;
				} // end 2nd loop over atoms


				// ------------------------------------------------------------------------ //

				atomJ_count += 1;

				// de-allocate memory
				for(k=0;k<nzVloc;k++)
				{
					for(j=0;j<nyVloc;j++)
					{
						delete [] VlocJ[k][j];
						delete [] bJ[k][j];
					}
					delete [] VlocJ[k];
					delete [] bJ[k];
				}
				delete [] VlocJ;
				delete [] bJ;	

				for(k=0;k<nzVloc-4*pSQ->FDn;k++)
				{
					for(j=0;j<nyVloc-4*pSQ->FDn;j++)
					{
						delete [] gradbx[k][j];
						delete [] gradby[k][j];
						delete [] gradbz[k][j];
					}
					delete [] gradbx[k];
					delete [] gradby[k];
					delete [] gradbz[k];
				}
				delete [] gradbx;
				delete [] gradby;
				delete [] gradbz;


			} // end for loop over replica atoms of current main atom (includes the main domain atom too)
			atom_count += 1;
		} // end for loop over main atoms of current type
		delete [] DVloc;
	} // end for loop over atom types


	for(k=0;k<pSQ->n_atm;k++)
	{
		pSQ->fcorrx[k] = fx_zz_exact[k] - fx_zz_reform[k];
		pSQ->fcorry[k] = fy_zz_exact[k] - fy_zz_reform[k];
		pSQ->fcorrz[k] = fz_zz_exact[k] - fz_zz_reform[k];
	}


	delete [] fx_zz_reform;
	delete [] fy_zz_reform;
	delete [] fz_zz_reform;
	delete [] fx_zz_exact;
	delete [] fy_zz_exact;
	delete [] fz_zz_exact;
}

