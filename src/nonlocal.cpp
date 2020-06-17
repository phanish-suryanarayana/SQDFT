/** \file nonlocal.cpp
  \brief This file contains functions to compute nonlocal component calculations in SQ.


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

// function to compute sub-nonlocal matrix times a vector in a matrix free way
// pSQ is the pointer to the SQ object of the DS_SQ data structure. p0,q0,r0 are indices of current sub-Hamiltonian's node w.r.t proc+Rcut domain
void VnlocsubTimesVec(DS_SQ* pSQ,double ***vec,int p0,int q0,int r0,double ***Vv)    // pSQ is the pointer to the SQ object of the DS_SQ data structure. p0,q0,r0 are indices of current sub-Hamiltonian's node w.r.t proc+Rcut domain
{
	int rank,ii,jj,kk,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int inn,jnn,knn,nlocc,ctr,l,m,aa,bb,cc,pp,qq,rr,cnt,count;
	double Atm_x,Atm_y,Atm_z;

	int xl,yl,zl,xr,yr,zr;  // rc-domain limits (l=left), (r=right)
	int xs_rcut,ys_rcut,zs_rcut,xe_rcut,ye_rcut,ze_rcut;  // Rcut domain limits (s=start), (e=end)

	int xstart=-1-pSQ->nloc,xend=-1-pSQ->nloc,ystart=-1-pSQ->nloc,yend=-1-pSQ->nloc,zstart=-1-pSQ->nloc,zend=-1-pSQ->nloc;// limits of intersection domain (rc & Rcut), if there is intersection then all the indices should be greater than or equal to (-pSQ->nloc), since all of them are in main+Rcut domain only

	double Denom,temp_scalar; 
	int cnt_tmp=0;
	double tt0,tt1,tt1_temp,tt2,tt2_temp=0.0,time1=0.0,time2=0.0;
	tt0 = MPI_Wtime();
	tt1_temp = tt0;

	// Loop over all the atoms in main+rb and compute non-local effect from those atoms whose rc-domain intersects proc domain

	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{
		nlocc = ceil(pSQ->Psd[JJ_typ].rc/pSQ->delta);

		for(int JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) // loop over main atoms of current type
		{	  
			for(int JJr=0;JJr<pSQ->Atm[JJ_typ].Rx[JJ].n_replica_Rcut;JJr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rc-domain intersects the processor+Rcut domain
			{

				Atm_x = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].coord; // x-coordinate of replica atom
				Atm_y = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtmRcut[JJr].coord; // y-coordinate of replica atom
				Atm_z = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtmRcut[JJr].coord; // z-coordinate of replica atom

				inn = round(Atm_x/pSQ->delta); // nearest node index of replica atom w.r.t main domain (zero at origin)
				jnn = round(Atm_y/pSQ->delta);
				knn = round(Atm_z/pSQ->delta);

				// If an atom is almost at the center of two nodes then assign the nearest node to the left of the atom
				if(fabs(fabs(inn*pSQ->delta-Atm_x)-pSQ->delta*0.5)<TEMP_TOL)
					inn = floor((Atm_x/pSQ->delta) - 10*TEMP_TOL);
				if(fabs(fabs(jnn*pSQ->delta-Atm_y)-pSQ->delta*0.5)<TEMP_TOL)
					jnn = floor((Atm_y/pSQ->delta) - 10*TEMP_TOL);
				if(fabs(fabs(knn*pSQ->delta-Atm_z)-pSQ->delta*0.5)<TEMP_TOL)
					knn = floor((Atm_z/pSQ->delta) - 10*TEMP_TOL);

				// Find out if the rc-domain of this atom intersects Rcut domain of current node p0,q0,r0. 
				// find start & end nodes (w.r.t main?? domain) of intersection of (rc)-domain around nearest node of current atom and the Rcut domain of p0,q0,r0 node
				// limits of rc-domain of current atom, w.r.t main domain indexing
				xl = inn-nlocc; xr = inn+nlocc;
				yl = jnn-nlocc; yr = jnn+nlocc;
				zl = knn-nlocc; zr = knn+nlocc;

				// limits of Rcut domain of p0,q0,r0 node w.r.t main domain indexing
				xs_rcut = pSQ->pnode_s[0]+(p0-pSQ->nloc)-pSQ->nloc; xe_rcut = pSQ->pnode_s[0]+(p0-pSQ->nloc)+pSQ->nloc; 
				ys_rcut = pSQ->pnode_s[1]+(q0-pSQ->nloc)-pSQ->nloc; ye_rcut = pSQ->pnode_s[1]+(q0-pSQ->nloc)+pSQ->nloc; 
				zs_rcut = pSQ->pnode_s[2]+(r0-pSQ->nloc)-pSQ->nloc; ze_rcut = pSQ->pnode_s[2]+(r0-pSQ->nloc)+pSQ->nloc;

				xstart=-1-pSQ->nloc;xend=-1-pSQ->nloc;ystart=-1-pSQ->nloc;yend=-1-pSQ->nloc;zstart=-1-pSQ->nloc;zend=-1-pSQ->nloc; // start & end nodes of intersection region

				if(xl >= xs_rcut && xl <= xe_rcut)
					xstart=xl;
				else if(xs_rcut >= xl && xs_rcut <= xr)
					xstart=xs_rcut;

				if(xr >= xs_rcut && xr <= xe_rcut)
					xend=xr;
				else if(xe_rcut >= xl && xe_rcut <= xr)
					xend=xe_rcut;

				if(yl >= ys_rcut && yl <= ye_rcut)
					ystart=yl;
				else if(ys_rcut >= yl && ys_rcut <= yr)
					ystart=ys_rcut;

				if(yr >= ys_rcut && yr <= ye_rcut)
					yend=yr;
				else if(ye_rcut >= yl && ye_rcut <= yr)
					yend=ye_rcut;

				if(zl >= zs_rcut && zl <= ze_rcut)
					zstart=zl;
				else if(zs_rcut >= zl && zs_rcut <= zr)
					zstart=zs_rcut;

				if(zr >= zs_rcut && zr <= ze_rcut)
					zend=zr;
				else if(ze_rcut >= zl && ze_rcut <= zr)
					zend=ze_rcut;



				tt1 = MPI_Wtime();
				time1 += (tt1-tt1_temp-tt2_temp);
				tt1_temp =tt1;
				tt2_temp=0.0;


				if((xstart!=-1-pSQ->nloc)&&(xend!=-1-pSQ->nloc)&&(ystart!=-1-pSQ->nloc)&&(yend!=-1-pSQ->nloc)&&(zstart!=-1-pSQ->nloc)&&(zend!=-1-pSQ->nloc))
				{// intersection between rc-domain and Rcut-domain exists, so this atom does contribute to Vnlocsub of this node p0,q0,r0

					count=0;
					for(l=0;l<=pSQ->Psd[JJ_typ].lmax;l++) // loop over quantum number l
					{

						if(l!=pSQ->Atm[JJ_typ].lloc)
						{

							if(l==0)
							{
								Denom = pSQ->Psd[JJ_typ].Denom_s;
							}else if(l==1)
							{
								Denom = pSQ->Psd[JJ_typ].Denom_p;
							}else if(l==2)
							{
								Denom = pSQ->Psd[JJ_typ].Denom_d;
							}else if(l==3)
							{
								Denom = pSQ->Psd[JJ_typ].Denom_f;
							}    

							for(m=-l;m<=l;m++) // loop over quantum number m
							{


								// Find Vnlocsub x vec

								// find dot product of projector and the vector
								temp_scalar=0.0;
								for(k=zstart;k<=zend;k++)
								{
									for(j=ystart;j<=yend;j++)
									{
										for(i=xstart;i<=xend;i++)
										{
											// find index of i,j,k w.r.t Rcut domain around p0,q0,r0 -- as pp,qq,rr
											pp = i-pSQ->pnode_s[0]-(p0-pSQ->nloc)+pSQ->nloc;
											qq = j-pSQ->pnode_s[1]-(q0-pSQ->nloc)+pSQ->nloc;
											rr = k-pSQ->pnode_s[2]-(r0-pSQ->nloc)+pSQ->nloc;

											// find ctr (index) corresponding to i,j,k in the local coords of rc-domain around the atom
											ctr = ((k-knn)+nlocc)*(2*nlocc+1)*(2*nlocc+1) + ((j-jnn)+nlocc)*(2*nlocc+1)+((i-inn)+nlocc); // +nlocc because indexing on 0 to 2*nlocc domain around atom
											if(Denom==0){printf("Error: Denom=0.0 in Nonlocal. \n");exit(1);}

											temp_scalar =  temp_scalar + vec[rr][qq][pp]*pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[count][ctr]/Denom;

											cnt_tmp += 1;

										}
									}
								}

								for(kk=zstart;kk<=zend;kk++) // indices are w.r.t main
								{
									for(jj=ystart;jj<=yend;jj++)
									{
										for(ii=xstart;ii<=xend;ii++)
										{
											// find index of ii,jj,kk w.r.t Rcut domain around p0,q0,r0 -- as aa,bb,cc
											aa = ii-pSQ->pnode_s[0]-(p0-pSQ->nloc)+pSQ->nloc;
											bb = jj-pSQ->pnode_s[1]-(q0-pSQ->nloc)+pSQ->nloc;
											cc = kk-pSQ->pnode_s[2]-(r0-pSQ->nloc)+pSQ->nloc;

											// find cnt (index) corresponding to ii,jj,kk in the local coords of rc-domain around the atom
											// first find index of ii,jj,kk w.r.t rc-domain of atom (0 to 2*nlocc). For this need to find ii,jj,kk relative to nearest node to the atom
											cnt = ((kk-knn)+nlocc)*(2*nlocc+1)*(2*nlocc+1) + ((jj-jnn)+nlocc)*(2*nlocc+1)+((ii-inn)+nlocc); // +nlocc because indexing on 0 to 2*nlocc domain around atom
											Vv[cc][bb][aa] =  Vv[cc][bb][aa] + temp_scalar*pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[count][cnt];

										}
									}
								}


								tt2 = MPI_Wtime();
								time2 += (tt2-tt1);
								tt2_temp =tt2-tt1;


								count += 1;
							} // end loop over m

						}

					} // end loop over l


				} // end if condition for overlap/intersection 

			} // end for loop over replica atoms of current main atom (includes the main domain atom too)

		} // end for loop over main atoms of current type

	} // end for loop over atom types

}


double SphericalHarmonics(double x,double y, double z, double r, int l, int m) // function to evaluate spherical harmonics at (x,y,z) relative to atom position
{

	// constants
	// l=0
	double C00 = 0.282094791773878; // 0.5*sqrt(1/pi)
	// l=1
	double C1m1 = 0.488602511902920; // sqrt(3/(4*pi))
	double C10 = 0.488602511902920; // sqrt(3/(4*pi))
	double C1p1 = 0.488602511902920; // sqrt(3/(4*pi))
	// l=2
	double C2m2 = 1.092548430592079; // 0.5*sqrt(15/pi)
	double C2m1 = 1.092548430592079; // 0.5*sqrt(15/pi)  
	double C20 =  0.315391565252520; // 0.25*sqrt(5/pi)
	double C2p1 = 1.092548430592079; // 0.5*sqrt(15/pi)  
	double C2p2 =  0.546274215296040; // 0.25*sqrt(15/pi)
	// l=3
	double C3m3 =  0.590043589926644; // 0.25*sqrt(35/(2*pi))   
	double C3m2 = 2.890611442640554; // 0.5*sqrt(105/(pi))
	double C3m1 = 0.457045799464466; // 0.25*sqrt(21/(2*pi))
	double C30 =  0.373176332590115; // 0.25*sqrt(7/pi)
	double C3p1 =  0.457045799464466; // 0.25*sqrt(21/(2*pi))
	double C3p2 = 1.445305721320277; //  0.25*sqrt(105/(pi))
	double C3p3 = 0.590043589926644; //  0.25*sqrt(35/(2*pi))

	//double r = sqrt(x*x + y*y + z*z);                                    


	double SH=0.0;

	if(l==0)
		SH = C00;

	if(fabs(r)>TEMP_TOL) //(r!=0)
	{
		if(l==0)
			SH = C00;
		else if(l==1)
		{
			if(m==-1)
				SH = C1m1*(y/r);
			else if(m==0)
				SH = C10*(z/r);
			else if(m==1)
				SH = C1p1*(x/r);
			else{
				printf("incorrect l: %d,m: %d\n",l,m);  
				exit(1);
			}       
		}
		else if(l==2)
		{
			if(m==-2)
				SH = C2m2*(x*y)/(r*r);
			else if(m==-1)
				SH = C2m1*(y*z)/(r*r);
			else if(m==0)
				SH = C20*(-x*x - y*y + 2*z*z)/(r*r);
			else if(m==1)
				SH = C2p1*(z*x)/(r*r);
			else if(m==2)
				SH = C2p2*(x*x - y*y)/(r*r);
			else{
				printf("incorrect l: %d,m: %d\n",l,m);  
				exit(1);
			}		                     
		}
		else if(l==3)
		{
			if(m==-3)
				SH = C3m3*(3*x*x - y*y)*y/(r*r*r);
			else if(m==-2)
				SH = C3m2*(x*y*z)/(r*r*r);
			else if(m==-1)
				SH = C3m1*y*(4*z*z - x*x - y*y)/(r*r*r);
			else if(m==0)
				SH = C30*z*(2*z*z-3*x*x-3*y*y)/(r*r*r);
			else if(m==1)
				SH = C3p1*x*(4*z*z - x*x - y*y)/(r*r*r);
			else if(m==2)
				SH = C3p2*z*(x*x - y*y)/(r*r*r);
			else if(m==3)
				SH = C3p3*x*(x*x-3*y*y)/(r*r*r);
			else{
				printf("incorrect l: %d,m: %d\n",l,m);  
				exit(1);
			}		   
		}else
		{
			printf("Only l=0,1,2,3,supported. Input l:%d\n",l);
			exit(1);
		}

	}// end if(r!=0)

	return SH;
}

// function to pre-compute nonlocal projectors UdVtm for all atoms associated with proc+Rcut domain (any atoms whose rc-domain intersects the proc+Rcut domain), for all l & m
void NonlocalProjectors(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure. 
{
	int rank,ii,jj,kk;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int inn,jnn,knn,nlocc,nnd,ctr,l,m,cnt;
	double Atm_x,Atm_y,Atm_z,*rad_dist,*x_dist,*y_dist,*z_dist,xrc,yrc,zrc;

	double r; // radial distance
	double VlocJ,VlJ,UlJ,SHlmJ; 

	// Loop over all the atoms in main+rb and compute non-local effect from those atoms whose rc-domain intersects proc domain
	// Need to compute VlocJ on intersected+FDn domain to apply FD stencil

	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{
		nlocc = ceil(pSQ->Psd[JJ_typ].rc/pSQ->delta);

		// Find no. of nodes in this intersection region and allocate rad_dist and UdVtm arrays
		nnd = (2*nlocc+1)*(2*nlocc+1)*(2*nlocc+1);
		rad_dist = new double [nnd]();
		x_dist = new double [nnd]();
		y_dist = new double [nnd]();
		z_dist = new double [nnd]();

		for(int JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) // loop over main atoms of current type
		{	
			for(int JJr=0;JJr<pSQ->Atm[JJ_typ].Rx[JJ].n_replica_Rcut;JJr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rb-domain intersects the processor domain
			{
				Atm_x = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].coord; // x-coordinate of replica atom
				Atm_y = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtmRcut[JJr].coord; // y-coordinate of replica atom
				Atm_z = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtmRcut[JJr].coord; // z-coordinate of replica atom

				inn = round(Atm_x/pSQ->delta); // nearest node index of replica atom w.r.t main domain (zero at origin)
				jnn = round(Atm_y/pSQ->delta);
				knn = round(Atm_z/pSQ->delta);

				// If an atom is almost at the center of two nodes then assign the nearest node to the left of the atom
				if(fabs(fabs(inn*pSQ->delta-Atm_x)-pSQ->delta*0.5)<TEMP_TOL)
					inn = floor((Atm_x/pSQ->delta) - 10*TEMP_TOL);
				if(fabs(fabs(jnn*pSQ->delta-Atm_y)-pSQ->delta*0.5)<TEMP_TOL)
					jnn = floor((Atm_y/pSQ->delta) - 10*TEMP_TOL);
				if(fabs(fabs(knn*pSQ->delta-Atm_z)-pSQ->delta*0.5)<TEMP_TOL)
					knn = floor((Atm_z/pSQ->delta) - 10*TEMP_TOL);

				// Compute UdVtm at all nodes in the rc-domain of this atom, for all l & m

				cnt = 0;
				for(l=0;l<=pSQ->Psd[JJ_typ].lmax;l++) // loop over quantum number l
				{
					if(l!=pSQ->Atm[JJ_typ].lloc)
					{  
						for(m=-l;m<=l;m++) // loop over quantum number m
						{
							cnt += 1;
						}
					}
				}

				pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm = new double* [cnt]();
				for(kk=0;kk<cnt;kk++)
				{
					pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[kk] = new double [nnd]; 
				}

				// Loop over nodes in rc-region and find their radial distances from current replica atom
				ctr=0;
				for(kk=-nlocc;kk<=nlocc;kk++)
				{
					for(jj=-nlocc;jj<=nlocc;jj++)
					{
						for(ii=-nlocc;ii<=nlocc;ii++)
						{
							xrc = (inn+ii)*pSQ->delta;
							yrc = (jnn+jj)*pSQ->delta;
							zrc = (knn+kk)*pSQ->delta;
							r = sqrt((xrc-Atm_x)*(xrc-Atm_x)  +  (yrc-Atm_y)*(yrc-Atm_y)  +  (zrc-Atm_z)*(zrc-Atm_z));
							rad_dist[ctr]=r;
							x_dist[ctr] = xrc-Atm_x;
							y_dist[ctr] = yrc-Atm_y;
							z_dist[ctr] = zrc-Atm_z;

							ctr += 1;
						}
					}
				}

				cnt = 0;
				for(l=0;l<=pSQ->Psd[JJ_typ].lmax;l++) // loop over quantum number l
				{

					if(l!=pSQ->Atm[JJ_typ].lloc)
					{  

						for(m=-l;m<=l;m++) // loop over quantum number m
						{

							for(ctr=0;ctr<nnd;ctr++) // loop over nodes in intersection region to find UdVtm
							{

								// find VlocJ, VlJ, UlJ at rad_dist[ctr]
								ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&VlocJ,1,pSQ->Psd[JJ_typ].DVloc);

								if(fabs(rad_dist[ctr]-0.0)<TEMP_TOL)
								{
									VlocJ = pSQ->Psd[JJ_typ].rVloc[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
								}else
								{
									VlocJ = VlocJ/rad_dist[ctr];
								}

								if(l==0)
								{
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVs, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&VlJ,1,pSQ->Psd[JJ_typ].DVsJ);
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rUs, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&UlJ,1,pSQ->Psd[JJ_typ].DUsJ);

									if(fabs(rad_dist[ctr]-0.0)<TEMP_TOL)
									{
										VlJ = pSQ->Psd[JJ_typ].rVs[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
										UlJ = pSQ->Psd[JJ_typ].rUs[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
									}else
									{
										VlJ = VlJ/rad_dist[ctr];
										UlJ = UlJ/rad_dist[ctr];
									}

								}else if(l==1)
								{
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVp, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&VlJ,1,pSQ->Psd[JJ_typ].DVpJ);
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rUp, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&UlJ,1,pSQ->Psd[JJ_typ].DUpJ);


									if(fabs(rad_dist[ctr]-0.0)<TEMP_TOL)
									{
										VlJ = pSQ->Psd[JJ_typ].rVp[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
										UlJ = pSQ->Psd[JJ_typ].rUp[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
									}else
									{
										VlJ = VlJ/rad_dist[ctr];
										UlJ = UlJ/rad_dist[ctr];
									}

								}else if(l==2)
								{
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVd, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&VlJ,1,pSQ->Psd[JJ_typ].DVdJ);
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rUd, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&UlJ,1,pSQ->Psd[JJ_typ].DUdJ);


									if(fabs(rad_dist[ctr]-0.0)<TEMP_TOL)
									{
										VlJ = pSQ->Psd[JJ_typ].rVd[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
										UlJ = pSQ->Psd[JJ_typ].rUd[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
									}else
									{
										VlJ = VlJ/rad_dist[ctr];
										UlJ = UlJ/rad_dist[ctr];
									}

								}else if(l==3)
								{
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVf, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&VlJ,1,pSQ->Psd[JJ_typ].DVfJ);
									ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rUf, pSQ->Psd[JJ_typ].size, &rad_dist[ctr],&UlJ,1,pSQ->Psd[JJ_typ].DUfJ);


									if(fabs(rad_dist[ctr]-0.0)<TEMP_TOL)
									{
										VlJ = pSQ->Psd[JJ_typ].rVf[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
										UlJ = pSQ->Psd[JJ_typ].rUf[1]/pSQ->Psd[JJ_typ].RadialGrid[1];
									}else
									{
										VlJ = VlJ/rad_dist[ctr];
										UlJ = UlJ/rad_dist[ctr];
									}

								}

								// Evaluate Spherical Harmonics at rad_dist[ctr] for current l,m
								SHlmJ = SphericalHarmonics(x_dist[ctr],y_dist[ctr],z_dist[ctr],rad_dist[ctr],l,m);
								pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[cnt][ctr] = UlJ*SHlmJ*(VlJ-VlocJ); // U*del(V)
							}


							cnt += 1;

						} // end loop over m

					}

				} // end loop over l

			} // end for loop over replica atoms of current main atom (includes the main domain atom too)

		} // end for loop over main atoms of current type

		delete [] rad_dist;
		delete [] x_dist;
		delete [] y_dist;
		delete [] z_dist;

	} // end for loop over atom types

}

// function to compute sub-nonlocal matrix times a vector in a matrix free way
// pSQ is the pointer to the SQ object of the DS_SQ data structure. p0,q0,r0 are indices of current sub-Hamiltonian's node w.r.t proc+Rcut domain
void VnlocsubTimesVec_J(DS_SQ* pSQ,double ***vec,int p0,int q0,int r0,double ***Vv, int j_atm, int j_typ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure. p0,q0,r0 are indices of current sub-Hamiltonian's node w.r.t proc+Rcut domain
{
	int rank,ii,jj,kk,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int inn,jnn,knn,nlocc,ctr,l,m,aa,bb,cc,pp,qq,rr,cnt,count;
	double Atm_x,Atm_y,Atm_z;

	int xl,yl,zl,xr,yr,zr;  // rc-domain limits (l=left), (r=right)
	int xs_rcut,ys_rcut,zs_rcut,xe_rcut,ye_rcut,ze_rcut;  // Rcut domain limits (s=start), (e=end)

	int xstart=-1-pSQ->nloc,xend=-1-pSQ->nloc,ystart=-1-pSQ->nloc,yend=-1-pSQ->nloc,zstart=-1-pSQ->nloc,zend=-1-pSQ->nloc;// limits of intersection domain (rc & Rcut), if there is intersection then all the indices should be greater than or equal to (-pSQ->nloc), since all of them are in main+Rcut domain only

	double Denom,temp_scalar; 


	int cnt_tmp=0;
	double tt0,tt1,tt1_temp,tt2,tt2_temp=0.0,time1=0.0,time2=0.0;
	tt0 = MPI_Wtime();
	tt1_temp = tt0;

	for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
			{    
				Vv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc] = 0.0;
			}
		}
	}


	int JJ_typ=j_typ;
	int JJ=j_atm;
	nlocc = ceil(pSQ->Psd[JJ_typ].rc/pSQ->delta);


	for(int JJr=0;JJr<pSQ->Atm[JJ_typ].Rx[JJ].n_replica_Rcut;JJr++) // loop over replica atoms of current main atom (includes the main domain atom too)--> only atoms whose rc-domain intersects the processor+Rcut domain
	{

		Atm_x = pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].coord; // x-coordinate of replica atom
		Atm_y = pSQ->Atm[JJ_typ].Ry[JJ].ProcAtmRcut[JJr].coord; // y-coordinate of replica atom
		Atm_z = pSQ->Atm[JJ_typ].Rz[JJ].ProcAtmRcut[JJr].coord; // z-coordinate of replica atom

		inn = round(Atm_x/pSQ->delta); // nearest node index of replica atom w.r.t main domain (zero at origin)
		jnn = round(Atm_y/pSQ->delta);
		knn = round(Atm_z/pSQ->delta);

		// If an atom is almost at the center of two nodes then assign the nearest node to the left of the atom
		if(fabs(fabs(inn*pSQ->delta-Atm_x)-pSQ->delta*0.5)<TEMP_TOL)
			inn = floor((Atm_x/pSQ->delta) - 10*TEMP_TOL);
		if(fabs(fabs(jnn*pSQ->delta-Atm_y)-pSQ->delta*0.5)<TEMP_TOL)
			jnn = floor((Atm_y/pSQ->delta) - 10*TEMP_TOL);
		if(fabs(fabs(knn*pSQ->delta-Atm_z)-pSQ->delta*0.5)<TEMP_TOL)
			knn = floor((Atm_z/pSQ->delta) - 10*TEMP_TOL);

		// Find out if the rc-domain of this atom intersects Rcut domain of current node p0,q0,r0. 
		// find start & end nodes (w.r.t main?? domain) of intersection of (rc)-domain around nearest node of current atom and the Rcut domain of p0,q0,r0 node
		// limits of rc-domain of current atom, w.r.t main domain indexing
		xl = inn-nlocc; xr = inn+nlocc;
		yl = jnn-nlocc; yr = jnn+nlocc;
		zl = knn-nlocc; zr = knn+nlocc;

		// limits of Rcut domain of p0,q0,r0 node w.r.t main domain indexing
		xs_rcut = pSQ->pnode_s[0]+(p0-pSQ->nloc)-pSQ->nloc; xe_rcut = pSQ->pnode_s[0]+(p0-pSQ->nloc)+pSQ->nloc; 
		ys_rcut = pSQ->pnode_s[1]+(q0-pSQ->nloc)-pSQ->nloc; ye_rcut = pSQ->pnode_s[1]+(q0-pSQ->nloc)+pSQ->nloc; 
		zs_rcut = pSQ->pnode_s[2]+(r0-pSQ->nloc)-pSQ->nloc; ze_rcut = pSQ->pnode_s[2]+(r0-pSQ->nloc)+pSQ->nloc;

		xstart=-1-pSQ->nloc;xend=-1-pSQ->nloc;ystart=-1-pSQ->nloc;yend=-1-pSQ->nloc;zstart=-1-pSQ->nloc;zend=-1-pSQ->nloc; // start & end nodes of intersection region

		if(xl >= xs_rcut && xl <= xe_rcut)
			xstart=xl;
		else if(xs_rcut >= xl && xs_rcut <= xr)
			xstart=xs_rcut;

		if(xr >= xs_rcut && xr <= xe_rcut)
			xend=xr;
		else if(xe_rcut >= xl && xe_rcut <= xr)
			xend=xe_rcut;

		if(yl >= ys_rcut && yl <= ye_rcut)
			ystart=yl;
		else if(ys_rcut >= yl && ys_rcut <= yr)
			ystart=ys_rcut;

		if(yr >= ys_rcut && yr <= ye_rcut)
			yend=yr;
		else if(ye_rcut >= yl && ye_rcut <= yr)
			yend=ye_rcut;

		if(zl >= zs_rcut && zl <= ze_rcut)
			zstart=zl;
		else if(zs_rcut >= zl && zs_rcut <= zr)
			zstart=zs_rcut;

		if(zr >= zs_rcut && zr <= ze_rcut)
			zend=zr;
		else if(ze_rcut >= zl && ze_rcut <= zr)
			zend=ze_rcut;



		tt1 = MPI_Wtime();
		time1 += (tt1-tt1_temp-tt2_temp);
		tt1_temp =tt1;
		tt2_temp=0.0;

		if((xstart!=-1-pSQ->nloc)&&(xend!=-1-pSQ->nloc)&&(ystart!=-1-pSQ->nloc)&&(yend!=-1-pSQ->nloc)&&(zstart!=-1-pSQ->nloc)&&(zend!=-1-pSQ->nloc))
		{// intersection between rc-domain and Rcut-domain exists, so this atom does contribute to Vnlocsub of this node p0,q0,r0

			count=0;
			for(l=0;l<=pSQ->Psd[JJ_typ].lmax;l++) // loop over quantum number l
			{

				if(l!=pSQ->Atm[JJ_typ].lloc)
				{

					if(l==0)
					{
						Denom = pSQ->Psd[JJ_typ].Denom_s;
					}else if(l==1)
					{
						Denom = pSQ->Psd[JJ_typ].Denom_p;
					}else if(l==2)
					{
						Denom = pSQ->Psd[JJ_typ].Denom_d;
					}else if(l==3)
					{
						Denom = pSQ->Psd[JJ_typ].Denom_f;
					}    

					for(m=-l;m<=l;m++) // loop over quantum number m
					{
						// Find Vnlocsub x vec

						// find dot product of projector and the vector
						temp_scalar=0.0;
						for(k=zstart;k<=zend;k++)
						{
							for(j=ystart;j<=yend;j++)
							{
								for(i=xstart;i<=xend;i++)
								{
									// find index of i,j,k w.r.t Rcut domain around p0,q0,r0 -- as pp,qq,rr
									pp = i-pSQ->pnode_s[0]-(p0-pSQ->nloc)+pSQ->nloc;
									qq = j-pSQ->pnode_s[1]-(q0-pSQ->nloc)+pSQ->nloc;
									rr = k-pSQ->pnode_s[2]-(r0-pSQ->nloc)+pSQ->nloc;

									// find ctr (index) corresponding to i,j,k in the local coords of rc-domain around the atom
									ctr = ((k-knn)+nlocc)*(2*nlocc+1)*(2*nlocc+1) + ((j-jnn)+nlocc)*(2*nlocc+1)+((i-inn)+nlocc); // +nlocc because indexing on 0 to 2*nlocc domain around atom
									if(Denom==0){printf("Error: Denom=0.0 in Nonlocal. \n");exit(1);}

									temp_scalar =  temp_scalar + vec[rr][qq][pp]*pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[count][ctr]/Denom;

									cnt_tmp += 1;

								}
							}
						}

						for(kk=zstart;kk<=zend;kk++) // indices are w.r.t main
						{
							for(jj=ystart;jj<=yend;jj++)
							{
								for(ii=xstart;ii<=xend;ii++)
								{
									// find index of ii,jj,kk w.r.t Rcut domain around p0,q0,r0 -- as aa,bb,cc
									aa = ii-pSQ->pnode_s[0]-(p0-pSQ->nloc)+pSQ->nloc;
									bb = jj-pSQ->pnode_s[1]-(q0-pSQ->nloc)+pSQ->nloc;
									cc = kk-pSQ->pnode_s[2]-(r0-pSQ->nloc)+pSQ->nloc;

									// find cnt (index) corresponding to ii,jj,kk in the local coords of rc-domain around the atom
									// first find index of ii,jj,kk w.r.t rc-domain of atom (0 to 2*nlocc). For this need to find ii,jj,kk relative to nearest node to the atom
									cnt = ((kk-knn)+nlocc)*(2*nlocc+1)*(2*nlocc+1) + ((jj-jnn)+nlocc)*(2*nlocc+1)+((ii-inn)+nlocc); // +nlocc because indexing on 0 to 2*nlocc domain around atom


									Vv[cc][bb][aa] =  Vv[cc][bb][aa] + temp_scalar*pSQ->Atm[JJ_typ].Rx[JJ].ProcAtmRcut[JJr].UdVtm[count][cnt];


								}
							}
						}


						tt2 = MPI_Wtime();
						time2 += (tt2-tt1);
						tt2_temp =tt2-tt1;


						count += 1;
					} // end loop over m

				}

			} // end loop over l


		} // end if condition for overlap/intersection 

	} // end for loop over replica atoms of current main atom (includes the main domain atom too)

}

