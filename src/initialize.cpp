/** \file initialize.cpp
  \brief This file contains functions that initialize the SQDFT calculation. 

  This involves storing replica atoms, creating MPI communication topologies, computing nuclear charge density, inital guess electron density and inital guess electrostatic potential.


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

// function to compute coordinates of replica atoms file
void Replica_atoms(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int px,py,pz,mx,my,mz; // no. of integer multiples in +ve or -ve x,y,z directions
	int count;
	int inn,jnn,knn;
	double Atm_x,Atm_y,Atm_z;
	int xl,yl,zl,xr,yr,zr;  // rb-domain limits (l=left), (r=right)
	int xs_main,ys_main,zs_main,xe_main,ye_main,ze_main;  // processor domain limits (s=start), (e=end)
	int xstart=-1,xend=-1,ystart=-1,yend=-1,zstart=-1,zend=-1; // limits of intersection domain, if there is intersection then all the indices should be greater than or equal to zero, since all of them are in main domain only

	// compute coordinates of replica atoms in main+rb domain
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{

		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			// Loop over all replicas in main+rb and find out if the processor domain overlaps with the rb-domain of each of those replicas
			px = floor(((pSQ->domain[0]+pSQ->Atm[ii].rb - pSQ->Atm[ii].Rx[jj].main_atm)/pSQ->domain[0])+ TEMP_TOL); // no. of positive integer multiples in x-direction
			py = floor(((pSQ->domain[1]+pSQ->Atm[ii].rb - pSQ->Atm[ii].Ry[jj].main_atm)/pSQ->domain[1])+ TEMP_TOL); // no. of positive integer multiples in y-direction
			pz = floor(((pSQ->domain[2]+pSQ->Atm[ii].rb - pSQ->Atm[ii].Rz[jj].main_atm)/pSQ->domain[2])+ TEMP_TOL); // no. of positive integer multiples in z-direction
			mx = ceil(((-pSQ->Atm[ii].rb - pSQ->Atm[ii].Rx[jj].main_atm)/pSQ->domain[0])- TEMP_TOL); // no. of negative integer multiples in x-direction
			my = ceil(((-pSQ->Atm[ii].rb - pSQ->Atm[ii].Ry[jj].main_atm)/pSQ->domain[1])- TEMP_TOL); // no. of negative integer multiples in y-direction
			mz = ceil(((-pSQ->Atm[ii].rb - pSQ->Atm[ii].Rz[jj].main_atm)/pSQ->domain[2])- TEMP_TOL); // no. of negative integer multiples in z-direction
			count=0;
			for(int rr=mz;rr<=pz;rr++)
			{
				for(int qq=my;qq<=py;qq++)
				{
					for(int pp=mx;pp<=px;pp++)
					{

						Atm_x = pSQ->Atm[ii].Rx[jj].main_atm + pp*pSQ->domain[0]; // x-coordinate of replica atom
						Atm_y = pSQ->Atm[ii].Ry[jj].main_atm + qq*pSQ->domain[1]; // y-coordinate of replica atom
						Atm_z = pSQ->Atm[ii].Rz[jj].main_atm + rr*pSQ->domain[2]; // z-coordinate of replica atom

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

						//--------- If the rb-domain around the nearest node of this atom intersects the processor domain, then store overlap info---------

						// find start & end nodes (w.r.t main domain) of intersection of (rb)-domain around nearest node and the processor domain
						// limits of rb-domain
						xl = inn-round(pSQ->Atm[ii].rb/pSQ->delta); xr = inn+round(pSQ->Atm[ii].rb/pSQ->delta);
						yl = jnn-round(pSQ->Atm[ii].rb/pSQ->delta); yr = jnn+round(pSQ->Atm[ii].rb/pSQ->delta);
						zl = knn-round(pSQ->Atm[ii].rb/pSQ->delta); zr = knn+round(pSQ->Atm[ii].rb/pSQ->delta);

						// limits of processor domain
						xs_main = pSQ->pnode_s[0]; xe_main = pSQ->pnode_e[0]; 
						ys_main = pSQ->pnode_s[1]; ye_main = pSQ->pnode_e[1]; 
						zs_main = pSQ->pnode_s[2]; ze_main = pSQ->pnode_e[2];

						xstart=-1;xend=-1;ystart=-1;yend=-1;zstart=-1;zend=-1;

						if(xl >= xs_main && xl <= xe_main)
							xstart=xl;
						else if(xs_main >= xl && xs_main <= xr)
							xstart=xs_main;

						if(xr >= xs_main && xr <= xe_main)
							xend=xr;
						else if(xe_main >= xl && xe_main <= xr)
							xend=xe_main;

						if(yl >= ys_main && yl <= ye_main)
							ystart=yl;
						else if(ys_main >= yl && ys_main <= yr)
							ystart=ys_main;

						if(yr >= ys_main && yr <= ye_main)
							yend=yr;
						else if(ye_main >= yl && ye_main <= yr)
							yend=ye_main;

						if(zl >= zs_main && zl <= ze_main)
							zstart=zl;
						else if(zs_main >= zl && zs_main <= zr)
							zstart=zs_main;

						if(zr >= zs_main && zr <= ze_main)
							zend=zr;
						else if(ze_main >= zl && ze_main <= zr)
							zend=ze_main;

						if((xstart!=-1)&&(xend!=-1)&&(ystart!=-1)&&(yend!=-1)&&(zstart!=-1)&&(zend!=-1))
						{// intersection between rb-domain and processor-domain exists
							count+=1;
						}

					}
				}
			}

			// Now allocate the replica atoms array and store the overlap info
			pSQ->Atm[ii].Rx[jj].n_replica=count; // no. of replica's of jjth atoms (including main domain's jjth atom itself) whose rb-domain intersects processor domain
			pSQ->Atm[ii].Rx[jj].ProcAtm = new DS_ProcAtm[pSQ->Atm[ii].Rx[jj].n_replica]; // need to de-allocate later
			pSQ->Atm[ii].Ry[jj].ProcAtm = new DS_ProcAtm[pSQ->Atm[ii].Rx[jj].n_replica]; // need to de-allocate later
			pSQ->Atm[ii].Rz[jj].ProcAtm = new DS_ProcAtm[pSQ->Atm[ii].Rx[jj].n_replica]; // need to de-allocate later
			count=0;
			for(int rr=mz;rr<=pz;rr++)
			{
				for(int qq=my;qq<=py;qq++)
				{
					for(int pp=mx;pp<=px;pp++)
					{

						Atm_x = pSQ->Atm[ii].Rx[jj].main_atm + pp*pSQ->domain[0]; // x-coordinate of replica atom
						Atm_y = pSQ->Atm[ii].Ry[jj].main_atm + qq*pSQ->domain[1]; // y-coordinate of replica atom
						Atm_z = pSQ->Atm[ii].Rz[jj].main_atm + rr*pSQ->domain[2]; // z-coordinate of replica atom

						inn = round((Atm_x/pSQ->delta)); // nearest node index of replica atom w.r.t main domain (zero at origin)
						jnn = round((Atm_y/pSQ->delta));
						knn = round((Atm_z/pSQ->delta));


						// If an atom is almost at the center of two nodes then assign the nearest node to the left of the atom
						if(fabs(fabs(inn*pSQ->delta-Atm_x)-pSQ->delta*0.5)<TEMP_TOL)
							inn = floor((Atm_x/pSQ->delta) - 10*TEMP_TOL);
						if(fabs(fabs(jnn*pSQ->delta-Atm_y)-pSQ->delta*0.5)<TEMP_TOL)
							jnn = floor((Atm_y/pSQ->delta) - 10*TEMP_TOL);
						if(fabs(fabs(knn*pSQ->delta-Atm_z)-pSQ->delta*0.5)<TEMP_TOL)
							knn = floor((Atm_z/pSQ->delta) - 10*TEMP_TOL);

						//--------- If the rb-domain around the nearest node of this atom intersects the processor domain, then store overlap info---------

						// find start & end nodes (w.r.t main domain) of intersection of (rb)-domain around nearest node and the processor domain
						// limits of rb-domain
						xl = inn-round(pSQ->Atm[ii].rb/pSQ->delta); xr = inn+round(pSQ->Atm[ii].rb/pSQ->delta);
						yl = jnn-round(pSQ->Atm[ii].rb/pSQ->delta); yr = jnn+round(pSQ->Atm[ii].rb/pSQ->delta);
						zl = knn-round(pSQ->Atm[ii].rb/pSQ->delta); zr = knn+round(pSQ->Atm[ii].rb/pSQ->delta);

						// limits of processor domain
						xs_main = pSQ->pnode_s[0]; xe_main = pSQ->pnode_e[0]; 
						ys_main = pSQ->pnode_s[1]; ye_main = pSQ->pnode_e[1]; 
						zs_main = pSQ->pnode_s[2]; ze_main = pSQ->pnode_e[2];

						xstart=-1;xend=-1;ystart=-1;yend=-1;zstart=-1;zend=-1;

						if(xl >= xs_main && xl <= xe_main)
							xstart=xl;
						else if(xs_main >= xl && xs_main <= xr)
							xstart=xs_main;

						if(xr >= xs_main && xr <= xe_main)
							xend=xr;
						else if(xe_main >= xl && xe_main <= xr)
							xend=xe_main;

						if(yl >= ys_main && yl <= ye_main)
							ystart=yl;
						else if(ys_main >= yl && ys_main <= yr)
							ystart=ys_main;

						if(yr >= ys_main && yr <= ye_main)
							yend=yr;
						else if(ye_main >= yl && ye_main <= yr)
							yend=ye_main;

						if(zl >= zs_main && zl <= ze_main)
							zstart=zl;
						else if(zs_main >= zl && zs_main <= zr)
							zstart=zs_main;

						if(zr >= zs_main && zr <= ze_main)
							zend=zr;
						else if(ze_main >= zl && ze_main <= zr)
							zend=ze_main;

						if((xstart!=-1)&&(xend!=-1)&&(ystart!=-1)&&(yend!=-1)&&(zstart!=-1)&&(zend!=-1))
						{// intersection between rb-domain and main-domain exists		      		      
							pSQ->Atm[ii].Rx[jj].ProcAtm[count].coord = Atm_x;
							pSQ->Atm[ii].Ry[jj].ProcAtm[count].coord = Atm_y;
							pSQ->Atm[ii].Rz[jj].ProcAtm[count].coord = Atm_z;
							pSQ->Atm[ii].Rx[jj].ProcAtm[count].start = xstart; pSQ->Atm[ii].Rx[jj].ProcAtm[count].end = xend;
							pSQ->Atm[ii].Ry[jj].ProcAtm[count].start = ystart; pSQ->Atm[ii].Ry[jj].ProcAtm[count].end = yend;
							pSQ->Atm[ii].Rz[jj].ProcAtm[count].start = zstart; pSQ->Atm[ii].Rz[jj].ProcAtm[count].end = zend;

							count+=1;
						}
					}
				}
			}
		} //end for over main atoms of each type
	} // end for over types

}

// function to compute end nodes of processor domain
void Processor_domain(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,count,countx,county,countz;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// number of nodes in each direction of processor domain //(different for the end processor in that direction)
	int ndpx,ndpy,ndpz,ndpx_curr,ndpy_curr,ndpz_curr;
	ndpx = (pSQ->n_int[0]/pSQ->nprocx); //ceil((double)pSQ->n_int[0]/pSQ->nprocx);
	ndpy = (pSQ->n_int[1]/pSQ->nprocy); //ceil((double)pSQ->n_int[1]/pSQ->nprocy);
	ndpz = (pSQ->n_int[2]/pSQ->nprocz); //ceil((double)pSQ->n_int[2]/pSQ->nprocz);

	// Based on the current processor's rank, compute the processor domain end nodes, ordering the processors as z,y,x
	count=0;
	for(countx=0;countx<pSQ->nprocx;countx++) // loop over all processors count
	{
		for(county=0;county<pSQ->nprocy;county++)
		{
			for(countz=0;countz<pSQ->nprocz;countz++)
			{
				if(rank==count) // current processor
				{		 
					// no. of nodes in each direction of current processor domain
					ndpx_curr = ndpx;
					ndpy_curr = ndpy;
					ndpz_curr = ndpz;

					// update the no. of nodes for processor domains that are at the right ends of main domain
					if(countx==pSQ->nprocx-1)
						ndpx_curr = pSQ->n_int[0]-(pSQ->nprocx-1)*ceil((double)pSQ->n_int[0]/pSQ->nprocx);
					if(county==pSQ->nprocy-1)
						ndpy_curr = pSQ->n_int[1]-(pSQ->nprocy-1)*ceil((double)pSQ->n_int[1]/pSQ->nprocy);
					if(countz==pSQ->nprocz-1)
						ndpz_curr = pSQ->n_int[2]-(pSQ->nprocz-1)*ceil((double)pSQ->n_int[2]/pSQ->nprocz);

					// Start and End node indices of processor domain, w.r.t main domain. Assuming main domain indexing starts from zero. 
					pSQ->pnode_s[0] = countx*ndpx; pSQ->pnode_e[0] = pSQ->pnode_s[0]+ndpx_curr-1;
					pSQ->pnode_s[1] = county*ndpy; pSQ->pnode_e[1] = pSQ->pnode_s[1]+ndpy_curr-1;
					pSQ->pnode_s[2] = countz*ndpz; pSQ->pnode_e[2] = pSQ->pnode_s[2]+ndpz_curr-1;
					// printf("Processor rank:%u, xs=%u,xe=%u,ys=%u,ye=%u,zs=%u,ze=%u  \n",count,pSQ->pnode_s[0],pSQ->pnode_e[0],pSQ->pnode_s[1],pSQ->pnode_e[1],pSQ->pnode_s[2],pSQ->pnode_e[2]);

					pSQ->np_x = pSQ->pnode_e[0] - pSQ->pnode_s[0] + 1; // no. of nodes in x-direction of processor domain
					pSQ->np_y = pSQ->pnode_e[1] - pSQ->pnode_s[1] + 1; // no. of nodes in y-direction of processor domain
					pSQ->np_z = pSQ->pnode_e[2] - pSQ->pnode_s[2] + 1; // no. of nodes in z-direction of processor domain
				}
				count=count+1;
			}
		}      
	} // end for over all processors count

	// Cubical domain assumption
	if(pSQ->n_int[0]!=pSQ->n_int[1] || pSQ->n_int[2]!=pSQ->n_int[1] || pSQ->n_int[0]!=pSQ->n_int[2] )
	{
		printf("Error:Domain is not cubical. Current code can only perform the calculations for a cubical domain. \n");
		exit(0);
	}
	// All equivalent procs assumption
	if((pSQ->n_int[0] % pSQ->nprocx)!=0 || (pSQ->n_int[1] % pSQ->nprocy)!=0 || (pSQ->n_int[2] % pSQ->nprocz)!=0)
	{

		printf("Error: Cannot divide the nodes equally among processors. Current code needs equal nodes in all procs. \n");
		exit(0);
	}

}

// function to create communicator topologies
void Comm_topologies(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tcomm1,tcomm2;
	int reorder=0,j,i,k;
	int pdims[3]={pSQ->nprocz,pSQ->nprocy,pSQ->nprocx}; // no. of processors in each direction --> order is z,y,x
	int periodicity[3]={1,1,1}; // periodic BC's in each direction
	MPI_Comm topocomm; // declare a new communicator to give topology attribute
	tcomm1 = MPI_Wtime();
	MPI_Cart_create(MPI_COMM_WORLD,3,pdims,periodicity,reorder,&topocomm); // create Cartesian topology

	int pcoords[3];
	MPI_Cart_coords(topocomm,rank,3,pcoords); // local coordinate indices of processors
	int rank_chk;
	MPI_Cart_rank(topocomm,pcoords,&rank_chk); // proc rank corresponding to coords
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;

	//////////////////////////////////////////////////////
	//// Cartesian Topology///////////////////////////////
	// This topology has atleast 6 nearest neighbors for communication in 3D
	//////////////////////////////////////////////////////

	int nneigh = 6*ceil((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x); // total number of neighbors
	int *neighs,count=0;
	neighs = new int [nneigh](); // need to de-allocate later
	pSQ->neighs_lap = new int [nneigh](); // need to de-allocate later

	int proc_l,proc_r,proc_u,proc_d,proc_f,proc_b; // procs on left,right, up,down, front,back of the current proc
	for(j=0;j<ceil((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x);j++) // no. of layers of nearest neighbors required for FD stencil
	{
		tcomm1 = MPI_Wtime();
		MPI_Cart_shift(topocomm,0,j+1,&proc_l,&proc_r);  // x-direction
		MPI_Cart_shift(topocomm,1,j+1,&proc_b,&proc_f);  // y-direction
		MPI_Cart_shift(topocomm,2,j+1,&proc_d,&proc_u);  // z-direction
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;

		neighs[count]=proc_l; neighs[count+1]=proc_r;
		neighs[count+2]=proc_b; neighs[count+3]=proc_f;
		neighs[count+4]=proc_d; neighs[count+5]=proc_u;

		pSQ->neighs_lap[count]=proc_l; pSQ->neighs_lap[count+1]=proc_r;
		pSQ->neighs_lap[count+2]=proc_b; pSQ->neighs_lap[count+3]=proc_f;
		pSQ->neighs_lap[count+4]=proc_d; pSQ->neighs_lap[count+5]=proc_u;
		count = count+5 + 1;
	}
	MPI_Comm comm_dist_graph_cart; // communicator with cartesian distributed graph topology
	tcomm1 = MPI_Wtime();
	MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,nneigh,neighs,(int *)MPI_UNWEIGHTED,nneigh,neighs,(int *)MPI_UNWEIGHTED,MPI_INFO_NULL,reorder,&comm_dist_graph_cart); // creates a distributed graph topology (adjacent, cartesian)
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	pSQ->comm_laplacian=comm_dist_graph_cart;

	// de-allocate neighs
	delete [] neighs;

	//////////////////////////////////////////////////////
	//// Cubical Topology///////////////////////////////
	// This topology has atleast 26 nearest neighbors for communication in 3D 
	// (includes 6 faces, 12 edges and 8 corners)
	// Needed for coomunication in Spectral Quadrature to find Veff outside proc domain within Rcut distance which inturn is needed do sub Hamiltonian times vector
	//////////////////////////////////////////////////////
	int ncoords[3],ncoords_mapped[3]; // coordinates of neigh proc
	int nnproc = ceil((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x); // no. of neigh procs in each direction on one side
	nneigh = pow(2*nnproc+1,3); // total number of neighbors (2n+1)^3 -1, n = no. of neigh layers
	neighs = new int [nneigh](); // need to de-allocate later
	pSQ->neighs_sq = new int [nneigh](); // need to de-allocate later

	// Go over all the (2n+1)^3-1 neigh processors (-1 because excluding the current proc)
	// For each neigh proc, find the corresponding proc coordinate in the main domain (using periodic bc)
	// Convert the so obtained proc coordinates to proc rank using MPI_Cart_rank
	count=0;
	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{
				ncoords[0] = pcoords[0]+i;
				ncoords[1] = pcoords[1]+j;
				ncoords[2] = pcoords[2]+k;

				ncoords_mapped[0] = ncoords[0];
				ncoords_mapped[1] = ncoords[1];
				ncoords_mapped[2] = ncoords[2];

				// Map the coordinates of the proc back into the main domain, if they are outside
				// x-dir
				if(ncoords[0]<0) // neigh proc is outside main domain
				{
					if(((0-ncoords[0])%pSQ->nprocx)==0)
					{
						ncoords_mapped[0]=0;
					}else
					{
						ncoords_mapped[0]=pSQ->nprocx-((0-ncoords[0])%pSQ->nprocx);
					}

				}
				if(ncoords[0]>pSQ->nprocx-1) // neigh proc is outside main domain
				{
					if(((ncoords[0]-(pSQ->nprocx-1))%pSQ->nprocx)==0)
					{
						ncoords_mapped[0]=pSQ->nprocx-1;
					}else
					{
						ncoords_mapped[0]=((ncoords[0]-(pSQ->nprocx-1))%pSQ->nprocx)-1;
					}

				}

				// y-dir
				if(ncoords[1]<0) // neigh proc is outside main domain
				{
					if(((0-ncoords[1])%pSQ->nprocx)==0)
					{
						ncoords_mapped[1]=0;
					}else
					{
						ncoords_mapped[1]=pSQ->nprocx-((0-ncoords[1])%pSQ->nprocx);
					}
				}
				if(ncoords[1]>pSQ->nprocx-1) // neigh proc is outside main domain
				{
					if(((ncoords[1]-(pSQ->nprocx-1))%pSQ->nprocx)==0)
					{
						ncoords_mapped[1]=pSQ->nprocx-1;
					}else
					{
						ncoords_mapped[1]=((ncoords[1]-(pSQ->nprocx-1))%pSQ->nprocx)-1;
					}
				}

				// z-dir
				if(ncoords[2]<0) // neigh proc is outside main domain
				{
					if(((0-ncoords[2])%pSQ->nprocx)==0)
					{
						ncoords_mapped[2]=0;
					}else
					{
						ncoords_mapped[2]=pSQ->nprocx-((0-ncoords[2])%pSQ->nprocx);
					}
				}
				if(ncoords[2]>pSQ->nprocx-1) // neigh proc is outside main domain
				{
					if(((ncoords[2]-(pSQ->nprocx-1))%pSQ->nprocx)==0)
					{
						ncoords_mapped[2]=pSQ->nprocx-1;
					}else
					{
						ncoords_mapped[2]=((ncoords[2]-(pSQ->nprocx-1))%pSQ->nprocx)-1;
					}
				}

				tcomm1 = MPI_Wtime();
				MPI_Cart_rank(topocomm,ncoords_mapped,&rank_chk); // proc rank corresponding to ncoords_mapped
				tcomm2 = MPI_Wtime();
				pSQ->mpi_time+=tcomm2-tcomm1;
				neighs[count]=rank_chk;
				pSQ->neighs_sq[count] = neighs[count];

				count = count + 1;

			}
		}
	}

	MPI_Comm comm_dist_graph_sq; // communicator with cartesian distributed graph topology for SQ (cubical)

	tcomm1 = MPI_Wtime();
	MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,nneigh,neighs,(int *)MPI_UNWEIGHTED,nneigh,neighs,(int *)MPI_UNWEIGHTED,MPI_INFO_NULL,reorder,&comm_dist_graph_sq); // creates a distributed graph topology (adjacent, cartesian cubical)
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	pSQ->comm_sq=comm_dist_graph_sq;  
	pSQ->nneigh_sq = nneigh;

	// de-allocate neighs
	delete [] neighs;


}

// function to compute b (charge density)
void ChargeDensity(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tcomm1,tcomm2;
	double ***VlocJ,***VlocJ_tilda; // can modify the code to do it without storing VlocJ 3d array
	int xstart=-1,xend=-1,ystart=-1,yend=-1,zstart=-1,zend=-1;// limits of intersection domain, if there is intersection then all the indices should be greater than or equal to zero, since all of them are in main domain only
	int nxVloc,nyVloc,nzVloc; // no. of nodes in each direction of intersection+FDn region
	int i,j,k,i0,j0,k0,a;
	double r; // radial distance
	double *DVloc,*Duu; // derivative of pseudopotential
	double b_temp,rho_temp;
	pSQ->rc_tilda=0.5;
	double rc_tilda = pSQ->rc_tilda;  
	double Eself_tilda=0.0,Ncharge_tilda=0.0,b_tilda_temp,Vcon=0.0; 
	pSQ->Eself=0.0;
	pSQ->Ncharge=0.0; 
	pSQ->Nelectron=0.0; 

	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->b[k][j][i]=0.0;
				pSQ->b_tilda[k][j][i]=0.0;
				pSQ->Vc[k][j][i]=0.0;
			}
		}
	}   

	// Loop over all the atoms in main+rb and compute VlocJ for those atoms whose rb-domain intersects
	// Need to compute VlocJ on intersected+FDn domain to apply FD stencil
	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) // loop over atom types
	{
		DVloc = new double [pSQ->Psd[JJ_typ].size](); // need to de-allocate later
		getYD_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].rVloc,DVloc,pSQ->Psd[JJ_typ].size); //derivatives for spline

		Duu = new double [pSQ->Psd[JJ_typ].size](); // need to de-allocate later
		getYD_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].uu,Duu,pSQ->Psd[JJ_typ].size); //derivatives for spline

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
				VlocJ_tilda = new double** [nzVloc](); // need to de-allocate later
				VlocJ = new double** [nzVloc](); // need to de-allocate later
				if(VlocJ == NULL)
				{
					if(rank==0)
						cout << "Memory allocation failed in VlocJ"<< endl;
					exit(1);
				}
				for(k=0;k<nzVloc;k++)
				{
					VlocJ_tilda[k] = new double* [nyVloc](); // need to de-allocate later
					VlocJ[k] = new double* [nyVloc](); // need to de-allocate later
					if(VlocJ[k] == NULL)
					{
						if(rank==0)
							cout << "Memory allocation failed in VlocJ[k]"<< endl;
						exit(1);
					}

					for(j=0;j<nyVloc;j++)
					{
						VlocJ_tilda[k][j] = new double [nxVloc](); // need to de-allocate later
						VlocJ[k][j] = new double [nxVloc](); // need to de-allocate later
						if(VlocJ == NULL)
						{
							if(rank==0)
								cout << "Memory allocation failed in VlocJ[k][j]"<< endl;
							exit(1);
						}
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
							if(i>=pSQ->FDn && i<(nxVloc-pSQ->FDn) && j>=pSQ->FDn && j<(nyVloc-pSQ->FDn)  &&  k>=pSQ->FDn && k<(nzVloc-pSQ->FDn))
							{
								ispline_gen(pSQ->Psd[JJ_typ].RadialGrid, pSQ->Psd[JJ_typ].uu, pSQ->Psd[JJ_typ].size, &r,&rho_temp,1,Duu);
								pSQ->rho_at[k+k0-pSQ->pnode_s[2]][j+j0-pSQ->pnode_s[1]][i+i0-pSQ->pnode_s[0]]+=rho_temp;
								pSQ->Nelectron += rho_temp*pSQ->delta*pSQ->delta*pSQ->delta; // integral(rho)
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

				// Now compute b = -(1/4pi)*Laplacian(VlocJ) by using FD stencil on points in intersection region, b[k][j][i]
				for(k=zstart;k<=zend;k++)
				{
					for(j=ystart;j<=yend;j++)
					{
						for(i=xstart;i<=xend;i++)
						{

							b_temp = VlocJ[k-k0][j-j0][i-i0]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));
							b_tilda_temp = VlocJ_tilda[k-k0][j-j0][i-i0]*3*pSQ->coeff_lap[0]*(-1/(4*M_PI));
							for(a=1;a<=pSQ->FDn;a++)
							{
								b_temp += (VlocJ[k-k0][j-j0][i-i0-a] + VlocJ[k-k0][j-j0][i-i0+a] + VlocJ[k-k0][j-j0-a][i-i0] + VlocJ[k-k0][j-j0+a][i-i0] + VlocJ[k-k0-a][j-j0][i-i0] + VlocJ[k-k0+a][j-j0][i-i0])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 	
								b_tilda_temp += (VlocJ_tilda[k-k0][j-j0][i-i0-a] + VlocJ_tilda[k-k0][j-j0][i-i0+a] + VlocJ_tilda[k-k0][j-j0-a][i-i0] + VlocJ_tilda[k-k0][j-j0+a][i-i0] + VlocJ_tilda[k-k0-a][j-j0][i-i0] + VlocJ_tilda[k-k0+a][j-j0][i-i0])*pSQ->coeff_lap[a]*(-1/(4*M_PI)); 			 
							}


							pSQ->b[k-pSQ->pnode_s[2]][j-pSQ->pnode_s[1]][i-pSQ->pnode_s[0]] += b_temp; // bJ(x), i,j,k indices are w.r.t main domain hence subtract processor domains starting node to get local indexing
							pSQ->Eself += 0.5*b_temp*VlocJ[k-k0][j-j0][i-i0]*pSQ->delta*pSQ->delta*pSQ->delta; //0.5*bJ*VJ
							pSQ->Ncharge += b_temp*pSQ->delta*pSQ->delta*pSQ->delta; // integral(b)

							pSQ->b_tilda[k-pSQ->pnode_s[2]][j-pSQ->pnode_s[1]][i-pSQ->pnode_s[0]] += b_tilda_temp; // bJ(x), i,j,k indices are w.r.t main domain hence subtract processor domains starting node to get local indexing
							Eself_tilda += 0.5*b_tilda_temp*VlocJ_tilda[k-k0][j-j0][i-i0]*pSQ->delta*pSQ->delta*pSQ->delta; //0.5*bJ*VJ
							Ncharge_tilda += b_tilda_temp*pSQ->delta*pSQ->delta*pSQ->delta; // integral(b)


							pSQ->Vc[k-pSQ->pnode_s[2]][j-pSQ->pnode_s[1]][i-pSQ->pnode_s[0]] += VlocJ_tilda[k-k0][j-j0][i-i0] - VlocJ[k-k0][j-j0][i-i0];

						}
					}
				}

				// find Vcon at 0,0,0 node of main domain
				if(xstart<=0 && xend>=0 && ystart<=0 && yend>=0 && zstart<=0 && zend>=0)
				{
					Vcon += VlocJ_tilda[0-k0][0-j0][0-i0] - VlocJ[0-k0][0-j0][0-i0];
				}


				// de-allocate memory
				for(k=0;k<nzVloc;k++)
				{
					for(j=0;j<nyVloc;j++)
					{
						delete [] VlocJ[k][j];
						delete [] VlocJ_tilda[k][j];
					}
					delete [] VlocJ[k];
					delete [] VlocJ_tilda[k];
				}
				delete [] VlocJ;
				delete [] VlocJ_tilda;



			} // end for loop over replica atoms of current main atom (includes the main domain atom too)	  
		} // end for loop over main atoms of current type
		delete [] DVloc;
		delete [] Duu;
	} // end for loop over atom types

	if(pSQ->MDstepCount==1) // set rho = rho_guess = rho_at only for first MD step, for 2nd and 3rd, rho_guess will be previous rho, beyond 3rd, we do extrapolation (if pSQ->ChgExtrap==1)
	{
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					pSQ->rho[k][j][i]=pSQ->rho_at[k][j][i];
				}
			}
		}   
	}

	if(pSQ->restart_scf==1) // SCF restart using rho(x) from .restart file
	{
		RestartSCF(pSQ); 
	}

	// Charge extrapolation (at the end of 3rd MD step, we would have first extrapolation which can be used in 4th MD step) --> to get better rho_guess
	if(pSQ->MDstepCount>3 && pSQ->ChgExtrap==1) 
	{
		if(rank==0)
			printf("Using Charge extrapolation for rho_guess! \n");

		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					pSQ->rho[k][j][i]=pSQ->rho_at[k][j][i] + (pSQ->drho_new[k][j][i]); 
					if(pSQ->rho[k][j][i] < 0.0)
						pSQ->rho[k][j][i] = pow(10.0,-6.0) ;
				}
			}
		}  
	}

	pSQ->Nelectron=0.0;  
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->Nelectron += pSQ->rho[k][j][i]*pSQ->delta*pSQ->delta*pSQ->delta;
			}
		}
	}   

	double reduce_temp[5];
	reduce_temp[0]=pSQ->Eself;      reduce_temp[1]=pSQ->Ncharge;     reduce_temp[2]=pSQ->Nelectron;   
	reduce_temp[3]=Ncharge_tilda;   reduce_temp[4]=Eself_tilda;
	tcomm1 = MPI_Wtime();
	MPI_Allreduce(MPI_IN_PLACE, reduce_temp, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	pSQ->Eself=reduce_temp[0];      pSQ->Ncharge=reduce_temp[1];     pSQ->Nelectron=reduce_temp[2];   
	Ncharge_tilda=reduce_temp[3];   Eself_tilda=reduce_temp[4];

	// Scale rho so that integral(rho) matches integral(b)
	double Nelectron_chk=0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->rho[k][j][i] = fabs(pSQ->Ncharge/pSQ->Nelectron)*pSQ->rho[k][j][i];
				Nelectron_chk += pSQ->rho[k][j][i]*pSQ->delta*pSQ->delta*pSQ->delta;

			}
		}
	}

	pSQ->Nelectron=-pSQ->Ncharge;

	// scale b_tilda
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->b_tilda[k][j][i] = fabs(pSQ->Ncharge/Ncharge_tilda)*pSQ->b_tilda[k][j][i];
			}
		}
	}
	// shift Vc
	tcomm1 = MPI_Wtime();
	double V000 = pSQ->Vc[0][0][0]; //Vc[0][0][0];
	MPI_Bcast(&V000,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Vcon,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	double Ecorr=0.0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				if(pSQ->MDstepCount==1) // set phi_guess only for first MD step
				{pSQ->phi_guess[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=1;}
				pSQ->rhs[k][j][i]=pSQ->b[k][j][i]+pSQ->rho[k][j][i];
				pSQ->Vc[k][j][i] += -V000+Vcon;
				Ecorr += 0.5*(pSQ->b_tilda[k][j][i]+pSQ->b[k][j][i])*pSQ->Vc[k][j][i]*pow(pSQ->delta,3);
			}
		}
	}

	double Ecorr_global;
	tcomm1 = MPI_Wtime();
	MPI_Allreduce(&Ecorr, &Ecorr_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	pSQ->Ecorr = Ecorr_global + pSQ->Eself - Eself_tilda;
	if(pSQ->Correction == 1)
		if(rank==0)
		{
			printf("Energy Correction (Ha/atom) : %.14f \n",pSQ->Ecorr/pSQ->n_atm);
		}
}

// function to generate initial guess for Poisson solver
void AllocateArrays(DS_SQ* pSQ)
{
	int i,j,k;//,ctr=0;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  


	//////////////////////////////////////////////////////
	// allocate memory to store b(x) in domain
	pSQ->b_tilda = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->Vc = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->rho = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->rho_at = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->drho_new = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->drho = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->drho_dt = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->drho_2dt = new double** [pSQ->np_z](); // need to de-allocate later
	pSQ->b = new double** [pSQ->np_z](); // need to de-allocate later
	if(pSQ->b == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in pSQ_>b"<< endl;
		exit(1);
	}
	for(k=0;k<pSQ->np_z;k++)
	{
		pSQ->b_tilda[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->Vc[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->rho[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->rho_at[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->drho_new[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->drho[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->drho_dt[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->drho_2dt[k] = new double* [pSQ->np_y](); // need to de-allocate later
		pSQ->b[k] = new double* [pSQ->np_y](); // need to de-allocate later
		if(pSQ->b[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->b[k]"<< endl;
			exit(1);
		}

		for(j=0;j<pSQ->np_y;j++)
		{
			pSQ->b_tilda[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->Vc[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->rho[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->rho_at[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->drho_new[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->drho[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->drho_dt[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->drho_2dt[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			pSQ->b[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			if(pSQ->b[k][j] == NULL)
			{
				if(rank==0)
					cout << "Memory allocation failed in pSQ->b[k][j]"<< endl;
				exit(1);
			}
		}
	}
	//////////////////////////////////////////////////////////

	// allocate memory to store phi(x) in domain
	pSQ->phi_guess = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later  
	pSQ->phi = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later
	if(pSQ->phi_guess == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in pSQ->phi_guess"<< endl;
		exit(1);
	}
	for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++)
	{
		pSQ->phi_guess[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later
		pSQ->phi[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later
		if(pSQ->phi_guess[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->phi[k]"<< endl;
			exit(1);
		}

		for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++)
		{
			pSQ->phi_guess[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			pSQ->phi[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			if(pSQ->phi_guess[k][j] == NULL)
			{
				if(rank==0)
					cout << "Memory allocation failed in pSQ->phi[k][j]"<< endl;
				exit(1);
			}
		}
	}


	pSQ->rhs = new double** [pSQ->np_z](); // need to de-allocate later 
	for(k=0;k<pSQ->np_z;k++)
	{
		pSQ->rhs[k] = new double* [pSQ->np_y](); // need to de-allocate 
		for(j=0;j<pSQ->np_y;j++)
		{
			pSQ->rhs[k][j] = new double [pSQ->np_x](); // need to de-allocate later
		}
	}


	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				pSQ->phi_guess[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=1;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Initialize Vxc & Veff

	// allocate memory to store Vxc(rho(x)) in domain
	pSQ->Vxc = new double** [pSQ->np_z](); // need to de-allocate later  
	if(pSQ->Vxc == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in pSQ->Vxc"<< endl;
		exit(1);
	}
	for(k=0;k<pSQ->np_z;k++)
	{
		pSQ->Vxc[k] = new double* [pSQ->np_y](); // need to de-allocate later
		if(pSQ->Vxc[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Vxc[k]"<< endl;
			exit(1);
		}

		for(j=0;j<pSQ->np_y;j++)
		{
			pSQ->Vxc[k][j] = new double [pSQ->np_x](); // need to de-allocate later
			if(pSQ->Vxc[k][j] == NULL)
			{
				if(rank==0)
					cout << "Memory allocation failed in pSQ->Vxc[k][j]"<< endl;
				exit(1);
			}
		}
	}

	// allocate memory to store Vxc(rho(x)) in domain
	pSQ->Veff = new double** [pSQ->np_z+2*pSQ->nloc](); // need to de-allocate later  
	if(pSQ->Veff == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in pSQ->phi"<< endl;
		exit(1);
	}
	for(k=0;k<pSQ->np_z+2*pSQ->nloc;k++)
	{
		pSQ->Veff[k] = new double* [pSQ->np_y+2*pSQ->nloc](); // need to de-allocate later
		if(pSQ->Veff[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->phi[k]"<< endl;
			exit(1);
		}

		for(j=0;j<pSQ->np_y+2*pSQ->nloc;j++)
		{
			pSQ->Veff[k][j] = new double [pSQ->np_x+2*pSQ->nloc](); // need to de-allocate later
			if(pSQ->Veff[k][j] == NULL)
			{
				if(rank==0)
					cout << "Memory allocation failed in pSQ->phi[k][j]"<< endl;
				exit(1);
			}
		}
	}

	// SQ stuff (scf.cpp)
	int nnode = pSQ->np_z*pSQ->np_y*pSQ->np_x; 
	pSQ->chi = new double [nnode]();
	pSQ->zee = new double [nnode]();
	pSQ->Ci  = new double* [nnode]();
	pSQ->rho_pj = new double* [nnode]();
	for(k=0;k<nnode;k++)
	{
		pSQ->rho_pj[k]=new double [pSQ->npl+1]();
		pSQ->Ci[k]=new double [pSQ->npl+1]();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// forces stuff (forces.cpp)
	pSQ->fx = new double [pSQ->n_atm]();
	pSQ->fy = new double [pSQ->n_atm]();
	pSQ->fz = new double [pSQ->n_atm]();
	pSQ->forces = new double [3*pSQ->n_atm]();

	pSQ->fnlocx = new double [pSQ->n_atm]();
	pSQ->fnlocy = new double [pSQ->n_atm]();
	pSQ->fnlocz = new double [pSQ->n_atm]();

	pSQ->fcorrx = new double [pSQ->n_atm]();
	pSQ->fcorry = new double [pSQ->n_atm]();
	pSQ->fcorrz = new double [pSQ->n_atm]();

	pSQ->flocx = new double [pSQ->n_atm]();
	pSQ->flocy = new double [pSQ->n_atm]();
	pSQ->flocz = new double [pSQ->n_atm]();
}

// function to compute coordinates of replica atoms whose rc-domain intersects with the proc+Rcut domain
void Replica_atoms_Rcut(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int px,py,pz,mx,my,mz; // no. of integer multiples in +ve or -ve x,y,z directions
	int count;
	int inn,jnn,knn;
	double Atm_x,Atm_y,Atm_z;
	int xl,yl,zl,xr,yr,zr;  // rb-domain limits (l=left), (r=right)
	int xs_main,ys_main,zs_main,xe_main,ye_main,ze_main;  // processor domain limits (s=start), (e=end)
	int xstart=-1-pSQ->nloc,xend=-1-pSQ->nloc,ystart=-1-pSQ->nloc,yend=-1-pSQ->nloc,zstart=-1-pSQ->nloc,zend=-1-pSQ->nloc; // limits of intersection domain, if there is intersection then all the indices should be greater than or equal to -pSQ->nloc, since all of them are in main+Rcut domain only

	// compute coordinates of replica atoms in main+rb domain
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{

		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			// Loop over all replicas in main+Rcut+rc and find out if the processor+Rcut domain overlaps/intersects with the rc-domain of each of those replicas
			px = floor(((pSQ->domain[0]+pSQ->Rcut+pSQ->Psd[ii].rc - pSQ->Atm[ii].Rx[jj].main_atm)/pSQ->domain[0])+ TEMP_TOL); // no. of positive integer multiples in x-direction
			py = floor(((pSQ->domain[1]+pSQ->Rcut+pSQ->Psd[ii].rc - pSQ->Atm[ii].Ry[jj].main_atm)/pSQ->domain[1])+ TEMP_TOL); // no. of positive integer multiples in y-direction
			pz = floor(((pSQ->domain[2]+pSQ->Rcut+pSQ->Psd[ii].rc - pSQ->Atm[ii].Rz[jj].main_atm)/pSQ->domain[2])+ TEMP_TOL); // no. of positive integer multiples in z-direction
			mx = ceil(((-(pSQ->Rcut+pSQ->Psd[ii].rc) - pSQ->Atm[ii].Rx[jj].main_atm)/pSQ->domain[0])- TEMP_TOL); // no. of negative integer multiples in x-direction
			my = ceil(((-(pSQ->Rcut+pSQ->Psd[ii].rc) - pSQ->Atm[ii].Ry[jj].main_atm)/pSQ->domain[1])- TEMP_TOL); // no. of negative integer multiples in y-direction
			mz = ceil(((-(pSQ->Rcut+pSQ->Psd[ii].rc) - pSQ->Atm[ii].Rz[jj].main_atm)/pSQ->domain[2])- TEMP_TOL); // no. of negative integer multiples in z-direction

			count=0;
			for(int rr=mz;rr<=pz;rr++)
			{
				for(int qq=my;qq<=py;qq++)
				{
					for(int pp=mx;pp<=px;pp++)
					{

						Atm_x = pSQ->Atm[ii].Rx[jj].main_atm + pp*pSQ->domain[0]; // x-coordinate of replica atom
						Atm_y = pSQ->Atm[ii].Ry[jj].main_atm + qq*pSQ->domain[1]; // y-coordinate of replica atom
						Atm_z = pSQ->Atm[ii].Rz[jj].main_atm + rr*pSQ->domain[2]; // z-coordinate of replica atom

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

						//--------- If the rc-domain around the nearest node of this atom intersects the processor+Rcut domain, then store overlap info---------

						// find start & end nodes (w.r.t main domain) of intersection of (rc)-domain around nearest node and the processor+Rcut domain
						// limits of rc-domain
						xl = inn-ceil(pSQ->Psd[ii].rc/pSQ->delta); xr = inn+ceil(pSQ->Psd[ii].rc/pSQ->delta);
						yl = jnn-ceil(pSQ->Psd[ii].rc/pSQ->delta); yr = jnn+ceil(pSQ->Psd[ii].rc/pSQ->delta);
						zl = knn-ceil(pSQ->Psd[ii].rc/pSQ->delta); zr = knn+ceil(pSQ->Psd[ii].rc/pSQ->delta);

						// limits of processor+Rcut domain
						xs_main = pSQ->pnode_s[0]-pSQ->nloc; xe_main = pSQ->pnode_e[0]+pSQ->nloc; 
						ys_main = pSQ->pnode_s[1]-pSQ->nloc; ye_main = pSQ->pnode_e[1]+pSQ->nloc; 
						zs_main = pSQ->pnode_s[2]-pSQ->nloc; ze_main = pSQ->pnode_e[2]+pSQ->nloc;

						xstart=-1-pSQ->nloc;xend=-1-pSQ->nloc;ystart=-1-pSQ->nloc;yend=-1-pSQ->nloc;zstart=-1;zend=-1-pSQ->nloc;

						if(xl >= xs_main && xl <= xe_main)
							xstart=xl;
						else if(xs_main >= xl && xs_main <= xr)
							xstart=xs_main;

						if(xr >= xs_main && xr <= xe_main)
							xend=xr;
						else if(xe_main >= xl && xe_main <= xr)
							xend=xe_main;

						if(yl >= ys_main && yl <= ye_main)
							ystart=yl;
						else if(ys_main >= yl && ys_main <= yr)
							ystart=ys_main;

						if(yr >= ys_main && yr <= ye_main)
							yend=yr;
						else if(ye_main >= yl && ye_main <= yr)
							yend=ye_main;

						if(zl >= zs_main && zl <= ze_main)
							zstart=zl;
						else if(zs_main >= zl && zs_main <= zr)
							zstart=zs_main;

						if(zr >= zs_main && zr <= ze_main)
							zend=zr;
						else if(ze_main >= zl && ze_main <= zr)
							zend=ze_main;

						if((xstart!=-1-pSQ->nloc)&&(xend!=-1-pSQ->nloc)&&(ystart!=-1-pSQ->nloc)&&(yend!=-1-pSQ->nloc)&&(zstart!=-1-pSQ->nloc)&&(zend!=-1-pSQ->nloc))
						{// intersection between rb-domain and processor-domain exists
							count+=1;
						}

					}
				}
			}

			// Now allocate the replica atoms array and store the overlap info
			pSQ->Atm[ii].Rx[jj].n_replica_Rcut=count; // no. of replica's of jjth atoms (including main domain's jjth atom itself) whose rc-domain intersects processor+Rcut domain
			pSQ->Atm[ii].Rx[jj].ProcAtmRcut = new DS_ProcAtmRcut[pSQ->Atm[ii].Rx[jj].n_replica_Rcut]; // need to de-allocate later
			pSQ->Atm[ii].Ry[jj].ProcAtmRcut = new DS_ProcAtmRcut[pSQ->Atm[ii].Rx[jj].n_replica_Rcut]; // need to de-allocate later
			pSQ->Atm[ii].Rz[jj].ProcAtmRcut = new DS_ProcAtmRcut[pSQ->Atm[ii].Rx[jj].n_replica_Rcut]; // need to de-allocate later
			count=0;
			for(int rr=mz;rr<=pz;rr++)
			{
				for(int qq=my;qq<=py;qq++)
				{
					for(int pp=mx;pp<=px;pp++)
					{

						Atm_x = pSQ->Atm[ii].Rx[jj].main_atm + pp*pSQ->domain[0]; // x-coordinate of replica atom
						Atm_y = pSQ->Atm[ii].Ry[jj].main_atm + qq*pSQ->domain[1]; // y-coordinate of replica atom
						Atm_z = pSQ->Atm[ii].Rz[jj].main_atm + rr*pSQ->domain[2]; // z-coordinate of replica atom

						inn = round((Atm_x/pSQ->delta)); // nearest node index of replica atom w.r.t main domain (zero at origin)
						jnn = round((Atm_y/pSQ->delta));
						knn = round((Atm_z/pSQ->delta));


						// If an atom is almost at the center of two nodes then assign the nearest node to the left of the atom
						if(fabs(fabs(inn*pSQ->delta-Atm_x)-pSQ->delta*0.5)<TEMP_TOL)
							inn = floor((Atm_x/pSQ->delta) - 10*TEMP_TOL);
						if(fabs(fabs(jnn*pSQ->delta-Atm_y)-pSQ->delta*0.5)<TEMP_TOL)
							jnn = floor((Atm_y/pSQ->delta) - 10*TEMP_TOL);
						if(fabs(fabs(knn*pSQ->delta-Atm_z)-pSQ->delta*0.5)<TEMP_TOL)
							knn = floor((Atm_z/pSQ->delta) - 10*TEMP_TOL);

						//--------- If the rc-domain around the nearest node of this atom intersects the processor+Rcut domain, then store overlap info---------

						// find start & end nodes (w.r.t main domain) of intersection of (rc)-domain around nearest node and the processor+Rcut domain
						// limits of rc-domain
						xl = inn-ceil(pSQ->Psd[ii].rc/pSQ->delta); xr = inn+ceil(pSQ->Psd[ii].rc/pSQ->delta);
						yl = jnn-ceil(pSQ->Psd[ii].rc/pSQ->delta); yr = jnn+ceil(pSQ->Psd[ii].rc/pSQ->delta);
						zl = knn-ceil(pSQ->Psd[ii].rc/pSQ->delta); zr = knn+ceil(pSQ->Psd[ii].rc/pSQ->delta);

						// limits of processor+Rcut domain
						xs_main = pSQ->pnode_s[0]-pSQ->nloc; xe_main = pSQ->pnode_e[0]+pSQ->nloc; 
						ys_main = pSQ->pnode_s[1]-pSQ->nloc; ye_main = pSQ->pnode_e[1]+pSQ->nloc; 
						zs_main = pSQ->pnode_s[2]-pSQ->nloc; ze_main = pSQ->pnode_e[2]+pSQ->nloc;

						xstart=-1-pSQ->nloc;xend=-1-pSQ->nloc;ystart=-1-pSQ->nloc;yend=-1-pSQ->nloc;zstart=-1;zend=-1-pSQ->nloc;

						if(xl >= xs_main && xl <= xe_main)
							xstart=xl;
						else if(xs_main >= xl && xs_main <= xr)
							xstart=xs_main;

						if(xr >= xs_main && xr <= xe_main)
							xend=xr;
						else if(xe_main >= xl && xe_main <= xr)
							xend=xe_main;

						if(yl >= ys_main && yl <= ye_main)
							ystart=yl;
						else if(ys_main >= yl && ys_main <= yr)
							ystart=ys_main;

						if(yr >= ys_main && yr <= ye_main)
							yend=yr;
						else if(ye_main >= yl && ye_main <= yr)
							yend=ye_main;

						if(zl >= zs_main && zl <= ze_main)
							zstart=zl;
						else if(zs_main >= zl && zs_main <= zr)
							zstart=zs_main;

						if(zr >= zs_main && zr <= ze_main)
							zend=zr;
						else if(ze_main >= zl && ze_main <= zr)
							zend=ze_main;

						if((xstart!=-1-pSQ->nloc)&&(xend!=-1-pSQ->nloc)&&(ystart!=-1-pSQ->nloc)&&(yend!=-1-pSQ->nloc)&&(zstart!=-1-pSQ->nloc)&&(zend!=-1-pSQ->nloc))
						{// intersection between rb-domain and main-domain exists		      		      
							pSQ->Atm[ii].Rx[jj].ProcAtmRcut[count].coord = Atm_x;
							pSQ->Atm[ii].Ry[jj].ProcAtmRcut[count].coord = Atm_y;
							pSQ->Atm[ii].Rz[jj].ProcAtmRcut[count].coord = Atm_z;

							count+=1;
						}
					}
				}
			}

		} //end for over main atoms of each type
	} // end for over types
}



