/** \file sq.cpp
  \brief This file contains functions related to computing spectral quadrature components.

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

// function to compute sq comm indices
void SQ_Comm_Indices(DS_SQ* pSQ)
{
	int nneigh = pSQ->nneigh_sq,k; 
	// allocate memory to store arrays for MPI communication
	pSQ->SqInd.displs_send = new int [nneigh]();
	pSQ->SqInd.ncounts_send = new int [nneigh]();
	pSQ->SqInd.displs_recv = new int [nneigh]();
	pSQ->SqInd.ncounts_recv = new int [nneigh]();
	pSQ->SqInd.eout_s = new int* [3](); // need to de-allocate later
	pSQ->SqInd.eout_e = new int* [3](); // need to de-allocate later
	pSQ->SqInd.ein_s = new int* [3](); // need to de-allocate later
	pSQ->SqInd.ein_e = new int* [3](); // need to de-allocate later
	for(k=0;k<3;k++)
	{
		pSQ->SqInd.eout_s[k] = new int [nneigh](); // need to de-allocate later
		pSQ->SqInd.eout_e[k] = new int [nneigh](); // need to de-allocate later
		pSQ->SqInd.ein_s[k] = new int [nneigh](); // need to de-allocate later
		pSQ->SqInd.ein_e[k] = new int [nneigh](); // need to de-allocate later      
	}

	if(pSQ->SqInd.ein_e == NULL)
	{
		cout << "Memory allocation failed in pSQ->phi"<< endl;
		exit(1);
	}

	// Compute edge indices for MPI communication
	RcutRegionIndicesForSQ(pSQ, pSQ->SqInd.eout_s,pSQ->SqInd.eout_e, pSQ->SqInd.ein_s,pSQ->SqInd.ein_e,pSQ->SqInd.displs_send,pSQ->SqInd.displs_recv,pSQ->SqInd.ncounts_send,pSQ->SqInd.ncounts_recv) ; // MOVE THIS FUNCTION TO INITIALIZATION PART (outside MD)

}




// function to compute Rcut region indices
void RcutRegionIndicesForSQ(DS_SQ* pSQ, int **eout_s,int **eout_e, int **ein_s,int **ein_e,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv)  
{
	// eout: start(s) and end(e) node indices w.r.t proc domain of the outgoing buffer indices
	// ein : start(s) and end(e) node indices w.r.t proc domain of the incoming buffer indices

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//if(rank==0)
	//printf("HAVE TO MOVE THIS FUNCTION TO INITIALIZATION PART (outside MD) \n");

	int i,j,k,ii,jj,kk;
	int count=0;
	int edge_s[3],edge_e[3]; // start and end nodes of the edge region in 3 directions, indicies w.r.t local processor domain
	int nnproc=ceil((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x);
	int nneigh=pSQ->nneigh_sq;
	int **ein_s_temp,**ein_e_temp;
	ein_s_temp = new int* [3](); // need to de-allocate later
	ein_e_temp = new int* [3](); // need to de-allocate later
	for(k=0;k<3;k++)
	{
		ein_s_temp[k] = new int [nneigh](); // need to de-allocate later
		ein_e_temp[k] = new int [nneigh](); // need to de-allocate later      
	}


	// compute the start/end indices of the regions for outgoing buffers. Indices w.r.t proc domain
	count=0;
	displs_send[0]=0;

	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{

				if(i==0 && j==0 && k==0) // exclude current processor
				{
					ncounts_send[count]=0;
					if(count>0)
						displs_send[count]=displs_send[count-1]+ncounts_send[count-1]; // relative displacement of index in out going array 
				}else
				{

					// store neighbor information
					ncounts_send[count]=pSQ->np_x*pSQ->np_y*pSQ->np_z; // no. of nodes to be communicated to this neighbor
					//displs[edge_count]=edge_count*ncounts[edge_count]; // relative displacement of index in out going array 
					if(count>0)
						displs_send[count]=displs_send[count-1]+ncounts_send[count-1]; // relative displacement of index in out going array 

					// for each neigh proc, compute start and end nodes of the communication region of size nloc x np x np
					edge_s[0]=0;edge_s[1]=0;edge_s[2]=0; // initialize all start nodes to start node of proc domain i.e. zero
					edge_e[0]=pSQ->np_x-1;edge_e[1]=pSQ->np_x-1;edge_e[2]=pSQ->np_x-1; // initialize all end nodes to end node of proc domain i.e. np-1


					// Update the index info for outermost neigh layer regions where only partial regions of proc domain are needed to be communicated
					if(abs(i)==nnproc || abs(j)==nnproc || abs(k)==nnproc) // for outermost neigh layer, need only part of the domains
					{			

						// if an index is equal to nnproc then that direction need to be updated. Whether start or end node to be updated depends on whether the index is negative or positive



						if(i==nnproc) // for right neigh proc, x-dir
							edge_s[0]=pSQ->np_x-(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)+1-1; // update the start node for the edge region in the direction of proc dir which is only nloc width
						if(j==nnproc) // for right neigh proc, y-dir
							edge_s[1]=pSQ->np_x-(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)+1-1; // update the start node for the edge region in the direction of proc dir which is only nloc width
						if(k==nnproc) // for right neigh proc, z-dir
							edge_s[2]=pSQ->np_x-(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)+1-1; // update the start node for the edge region in the direction of proc dir which is only nloc width


						if(i==-nnproc) // for left neigh proc
							edge_e[0]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width
						if(j==-nnproc) // for left neigh proc
							edge_e[1]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width
						if(k==-nnproc) // for left neigh proc
							edge_e[2]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width



						ncounts_send[count]= (edge_e[0]-edge_s[0]+1)*(edge_e[1]-edge_s[1]+1)*(edge_e[2]-edge_s[2]+1); // update no. of nodes (only partial required)	     
						if(count>0)
							displs_send[count]=displs_send[count-1]+ncounts_send[count-1]; // relative displacement of index in out going array 

					}


					eout_s[0][count]=edge_s[0]; eout_e[0][count]=edge_e[0];
					eout_s[1][count]=edge_s[1]; eout_e[1][count]=edge_e[1];
					eout_s[2][count]=edge_s[2]; eout_e[2][count]=edge_e[2];


				}
				count += 1;
			}
		}
	}


	// Compute the ncount_recv and displs_recv arrays for the current processor. For "count" index, ncounts_recv[recv_ind]=ncounts_send[count]. Need to find "recv_ind" index relative to [i,j,k] proc, of the current processor[0,0,0]. ("count" is the index of [i,j,k] relative to [0,0,0], so we want to find the other way around, i.e. [0,0,0] w.r.t [i,j,k], i.e. [-i,-j,-k] w.r.t [0,0,0])
	int *mark_off,*send_ind_arry; // array used to assign "1" as the replica procs are assinged in recv buff, so that in the next round that marked off proc will not be counted to get the min distance from current proc
	mark_off = new int [pSQ->nneigh_sq]();
	send_ind_arry = new int [pSQ->nneigh_sq]();
	int send_ind,ccnt,rep_dist,ctr,rep_dist_old=0,rep_ind,rep_ind_old=0;
	count=0;
	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{

				if(i==0 && j==0 && k==0) // exclude current processor
				{
					ncounts_recv[count]=0; // do this if we are including current processor in the comm topology
				}else
				{

					// Need to find "send_ind" corresponding to the index of replica atom of current neigh[count] with closest relative distance from current proc [0,0,0]. Do this by going over all (2*nloc+1)^3 procs and finding the replica procs and then by finding the closest among them to [0,0,0] proc.




					ccnt=0;ctr=0;
					for(kk=-nnproc;kk<=nnproc;kk++) // no. of neigh procs in each direction
					{
						for(jj=-nnproc;jj<=nnproc;jj++) // no. of neigh procs in each direction
						{
							for(ii=-nnproc;ii<=nnproc;ii++) // no. of neigh procs in each direction
							{

								if(ii==0 && jj==0 && kk==0)
								{
								}else
								{
									if(pSQ->neighs_sq[ccnt]==pSQ->neighs_sq[count] && mark_off[ccnt]==0) // if we match current neighbor proc in the outer for loops, that means we are at a replica (or the original neigh proc) and the replica proc has not been used yet since its mark=0
									{
										rep_dist = (nnproc - ii) + (nnproc - jj)*(2*nnproc+1) + (nnproc - kk)*(2*nnproc+1)*(2*nnproc+1); // index current proc relative to replica proc i.e. [0,0,0] w.r.t [ii,jj,kk], i.e. [-ii,-jj,-kk] w.r.t [0,0,0]). Smallest rep_dist w.r.t current proc will be the first to be received by current proc. Hence the index of the smallest over rep_dist should go into ncounts_recv buffer
										rep_ind = ccnt;

										if(ctr==0)
										{
											rep_dist_old=rep_dist;
											rep_ind_old = rep_ind;
										}

										if(ctr>0)
										{
											if(rep_dist<rep_dist_old)
											{
												rep_dist_old=rep_dist;
												rep_ind_old=rep_ind;
											}
										}
										ctr +=1;				  
									}	      
								}			      

								ccnt +=1;
							}
						}
					}

					send_ind = rep_ind_old;
					mark_off[send_ind] = 1; // min distance proc marked off. Will not be considered next time.
					ncounts_recv[count]=ncounts_send[send_ind];
					send_ind_arry[count]=send_ind; // stores the proc index to which buffer will be sent from current proc in the order of "count". So when count=0, buffer of size ncounts_recv[0] will first be received from neigh proc send_ind_arry[0] whose corresponding proc rank is pSQ->neighs_sq[send_ind_arry[0]].

				}
				count += 1;
			}
		}
	}

	delete [] mark_off;



	// find displs_recv
	count=0;
	displs_recv[0]=0;
	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{

				if(count>0)
					displs_recv[count]=displs_recv[count-1]+ncounts_recv[count-1]; // relative displacement of index in out going array 
				count += 1;
			}
		}
	}


	// compute the start/end indices of the regions for incoming buffers. Indices w.r.t proc+nloc domain
	count=0;	
	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{

				if(i==0 && j==0 && k==0) // exclude current processor
				{
				}else
				{

					// Need to compute start/end nodes based on the recv_buffer indexing, WE RECEIVE DATA CORRESPONDING TO THE "send_ind" proc 
					// CAN USE START/END INDICES FROM THE OUTGOING BUFFER (of the "send_ind" proc)
					// For this we will first find start/end indices in the order of send buffer. After that loop again to arrange these indices according to recv buffer

					// for each neigh proc, compute start and end nodes of the communication region of size nloc x np x np
					edge_s[0]=pSQ->nloc-1+1;edge_s[1]=pSQ->nloc-1+1;edge_s[2]=pSQ->nloc-1+1; // initialize all start nodes to start node of proc domain i.e. nloc-1+1 (w.r.t proc+nloc)
					edge_e[0]=pSQ->np_x-1+pSQ->nloc;edge_e[1]=pSQ->np_x-1+pSQ->nloc;edge_e[2]=pSQ->np_x-1+pSQ->nloc; // initialize all end nodes to end node of proc domain i.e. np-1+nloc


					// for left proc (i.e i<0) need to do -abs(i) so use +i below
					edge_s[0]=edge_s[0]+i*(pSQ->np_x); // update the start node for the edge region in the direction of proc dir which is only FDn width
					edge_e[0]=edge_e[0]+i*(pSQ->np_x);

					edge_s[1]=edge_s[1]+j*(pSQ->np_x); // update the start node for the edge region in the direction of proc dir which is only FDn width
					edge_e[1]=edge_e[1]+j*(pSQ->np_x);

					edge_s[2]=edge_s[2]+k*(pSQ->np_x); // update the start node for the edge region in the direction of proc dir which is only FDn width
					edge_e[2]=edge_e[2]+k*(pSQ->np_x);

					// Update the index info for outermost neigh layer regions where only partial regions of proc domain are needed to be communicated
					if(abs(i)==nnproc || abs(j)==nnproc || abs(k)==nnproc) // for outermost neigh layer, need only part of the domains
					{			

						if(i==nnproc) // for right neigh proc, x-dir 
						{
							edge_s[0]=pSQ->np_x+pSQ->nloc+1-1+floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x; // update the start node for the edge region in the direction of proc dir which is only nloc width
							edge_e[0]=pSQ->np_x+2*pSQ->nloc-1;
						}
						if(j==nnproc) // for right neigh proc, y-dir
						{
							edge_s[1]=pSQ->np_x+pSQ->nloc+1-1+floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x; // update the start node for the edge region in the direction of proc dir which is only nloc width
							edge_e[1]=pSQ->np_x+2*pSQ->nloc-1;
						}
						if(k==nnproc) // for right neigh proc, z-dir
						{
							edge_s[2]=pSQ->np_x+pSQ->nloc+1-1+floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x; // update the start node for the edge region in the direction of proc dir which is only nloc width
							edge_e[2]=pSQ->np_x+2*pSQ->nloc-1;
						}

						if(i==-nnproc) // for left neigh proc
						{
							edge_e[0]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width
							edge_s[0]=0;
						}
						if(j==-nnproc) // for left neigh proc
						{
							edge_e[1]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width
							edge_s[1]=0;
						}
						if(k==-nnproc) // for left neigh proc
						{
							edge_e[2]=(pSQ->nloc-floor((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only nloc width
							edge_s[2]=0;
						}
					}

					ein_s_temp[0][count]=edge_s[0]; ein_e_temp[0][count]=edge_e[0];
					ein_s_temp[1][count]=edge_s[1]; ein_e_temp[1][count]=edge_e[1];
					ein_s_temp[2][count]=edge_s[2]; ein_e_temp[2][count]=edge_e[2];

				}
				count += 1;
			}
		}
	}

	// Now loop again and find the proper order of the indices for the recv buffer based on send_ind_arry that was used to set up ncounts_recv
	// compute the start/end indices of the regions for incoming buffers. Indices w.r.t proc+nloc domain
	count=0;	
	for(k=-nnproc;k<=nnproc;k++) // no. of neigh procs in each direction
	{
		for(j=-nnproc;j<=nnproc;j++) // no. of neigh procs in each direction
		{
			for(i=-nnproc;i<=nnproc;i++) // no. of neigh procs in each direction
			{

				if(i==0 && j==0 && k==0) // exclude current processor
				{
				}else
				{

					// Need to compute start/end nodes based on the recv_buffer indexing, WE RECEIVE DATA CORRESPONDING TO THE "send_ind" proc 
					// CAN USE START/END INDICES FROM THE OUTGOING BUFFER (of the "send_ind" proc)
					// For this we will first find start/end indices in the order of send buffer. After that loop again to arrange these indices according to recv buffer


					ein_s[0][count]=ein_s_temp[0][send_ind_arry[count]]; ein_e[0][count]=ein_e_temp[0][send_ind_arry[count]];
					ein_s[1][count]=ein_s_temp[1][send_ind_arry[count]]; ein_e[1][count]=ein_e_temp[1][send_ind_arry[count]];
					ein_s[2][count]=ein_s_temp[2][send_ind_arry[count]]; ein_e[2][count]=ein_e_temp[2][send_ind_arry[count]];

				}
				count += 1;
			}
		}
	}



	// de-allocate memory
	for(k=0;k<3;k++)
	{
		delete [] ein_s_temp[k];
		delete [] ein_e_temp[k];
	}
	delete [] ein_s_temp;
	delete [] ein_e_temp;

	delete [] send_ind_arry;

}

// function clenshaw-curtis components
void ClenshawCurtisSpectralQuadrature(DS_SQ* pSQ,int scf_iter)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k,ii,jj,kk,nq=0,nd=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	pSQ->mpi_time=0.0; 
	double tcomm1,tcomm2;

	double time4,time0,tcomm_mpi=0.0,tferm_mpi=0.0,teng_mpi=0.0,time1,time2,time3;
	time0 = MPI_Wtime();

	MPI_Comm comm_dist_graph_cart = pSQ->comm_sq; // communicator with cartesian distributed graph topology

	int np_edge = pow(2*pSQ->nloc+pSQ->np_x,3)-pSQ->np_x*pSQ->np_x*pSQ->np_x; // no. of nodes in the Rcut communication region 
	double *V_edge_in,*V_edge_out; // linear array to store input and output communication data/buffer to and from the processor
	V_edge_in = new double [np_edge](); // need to de-allocate later
	V_edge_out = new double [np_edge](); // need to de-allocate later

	int neigh_count,edge_count;
	int nnproc=ceil((double)(pSQ->nloc-TEMP_TOL)/pSQ->np_x);

	/// Communicate with neighboring processors to find out Veff in outside Rcut region of the processor domain

	// Store the edge data into outgoing buffer
	neigh_count=0;
	edge_count=0;
	for(kk=-nnproc;kk<=nnproc;kk++) // no. of neigh procs in each direction
	{
		for(jj=-nnproc;jj<=nnproc;jj++) // no. of neigh procs in each direction
		{
			for(ii=-nnproc;ii<=nnproc;ii++) // no. of neigh procs in each direction
			{


				if(ii==0 && jj==0 && kk==0) // exclude current processor
				{
				}else
				{

					// i,j,k are w.r.t proc domain
					// Now read Veff values in the edge region into the array
					for(k=pSQ->SqInd.eout_s[2][neigh_count];k<=pSQ->SqInd.eout_e[2][neigh_count];k++) // z-direction of edge region
					{
						for(j=pSQ->SqInd.eout_s[1][neigh_count];j<=pSQ->SqInd.eout_e[1][neigh_count];j++) // y-direction of edge region
						{
							for(i=pSQ->SqInd.eout_s[0][neigh_count];i<=pSQ->SqInd.eout_e[0][neigh_count];i++) // x-direction of edge region
							{
								V_edge_out[edge_count]=pSQ->Veff[k+pSQ->nloc][j+pSQ->nloc][i+pSQ->nloc];
								edge_count = edge_count + 1;      

							}
						}
					}

				}

				neigh_count = neigh_count + 1;
			}
		}
	}

	tcomm1 = MPI_Wtime();  
	// Using Neighborhood collective for local communication between processors
	MPI_Request request;
	if(pSQ->non_blocking==0)
	{
		// MPI_Neighbor_alltoall(phi_edge_out,np_edge,MPI_DOUBLE,phi_edge_in,np_edge,MPI_DOUBLE,comm_dist_graph_cart);
		MPI_Neighbor_alltoallv(V_edge_out,pSQ->SqInd.ncounts_send,pSQ->SqInd.displs_send,MPI_DOUBLE,V_edge_in,pSQ->SqInd.ncounts_recv,pSQ->SqInd.displs_recv,MPI_DOUBLE,comm_dist_graph_cart);
	}
	else
	{	    
		//MPI_Ineighbor_alltoall(phi_edge_out,np_edge,MPI_DOUBLE,phi_edge_in,np_edge,MPI_DOUBLE,comm_dist_graph_cart,&request); // non-blocking
		MPI_Ineighbor_alltoallv(V_edge_out,pSQ->SqInd.ncounts_send,pSQ->SqInd.displs_send,MPI_DOUBLE,V_edge_in,pSQ->SqInd.ncounts_recv,pSQ->SqInd.displs_recv,MPI_DOUBLE,comm_dist_graph_cart,&request); // non-blocking
	}


	// Make sure communication has finished and then proceed to next task. This is to be done when using non-blocking routine.
	if(pSQ->non_blocking==1)
		MPI_Wait(&request,MPI_STATUS_IGNORE);

	tcomm2 = MPI_Wtime();
	tcomm_mpi+=tcomm2-tcomm1;
	pSQ->mpi_time+=tcomm_mpi;

	// Store the incoming buffer into the outer Rcut region of Veff array
	neigh_count=0;
	edge_count=0;
	for(kk=-nnproc;kk<=nnproc;kk++) // no. of neigh procs in each direction
	{
		for(jj=-nnproc;jj<=nnproc;jj++) // no. of neigh procs in each direction
		{
			for(ii=-nnproc;ii<=nnproc;ii++) // no. of neigh procs in each direction
			{

				if(ii==0 && jj==0 && kk==0) // exclude current processor
				{
				}else
				{
					// i,j,k are w.r.t proc+nloc domain
					// Now read Veff values in the edge region from the array
					for(k=pSQ->SqInd.ein_s[2][neigh_count];k<=pSQ->SqInd.ein_e[2][neigh_count];k++) // z-direction of edge region
					{
						for(j=pSQ->SqInd.ein_s[1][neigh_count];j<=pSQ->SqInd.ein_e[1][neigh_count];j++) // y-direction of edge region
						{
							for(i=pSQ->SqInd.ein_s[0][neigh_count];i<=pSQ->SqInd.ein_e[0][neigh_count];i++) // x-direction of edge region
							{
								pSQ->Veff[k][j][i]=V_edge_in[edge_count];
								edge_count = edge_count + 1;      
							}
						}
					}
				}
				neigh_count = neigh_count + 1;
			}
		}
	}

	time1 = MPI_Wtime();

	double ***vec,***Hv,***t0,***t1,***t2,lambda_min, lambda_max;
	vec = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	Hv = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t0 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t1 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	t2 = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	for(k=0;k<2*pSQ->nloc+1;k++)
	{
		vec[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		Hv[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t0[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t1[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		t2[k] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		for(j=0;j<2*pSQ->nloc+1;j++)
		{
			vec[k][j]=new double [2*pSQ->nloc+1](); // need to de-allocate later
			Hv[k][j]=new double [2*pSQ->nloc+1](); // need to de-allocate later
			t0[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			t1[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
			t2[k][j] = new double [2*pSQ->nloc+1](); // need to de-allocate later
		}
	}

	/// loop over finite difference nodes in the processor domain
	/// Within the loop: First find min & max eigenvalues of Hsub for each finite difference node using Lanczos, Use the HsubTimesVec function for this.
	double lambda_min_MIN,lambda_max_MAX;

	if(fabs((double)0.25*pSQ->npl-round((double)0.25*pSQ->npl))>TEMP_TOL && rank==0)
	{printf("Error: npl should be a multiple of 4. \n");exit(1);}

	double time_l1,time_l2,tt_lanczos=0.0;
	double dotprod1,dotprod2,htemp;

	nd=0; // node count
	for(k=pSQ->nloc+1-1;k<=pSQ->np_z+pSQ->nloc-1;k++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(j=pSQ->nloc+1-1;j<=pSQ->np_y+pSQ->nloc-1;j++) // no. of proc nodes in each direction
		{
			for(i=pSQ->nloc+1-1;i<=pSQ->np_x+pSQ->nloc-1;i++) // no. of proc nodes in each direction
			{
				// Set-up initial guess vector for Lanczos
				for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
				{
					for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
					{
						for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
						{
							vec[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=1.0;

						}
					}
				}

				time_l1 = MPI_Wtime();
				LanczosAlgorithm(pSQ,vec,i,j,k,&lambda_min,&lambda_max,nd); // vec is the initial guess vector and i,j,k are the indices of node in proc domain to which the Hsub corresponds to and the indices are w.r.t proc+Rcut domain. lambda_min and lambda_max are the extremal eigenvalues of Hsub.

				time_l2 = MPI_Wtime();
				tt_lanczos+=(time_l2-time_l1);

				if(nd==0)
				{lambda_min_MIN=lambda_min; lambda_max_MAX=lambda_max;
				}else
				{
					if(lambda_min<lambda_min_MIN)
						lambda_min_MIN=lambda_min;
					if(lambda_max>lambda_max_MAX)
						lambda_max_MAX=lambda_max;
				}
				pSQ->chi[nd] = (lambda_max+lambda_min)/2;
				pSQ->zee[nd] = (lambda_max-lambda_min)/2;

				/// For each FD node, after getting min/max eigenvalues, loop over quadrature order to find Chebyshev expansion components, Use the HsubTimesVec function this.
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
						}
					}
				}
				htemp = t1[pSQ->nloc][pSQ->nloc][pSQ->nloc];

				// loop over quadrature order 0 to npl/4
				for(nq=0;nq<=round(0.25*pSQ->npl);nq++)
				{
					if(nq==0)
					{
						pSQ->rho_pj[nd][nq] = (double)(t0[pSQ->nloc][pSQ->nloc][pSQ->nloc]/(pow(pSQ->delta,3)));
					}else if(nq==1)
					{
						pSQ->rho_pj[nd][nq] = (double)(t1[pSQ->nloc][pSQ->nloc][pSQ->nloc]/(pow(pSQ->delta,3)));
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

									t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];
									t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];

								}
							}
						}

						pSQ->rho_pj[nd][nq] = (double)(t2[pSQ->nloc][pSQ->nloc][pSQ->nloc]/(pow(pSQ->delta,3)));

					}
				} // end loop over nq



				// loop over quadrature order 0.25*npl+1 to npl/2
				for(nq=round(0.25*pSQ->npl)+1;nq<=round(0.5*pSQ->npl);nq++)
				{
					dotprod1=0.0; dotprod2=0.0;
					HsubTimesVec(pSQ,t1,i,j,k,Hv); // Hv=Hsub*t1
					for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
					{
						for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
						{
							for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
							{
								t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=(2*(Hv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]-pSQ->chi[nd]*t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc])/pSQ->zee[nd]) - t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];		      

								dotprod1 += t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]*t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]; // T_2n
								dotprod2 += t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]*t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]; // T_2n-1

								t0[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];
								t1[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=t2[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];

							}
						}
					}

					pSQ->rho_pj[nd][nq] = (double)(t2[pSQ->nloc][pSQ->nloc][pSQ->nloc]/(pow(pSQ->delta,3)));

					pSQ->rho_pj[nd][2*nq] = (double)((2*dotprod1-1)/(pow(pSQ->delta,3)));
					pSQ->rho_pj[nd][2*nq-1] = (double)((2*dotprod2-htemp)/(pow(pSQ->delta,3)));

				} // end loop over nq

				nd = nd+1;

			}
		}
	}

	time2 = MPI_Wtime();

	///////////////////////////////////////////////////////////////////////////////////////////////
	//------------------------ Find out Fermi energy (lambda_f)----------------------------------//
	///////////////////////////////////////////////////////////////////////////////////////////////

	// communicate across procs to find global lambda_min and max to give as input guesses for Brent's algorithm

	int nproc=pSQ->nprocx*pSQ->nprocy*pSQ->nprocz,n_count=1;
	double *recv_data;
	recv_data = new double [nproc](); // need to de-allocate later 


	if(scf_iter<=1) // /////////////////////////
	{
		tcomm1 = MPI_Wtime();
		MPI_Allgather(&lambda_min_MIN,n_count,MPI_DOUBLE,recv_data,n_count,MPI_DOUBLE,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		tferm_mpi+=tcomm2-tcomm1;

		for(k=0;k<nproc;k++)
		{
			if(k==0)
			{
				lambda_min = recv_data[0];
			}else
			{
				if(recv_data[k]<lambda_min)
					lambda_min=recv_data[k];
			}

		}
		tcomm1 = MPI_Wtime();
		MPI_Allgather(&lambda_max_MAX,n_count,MPI_DOUBLE,recv_data,n_count,MPI_DOUBLE,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		tferm_mpi+=tcomm2-tcomm1;

		for(k=0;k<nproc;k++)
		{
			if(k==0)
			{
				lambda_max = recv_data[0];
			}else
			{
				if(recv_data[k]>lambda_max)
					lambda_max=recv_data[k];
			}

		}

	} // end if on scf_iter

	double lambda_f,*Ci,*di;
	Ci = new double [pSQ->npl+1](); // need to de-allocate later
	di = new double [pSQ->npl+1](); // need to de-allocate later

	if(pSQ->fermi_alg == 1)
		lambda_f = BrentsAlgorithm(pSQ,lambda_min-1,lambda_max+1,Ci,di,scf_iter,&tferm_mpi); // ///////////////////////// NOTE: USE THIS WHEN USING LOW TEMPERATURES and comment the if loop below
	else if(pSQ->fermi_alg == 2)
	{if(scf_iter<=1) // /////////////////////////
		{
			pSQ->lambda_f=0.2; 
			lambda_f = NewtonRaphson(pSQ,Ci,di,&tferm_mpi); 
		}else
		{
			lambda_f = NewtonRaphson(pSQ,Ci,di,&tferm_mpi);
		}
	}
	pSQ->mpi_time+=tferm_mpi;

	pSQ->lambda_f = lambda_f;
	time3 = MPI_Wtime();

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//--------------- Compute Electron density, Band structure energy and Entropy energy -----------------//
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	double bet0=pSQ->bet; //(27.211384523/((8.617343e-5)*pSQ->T)); // inverse of smearing (in 1/Ha)
	double lambda_fp,bet,rho_pc,ebs,ent;

	// loop over all the nodes in the proc domain
	pSQ->Ebs=0.0;
	pSQ->Eent=0.0;
	nd=0; // node count
	for(k=0;k<pSQ->np_z;k++) // no. of proc nodes in each direction, indices w.r.t proc domain
	{
		for(j=0;j<pSQ->np_y;j++) // no. of proc nodes in each direction
		{
			for(i=0;i<pSQ->np_x;i++) // no. of proc nodes in each direction
			{
				bet = bet0*pSQ->zee[nd];
				lambda_fp = (lambda_f-pSQ->chi[nd])/pSQ->zee[nd];

				// electron density
				ChebyshevCoeff(pSQ->npl,lambda_fp,bet,pSQ->Ci[nd],di,FermiDirac); // Need to save coeff in a separate global variable
				rho_pc = 0.0;
				for(jj=0;jj<pSQ->npl+1;jj++)
				{
					rho_pc = rho_pc + 2*pSQ->rho_pj[nd][jj]*pSQ->Ci[nd][jj];
				}
				pSQ->rho[k][j][i]=fabs(rho_pc);
				if(rho_pc<0 && fabs(rho_pc) > 1e-6)
				{printf("Error: Element of rho is negative at proc rank=%d and [i,j,k]=[%d %d %d] \n",rank,i,j,k);exit(1);}

				// band structure energy
				ChebyshevCoeff(pSQ->npl,lambda_fp,bet,Ci,di,UbandFunc); 
				ebs = 0.0;
				for(jj=0;jj<pSQ->npl+1;jj++)
				{
					ebs = ebs + 2*pSQ->rho_pj[nd][jj]*Ci[jj];
				}
				pSQ->Ebs = pSQ->Ebs + (pSQ->zee[nd]*ebs + pSQ->chi[nd]*rho_pc)*pow(pSQ->delta,3);

				// entropy energy
				ChebyshevCoeff(pSQ->npl,lambda_fp,bet,Ci,di,EentFunc); 
				ent = 0.0;
				for(jj=0;jj<pSQ->npl+1;jj++)
				{
					ent = ent + 2*pSQ->rho_pj[nd][jj]*Ci[jj];
				}
				pSQ->Eent = pSQ->Eent + (1/bet0)*ent*pow(pSQ->delta,3);

				nd=nd+1;
			}
		}
	}

	double engy_temp[2]={pSQ->Ebs,pSQ->Eent};
	tcomm1 = MPI_Wtime();
	MPI_Allreduce(MPI_IN_PLACE, engy_temp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	teng_mpi+=tcomm2-tcomm1;
	pSQ->mpi_time+=teng_mpi;
	pSQ->Ebs = engy_temp[0];
	pSQ->Eent = engy_temp[1];

	time4 = MPI_Wtime();

	if(rank==0)
		printf("Time stats (sec): MPI_Comm = %.2f (%.2f), Lanczos = %.2f, Cheb_Comps = %.2f, Fermi_energy calc = %.2f (%.2f), Rho_Ebs_Eent = %.2f (%.2f) \n",time1-time0,tcomm_mpi,tt_lanczos,time2-time1-tt_lanczos,time3-time2,tferm_mpi,time4-time3,teng_mpi);

	pSQ->sq_time +=time4-time0;

	delete [] V_edge_in; //de-allocate memory
	delete [] V_edge_out; //de-allocate memory

	for(k=0;k<2*pSQ->nloc+1;k++)
	{
		for(j=0;j<2*pSQ->nloc+1;j++)
		{
			delete [] vec[k][j];
			delete [] Hv[k][j];
			delete [] t0[k][j];
			delete [] t1[k][j];
			delete [] t2[k][j];
		}
		delete [] vec[k];
		delete [] Hv[k];
		delete [] t0[k];
		delete [] t1[k];
		delete [] t2[k];
	}
	delete [] vec;
	delete [] Hv;
	delete [] t0;
	delete [] t1;
	delete [] t2;
	delete [] recv_data;
	delete [] Ci;
	delete [] di;

}


// function to compute Hsub times a vector, without explicitly using Hsub matrix
void HsubTimesVec(DS_SQ* pSQ,double ***vec,int i,int j,int k,double ***Hv)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,ii,jj,kk,a,p,q,r;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);


	double ***vec_fd,temp;
	vec_fd = new double** [2*pSQ->FDn+2*pSQ->nloc+1](); // need to de-allocate later

	for(kk=0;kk<2*pSQ->FDn+2*pSQ->nloc+1;kk++)
	{
		vec_fd[kk] = new double* [2*pSQ->FDn+2*pSQ->nloc+1](); // need to de-allocate later
		for(jj=0;jj<2*pSQ->FDn+2*pSQ->nloc+1;jj++)
		{
			vec_fd[kk][jj]=new double [2*pSQ->FDn+2*pSQ->nloc+1](); // need to de-allocate later
		}
	}

	for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
			{
				p=kk+pSQ->nloc+pSQ->FDn;
				q=jj+pSQ->nloc+pSQ->FDn;
				r=ii+pSQ->nloc+pSQ->FDn;
				vec_fd[p][q][r] = vec[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];
			}
		}
	}


	double ***Vv;
	Vv = new double** [2*pSQ->nloc+1](); // need to de-allocate later

	for(kk=0;kk<2*pSQ->nloc+1;kk++)
	{
		Vv[kk] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		for(jj=0;jj<2*pSQ->nloc+1;jj++)
		{
			Vv[kk][jj]=new double [2*pSQ->nloc+1](); // need to de-allocate later
		}
	}

	VnlocsubTimesVec(pSQ,vec,i,j,k,Vv); //Vv = Vnclosub x vec

	for(kk=-pSQ->nloc;kk<=pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=-pSQ->nloc;jj<=pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=-pSQ->nloc;ii<=pSQ->nloc;ii++) // no. of proc nodes in each direction
			{
				p=kk+pSQ->nloc+pSQ->FDn;
				q=jj+pSQ->nloc+pSQ->FDn;
				r=ii+pSQ->nloc+pSQ->FDn;
				temp = vec_fd[p][q][r]*3*pSQ->coeff_lap[0];
				for(a=1;a<=pSQ->FDn;a++)
				{
					temp += (vec_fd[p][q][r-a] + vec_fd[p][q][r+a] + vec_fd[p][q-a][r] + vec_fd[p][q+a][r] + vec_fd[p-a][q][r] + vec_fd[p+a][q][r])*pSQ->coeff_lap[a]; 				 
				}

				Hv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]=-0.5*temp + pSQ->Veff[k+kk][j+jj][i+ii]*vec[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc]+ Vv[kk+pSQ->nloc][jj+pSQ->nloc][ii+pSQ->nloc];//-0.5Lapl+Veff+Vnloc

			}
		}
	}




	// de-allocate memory
	for(kk=0;kk<2*pSQ->FDn+2*pSQ->nloc+1;kk++)
	{
		for(jj=0;jj<2*pSQ->FDn+2*pSQ->nloc+1;jj++)
		{
			delete [] vec_fd[kk][jj];
		}
		delete [] vec_fd[kk];
	}
	delete [] vec_fd;
	for(kk=0;kk<2*pSQ->nloc+1;kk++)
	{
		for(jj=0;jj<2*pSQ->nloc+1;jj++)
		{
			delete [] Vv[kk][jj];
		}
		delete [] Vv[kk];
	}
	delete [] Vv;

}

void Vector2Norm_local(double ***vec,int nx,double *val)
{
	int ii,jj,kk;
	*val=0.0;
	for(kk=0;kk<=2*nx;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=0;jj<=2*nx;jj++) // no. of proc nodes in each direction
		{
			for(ii=0;ii<=2*nx;ii++) // no. of proc nodes in each direction
			{
				*val = *val + vec[kk][jj][ii]*vec[kk][jj][ii];
			}
		}
	}
	*val = sqrt(*val);
}

void VectorDotProduct_local(double ***v1,double ***v2,int nx,double *val)
{
	int ii,jj,kk;
	*val=0.0;
	for(kk=0;kk<=2*nx;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=0;jj<=2*nx;jj++) // no. of proc nodes in each direction
		{
			for(ii=0;ii<=2*nx;ii++) // no. of proc nodes in each direction
			{
				*val = *val + v1[kk][jj][ii]*v2[kk][jj][ii];
			}
		}
	}

}

void LanczosAlgorithm(DS_SQ* pSQ,double ***vkm1,int i,int j,int k,double *lambda_min,double *lambda_max,int nd)
{
	int rank,ii,jj,kk,max_iter=500;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double ***vk,***vkp1,val,*aa,*bb;
	int nn=pow(2*pSQ->nloc+1,3);
	aa = new double [nn](); // need to de-allocate later
	bb = new double [nn](); // need to de-allocate later
	vk = new double** [2*pSQ->nloc+1](); // need to de-allocate later
	vkp1 = new double** [2*pSQ->nloc+1](); // need to de-allocate later  
	for(kk=0;kk<2*pSQ->nloc+1;kk++)
	{
		vk[kk] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		vkp1[kk] = new double* [2*pSQ->nloc+1](); // need to de-allocate later
		for(jj=0;jj<2*pSQ->nloc+1;jj++)
		{
			vk[kk][jj]=new double [2*pSQ->nloc+1](); // need to de-allocate later
			vkp1[kk][jj]=new double [2*pSQ->nloc+1](); // need to de-allocate later
		}
	}

	int nx = pSQ->nloc;
	Vector2Norm_local(vkm1,nx,&val);
	for(kk=0;kk<=2*pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=0;jj<=2*pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=0;ii<=2*pSQ->nloc;ii++) // no. of proc nodes in each direction
			{
				vkm1[kk][jj][ii]=vkm1[kk][jj][ii]/val;
			}
		}
	}
	HsubTimesVec(pSQ,vkm1,i,j,k,vk); // vk=Hsub*vkm1. Here vkm1 is the vector and i,j,k are the indices of node in proc domain to which the Hsub corresponds to and the indices are w.r.t proc+Rcut domain
	VectorDotProduct_local(vkm1,vk,nx,&val);
	aa[0]=val;
	for(kk=0;kk<=2*pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=0;jj<=2*pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=0;ii<=2*pSQ->nloc;ii++) // no. of proc nodes in each direction
			{
				vk[kk][jj][ii]=vk[kk][jj][ii]-aa[0]*vkm1[kk][jj][ii];
			}
		}
	}
	Vector2Norm_local(vk,nx,&val);
	bb[0]=val;
	for(kk=0;kk<=2*pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
	{
		for(jj=0;jj<=2*pSQ->nloc;jj++) // no. of proc nodes in each direction
		{
			for(ii=0;ii<=2*pSQ->nloc;ii++) // no. of proc nodes in each direction
			{
				vk[kk][jj][ii]=vk[kk][jj][ii]/bb[0];
			}
		}
	}

	int count=0;
	double lmin_prev=0.0,lmax_prev=0.0,dl=1.0,dm=1.0;
	while(count<max_iter && (dl>pSQ->lanczos_tol || dm>pSQ->lanczos_tol))
	{
		HsubTimesVec(pSQ,vk,i,j,k,vkp1); // vkp1=Hsub*vk
		VectorDotProduct_local(vk,vkp1,nx,&val); // val=vk'*vkp1
		aa[count+1]=val;
		for(kk=0;kk<=2*pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
		{
			for(jj=0;jj<=2*pSQ->nloc;jj++) // no. of proc nodes in each direction
			{
				for(ii=0;ii<=2*pSQ->nloc;ii++) // no. of proc nodes in each direction
				{
					vkp1[kk][jj][ii]=vkp1[kk][jj][ii]-aa[count+1]*vk[kk][jj][ii]-bb[count]*vkm1[kk][jj][ii];
				}
			}
		}
		Vector2Norm_local(vkp1,nx,&val);
		bb[count+1]=val;
		for(kk=0;kk<=2*pSQ->nloc;kk++) // no. of proc nodes in each direction, indices w.r.t proc+Rcut domain
		{
			for(jj=0;jj<=2*pSQ->nloc;jj++) // no. of proc nodes in each direction
			{
				for(ii=0;ii<=2*pSQ->nloc;ii++) // no. of proc nodes in each direction
				{
					vkm1[kk][jj][ii] = vk[kk][jj][ii];
					vk[kk][jj][ii] = vkp1[kk][jj][ii]/bb[count+1];
				}
			}
		}

		// Eigendecompose the tridiagonal matrix

		TridiagEigenSolve(aa,bb,count+2,lambda_min,lambda_max);

		dl = fabs((*lambda_min)-lmin_prev);
		dm = fabs((*lambda_max)-lmax_prev);
		lmin_prev = *lambda_min;
		lmax_prev = *lambda_max;          
		count = count + 1;
	}

	*lambda_min -= pSQ->lanczos_tol; 
	*lambda_max += pSQ->lanczos_tol; 

	if(rank==0 && nd==0)
		printf("Lanczos took %d iterations. \n",count);

	if(count==max_iter)
		printf("Error: Lanczos exceeded max_iter. count=%d, dl=%f,dm=%f \n",count,dl,dm);

	for(kk=0;kk<2*pSQ->nloc+1;kk++)
	{
		for(jj=0;jj<2*pSQ->nloc+1;jj++)
		{
			delete [] vk[kk][jj];
			delete [] vkp1[kk][jj];
		}
		delete [] vk[kk];
		delete [] vkp1[kk];
	}
	delete [] vk;
	delete [] vkp1;

	delete [] aa;
	delete [] bb;

}

void TridiagEigenSolve(double *diag,double *subdiag,int n,double *lambda_min,double *lambda_max)
{

	int m,l,iter,i;
	double s,r,p,g,f,dd,c,b;

	// allocate memory
	double *d,*e; // d has diagonal and e has subdiagonal
	d = new double [n](); // need to de-allocate later
	e = new double [n](); // need to de-allocate later

	//create copy of diag and subdiag in d and e
	for(i=0;i<n;i++)
	{
		d[i] = diag[i];
		e[i] = subdiag[i];
	}

	//e has the subdiagonal elements //ignore last element(n-1) of e, make it zero
	e[n-1] = 0.0;

	for(l=0;l<=n-1;l++)
	{
		iter=0;
		do{
			for(m=l;m<=n-2;m++)
			{
				dd = fabs(d[m])+fabs(d[m+1]);
				if((double)(fabs(e[m])+dd) == dd) break;
			}
			if(m != l) {
				if(iter++ == 200) {  printf("Too many iterations in Tridiagonal solver\n");exit(1);}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=sqrt(g*g+1.0); // pythag
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); //temp=SIGN(r,g)=r*sign(g)
				s=c=1.0;
				p=0.0;

				for(i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1] = (r=sqrt(g*g+f*f));
					if (r == 0.0) {
						d[i+1] -=p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s +2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;

				}
				if(r==0.0 && i>=l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		}while(m!=l);
	}

	// go over the array d to find the smallest and largest eigenvalue
	//initial
	*lambda_min = d[0];
	*lambda_max = d[0];

	for(i=1;i<n;i++)
	{
		if(d[i]>*lambda_max)
		{
			*lambda_max =d[i];
		}
		else if(d[i]<*lambda_min)
		{
			*lambda_min =d[i];
		}

	}

	// free memory
	delete [] d;
	delete [] e;

	return;
}


double NeConstraint(DS_SQ* pSQ,double lambda_f,double *Ci,double *di,double *tferm_mpi) // number of electrons constraint function that can be used to find Fermi energy
{
	int k,j,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double g=0.0,fval,lambda_fp,bet;//,*Ci;
	double bet0=pSQ->bet; //(27.211384523/((8.617343e-5)*pSQ->T)); // inverse of smearing (in 1/Ha)

	double tcomm1,tcomm2;
	double tt0,tt1,tt2,ts0=0.0,ts1=0.0;

	for(k=0;k<pSQ->np_x*pSQ->np_y*pSQ->np_z;k++) // loop over nodes in the proc domain
	{
		lambda_fp = (lambda_f-pSQ->chi[k])/pSQ->zee[k];
		bet = bet0*pSQ->zee[k];
		tt0 = MPI_Wtime();

		ChebyshevCoeff(pSQ->npl,lambda_fp,bet,Ci,di,FermiDirac);

		tt1 = MPI_Wtime();
		ts0 += tt1-tt0;

		for(j=0;j<pSQ->npl+1;j++)
		{
			g = g + 2*pSQ->rho_pj[k][j]*Ci[j];
		}

		tt2 = MPI_Wtime();
		ts1 += tt2-tt1;

	}
	g = g*pow(pSQ->delta,3);

	tt0 = MPI_Wtime();

	tcomm1 = MPI_Wtime();
	double g_global;
	MPI_Allreduce(&g,&g_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	g=g_global;
	tcomm2 = MPI_Wtime();
	*tferm_mpi+=tcomm2-tcomm1;

	tt1 = MPI_Wtime();
	fval = g-(-pSQ->Ncharge);
	return fval;

}

void ChebyshevCoeff(int npl,double lambda_fp,double bet,double *Ci,double *d,double (*func)(double,double,double)) // compute chenbyshev coefficients Ci (length npl+1) to fit the function "func"
{

	// NOTE: THIS FUNCTION CURRENTLY USES DISCRETE ORTHOGONALITY PROPERTY OF CHEBYSHEV POLYNOMIALS WHICH ONLY GIVE "approximate" COEFFICIENTS. SO IF THE FUNCTION IS NOT SMOOTH ENOUGH THE COEFFICIENTS WILL NOT BE ACCURATE ENOUGH. SO THIS FUNCTION HAS TO BE REPLACED WITH ITS CONTINUOUS COUNTER PART WHICH EVALUATES OSCILLATORY INTEGRALS AND USES FFT.


	int k,j;
	double y,fac,sum;

	for(k=0;k<npl+1;k++)
	{
		y = cos(M_PI*((double)((k-0.5+1)/(npl+1))));
		d[k]=(*func)(y,lambda_fp,bet);
	}

	fac = 2.0/(npl+1);

	for(j=0;j<npl+1;j++)
	{
		sum = 0.0;
		for(k=0;k<npl+1;k++)
		{
			sum = sum + d[k]*cos((M_PI*(j-1+1))*((double)((k-0.5+1)/(npl+1))));
		}
		Ci[j]=fac*sum;
	}

	Ci[0] = Ci[0]/2.0;

}

double FermiDirac(double t,double lambda_f,double bet)
{
	double v;
	v = 1/(1+exp(bet*(t-lambda_f)));
	return v;
}

double UbandFunc(double t,double lambda_f,double bet)
{
	double v;
	v = t*(1/(1+exp(bet*(t-lambda_f))));
	return v;
}

double EentFunc(double t,double lambda_f,double bet)
{
	double v,focc;
	focc = (1/(1+exp(bet*(t-lambda_f))));
	if(fabs(focc)<0.01*TEMP_TOL || fabs(focc-1.0)<0.01*TEMP_TOL)
	{
		v=0.0;
	}else
	{
		v=(focc*log(focc) + (1-focc)*log(1-focc));
	}
	return v;
}

double BrentsAlgorithm(DS_SQ *pSQ,double x1,double x2,double *Ci,double *di,int scf_iter,double *tferm_mpi) // finds the root of the function NeConstraint()
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double EPS=2.5e-16,tol=pSQ->fermi_tol; //1e-6;
	int ITMAXBRENTS=100;
	int iter;
	double tol1q,eq;
	double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;

	if(scf_iter>1)
	{
		a = pSQ->lambda_f-2; b = pSQ->lambda_f+2;
	}
	double fa=NeConstraint(pSQ,a,Ci,di,tferm_mpi),fb=NeConstraint(pSQ,b,Ci,di,tferm_mpi),fc,p,q,r,s,tol1,xm;

	if((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0))
	{
		printf("Rank=%d, Error: Root must be bracketed in Brent's method. fa=%f, fb=%f\n",rank,fa,fb);
		exit(1);
	}
	fc=fb;
	for(iter=1;iter<=ITMAXBRENTS;iter++)
	{
		if((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0))
		{
			c=a;
			fc=fa;
			e=d=b-a;
		}

		if(fabs((double)fc)<fabs((double)fb))
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs((double)b)+0.5*tol; // convergence check
		xm=0.5*(c-b);
		if(fabs((double)xm)<=tol1 || fb==0.0)
		{
			if(rank==0)
				printf("Brent's algorithm took %d iterations. \n", iter);
			return b;
		}

		if(fabs((double)e)>=tol1 && fabs((double)fa) > fabs((double)fb))
		{
			// attempt inverse quadratic interpolation

			s = fb/fa;
			if(a==c)
			{
				p = 2.0*xm*s;
				q = 1.0-s;
			}else
			{
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if(p>0.0)	
			{
				//check whether in bounds
				q = -q;
			}
			p = fabs((double)p);

			tol1q =tol1*q;

			min1 = 3.0*xm*q - fabs((double)tol1q);
			eq = e*q;
			min2 = fabs((double)eq);
			if(2.0*p < (min1<min2 ? min1:min2))
			{	//accept interpolation
				e = d;
				d = p/q;
			}else
			{  //Bounds decreasing too slowly, use bisection
				d = xm;
				e = d;
			}
		}else
		{
			d=xm;
			e=d;
		}

		//move last best guess to a

		a = b;
		fa = fb;
		if(fabs((double)d) > tol1)
		{

			//Evaluate new trial root
			b +=d;
		}else
		{
			b += SIGN(tol1,xm);

		}
		fb = NeConstraint(pSQ,b,Ci,di,tferm_mpi);
	}

	printf("Rank=%d, Maximum iterations exceeded in brents root finding method...exiting\n",rank);
	exit(1);
	return 0.0;

}


double NeConstraint_derv(DS_SQ* pSQ,double lambda_f,double *Ci,double *di,double *fa,double *dfa,double *tferm_mpi) // number of electrons constraint function that can be used to find Fermi energy
{
	int k,j,rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double g=0.0,dg=0.0,lambda_fp,bet;
	double bet0=pSQ->bet; //(27.211384523/((8.617343e-5)*pSQ->T)); // inverse of smearing (in 1/Ha)

	double tcomm1,tcomm2;

	for(k=0;k<pSQ->np_x*pSQ->np_y*pSQ->np_z;k++) // loop over nodes in the proc domain
	{
		lambda_fp = (lambda_f-pSQ->chi[k])/pSQ->zee[k];
		bet = bet0*pSQ->zee[k];

		ChebyshevCoeff(pSQ->npl,lambda_fp,bet,Ci,di,FermiDirac);

		for(j=0;j<pSQ->npl+1;j++)
		{
			g = g + 2*pSQ->rho_pj[k][j]*Ci[j];
		}

		ChebyshevCoeff(pSQ->npl,lambda_fp,bet,Ci,di,FermiDirac_derv);

		for(j=0;j<pSQ->npl+1;j++)
		{
			dg = dg + 2*pSQ->rho_pj[k][j]*Ci[j];
		}

	}
	dg=dg*bet0; // since derivative is d(NeConstraint)/d(lambda_f) = 2*Tr(bet0*g*(I-g)), where g is Density matrix
	dg = dg*pow(pSQ->delta,3);
	g = g*pow(pSQ->delta,3);

	double temp[2]={g,dg};
	tcomm1 = MPI_Wtime();
	MPI_Allreduce(MPI_IN_PLACE,temp,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	*tferm_mpi+=tcomm2-tcomm1;

	*fa=temp[0]-(-pSQ->Ncharge);
	*dfa=temp[1];
	return 0.0;
}

double FermiDirac_derv(double t,double lambda_f,double bet)
{
	double v;
	v = 1/(1+exp(bet*(t-lambda_f)));
	v = v*(1-v); // derivative of density matrix w.r.t lambda_f
	return v;
}

double NewtonRaphson(DS_SQ *pSQ,double *Ci,double *di,double *tferm_mpi) // finds the root of the function NeConstraint()
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double tol=pSQ->fermi_tol,lambda_f_new=0.0,lambda_f_old=pSQ->lambda_f,fa,dfa; //1e-6;
	int ITMAX=200,count=0;

	double err=1.0;
	while(err>tol && count<ITMAX)
	{
		NeConstraint_derv(pSQ,lambda_f_old,Ci,di,&fa,&dfa,tferm_mpi); // just uses one allreduce per iteration
		lambda_f_new = lambda_f_old - (fa/dfa);
		err = fabs(lambda_f_new-lambda_f_old);
		lambda_f_old = lambda_f_new;      
		count+=1;
	}


	if(count==ITMAX)
	{
		if(rank==0)
			printf("WARNING: Newton-Raphson exceeded max (%d) iterations with err=%g. \n",count,err);
	}else
	{
		if(rank==0)
			printf("Newton-Raphson took %d iterations. \n",count);
	}

	return lambda_f_new;
}
