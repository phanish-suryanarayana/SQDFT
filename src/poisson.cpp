/** \file poisson.cpp
  \brief This file contains some of the functions needed to solve the Poisson's equation using AAJ method.


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

// function to solve Poisson's equation using AAJ
void PoissonSolver_AAJ(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tcomm1,tcomm2;
	pSQ->mpi_time=0.0; 
	double t_poiss0 = MPI_Wtime(),t_poiss1;
	double beta_aaj=pSQ->beta_aaj; //0.5; 
	int m_aaj=pSQ->m_aaj; //3;
	int p_aaj = pSQ->p_aaj; //1;
	int i,j,k,m;
	double ***phi_new,***phi_old,***phi_res,***phi_temp,*phi_res_vec,*phi_old_vec,*Xold,*Fold,**DX,**DF,*am_vec; 
	MPI_Comm comm_dist_graph_cart = pSQ->comm_laplacian; // communicator with cartesian distributed graph topology

	// allocate memory to store phi(x) in domain
	phi_new = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later
	phi_old = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later
	phi_res = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later
	phi_temp = new double** [pSQ->np_z+2*pSQ->FDn](); // need to de-allocate later
	if(phi_new == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in pSQ->phi"<< endl;
		exit(1);
	}
	for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++)
	{
		phi_new[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later
		phi_old[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later
		phi_res[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later
		phi_temp[k] = new double* [pSQ->np_y+2*pSQ->FDn](); // need to de-allocate later

		if(phi_new[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->phi[k]"<< endl;
			exit(1);
		}

		for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++)
		{
			phi_new[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			phi_old[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			phi_res[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			phi_temp[k][j] = new double [pSQ->np_x+2*pSQ->FDn](); // need to de-allocate later
			if(phi_new[k][j] == NULL)
			{
				if(rank==0)
					cout << "Memory allocation failed in pSQ->phi[k][j]"<< endl;
				exit(1);
			}
		}
	}

	// allocate memory to store phi(x) in domain
	phi_res_vec = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	phi_old_vec = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Xold = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Fold = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	am_vec = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 

	// allocate memory to store DX, DF history matrices
	DX = new double* [m_aaj](); // need to de-allocate later
	DF = new double* [m_aaj](); // need to de-allocate later
	if(DF == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in DF"<< endl;
		exit(1);
	}
	for(k=0;k<m_aaj;k++)
	{
		DX[k] = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later
		DF[k] = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later
		if(DF[k] == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in DF[k]"<< endl;
			exit(1);
		}
	}

	// STEPS TO FOLLOW:
	// First initialize guess phi_old
	// Next create a contiguous vector of the stencil points to send to neigh procs and receive in a similar vector the ghost point data for the edge stencil points
	// Before the communication (for non-blocking code) finishes, compute lap_phi for nodes that are inside proc domain, which are stencil width away from domain boundary
	// After the coomunication is done, using the ghost point data, evaluate lap_phi for the points in the stencil region
	// Find phi_new from lap_phi and phi_old


	// Initialize phi_old from phi_guess
	int ctr=0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=pSQ->phi_guess[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]; 
				phi_temp[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];
				ctr+=1;
			}
		}
	}

	// calculate norm of right hand side (required for relative residual)
	double rhs_norm=1.0,*rhs_vec;
	rhs_vec = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	ctr=0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				rhs_vec[ctr]=pSQ->rhs[k][j][i]; //rhs=rho+b
				ctr=ctr+1;
			}
		}
	}
	Vector2Norm(rhs_vec,pSQ->np_x*pSQ->np_y*pSQ->np_z,&rhs_norm);
	rhs_norm*=4*M_PI*fabs(-1/(3*pSQ->coeff_lap[0])); // need inv(D)*4*pi*(rho+b) for relres
	delete [] rhs_vec;
	double tempp=0.0, tt1, tt2;
	int max_iter=pSQ->poisson_maxiter; //1000;
	int iter=1;
	double tol = pSQ->poisson_tol;
	double res = tol+1;

	// begin while loop
	while (res>tol && iter<=max_iter)
	{
		PoissonResidual(pSQ,phi_temp,phi_res,iter,comm_dist_graph_cart,pSQ->LapInd.eout_s,pSQ->LapInd.eout_e, pSQ->LapInd.ein_s,pSQ->LapInd.ein_e, pSQ->LapInd.ereg_s,pSQ->LapInd.ereg_e,pSQ->LapInd.stencil_sign,pSQ->LapInd.edge_ind, pSQ->LapInd.displs_send,pSQ->LapInd.displs_recv,pSQ->LapInd.ncounts_send,pSQ->LapInd.ncounts_recv);

		// -------------- Update phi --------------- //
		ctr=0;
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					phi_new[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]+beta_aaj*phi_res[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];
					phi_res_vec[ctr]=phi_res[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];
					phi_old_vec[ctr]=phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];
					ctr=ctr+1;
				}
			}
		}
		Vector2Norm(phi_res_vec,pSQ->np_x*pSQ->np_y*pSQ->np_z,&res);
		res = res/rhs_norm; // relative residual

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//----------Store Residual & Iterate History----------//
		if(iter>1)
		{
			m = ((iter-2) % m_aaj)+1-1; //-1 because the index starts from 0
			ctr=0;
			for(k=0;k<pSQ->np_z;k++)
			{
				for(j=0;j<pSQ->np_y;j++)
				{
					for(i=0;i<pSQ->np_x;i++)
					{
						DX[m][ctr]=phi_old_vec[ctr]-Xold[ctr];
						DF[m][ctr]=phi_res_vec[ctr]-Fold[ctr];
						ctr=ctr+1;
					}
				}
			}

		}

		ctr=0;
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					Xold[ctr] = phi_old_vec[ctr];
					Fold[ctr] = phi_res_vec[ctr];		    
					ctr=ctr+1;
				}
			}
		}

		//----------Anderson update-----------//
		if(iter % p_aaj == 0 && iter>1)
		{
			tt1 = MPI_Wtime();
			AndersonExtrapolation(DX,DF,phi_res_vec,beta_aaj,m_aaj,pSQ->np_z*pSQ->np_y*pSQ->np_x,am_vec,&pSQ->mpi_time);
			tt2 = MPI_Wtime();
			tempp += tt2-tt1;	    

			ctr=0;
			for(k=0;k<pSQ->np_z;k++)
			{
				for(j=0;j<pSQ->np_y;j++)
				{
					for(i=0;i<pSQ->np_x;i++)
					{
						// phi_k+1 = phi_k+1 - dp (AAJ Update)
						phi_new[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=phi_new[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]-am_vec[ctr];
						ctr=ctr+1;
					}
				}
			}

		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// set phi_old = phi_new
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					if(res<=tol || iter==max_iter)
					{
						pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]; // store the solution of Poisson's equation
					}
					phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=phi_new[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];		    
				}
			}
		}

		for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++) // z-direction of interior region
		{
			for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++) // y-direction of interior region
			{
				for(i=0;i<pSQ->np_x+2*pSQ->FDn;i++) // x-direction of interior region
				{//i,j,k indices are w.r.t proc+FDn domain
					phi_temp[k][j][i]=phi_old[k][j][i];
				}
			}
		}

		iter=iter+1;
	} // end while loop

	t_poiss1 = MPI_Wtime();
	pSQ->poiss_time+=t_poiss1-t_poiss0;
	if(rank==0)
	{
		if(iter<max_iter && res<=tol)
			printf("Poisson solver (AAJ) converged!:  Iterations = %d, Relative Residual = %g, Time = %.4f (%.4f) sec\n",iter-1,res,t_poiss1-t_poiss0,pSQ->mpi_time);
		if(iter>=max_iter)
		{printf("WARNING: AAJ exceeded maximum iterations.\n");printf("Poisson solve (AAJ):  Iterations = %d, Residual = %g, Time = %.4f sec \n",iter-1,res,t_poiss1-t_poiss0);}
	}

	// shift phi since it can differ by a constant
	double phi_fix=0.0;
	double phi000=pSQ->phi[pSQ->FDn][pSQ->FDn][pSQ->FDn];
	tcomm1 = MPI_Wtime();
	MPI_Bcast(&phi000,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;
	pSQ->tpoiss_mpi+=pSQ->mpi_time;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn] = pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]+phi_fix-phi000;
			}
		}
	}

	// de-allocate memory
	for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++)
	{
		for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++)
		{
			delete [] phi_new[k][j];
			delete [] phi_old[k][j];
			delete [] phi_res[k][j];
			delete [] phi_temp[k][j];
		}
		delete [] phi_new[k];
		delete [] phi_old[k];
		delete [] phi_res[k];
		delete [] phi_temp[k];
	}
	delete [] phi_new;
	delete [] phi_old;
	delete [] phi_res;
	delete [] phi_temp;

	delete [] phi_res_vec;
	delete [] phi_old_vec;
	delete [] Xold;
	delete [] Fold;
	delete [] am_vec;

	for(k=0;k<m_aaj;k++)
	{
		delete [] DX[k];
		delete [] DF[k];
	}
	delete [] DX;
	delete [] DF;
}

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// function to compute comm indices and arrays
void Laplacian_Comm_Indices(DS_SQ* pSQ)
{
	// compute edge index information to compute residual to be used in the solver
	int nneigh = 6*ceil((double)pSQ->FDn-TEMP_TOL/pSQ->np_x),k;
	// allocate memory to store arrays for MPI communication
	pSQ->LapInd.displs_send = new int [nneigh]();
	pSQ->LapInd.ncounts_send = new int [nneigh]();
	pSQ->LapInd.displs_recv = new int [nneigh]();
	pSQ->LapInd.ncounts_recv = new int [nneigh]();
	pSQ->LapInd.eout_s = new int* [3](); // need to de-allocate later
	pSQ->LapInd.eout_e = new int* [3](); // need to de-allocate later
	pSQ->LapInd.ein_s = new int* [3](); // need to de-allocate later
	pSQ->LapInd.ein_e = new int* [3](); // need to de-allocate later
	pSQ->LapInd.stencil_sign = new int* [3](); // need to de-allocate later
	pSQ->LapInd.edge_ind = new int* [3](); // need to de-allocate later
	for(k=0;k<3;k++)
	{
		pSQ->LapInd.eout_s[k] = new int [nneigh](); // need to de-allocate later
		pSQ->LapInd.eout_e[k] = new int [nneigh](); // need to de-allocate later
		pSQ->LapInd.ein_s[k] = new int [nneigh](); // need to de-allocate later
		pSQ->LapInd.ein_e[k] = new int [nneigh](); // need to de-allocate later
		pSQ->LapInd.stencil_sign[k] = new int [pSQ->FDn*pSQ->FDn*pSQ->FDn*8 + pSQ->FDn*pSQ->FDn*(pSQ->np_x-2*pSQ->FDn)*12 + pSQ->FDn*(pSQ->np_x-2*pSQ->FDn)*(pSQ->np_x-2*pSQ->FDn)*6](); // need to de-allocate later
		pSQ->LapInd.edge_ind[k] = new int [pSQ->FDn*pSQ->FDn*pSQ->FDn*8 + pSQ->FDn*pSQ->FDn*(pSQ->np_x-2*pSQ->FDn)*12 + pSQ->FDn*(pSQ->np_x-2*pSQ->FDn)*(pSQ->np_x-2*pSQ->FDn)*6](); // need to de-allocate later
	}

	EdgeIndicesForPoisson(pSQ, pSQ->LapInd.eout_s,pSQ->LapInd.eout_e, pSQ->LapInd.ein_s,pSQ->LapInd.ein_e, pSQ->LapInd.ereg_s,pSQ->LapInd.ereg_e, pSQ->LapInd.stencil_sign,pSQ->LapInd.edge_ind, pSQ->LapInd.displs_send,pSQ->LapInd.displs_recv,pSQ->LapInd.ncounts_send,pSQ->LapInd.ncounts_recv); // MOVE THIS FUNCTION TO INITIALIZATION PART (outside MD)

}

// function to compute edge region indices
void EdgeIndicesForPoisson(DS_SQ* pSQ, int **eout_s,int **eout_e, int **ein_s,int **ein_e, int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv)  
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int i,j,k,ii,jj,kk;
	int np = pSQ->np_x; // no. of nodes in each direction of processor domain
	int edge_count=0,neigh_count;
	int proc_dir,proc_lr;
	int edge_s[3],edge_e[3]; // start and end nodes of the edge region in 3 directions, indicies w.r.t local processor domain
	int neigh_level,nneigh=ceil((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x);

	// Setup the outgoing array phi_edge_out with the phi values from edges of the proc domain. Order: First loop over procs x,y,z direc and then for each proc, x,y,z over nodes in the edge region
	edge_count=0;
	displs_send[0]=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{
				// store neighbor information
				ncounts_send[edge_count]=pSQ->np_x*pSQ->np_y*pSQ->np_z; // no. of nodes to be communicated to this neighbor
				//displs[edge_count]=edge_count*ncounts[edge_count]; // relative displacement of index in out going array 
				if(edge_count>0)
					displs_send[edge_count]=displs_send[edge_count-1]+ncounts_send[edge_count-1]; // relative displacement of index in out going array 

				// for each neigh proc, compute start and end nodes of the communication region of size FDn x np x np
				edge_s[0]=0;edge_s[1]=0;edge_s[2]=0; // initialize all start nodes to start node of proc domain i.e. zero
				edge_e[0]=np-1;edge_e[1]=np-1;edge_e[2]=np-1; // initialize all end nodes to end node of proc domain i.e. np-1

				if(neigh_level+1==nneigh) // for outermost neigh layer, need only part of the domains
				{			
					ncounts_send[edge_count]=pSQ->np_x*pSQ->np_x*(pSQ->FDn-floor((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x)*pSQ->np_x); // update no. of nodes (only partial required)
					if(edge_count>0)
						displs_send[edge_count]=displs_send[edge_count-1]+ncounts_send[edge_count-1]; // relative displacement of index in out going array 

					if(proc_lr==1) // for right neigh proc
						edge_s[proc_dir]=np-(pSQ->FDn-floor((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)+1-1; // update the start node for the edge region in the direction of proc dir which is only FDn width
					if(proc_lr==0) // for left neigh proc
						edge_e[proc_dir]=(pSQ->FDn-floor((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x)*pSQ->np_x)-1; // update the end node for the edge region in the direction of proc dir which is only FDn width
				}
				eout_s[0][edge_count]=edge_s[0]; eout_e[0][edge_count]=edge_e[0];
				eout_s[1][edge_count]=edge_s[1]; eout_e[1][edge_count]=edge_e[1];
				eout_s[2][edge_count]=edge_s[2]; eout_e[2][edge_count]=edge_e[2];

				edge_count += 1;
			}
		}
	}

	int **ein_s_temp,**ein_e_temp;
	ein_s_temp = new int* [3](); // need to de-allocate later
	ein_e_temp = new int* [3](); // need to de-allocate later
	for(k=0;k<3;k++)
	{
		ein_s_temp[k] = new int [6*nneigh](); // need to de-allocate later
		ein_e_temp[k] = new int [6*nneigh](); // need to de-allocate later      
	}


	int *mark_off,*send_ind_arry; // array used to assign "1" as the replica procs are assinged in recv buff, so that in the next round that marked off proc will not be counted to get the min distance from current proc
	mark_off = new int [6*nneigh]();
	send_ind_arry = new int [6*nneigh]();
	int send_ind,ccnt,rep_dist,ctr,rep_dist_old=0,rep_ind,rep_ind_old=0;


	// Store the incoming buffer data from phi_edge_in into the outer stencil regions of phi_old array
	edge_count=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{
				////////////////////////////////////////
				// Need to find "send_ind" corresponding to the index of replica atom of current neigh[count] with closest relative distance from current proc [0,0,0]. Do this by going over all (2*nloc+1)^3 procs and finding the replica procs and then by finding the closest among them to [0,0,0] proc.

				ccnt=0;ctr=0;
				for(kk=0;kk<nneigh;kk++) // neigh_level
				{
					for(jj=0;jj<3;jj++) // proc_dir
					{
						for(ii=0;ii<2;ii++) // proc_lr
						{


							if(pSQ->neighs_lap[ccnt]==pSQ->neighs_lap[edge_count] && mark_off[ccnt]==0) // if we match current neighbor proc in the outer for loops, that means we are at a replica (or the original neigh proc) and the replica proc has not been used yet since its mark=0
							{
								// This is for SQ ----> //rep_dist = (nnproc - ii) + (nnproc - jj)*(2*nnproc+1) + (nnproc - kk)*(2*nnproc+1)*(2*nnproc+1); // index current proc relative to replica proc i.e. [0,0,0] w.r.t [ii,jj,kk], i.e. [-ii,-jj,-kk] w.r.t [0,0,0]). Smallest rep_dist w.r.t current proc will be the first to be received by current proc. Hence the index of the smallest over rep_dist should go into ncounts_recv buffer

								rep_dist = kk*6 + jj*2 + (1-ii); // [ii,jj,kk] relative current proc would have been kk*6 + 2*jj + ii

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


							ccnt +=1;
						}
					}
				}

				send_ind = rep_ind_old;
				mark_off[send_ind] = 1; // min distance proc marked off. Will not be considered next time.
				ncounts_recv[edge_count]=ncounts_send[send_ind];
				send_ind_arry[edge_count]=send_ind; // stores the proc index to which buffer will be sent from current proc in the order of "count". So when count=0, buffer of size ncounts_recv[0] will first be received from neigh proc send_ind_arry[0] whose corresponding proc rank is pSQ->neighs_lap[send_ind_arry[0]].


				///////////////////////////////////////////

				edge_count += 1;
			}
		}
	}


	delete [] mark_off;


	// find displs_recv
	edge_count=0; displs_recv[0]=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{

				if(edge_count>0)
					displs_recv[edge_count]=displs_recv[edge_count-1]+ncounts_recv[edge_count-1]; // relative displacement of index in out going array 

				edge_count += 1;
			}
		}
	}

	////////////////////////////////////////
	// Store the incoming buffer data from phi_edge_in into the outer stencil regions of phi_old array
	edge_count=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{

				// for each neigh proc, compute start and end nodes of the communication region of size FDn x np x np
				edge_s[0]=pSQ->FDn-1+1;edge_s[1]=pSQ->FDn-1+1;edge_s[2]=pSQ->FDn-1+1; // initialize all start nodes to start node of proc domain i.e. FDn-1+1 (w.r.t proc+FDn)
				edge_e[0]=np-1+pSQ->FDn;edge_e[1]=np-1+pSQ->FDn;edge_e[2]=np-1+pSQ->FDn; // initialize all end nodes to end node of proc domain i.e. np-1+FDn


				if(neigh_level+1==nneigh) // for outermost neigh layer, need only part of the domain
				{
					if((proc_lr)==1) // for right neigh proc
					{
						edge_s[proc_dir]=np+pSQ->FDn+1-1+floor((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x)*pSQ->np_x; // update the start node for the edge region in the direction of proc dir which is only FDn width
						edge_e[proc_dir]=np+2*pSQ->FDn-1;
					}
					if((proc_lr)==0) // for left neigh proc
					{
						edge_e[proc_dir]=pSQ->FDn-1-floor((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x)*pSQ->np_x; // update the end node for the edge region in the direction of proc dir which is only FDn width
						edge_s[proc_dir]=0;
					}
				}else // inner neigh layers
				{
					if((proc_lr)==1) // for right neigh proc
					{
						edge_s[proc_dir]=edge_s[proc_dir]+(neigh_level+1)*(np); // update the start node for the edge region in the direction of proc dir which is only FDn width
						edge_e[proc_dir]=edge_e[proc_dir]+(neigh_level+1)*(np);
					}
					if((proc_lr)==0) // for left neigh proc
					{
						edge_e[proc_dir]=edge_e[proc_dir]-(neigh_level+1)*(np); // update the end node for the edge region in the direction of proc dir which is only FDn width
						edge_s[proc_dir]=edge_s[proc_dir]-(neigh_level+1)*(np);
					}
				}

				ein_s_temp[0][edge_count]=edge_s[0]; ein_e_temp[0][edge_count]=edge_e[0];
				ein_s_temp[1][edge_count]=edge_s[1]; ein_e_temp[1][edge_count]=edge_e[1];
				ein_s_temp[2][edge_count]=edge_s[2]; ein_e_temp[2][edge_count]=edge_e[2];
				edge_count += 1;
			}
		}
	}

	// Now loop again and find the proper order of the indices for the recv buffer based on send_ind_arry that was used to set up ncounts_recv
	// compute the start/end indices of the regions for incoming buffers. Indices w.r.t proc+FDn domain
	edge_count=0; 
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{

				// Need to compute start/end nodes based on the recv_buffer indexing, WE RECEIVE DATA CORRESPONDING TO THE "send_ind" proc 
				// CAN USE START/END INDICES FROM THE OUTGOING BUFFER (of the "send_ind" proc)
				// For this we will first find start/end indices in the order of send buffer. After that loop again to arrange these indices according to recv buffer


				ein_s[0][edge_count]=ein_s_temp[0][send_ind_arry[edge_count]]; ein_e[0][edge_count]=ein_e_temp[0][send_ind_arry[edge_count]];
				ein_s[1][edge_count]=ein_s_temp[1][send_ind_arry[edge_count]]; ein_e[1][edge_count]=ein_e_temp[1][send_ind_arry[edge_count]];
				ein_s[2][edge_count]=ein_s_temp[2][send_ind_arry[edge_count]]; ein_e[2][edge_count]=ein_e_temp[2][send_ind_arry[edge_count]];

				edge_count += 1;
			}
		}
	}


	//////////////////////////////////////

	if(pSQ->np_x > 2*pSQ->FDn)
	{
		// Indices for boundary regions to do remaining partial stencil correction -- 26 boundary regions = (6 faces,12 edges, 8 corners)

		edge_count=0;
		// loop over 6 faces = 3x2
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{
				edge_s[0]=0+pSQ->FDn;edge_s[1]=0+pSQ->FDn;edge_s[2]=0+pSQ->FDn; // initialize all start nodes to start node of proc domain-FDn i.e. zero+FDn
				edge_e[0]=np-1-pSQ->FDn;edge_e[1]=np-1-pSQ->FDn;edge_e[2]=np-1-pSQ->FDn; // initialize all end nodes to end node of proc domain-FDn i.e. np-1-FDn
				if(proc_lr==1) // for right neigh proc
				{
					edge_s[proc_dir]=np-pSQ->FDn+1-1; // update the start node for the edge region in the direction of proc dir which is only FDn width
					edge_e[proc_dir]=np-1;
				}

				if(proc_lr==0) // for left neigh proc
				{
					edge_e[proc_dir]=pSQ->FDn-1; // update the end node for the edge region in the direction of proc dir which is only FDn width
					edge_s[proc_dir]=0;
				}	  

				ereg_s[0][edge_count]=edge_s[0]; ereg_e[0][edge_count]=edge_e[0];
				ereg_s[1][edge_count]=edge_s[1]; ereg_e[1][edge_count]=edge_e[1];
				ereg_s[2][edge_count]=edge_s[2]; ereg_e[2][edge_count]=edge_e[2];

				edge_count += 1;
			}
		}

		// loop over 12 edges = 3x2x2 --> 3 long edge dirs (x,y,z) then 4 regs for each long edge dir --> these 4 regs split as left/right and up/down in order of x,y,z
		// long edge along x-direction
		ereg_s[0][6]=pSQ->FDn; ereg_e[0][6]=np-1-pSQ->FDn;	ereg_s[1][6]=0;           ereg_e[1][6]=pSQ->FDn-1;	ereg_s[2][6]=0;           ereg_e[2][6]=pSQ->FDn-1;
		ereg_s[0][7]=pSQ->FDn; ereg_e[0][7]=np-1-pSQ->FDn;	ereg_s[1][7]=np-pSQ->FDn; ereg_e[1][7]=np-1;	        ereg_s[2][7]=0;           ereg_e[2][7]=pSQ->FDn-1;
		ereg_s[0][8]=pSQ->FDn; ereg_e[0][8]=np-1-pSQ->FDn;	ereg_s[1][8]=0;           ereg_e[1][8]=pSQ->FDn-1;	ereg_s[2][8]=np-pSQ->FDn; ereg_e[2][8]=np-1;
		ereg_s[0][9]=pSQ->FDn; ereg_e[0][9]=np-1-pSQ->FDn;	ereg_s[1][9]=np-pSQ->FDn; ereg_e[1][9]=np-1;	        ereg_s[2][9]=np-pSQ->FDn; ereg_e[2][9]=np-1;
		// long edge along y-direction
		ereg_s[1][10]=pSQ->FDn; ereg_e[1][10]=np-1-pSQ->FDn;	ereg_s[0][10]=0;          ereg_e[0][10]=pSQ->FDn-1;	ereg_s[2][10]=0;          ereg_e[2][10]=pSQ->FDn-1;
		ereg_s[1][11]=pSQ->FDn; ereg_e[1][11]=np-1-pSQ->FDn;	ereg_s[0][11]=np-pSQ->FDn;ereg_e[0][11]=np-1;	        ereg_s[2][11]=0;          ereg_e[2][11]=pSQ->FDn-1;
		ereg_s[1][12]=pSQ->FDn; ereg_e[1][12]=np-1-pSQ->FDn;	ereg_s[0][12]=0;          ereg_e[0][12]=pSQ->FDn-1;	ereg_s[2][12]=np-pSQ->FDn;ereg_e[2][12]=np-1;
		ereg_s[1][13]=pSQ->FDn; ereg_e[1][13]=np-1-pSQ->FDn;	ereg_s[0][13]=np-pSQ->FDn;ereg_e[0][13]=np-1;	        ereg_s[2][13]=np-pSQ->FDn;ereg_e[2][13]=np-1;
		// long edge along z-direction
		ereg_s[2][14]=pSQ->FDn; ereg_e[2][14]=np-1-pSQ->FDn;	ereg_s[1][14]=0;           ereg_e[1][14]=pSQ->FDn-1;	ereg_s[0][14]=0;           ereg_e[0][14]=pSQ->FDn-1;
		ereg_s[2][15]=pSQ->FDn; ereg_e[2][15]=np-1-pSQ->FDn;	ereg_s[1][15]=np-pSQ->FDn; ereg_e[1][15]=np-1;	        ereg_s[0][15]=0;           ereg_e[0][15]=pSQ->FDn-1;
		ereg_s[2][16]=pSQ->FDn; ereg_e[2][16]=np-1-pSQ->FDn;	ereg_s[1][16]=0;           ereg_e[1][16]=pSQ->FDn-1;	ereg_s[0][16]=np-pSQ->FDn; ereg_e[0][16]=np-1;
		ereg_s[2][17]=pSQ->FDn; ereg_e[2][17]=np-1-pSQ->FDn;	ereg_s[1][17]=np-pSQ->FDn; ereg_e[1][17]=np-1;	        ereg_s[0][17]=np-pSQ->FDn; ereg_e[0][17]=np-1;


		// loop over 8 corners
		ereg_s[0][18]=0; ereg_e[0][18]=pSQ->FDn-1;	ereg_s[1][18]=0;           ereg_e[1][18]=pSQ->FDn-1;	ereg_s[2][18]=0;           ereg_e[2][18]=pSQ->FDn-1;
		ereg_s[0][19]=np-pSQ->FDn; ereg_e[0][19]=np-1;	ereg_s[1][19]=0;           ereg_e[1][19]=pSQ->FDn-1;	ereg_s[2][19]=0;           ereg_e[2][19]=pSQ->FDn-1;
		ereg_s[0][20]=0; ereg_e[0][20]=pSQ->FDn-1;	ereg_s[1][20]=np-pSQ->FDn; ereg_e[1][20]=np-1;     	ereg_s[2][20]=0;           ereg_e[2][20]=pSQ->FDn-1;
		ereg_s[0][21]=np-pSQ->FDn; ereg_e[0][21]=np-1;	ereg_s[1][21]=np-pSQ->FDn; ereg_e[1][21]=np-1;	        ereg_s[2][21]=0;           ereg_e[2][21]=pSQ->FDn-1;

		ereg_s[0][22]=0; ereg_e[0][22]=pSQ->FDn-1;	ereg_s[1][22]=0;           ereg_e[1][22]=pSQ->FDn-1;	ereg_s[2][22]=np-pSQ->FDn; ereg_e[2][22]=np-1;
		ereg_s[0][23]=np-pSQ->FDn; ereg_e[0][23]=np-1;	ereg_s[1][23]=0;           ereg_e[1][23]=pSQ->FDn-1;	ereg_s[2][23]=np-pSQ->FDn; ereg_e[2][23]=np-1;
		ereg_s[0][24]=0; ereg_e[0][24]=pSQ->FDn-1;	ereg_s[1][24]=np-pSQ->FDn; ereg_e[1][24]=np-1;     	ereg_s[2][24]=np-pSQ->FDn; ereg_e[2][24]=np-1;
		ereg_s[0][25]=np-pSQ->FDn; ereg_e[0][25]=np-1;	ereg_s[1][25]=np-pSQ->FDn; ereg_e[1][25]=np-1;	        ereg_s[2][25]=np-pSQ->FDn; ereg_e[2][25]=np-1;

		int temp=0;
		for(neigh_count=0;neigh_count<26;neigh_count++)
		{
			for(k=ereg_s[2][neigh_count];k<=ereg_e[2][neigh_count];k++) // z-direction of edge region
			{
				for(j=ereg_s[1][neigh_count];j<=ereg_e[1][neigh_count];j++) // y-direction of edge region
				{
					for(i=ereg_s[0][neigh_count];i<=ereg_e[0][neigh_count];i++) // x-direction of edge region
					{//i,j,k indices are w.r.t proc domain, but phi_old and new arrays are on proc+FDn domain

						if(i-0 <= np-1-i) // index i closer to left edge
						{
							stencil_sign[0][temp]=1; // full half stencil can be updated on right side
							edge_ind[0][temp]=min(pSQ->FDn,i-0);
						}else // index i closer to right edge
						{
							stencil_sign[0][temp]=-1; // full half stencil can be updated on left side
							edge_ind[0][temp]=min(pSQ->FDn,np-1-i);
						}

						if(j-0 <= np-1-j) // index i closer to left edge
						{
							stencil_sign[1][temp]=1; // full half stencil can be updated on right side
							edge_ind[1][temp]=min(pSQ->FDn,j-0);
						}else // index i closer to right edge
						{
							stencil_sign[1][temp]=-1; // full half stencil can be updated on left side
							edge_ind[1][temp]=min(pSQ->FDn,np-1-j);
						}

						if(k-0 <= np-1-k) // index i closer to left edge
						{
							stencil_sign[2][temp]=1; // full half stencil can be updated on right side
							edge_ind[2][temp]=min(pSQ->FDn,k-0);
						}else // index i closer to right edge
						{
							stencil_sign[2][temp]=-1; // full half stencil can be updated on left side
							edge_ind[2][temp]=min(pSQ->FDn,np-1-k);
						}
						temp+=1;
					}
				}
			}
		}


	} // end if condition (pSQ->np_x > 2*pSQ->FDn)



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

// function to compute residual of Jacobi iteration for Poisson equation
void PoissonResidual(DS_SQ* pSQ,double ***phi_old,double ***phi_res,int iter,MPI_Comm comm_dist_graph_cart, int **eout_s,int **eout_e, int **ein_s,int **ein_e, int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tcomm1,tcomm2;
	// STEPS TO FOLLOW:
	// Create a contiguous vector of the stencil points to send to neigh procs and receive in a similar vector the ghost point data for the edge stencil points
	// Before the communication (for non-blocking code) finishes, compute lap_phi for nodes that are inside proc domain, which are stencil width away from domain boundary
	// After the coomunication is done, using the ghost point data, evaluate lap_phi for the points in the stencil region
	// Find phi_new from lap_phi and phi_old

	// Assemble an array of phi at stencil points for communication.
	// This array should go over the list of neigh procs and accumulate the points, in the order x,y,z with left and rights in each direction
	int np = pSQ->np_x; // no. of nodes in each direction of processor domain
	//int Nd = pSQ->n_int[0]; // no. of nodes in each direction (same in all directions) of main domain
	int np_edge = pSQ->FDn*np*np; // no. of nodes in the stencil communication region across each face of processor domain. Six such regions communicate.
	int i,j,k;
	double *phi_edge_in,*phi_edge_out;//,***phi_temp; // linear array to store input and output communication data/buffer to and from the processor
	phi_edge_in = new double [6*np_edge](); // need to de-allocate later
	phi_edge_out = new double [6*np_edge](); // need to de-allocate later

	int edge_count,neigh_count;
	int proc_dir,proc_lr;
	double lap_phi_k;
	int a;
	int neigh_level,nneigh=ceil((double)(pSQ->FDn-TEMP_TOL)/pSQ->np_x);
	// Setup the outgoing array phi_edge_out with the phi values from edges of the proc domain. Order: First loop over procs x,y,z direc and then for each proc, x,y,z over nodes in the edge region
	edge_count=0;
	neigh_count=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{

				// Now read phi values in the edge region into the array
				for(k=eout_s[2][neigh_count];k<=eout_e[2][neigh_count];k++) // z-direction of edge region
				{
					for(j=eout_s[1][neigh_count];j<=eout_e[1][neigh_count];j++) // y-direction of edge region
					{
						for(i=eout_s[0][neigh_count];i<=eout_e[0][neigh_count];i++) // x-direction of edge region
						{
							phi_edge_out[edge_count]=phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];

							edge_count = edge_count + 1;      

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
		MPI_Neighbor_alltoallv(phi_edge_out,ncounts_send,displs_send,MPI_DOUBLE,phi_edge_in,ncounts_recv,displs_recv,MPI_DOUBLE,comm_dist_graph_cart);
	}
	else
	{	    
		//MPI_Ineighbor_alltoall(phi_edge_out,np_edge,MPI_DOUBLE,phi_edge_in,np_edge,MPI_DOUBLE,comm_dist_graph_cart,&request); // non-blocking
		MPI_Ineighbor_alltoallv(phi_edge_out,ncounts_send,displs_send,MPI_DOUBLE,phi_edge_in,ncounts_recv,displs_recv,MPI_DOUBLE,comm_dist_graph_cart,&request); // non-blocking
	}


	if(pSQ->np_x > 2*pSQ->FDn)
	{

		// Overlapping Computation with Communication when using non-blocking routines
		// Now find lap_phi and update phi using Jacobi iteration. Do this on interior and edge domains of proc domain separately so that we can use non-blocking communication

		// Update phi on the proc domain assuming it is zero outside the proc domain (will correct for this after stencil communication)
		for(k=0+pSQ->FDn;k<=np-1+pSQ->FDn;k++) // z-direction of interior region
		{
			for(j=0+pSQ->FDn;j<=np-1+pSQ->FDn;j++) // y-direction of interior region
			{
				for(i=0+pSQ->FDn;i<=np-1+pSQ->FDn;i++) // x-direction of interior region
				{//i,j,k indices are w.r.t proc+FDn domain
					lap_phi_k = phi_old[k][j][i]*3*pSQ->coeff_lap[0];
					for(a=1;a<=pSQ->FDn;a++)
					{
						lap_phi_k += (phi_old[k][j][i-a] + phi_old[k][j][i+a] + phi_old[k][j-a][i] + phi_old[k][j+a][i] + phi_old[k-a][j][i] + phi_old[k+a][j][i])*pSQ->coeff_lap[a]; 				 
					}

					phi_res[k][j][i] = -(((4*M_PI)*(pSQ->rhs[k-pSQ->FDn][j-pSQ->FDn][i-pSQ->FDn]) + lap_phi_k)/(3*pSQ->coeff_lap[0])); // Jacobi update for nodes in interior of proc domain
				}
			}
		}

	}

	// Make sure communication has finished and then proceed to next task. This is to be done when using non-blocking routine.
	if(pSQ->non_blocking==1)
		MPI_Wait(&request,MPI_STATUS_IGNORE);

	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;


	// Store the incoming buffer data from phi_edge_in into the outer stencil regions of phi_old array
	neigh_count=0;
	edge_count=0;
	for(neigh_level=0;neigh_level<nneigh;neigh_level++) // loop over layers of nearest neighbors
	{
		for(proc_dir=0;proc_dir<3;proc_dir++) // loop over directions (loop over neighbor procs) ---> 0,1,2=x,y,z direcs
		{
			for(proc_lr=0;proc_lr<2;proc_lr++) // loop over left & right (loop over neighbor procs) ----> 0,1=left,right
			{

				// Now read phi values in the edge region into the array
				for(k=ein_s[2][neigh_count];k<=ein_e[2][neigh_count];k++) // z-direction of edge region
				{
					for(j=ein_s[1][neigh_count];j<=ein_e[1][neigh_count];j++) // y-direction of edge region
					{
						for(i=ein_s[0][neigh_count];i<=ein_e[0][neigh_count];i++) // x-direction of edge region
						{
							phi_old[k][j][i]=phi_edge_in[edge_count];

							edge_count = edge_count + 1;

						}
					}
				}

				neigh_count = neigh_count + 1;

			}
		}

	}

	if(pSQ->np_x <= 2*pSQ->FDn)
	{


		// NOTE: IN THE CURRENT CODE WHEN TWICE STENCIL WIDTH IS EQUAL OR GREATER THAN PROC DOMAIN SIZE, WE DO NOT OVERALP COMPUTATION AND COMMUNICATION. THIS IS BECAUSE THE PARTIAL STENCIL UPDATE NEEDS TO DONE CAREFULLY AND DIFFERENTLY THAN THE OTHER CASE. SO WE JUST COMPUTE LAP_phi ON ENTIRE PROC DOMAIN AFTER THE COMMUNICATION IS COMPLETED.


		// Now find lap_phi and update phi using Jacobi iteration. Do this on interior and edge domains of proc domain separately so that we can use non-blocking communication

		// Update phi on the proc domain assuming it is zero outside the proc domain (will correct for this after stencil communication)
		for(k=0+pSQ->FDn;k<=np-1+pSQ->FDn;k++) // z-direction of interior region
		{
			for(j=0+pSQ->FDn;j<=np-1+pSQ->FDn;j++) // y-direction of interior region
			{
				for(i=0+pSQ->FDn;i<=np-1+pSQ->FDn;i++) // x-direction of interior region
				{//i,j,k indices are w.r.t proc+FDn domain
					lap_phi_k = phi_old[k][j][i]*3*pSQ->coeff_lap[0];
					for(a=1;a<=pSQ->FDn;a++)
					{
						lap_phi_k += (phi_old[k][j][i-a] + phi_old[k][j][i+a] + phi_old[k][j-a][i] + phi_old[k][j+a][i] + phi_old[k-a][j][i] + phi_old[k+a][j][i])*pSQ->coeff_lap[a]; 				 
					}

					phi_res[k][j][i] = -(((4*M_PI)*(pSQ->rhs[k-pSQ->FDn][j-pSQ->FDn][i-pSQ->FDn]) + lap_phi_k)/(3*pSQ->coeff_lap[0])); // Jacobi update for nodes in interior of proc domain
				}
			}
		}

	}

	if(pSQ->np_x > 2*pSQ->FDn)
	{
		int temp=0;
		for(neigh_count=0;neigh_count<26;neigh_count++)
		{
			for(k=ereg_s[2][neigh_count];k<=ereg_e[2][neigh_count];k++) // z-direction of edge region
			{
				for(j=ereg_s[1][neigh_count];j<=ereg_e[1][neigh_count];j++) // y-direction of edge region
				{
					for(i=ereg_s[0][neigh_count];i<=ereg_e[0][neigh_count];i++) // x-direction of edge region
					{//i,j,k indices are w.r.t proc domain, but phi_old and new arrays are on proc+FDn domain
						lap_phi_k = 0.0; //phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]*3*pSQ->coeff_lap[0];


						//if(edge_ind[0][temp] < pSQ->FDn)
						for(a=(edge_ind[0][temp]+1);a<=pSQ->FDn;a++) // x-direction
						{
							lap_phi_k += (phi_old[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn-stencil_sign[0][temp]*a])*pSQ->coeff_lap[a]; // part of the other half stencil inside domain
						}


						//if(edge_ind[1][temp] < pSQ->FDn)
						for(a=(edge_ind[1][temp]+1);a<=pSQ->FDn;a++) // y-direction
						{
							lap_phi_k += (phi_old[k+pSQ->FDn][j+pSQ->FDn-stencil_sign[1][temp]*a][i+pSQ->FDn])*pSQ->coeff_lap[a]; // part of the other half stencil inside domain
						}

						//if(edge_ind[2][temp] < pSQ->FDn)
						for(a=(edge_ind[2][temp]+1);a<=pSQ->FDn;a++) // z-direction
						{
							lap_phi_k += (phi_old[k+pSQ->FDn-stencil_sign[2][temp]*a][j+pSQ->FDn][i+pSQ->FDn])*pSQ->coeff_lap[a]; // part of the other half stencil inside domain
						}


						phi_res[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn] += -((lap_phi_k)/(3*pSQ->coeff_lap[0])); // Jacobi update for nodes in edge region of proc domain
						temp+=1;
					}
				}
			}
		}

	}

	delete [] phi_edge_in; //de-allocate memory
	delete [] phi_edge_out; //de-allocate memory

}



// function to compute 2-norm of a vector
void Vector2Norm(double* Vec, int len, double* ResVal)    // Vec is the pointer to the vector, len is the length of the vector Vec, and ResVal is the pointer to the residual value
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double res=0;
	for(int k=0;k<len;k++)
	{
		res = res + Vec[k]*Vec[k];
	}

	double res_global;
	MPI_Allreduce(&res, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	res_global = sqrt(res_global);

	*ResVal = res_global;

}


//function to compute time taken for Poisson solver
void AvgTime_Poisson(double dt, int nproc) // dt=time taken by current proc
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//int n_count=1;
	double *recv_data;
	recv_data = new double [nproc](); // need to de-allocate later
	recv_data[rank]=dt;
	MPI_Allreduce(MPI_IN_PLACE, recv_data, nproc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double ttemp_max=0.0,ttemp_min=0.0,ttemp_avg=0.0;
	if(rank==0)
	{
		for(int k=0;k<nproc;k++)
		{
			ttemp_avg = ttemp_avg + recv_data[k];
			if(recv_data[k]>ttemp_max)
				ttemp_max=recv_data[k];
		}
	}
	ttemp_avg=ttemp_avg/nproc;
	ttemp_min=ttemp_max;
	if(rank==0)
	{
		for(int k=0;k<nproc;k++)
		{
			if(recv_data[k]<ttemp_min)
				ttemp_min=recv_data[k];
		}
	}

	if(rank==0)
	{printf("Poisson solve time : Avg=%f sec, Max=%f sec, Min=%f sec \n",ttemp_avg,ttemp_max,ttemp_min); cout << " "<<endl;}

	if(rank==0)
	{printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
		cout << " "<<endl;}

	delete [] recv_data;

}





