/** \file deallocate.cpp
  \brief This file contains the function to de-allocate memory.


*/
#include "headers.h"
#include "func_sq.h"

// function to de-allocate memory
void Deallocate_memory(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int j,k,cnt,l,m;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		delete [] pSQ->Psd[ii].RadialGrid; 
		delete [] pSQ->Psd[ii].rVs;  
		delete [] pSQ->Psd[ii].rVp; 
		delete [] pSQ->Psd[ii].rVd; 
		delete [] pSQ->Psd[ii].rVf;
		delete [] pSQ->Psd[ii].rUs;      
		delete [] pSQ->Psd[ii].rUp; 
		delete [] pSQ->Psd[ii].rUd; 
		delete [] pSQ->Psd[ii].rUf;
		delete [] pSQ->Psd[ii].rVloc;
		delete [] pSQ->Psd[ii].uu;


		delete [] pSQ->Psd[ii].DVloc; 
		delete [] pSQ->Psd[ii].DVsJ;
		delete [] pSQ->Psd[ii].DVpJ;
		delete [] pSQ->Psd[ii].DVdJ;
		delete [] pSQ->Psd[ii].DVfJ; 
		delete [] pSQ->Psd[ii].DUsJ;
		delete [] pSQ->Psd[ii].DUpJ;
		delete [] pSQ->Psd[ii].DUdJ;
		delete [] pSQ->Psd[ii].DUfJ;    

		// need to delete replicas first
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			cnt = 0;
			for(l=0;l<=pSQ->Psd[ii].lmax;l++) // loop over quantum number l
			{
				if(l!=pSQ->Atm[ii].lloc)
				{  
					for(m=-l;m<=l;m++) // loop over quantum number m
					{
						cnt += 1;
					}
				}
			}

			// non-local Rcut stuff

			for(int JJr=0;JJr<pSQ->Atm[ii].Rx[jj].n_replica_Rcut;JJr++)
			{
				for(int kk=0;kk<cnt;kk++)
				{
					delete [] pSQ->Atm[ii].Rx[jj].ProcAtmRcut[JJr].UdVtm[kk];; 
				}
				delete [] pSQ->Atm[ii].Rx[jj].ProcAtmRcut[JJr].UdVtm;
			}


			delete [] pSQ->Atm[ii].Rx[jj].ProcAtmRcut;
			delete [] pSQ->Atm[ii].Ry[jj].ProcAtmRcut;
			delete [] pSQ->Atm[ii].Rz[jj].ProcAtmRcut;

			delete [] pSQ->Atm[ii].Rx[jj].ProcAtm;
			delete [] pSQ->Atm[ii].Ry[jj].ProcAtm;
			delete [] pSQ->Atm[ii].Rz[jj].ProcAtm;
		}

		delete [] pSQ->Atm[ii].Rx;  
		delete [] pSQ->Atm[ii].Ry;
		delete [] pSQ->Atm[ii].Rz;              
	}

	delete [] pSQ->Atm;
	delete [] pSQ->Psd;
	delete [] pSQ->coeff_lap;
	delete [] pSQ->coeff_grad;

	// de-allocate b & rho
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			delete [] pSQ->b[k][j];
			delete [] pSQ->b_tilda[k][j];
			delete [] pSQ->Vc[k][j];
			delete [] pSQ->rho[k][j];
			delete [] pSQ->rho_at[k][j];
			delete [] pSQ->drho_new[k][j];
			delete [] pSQ->drho[k][j];
			delete [] pSQ->drho_dt[k][j];
			delete [] pSQ->drho_2dt[k][j];
			delete [] pSQ->rhs[k][j];
		}
		delete [] pSQ->b[k];
		delete [] pSQ->b_tilda[k];
		delete [] pSQ->Vc[k];
		delete [] pSQ->rho[k];
		delete [] pSQ->rho_at[k];
		delete [] pSQ->drho_new[k];
		delete [] pSQ->drho[k];
		delete [] pSQ->drho_dt[k];
		delete [] pSQ->drho_2dt[k];
		delete [] pSQ->rhs[k];
	}
	delete [] pSQ->b;
	delete [] pSQ->b_tilda;
	delete [] pSQ->Vc;
	delete [] pSQ->rho;
	delete [] pSQ->rho_at;
	delete [] pSQ->drho_new;
	delete [] pSQ->drho_dt;
	delete [] pSQ->drho_2dt;
	delete [] pSQ->rhs;

	// de-allocate phi_guess
	for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++)
	{
		for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++)
		{
			delete [] pSQ->phi_guess[k][j];
		}
		delete [] pSQ->phi_guess[k];
	}
	delete [] pSQ->phi_guess;

	// de-allocate comm topologies
	delete [] pSQ->neighs_lap;
	delete [] pSQ->neighs_sq;
	MPI_Comm_free(&pSQ->comm_laplacian);
	MPI_Comm_free(&pSQ->comm_sq);

	// de-allocate Laplacian comm stuff
	for(k=0;k<3;k++)
	{
		delete [] pSQ->LapInd.eout_s[k];
		delete [] pSQ->LapInd.eout_e[k];
		delete [] pSQ->LapInd.ein_s[k];
		delete [] pSQ->LapInd.ein_e[k];
		delete [] pSQ->LapInd.stencil_sign[k];
		delete [] pSQ->LapInd.edge_ind[k];
	}
	delete [] pSQ->LapInd.eout_s;
	delete [] pSQ->LapInd.eout_e;
	delete [] pSQ->LapInd.ein_s;
	delete [] pSQ->LapInd.ein_e;
	delete [] pSQ->LapInd.stencil_sign;
	delete [] pSQ->LapInd.edge_ind;

	delete [] pSQ->LapInd.displs_send;
	delete [] pSQ->LapInd.ncounts_send;
	delete [] pSQ->LapInd.displs_recv;
	delete [] pSQ->LapInd.ncounts_recv;  

	// de-allocate phi
	for(k=0;k<pSQ->np_z+2*pSQ->FDn;k++)
	{
		for(j=0;j<pSQ->np_y+2*pSQ->FDn;j++)
		{
			delete [] pSQ->phi[k][j];
		}
		delete [] pSQ->phi[k];
	}
	delete [] pSQ->phi;





	// de-allocate Vxc
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			delete [] pSQ->Vxc[k][j];
		}
		delete [] pSQ->Vxc[k];
	}
	delete [] pSQ->Vxc;

	// de-allocate Veff
	for(k=0;k<pSQ->np_z+2*pSQ->nloc;k++)
	{
		for(j=0;j<pSQ->np_y+2*pSQ->nloc;j++)
		{
			delete [] pSQ->Veff[k][j];
		}
		delete [] pSQ->Veff[k];
	}
	delete [] pSQ->Veff;



	// de-allocate sq stuff

	for(k=0;k<3;k++)
	{
		delete [] pSQ->SqInd.eout_s[k];
		delete [] pSQ->SqInd.eout_e[k];
		delete [] pSQ->SqInd.ein_s[k];
		delete [] pSQ->SqInd.ein_e[k];
	}
	delete [] pSQ->SqInd.eout_s;
	delete [] pSQ->SqInd.eout_e;
	delete [] pSQ->SqInd.ein_s;
	delete [] pSQ->SqInd.ein_e;

	delete [] pSQ->SqInd.displs_send;
	delete [] pSQ->SqInd.ncounts_send;
	delete [] pSQ->SqInd.displs_recv;
	delete [] pSQ->SqInd.ncounts_recv;



	delete [] pSQ->chi;
	delete [] pSQ->zee;  
	for(k=0;k<pSQ->np_y*pSQ->np_z*pSQ->np_x;k++)
	{
		delete [] pSQ->rho_pj[k];
		delete [] pSQ->Ci[k];
	}
	delete [] pSQ->rho_pj;
	delete [] pSQ->Ci;


	// de-allocate forces

	delete [] pSQ->flocx;
	delete [] pSQ->flocy;
	delete [] pSQ->flocz;

	delete [] pSQ->fcorrx;
	delete [] pSQ->fcorry;
	delete [] pSQ->fcorrz;

	delete [] pSQ->fnlocx;
	delete [] pSQ->fnlocy;
	delete [] pSQ->fnlocz;

	delete [] pSQ->fx;
	delete [] pSQ->fy;
	delete [] pSQ->fz;
	delete [] pSQ->forces;

	// de-allocate MD stuff
	delete [] pSQ->vels;
	delete [] pSQ->vels_old;
	delete [] pSQ->accels;
}
