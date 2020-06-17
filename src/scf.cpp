/** \file scf.cpp
  \brief This file contains most of the functions required in the SCF iteration. The SQ related function are in a separate file.


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

// function to compute exchange correlation potential
void ExchangeCorrelationPotential(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double p,A,alpha1,beta1,beta2,beta3,beta4,C3,Vxci,rhoi;
	// Perdew-Wang (Ceperley-Alder)
	p = 1.0 ;
	A = 0.031091 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;

	C3 = 0.9847450218427;  //  Exchange potential parameter


	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{
				rhoi = pSQ->rho[k][j][i];


				if (rhoi==0)
				{
					Vxci = 0.0 ;
				}
				else
				{
					Vxci = pow((0.75/(M_PI*rhoi)),(1.0/3.0)) ;
					Vxci = (-2.0*A*(1.0+alpha1*Vxci))*log(1.0+1.0/(2.0*A*(beta1*pow(Vxci,0.5) + beta2*Vxci + beta3*pow(Vxci,1.5) + beta4*pow(Vxci,(p+1.0)))))
						- (Vxci/3.0)*(-2.0*A*alpha1*log(1.0+1.0/(2.0*A*( beta1*pow(Vxci,0.5) + beta2*Vxci + beta3*pow(Vxci,1.5) + beta4*pow(Vxci,(p+1.0)))))
								- ((-2.0*A*(1.0+alpha1*Vxci))*(A*( beta1*pow(Vxci,-0.5)+ 2.0*beta2 + 3.0*beta3*pow(Vxci,0.5) + 2.0*(p+1.0)*beta4*pow(Vxci,p) )))
								/((2.0*A*( beta1*pow(Vxci,0.5) + beta2*Vxci + beta3*pow(Vxci,1.5) + beta4*pow(Vxci,(p+1.0)) ) )*(2.0*A*( beta1*pow(Vxci,0.5) + beta2*Vxci + beta3*pow(Vxci,1.5) + beta4*pow(Vxci,(p+1.0)) ) )+(2.0*A*( beta1*pow(Vxci,0.5) + beta2*Vxci + beta3*pow(Vxci,1.5) + beta4*pow(Vxci,(p+1.0)) ) )) ) ;

				}
				Vxci = Vxci - C3*pow(rhoi,1.0/3.0) ;
				pSQ->Vxc[k][j][i] = Vxci;	      
			}
		}
	} 

}

// function to compute exchange correlation potential
void EvaluateEffectivePotential(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_z;j++)
		{
			for(i=0;i<pSQ->np_z;i++)
			{/// Veff defined on processor+Rcut domain
				pSQ->Veff[k+pSQ->nloc][j+pSQ->nloc][i+pSQ->nloc] = pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn] + pSQ->Vxc[k][j][i];
			}
		}
	} 
}

void SCF_iteration(DS_SQ *pSQ)
{

	int rank,count=1,MAX_SCF=pSQ->scf_maxiter,MIN_SCF=pSQ->scf_miniter,i,j,k,ctr,m;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tscf_mpi=0.0,tscf_mpi_full=0.0;
	double SCF_TOL=pSQ->scf_tol, err=SCF_TOL+1.0,temp; //phi_fix=0.0

	// Anderson stuff
	double beta_scf=pSQ->beta_scf;  
	int m_scf=pSQ->m_scf; 
	int p_scf = pSQ->p_scf; 
	double *Veff_new,*Veff_old,*Veff_res,*Xold,*Fold,**DX,**DF,*am_vec; 

	// allocate memory to store phi(x) in domain
	Veff_new = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Veff_old = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Veff_res = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Xold = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	Fold = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 
	am_vec = new double [pSQ->np_z*pSQ->np_y*pSQ->np_x](); // need to de-allocate later 

	// allocate memory to store DX, DF history matrices
	DX = new double* [m_scf](); // need to de-allocate later
	DF = new double* [m_scf](); // need to de-allocate later
	if(DF == NULL)
	{
		if(rank==0)
			cout << "Memory allocation failed in DF"<< endl;
		exit(1);
	}
	for(k=0;k<m_scf;k++)
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

	// initialize Veff_old
	ctr=0;
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				Veff_old[ctr] = pSQ->Veff[k+pSQ->nloc][j+pSQ->nloc][i+pSQ->nloc];
				ctr=ctr+1;
			}
		}
	}

	double tt0,tt1,tt2,tempp,tengy0,tengy1,t_end,t_begin;
	tt0 = MPI_Wtime();
	t_begin = MPI_Wtime();

	while(err > SCF_TOL && count<=MAX_SCF)
	{
		tscf_mpi=0.0;

		/// Compute electron density, Cheb coeff/components, scaling factors for Hsub, Fermi energy, band structure energy and entropy energy
		ClenshawCurtisSpectralQuadrature(pSQ,count);
		tscf_mpi+=pSQ->mpi_time;
		pSQ->tsq_mpi+=pSQ->mpi_time;

		/// Evaluate total energy
		tengy0 = MPI_Wtime();
		EvaluateTotalEnergy(pSQ);
		tengy1 = MPI_Wtime();
		pSQ->engy_time+=tengy1-tengy0;
		tscf_mpi+=pSQ->mpi_time;
		pSQ->tengy_mpi+=pSQ->mpi_time;

		if(pSQ->prnt_scf==1)
		{PrintSCF(pSQ);}

		/// Solve Poisson's equation to update phi
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_z;j++)
			{
				for(i=0;i<pSQ->np_z;i++)
				{
					pSQ->phi_guess[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn]=pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn];
					pSQ->rhs[k][j][i]=pSQ->b[k][j][i]+pSQ->rho[k][j][i];
				}
			}
		}
		PoissonSolver_AAJ(pSQ);
		tscf_mpi+=pSQ->mpi_time;

		/// Compute Vxc and Veff
		ExchangeCorrelationPotential(pSQ);
		/// Compute residual of Veff i.e Veff_res = Veff_new-Veff_old
		ctr=0;
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					Veff_new[ctr] = pSQ->phi[k+pSQ->FDn][j+pSQ->FDn][i+pSQ->FDn] + pSQ->Vxc[k][j][i];// Effective potential //pSQ->Veff[k+pSQ->nloc][j+pSQ->nloc][i+pSQ->nloc];
					Veff_res[ctr] = Veff_new[ctr] - Veff_old[ctr];
					ctr=ctr+1;
				}
			}
		}

		/// Compute error in this SCF iteration
		Vector2Norm(Veff_res,pSQ->np_x*pSQ->np_y*pSQ->np_z,&err);
		Vector2Norm(Veff_new,pSQ->np_x*pSQ->np_y*pSQ->np_z,&temp);
		err = err/temp;
		if(rank==0)
			printf("Iter:%d,SCF_err = %g, Fermi energy (Ha) : %.12f  and Total energy (Ha/atom) : %.12f \n",count,err,pSQ->lambda_f,pSQ->Etot/pSQ->n_atm);

		/// Simple mixing update first
		ctr=0;
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					Veff_new[ctr] = Veff_old[ctr] + beta_scf*Veff_res[ctr];
					ctr=ctr+1;
				}
			}
		}

		/// Store residual and iterate history
		if(count>1)
		{
			m = ((count-2) % m_scf)+1-1; //-1 because the index starts from 0
			ctr=0;
			for(k=0;k<pSQ->np_z;k++)
			{
				for(j=0;j<pSQ->np_y;j++)
				{
					for(i=0;i<pSQ->np_x;i++)
					{
						DX[m][ctr]=Veff_old[ctr]-Xold[ctr];
						DF[m][ctr]=Veff_res[ctr]-Fold[ctr];
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
					Xold[ctr] = Veff_old[ctr];
					Fold[ctr] = Veff_res[ctr];		    
					ctr=ctr+1;
				}
			}
		}

		/// Anderson update
		if(count % p_scf == 0 && count>1)
		{
			tt1 = MPI_Wtime();
			AndersonExtrapolation(DX,DF,Veff_res,beta_scf,m_scf,pSQ->np_z*pSQ->np_y*pSQ->np_x,am_vec,&tscf_mpi);
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
						Veff_new[ctr]=Veff_new[ctr]-am_vec[ctr];
						ctr=ctr+1;
					}
				}
			}

		}


		// swap variables
		ctr=0;
		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					Veff_old[ctr] = Veff_new[ctr];
					if(err > SCF_TOL && count < MAX_SCF)
					{pSQ->Veff[k+pSQ->nloc][j+pSQ->nloc][i+pSQ->nloc] = Veff_new[ctr]; // do not need alltoall comm after scf convergence
					}
					ctr=ctr+1;
				}
			}
		}

		tscf_mpi_full+=tscf_mpi;

		t_end = MPI_Wtime();
		if(rank==0)
			printf("Time taken for current SCF iteration = %.2f (%.2f) sec \n",t_end-t_begin,tscf_mpi);
		t_begin=t_end;

		count = count + 1;
	}
	tt1 = MPI_Wtime();

	if(rank==0)
	{
		if(count<=MAX_SCF && err <=SCF_TOL && count >=MIN_SCF)
			printf("SCF converged! \n");
		if(count>MAX_SCF)
			printf("WARNING: SCF exceeded maximum iterations. \n");
		if(count<MIN_SCF)
			printf("WARNING: SCF iterations might not converged. \n");
		printf("No. of  SCF iterations : %d \n",count-1);
		printf("Free energy (Ha/atom)    : %.12f \n",pSQ->Etot/pSQ->n_atm);
		printf("Total SCF time         : %.2f (%.2f) sec \n\n\n",tt1-tt0,tscf_mpi_full);
	}


	// de-allocate memory
	delete [] Veff_new;
	delete [] Veff_old;
	delete [] Veff_res;
	delete [] Xold;
	delete [] Fold;
	delete [] am_vec;

	for(k=0;k<m_scf;k++)
	{
		delete [] DX[k];
		delete [] DF[k];
	}
	delete [] DX;
	delete [] DF;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// function to perform SCF and compute DFT forces
void SQDFT_forces(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double t_pois,t_scf,t_forc,t_prnt,t_nloc,t_begin,t_end;

	// ---------------- Compute DFT Forces (accelerations) --------------------
	t_begin = MPI_Wtime();
	/// Set up ionic/electron density in space
	pSQ->mpi_time=0.0;   
	Replica_atoms(pSQ);
	ChargeDensity(pSQ);
	if(pSQ->Correction == 2)
		OverlapCorrection(pSQ);
	Replica_atoms_Rcut(pSQ);
	NonlocalProjectors(pSQ);
	t_nloc = MPI_Wtime();

	/// Solve Poisson's equation to find phi.    
	double t_end_p,t_begin_p;
	t_begin_p = MPI_Wtime();  
	PoissonSolver_AAJ(pSQ);
	t_end_p = MPI_Wtime();
	if(rank == 0) 
		printf("Rank=%d,Poisson Solver took %.4f seconds. \n",rank,t_end_p-t_begin_p);
	MPI_Barrier(MPI_COMM_WORLD);
	t_pois = MPI_Wtime();

	/// Compute Vxc and Veff
	ExchangeCorrelationPotential(pSQ);
	EvaluateEffectivePotential(pSQ);

	/// SCF iteration
	SCF_iteration(pSQ);
	t_scf = MPI_Wtime();

	/// Forces
	ForcesTotal(pSQ);
	t_forc = MPI_Wtime();

	t_prnt = t_forc;
	if(pSQ->prnt_atoms==1 && rank==0 && (pSQ->MaxMDsteps==0 || pSQ->MDstepCount==pSQ->MaxMDsteps))
	{
		PrintAtoms(pSQ);
	}

	t_end = MPI_Wtime();

	if(rank==0)
	{printf("Total wall time (comm time)   = %.4f  seconds. \n",t_end-t_begin);
		printf("Break down: \n");
		printf("MPI_Barrier after 1st Poisson = %.4f seconds. \n",t_pois-t_end_p);
		printf("SCF time                      = %.4f seconds. \n",t_scf-t_pois);
		printf("Forces time                   = %.4f seconds. \n",t_forc-t_scf);
	}

}


