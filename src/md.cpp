/** \file md.cpp
  \brief This file contains functions required for Molecular Dynamics (MD).


*/

#include "headers.h"
#include "ds_sq.h"
#include "func_sq.h"
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef TEMP_TOL
#define TEMP_TOL 1e-12
#endif

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// MD stuff ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
void datacollect(DS_SQ* pSQ) {
	ofstream myfile;
	if(pSQ->restart_md==1)
		myfile.open (pSQ->name,ios::out|ios::app);
	else
		myfile.open ("Position.txt",ios::out|ios::app);
	if(pSQ->MDstepCount == 1)
	{
		myfile << "Trajectory Data\n";
		myfile << "---------------.\n\n";
		myfile << pSQ->MaxMDsteps<<endl ;
		myfile << pSQ->domain[0] << " " << pSQ->domain[1] << " " << pSQ->domain[2] << endl ;
		myfile << pSQ->n_atm << endl << endl ; 
	}

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			myfile << " " << pSQ->Atm[ii].Rx[jj].main_atm << " ";
			myfile << " " << pSQ->Atm[ii].Ry[jj].main_atm << " ";
			myfile << " " << pSQ->Atm[ii].Rz[jj].main_atm << endl;
		}
	}
	myfile << endl ;
	myfile.close();
}

// function to perform NVE (microcanonical) MD
void NVE_MD(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Leapfrog step (part-1)
	Leapfrog_part1(pSQ);

	// Apply boundary condition (PBC) - wrap around --> do this before any replicas are stored
	WrapAround_PBC(pSQ);

	// Charge extrapolation (for better rho_guess)
	ChargeExtrapolation_MD(pSQ);

	// Compute DFT energy and forces from SCF
	SQDFT_forces(pSQ);

	// Leapfrog step (part-2)
	Leapfrog_part2(pSQ);

	// Compute MD energies (TE=KE+PE)
	MD_energies(pSQ);

	if(rank == 0)
		datacollect(pSQ); 

}


// function to update position using Leapfrog method
void Leapfrog_part1(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	/*
	// Leapfrog algorithm
	PART-I
	v_(t+dt/2) = v_t + 0.5*dt*a_t
	r_(t+dt) = r_t + dt*v_(t+dt/2)

	PART-II (after computing a_(t+dt) for current positions)
	v_(t+dt) = v_(t+dt/2) + 0.5*dt*a_(t+dt)
	*/

	int k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			// v_(t+dt/2) = v_t + 0.5*dt*a_t
			pSQ->vels[k] = pSQ->vels[k] + pSQ->time_step*0.5*pSQ->accels[k];
			pSQ->vels[k+pSQ->n_atm] = pSQ->vels[k+pSQ->n_atm] + pSQ->time_step*0.5*pSQ->accels[k+pSQ->n_atm];
			pSQ->vels[k+2*pSQ->n_atm] = pSQ->vels[k+2*pSQ->n_atm] + pSQ->time_step*0.5*pSQ->accels[k+2*pSQ->n_atm];

			if(pSQ->MDstepCount==1) // non-mapped atom coords required for charge extrapolation
			{
				pSQ->Atm[ii].Rx[jj].main_atm_nm = pSQ->Atm[ii].Rx[jj].main_atm;
				pSQ->Atm[ii].Ry[jj].main_atm_nm = pSQ->Atm[ii].Ry[jj].main_atm;
				pSQ->Atm[ii].Rz[jj].main_atm_nm = pSQ->Atm[ii].Rz[jj].main_atm;
			}

			// r_(t+dt) = r_t + dt*v_(t+dt/2)
			pSQ->Atm[ii].Rx[jj].main_atm = pSQ->Atm[ii].Rx[jj].main_atm + pSQ->time_step*pSQ->vels[k];
			pSQ->Atm[ii].Ry[jj].main_atm = pSQ->Atm[ii].Ry[jj].main_atm + pSQ->time_step*pSQ->vels[k+pSQ->n_atm];
			pSQ->Atm[ii].Rz[jj].main_atm = pSQ->Atm[ii].Rz[jj].main_atm + pSQ->time_step*pSQ->vels[k+2*pSQ->n_atm];

			pSQ->Atm[ii].Rx[jj].main_atm_nm = pSQ->Atm[ii].Rx[jj].main_atm_nm + pSQ->time_step*pSQ->vels[k];
			pSQ->Atm[ii].Ry[jj].main_atm_nm = pSQ->Atm[ii].Ry[jj].main_atm_nm + pSQ->time_step*pSQ->vels[k+pSQ->n_atm];
			pSQ->Atm[ii].Rz[jj].main_atm_nm = pSQ->Atm[ii].Rz[jj].main_atm_nm + pSQ->time_step*pSQ->vels[k+2*pSQ->n_atm];

			k+=1;
		}
	}

}

// function to update velocity using Leapfrog method
void Leapfrog_part2(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	/*
	// Leapfrog algorithm
	PART-I
	v_(t+dt/2) = v_t + 0.5*dt*a_t
	r_(t+dt) = r_t + dt*v_(t+dt/2)

	PART-II (after computing a_(t+dt) for current positions)
	v_(t+dt) = v_(t+dt/2) + 0.5*dt*a_(t+dt)
	*/

	// update accelerations
	int k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			// a = F/m in Bohr/(fs)^2
			pSQ->accels[k]=pSQ->forces[k]/pSQ->Atm[ii].mass; 
			pSQ->accels[k+pSQ->n_atm]=pSQ->forces[k+pSQ->n_atm]/pSQ->Atm[ii].mass;
			pSQ->accels[k+2*pSQ->n_atm]=pSQ->forces[k+2*pSQ->n_atm]/pSQ->Atm[ii].mass;

			k+=1;
		}
	}
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			// v_(t+dt) = v_(t+dt/2) + 0.5*dt*a_(t+dt)
			pSQ->vels[k] = pSQ->vels[k] + pSQ->time_step*0.5*pSQ->accels[k];
			pSQ->vels[k+pSQ->n_atm] = pSQ->vels[k+pSQ->n_atm] + pSQ->time_step*0.5*pSQ->accels[k+pSQ->n_atm];
			pSQ->vels[k+2*pSQ->n_atm] = pSQ->vels[k+2*pSQ->n_atm] + pSQ->time_step*0.5*pSQ->accels[k+2*pSQ->n_atm];

			k+=1;
		}
	}

}

// function perform charge extrapolation to provide better rho_guess for MD steps
void ChargeExtrapolation_MD(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank,i,j,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double alpha,beta;

	// drho's 
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->drho[k][j][i]=pSQ->rho[k][j][i]-pSQ->rho_at[k][j][i]; 
			}
		}
	}  

	if(pSQ->MDstepCount>3) // at the end of 3rd step, can do the first extrapolation
	{
		// first find the 2 x 2 matrix and 2 x 1 rhs through dot products
		double FtF[4],Ftf[2]; // matrix and rhs
		double svec[2];

		double temp=0.0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				temp += (pSQ->Atm[ii].Rx[jj].main_atm_0dt - pSQ->Atm[ii].Rx[jj].main_atm_1dt)*(pSQ->Atm[ii].Rx[jj].main_atm_0dt - pSQ->Atm[ii].Rx[jj].main_atm_1dt);
				temp += (pSQ->Atm[ii].Ry[jj].main_atm_0dt - pSQ->Atm[ii].Ry[jj].main_atm_1dt)*(pSQ->Atm[ii].Ry[jj].main_atm_0dt - pSQ->Atm[ii].Ry[jj].main_atm_1dt);
				temp += (pSQ->Atm[ii].Rz[jj].main_atm_0dt - pSQ->Atm[ii].Rz[jj].main_atm_1dt)*(pSQ->Atm[ii].Rz[jj].main_atm_0dt - pSQ->Atm[ii].Rz[jj].main_atm_1dt);

			}
		}
		FtF[0] = temp; // a_11

		temp=0.0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				temp += (pSQ->Atm[ii].Rx[jj].main_atm_0dt - pSQ->Atm[ii].Rx[jj].main_atm_1dt)*(pSQ->Atm[ii].Rx[jj].main_atm_1dt - pSQ->Atm[ii].Rx[jj].main_atm_2dt);
				temp += (pSQ->Atm[ii].Ry[jj].main_atm_0dt - pSQ->Atm[ii].Ry[jj].main_atm_1dt)*(pSQ->Atm[ii].Ry[jj].main_atm_1dt - pSQ->Atm[ii].Ry[jj].main_atm_2dt);
				temp += (pSQ->Atm[ii].Rz[jj].main_atm_0dt - pSQ->Atm[ii].Rz[jj].main_atm_1dt)*(pSQ->Atm[ii].Rz[jj].main_atm_1dt - pSQ->Atm[ii].Rz[jj].main_atm_2dt);

			}
		}
		FtF[1] = temp; // a_12
		FtF[2] = temp; // a_21

		temp=0.0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				temp += (pSQ->Atm[ii].Rx[jj].main_atm_1dt - pSQ->Atm[ii].Rx[jj].main_atm_2dt)*(pSQ->Atm[ii].Rx[jj].main_atm_1dt - pSQ->Atm[ii].Rx[jj].main_atm_2dt);
				temp += (pSQ->Atm[ii].Ry[jj].main_atm_1dt - pSQ->Atm[ii].Ry[jj].main_atm_2dt)*(pSQ->Atm[ii].Ry[jj].main_atm_1dt - pSQ->Atm[ii].Ry[jj].main_atm_2dt);
				temp += (pSQ->Atm[ii].Rz[jj].main_atm_1dt - pSQ->Atm[ii].Rz[jj].main_atm_2dt)*(pSQ->Atm[ii].Rz[jj].main_atm_1dt - pSQ->Atm[ii].Rz[jj].main_atm_2dt);

			}
		}
		FtF[3] = temp; // a_22

		temp=0.0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				temp += (pSQ->Atm[ii].Rx[jj].main_atm_0dt - pSQ->Atm[ii].Rx[jj].main_atm_1dt)*(pSQ->Atm[ii].Rx[jj].main_atm_nm - pSQ->Atm[ii].Rx[jj].main_atm_0dt);
				temp += (pSQ->Atm[ii].Ry[jj].main_atm_0dt - pSQ->Atm[ii].Ry[jj].main_atm_1dt)*(pSQ->Atm[ii].Ry[jj].main_atm_nm - pSQ->Atm[ii].Ry[jj].main_atm_0dt);
				temp += (pSQ->Atm[ii].Rz[jj].main_atm_0dt - pSQ->Atm[ii].Rz[jj].main_atm_1dt)*(pSQ->Atm[ii].Rz[jj].main_atm_nm - pSQ->Atm[ii].Rz[jj].main_atm_0dt);

			}
		}
		Ftf[0] = temp; // b_1

		temp=0.0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				temp += (pSQ->Atm[ii].Rx[jj].main_atm_1dt - pSQ->Atm[ii].Rx[jj].main_atm_2dt)*(pSQ->Atm[ii].Rx[jj].main_atm_nm - pSQ->Atm[ii].Rx[jj].main_atm_0dt);
				temp += (pSQ->Atm[ii].Ry[jj].main_atm_1dt - pSQ->Atm[ii].Ry[jj].main_atm_2dt)*(pSQ->Atm[ii].Ry[jj].main_atm_nm - pSQ->Atm[ii].Ry[jj].main_atm_0dt);
				temp += (pSQ->Atm[ii].Rz[jj].main_atm_1dt - pSQ->Atm[ii].Rz[jj].main_atm_2dt)*(pSQ->Atm[ii].Rz[jj].main_atm_nm - pSQ->Atm[ii].Rz[jj].main_atm_0dt);

			}
		}
		Ftf[1] = temp; // b_2

		// second find alpha and beta by taking pseudoinverse of the 2 x 2 matrix
		PseudoInverseTimesVec(FtF,Ftf,svec,2); // svec=inv(FtF)*Ftf
		alpha=svec[0];
		beta=svec[1];

		// Extrapolation: drho_new = drho + alpha*(drho - drho_dt) + beta*(drho_dt - drho_2dt) = (1+alpha)*drho + (beta-alpha)*drho_dt - beta*drho_2dt 

		for(k=0;k<pSQ->np_z;k++)
		{
			for(j=0;j<pSQ->np_y;j++)
			{
				for(i=0;i<pSQ->np_x;i++)
				{
					pSQ->drho_new[k][j][i] = (1+alpha)*pSQ->drho[k][j][i] + (beta-alpha)*pSQ->drho_dt[k][j][i] - beta*pSQ->drho_2dt[k][j][i];

				}
			}
		} 

	}

	// drho's 
	for(k=0;k<pSQ->np_z;k++)
	{
		for(j=0;j<pSQ->np_y;j++)
		{
			for(i=0;i<pSQ->np_x;i++)
			{
				pSQ->drho_2dt[k][j][i]=pSQ->drho_dt[k][j][i];
				pSQ->drho_dt[k][j][i]=pSQ->drho[k][j][i];
			}
		}
	}   

	// store atomic positions of 1*dt into 2*dt, 0*dt into 1*dt and current atomic positions into 0*dt (after extrapolation)
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->Atm[ii].Rx[jj].main_atm_2dt = pSQ->Atm[ii].Rx[jj].main_atm_1dt;
			pSQ->Atm[ii].Ry[jj].main_atm_2dt = pSQ->Atm[ii].Ry[jj].main_atm_1dt;
			pSQ->Atm[ii].Rz[jj].main_atm_2dt = pSQ->Atm[ii].Rz[jj].main_atm_1dt;

			pSQ->Atm[ii].Rx[jj].main_atm_1dt = pSQ->Atm[ii].Rx[jj].main_atm_0dt;
			pSQ->Atm[ii].Ry[jj].main_atm_1dt = pSQ->Atm[ii].Ry[jj].main_atm_0dt;
			pSQ->Atm[ii].Rz[jj].main_atm_1dt = pSQ->Atm[ii].Rz[jj].main_atm_0dt;

			pSQ->Atm[ii].Rx[jj].main_atm_0dt = pSQ->Atm[ii].Rx[jj].main_atm_nm;
			pSQ->Atm[ii].Ry[jj].main_atm_0dt = pSQ->Atm[ii].Ry[jj].main_atm_nm;
			pSQ->Atm[ii].Rz[jj].main_atm_0dt = pSQ->Atm[ii].Rz[jj].main_atm_nm;
		}
	}

}


// function to compute Total energy (TE) = Kinetic energy (KE) + Potential energy (PE)
void MD_energies(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	pSQ->PE = pSQ->Etot/pSQ->n_atm; // Ha/atom

	double vsum_x=0.0,vsum_y=0.0,vsum_z=0.0,fsum_x=0.0,fsum_y=0.0,fsum_z=0.0;

	pSQ->KE=0.0;
	int k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->KE+=0.5*pSQ->Atm[ii].mass*(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]); // 1/2 m v^2

			vsum_x += pSQ->vels[k];
			vsum_y += pSQ->vels[k+pSQ->n_atm];
			vsum_z += pSQ->vels[k+2*pSQ->n_atm];

			fsum_x += pSQ->forces[k];
			fsum_y += pSQ->forces[k+pSQ->n_atm];
			fsum_z += pSQ->forces[k+2*pSQ->n_atm];

			k+=1;
		}
	}

	pSQ->KE=pSQ->KE/pSQ->n_atm;

	pSQ->TE=pSQ->PE+pSQ->KE;

	pSQ->T_MD = 2*pSQ->KE*pSQ->n_atm/(pSQ->kB*pSQ->dof);//DOF check, Update temperature --> in NVE (temp can change drastically and this means that SQ parameters also need to change accordingly for accuracy), Need NVT for fixed temperature
	pSQ->T=pSQ->T_MD; // Need to uncomment this (if ionic and electronic temperature are the same)

	if(rank==0)
	{
		printf("Force sum (Ha/Bohr)          : %.12f, %.12f, %.12f \n",fsum_x,fsum_y,fsum_z);
		printf("Velocity sum (Bohr/fs)       : %.12f, %.12f, %.12f \n",vsum_x,vsum_y,vsum_z);
		printf("Temperature (K)              : %.12f \n",pSQ->T_MD);
		printf("Energy (TE, PE, KE) (Ha/atom): %.12f %.12f %.12f \n",pSQ->TE, pSQ->PE, pSQ->KE);

	}
}

// function to initialize velocities and accelerations for Molecular Dynamics (MD)
void Initialize_MD(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	pSQ->snose = 0.0 ;
	pSQ->xi_nose = 0.0 ;
	if(rank == 0)
	{
		if(pSQ->ensemble == 1)
			printf("Canonical Ensemble in progress\n");
		else if(pSQ->ensemble == 2)
			printf("Microcanonical Ensemble in progress\n");
		else
			printf("Invalid option");
	}
	if(rank==0)
		printf("Using MD time step of %f fs \n",pSQ->time_step);
	pSQ->vels = new double [3*pSQ->n_atm]();
	pSQ->accels = new double [3*pSQ->n_atm]();
	pSQ->T_MD = pSQ->T;
	pSQ->Treq = pSQ->T_MD; // needed only for NVT MD
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		pSQ->Atm[ii].mass *= 1.066572305; // Conversion of units: 1 amu = 1.066572305 Ha (fs)^2 / Bohr^2 --> units of mass
	}

	// Initial velocity magnitudes
	double *velMag; // initial magnitude of all velocities (in Bohr/fs) based on temperature
	velMag = new double [pSQ->n_atm]();

	double mass_sum=0.0;
	int k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			mass_sum = mass_sum + pSQ->Atm[ii].mass;           
			k+=1;
		}
	}

	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			velMag[k]=sqrt(pSQ->dof*pSQ->Treq*pSQ->kB/mass_sum);
			k+=1;
		}
	}

	// Initial velocity direction (uniformly distributed random unit vectors)

	double s,x,y,vsum_x=0.0,vsum_y=0.0,vsum_z=0.0;
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		srand(ii+pSQ->rand_seed); // random seed, NOTE: Cannot use "srand(time(NULL))" as we want the seed to be the same across all the processors. Instead can have an option to input a different integer (rand_seed) through .inpt file and do "srand(ii+rand_seed)". NOTE: Call srand() before the following loop, otherwise, different procs are generating different sequence
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{  
			if(pSQ->vds == 2)
			{
				s=2.0;
				while(s>1.0)
				{
					x=2.0*((double)(rand()) / (double)(RAND_MAX))-1; // random number between -1 and +1
					y=2.0*((double)(rand()) / (double)(RAND_MAX))-1;
					s = x*x + y*y;
				}

				pSQ->vels[k+2*pSQ->n_atm]=velMag[k]*(1.0-2.0*s); // z-component
				s = 2.0*sqrt(1.0-s);
				pSQ->vels[k]=s*x*velMag[k]; // x-component
				pSQ->vels[k+pSQ->n_atm]=s*y*velMag[k]; // y-component
			}
			else
			{
				pSQ->vels[k] = sqrt(pSQ->kB*pSQ->Treq/pSQ->Atm[ii].mass)*cos(2*M_PI*((double)(rand()) / (double)(RAND_MAX)));
				pSQ->vels[k] *= sqrt(-2.0*log((double)(rand()) / (double)(RAND_MAX)));
				pSQ->vels[k+pSQ->n_atm] = sqrt(pSQ->kB*pSQ->Treq/pSQ->Atm[ii].mass)*cos(2*M_PI*((double)(rand()) / (double)(RAND_MAX)));
				pSQ->vels[k+pSQ->n_atm] *= sqrt(-2.0*log((double)(rand()) / (double)(RAND_MAX)));
				pSQ->vels[k+2*pSQ->n_atm] = sqrt(pSQ->kB*pSQ->Treq/pSQ->Atm[ii].mass)*cos(2*M_PI*((double)(rand()) / (double)(RAND_MAX)));
				pSQ->vels[k+2*pSQ->n_atm] *= sqrt(-2.0*log((double)(rand()) / (double)(RAND_MAX)));
			}  
			vsum_x += pSQ->vels[k];
			vsum_y += pSQ->vels[k+pSQ->n_atm];
			vsum_z += pSQ->vels[k+2*pSQ->n_atm];

			pSQ->accels[k]=pSQ->forces[k]/pSQ->Atm[ii].mass; 
			pSQ->accels[k+pSQ->n_atm]=pSQ->forces[k+pSQ->n_atm]/pSQ->Atm[ii].mass;
			pSQ->accels[k+2*pSQ->n_atm]=pSQ->forces[k+2*pSQ->n_atm]/pSQ->Atm[ii].mass;

			k+=1;
		}
	}

	// Remove translation (rotation not possible because of PBC)
	double KE=0.0;
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->vels[k] = pSQ->vels[k] - (vsum_x/pSQ->n_atm);
			pSQ->vels[k+pSQ->n_atm] = pSQ->vels[k+pSQ->n_atm] - (vsum_y/pSQ->n_atm);
			pSQ->vels[k+2*pSQ->n_atm] = pSQ->vels[k+2*pSQ->n_atm] - (vsum_z/pSQ->n_atm);

			KE+=0.5*pSQ->Atm[ii].mass*(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]); // 1/2 m v^2

			k+=1;
		}
	}
	// rescale velocities   
	double vel_rescale ;
	vel_rescale = sqrt(pSQ->dof*pSQ->kB*pSQ->Treq/(2.0*KE)) ;
	k = 0;KE = 0.0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->vels[k] *= vel_rescale;
			pSQ->vels[k+pSQ->n_atm] *= vel_rescale;
			pSQ->vels[k+2*pSQ->n_atm] *= vel_rescale;
			KE += 0.5*pSQ->Atm[ii].mass*(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]);
			k += 1 ;
		}
	}
	pSQ->KE = KE/pSQ->n_atm ;
	delete [] velMag;
}

// function to wrap around atom positions that lie outside main domain (periodic boundary conditions)
void WrapAround_PBC(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//double temp;

	if(rank==0)
		printf("Atom positions:\n");

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{

			if(pSQ->Atm[ii].Rx[jj].main_atm >= pSQ->domain[0])
				pSQ->Atm[ii].Rx[jj].main_atm = pSQ->Atm[ii].Rx[jj].main_atm - pSQ->domain[0];
			if(pSQ->Atm[ii].Ry[jj].main_atm >= pSQ->domain[1])
				pSQ->Atm[ii].Ry[jj].main_atm = pSQ->Atm[ii].Ry[jj].main_atm - pSQ->domain[1];
			if(pSQ->Atm[ii].Rz[jj].main_atm >= pSQ->domain[2])
				pSQ->Atm[ii].Rz[jj].main_atm = pSQ->Atm[ii].Rz[jj].main_atm - pSQ->domain[2];

			if(pSQ->Atm[ii].Rx[jj].main_atm<0)
				pSQ->Atm[ii].Rx[jj].main_atm = pSQ->Atm[ii].Rx[jj].main_atm + pSQ->domain[0];
			if(pSQ->Atm[ii].Ry[jj].main_atm<0)
				pSQ->Atm[ii].Ry[jj].main_atm = pSQ->Atm[ii].Ry[jj].main_atm + pSQ->domain[1];
			if(pSQ->Atm[ii].Rz[jj].main_atm<0)
				pSQ->Atm[ii].Rz[jj].main_atm = pSQ->Atm[ii].Rz[jj].main_atm + pSQ->domain[2];

			if(rank==0 && pSQ->n_atm<=8)
				printf("  %.12f    %.12f    %.12f \n",pSQ->Atm[ii].Rx[jj].main_atm,pSQ->Atm[ii].Ry[jj].main_atm,pSQ->Atm[ii].Rz[jj].main_atm);

		}
	}


}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// NVT MD stuff ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

//Abhiraj: NVT NH ABINIT (ionmov = 8) -> an implicit scheme

void NVT_NH_mov8(DS_SQ* pSQ)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank) ;
	double TE_ext = 0.0, fsnose = 0.0;

	// First step velocity verlet
	vv_1(pSQ);

	// Update thermostat
	fsnose = (pSQ->v2nose - pSQ->dof*pSQ->kB*pSQ->Treq)/pSQ->Ms ;
	pSQ->snose += pSQ->time_step*(pSQ->xi_nose + pSQ->time_step*fsnose/2.0) ;
	pSQ->xi_nose += pSQ->time_step*fsnose/2.0 ;

	// Apply periodic boundary condition
	WrapAround_PBC(pSQ);

	//For better rho_guess to be used in SCF iteration
	ChargeExtrapolation_MD(pSQ);

	//Compute forces at next time step i.e. F(t+delta_t)
	SQDFT_forces(pSQ);

	// Second Step velocity verlet
	pSQ->v2nose = 0.0 ;int k = 0;

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->v2nose += pSQ->Atm[ii].mass*(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]);
			k += 1 ;
		}
	}

	// update accelerations
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->accels[k]=pSQ->forces[k]/pSQ->Atm[ii].mass; 
			pSQ->accels[k+pSQ->n_atm]=pSQ->forces[k+pSQ->n_atm]/pSQ->Atm[ii].mass;
			pSQ->accels[k+2*pSQ->n_atm]=pSQ->forces[k+2*pSQ->n_atm]/pSQ->Atm[ii].mass;
			k+=1;
		}
	}

	// Newton-Raphson Loop
	NR(pSQ);

	// Energy computation
	MD_energies(pSQ);  

	// Ouput data
	if(rank == 0)
		datacollect(pSQ);  

	// Extended System ( Ionic system + Thermostat) energy
	TE_ext += (0.5*pSQ->Ms*pSQ->xi_nose*pSQ->xi_nose + pSQ->dof*pSQ->kB*pSQ->Treq*pSQ->snose)/pSQ->n_atm ;
	TE_ext += pSQ->TE ;
	if(rank==0)
		printf("TE_extended (Ha/atom)        : %.12f \n",TE_ext);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void vv_1(DS_SQ* pSQ)
{
	int k=0;
	pSQ->v2nose = 0.0 ;

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->v2nose += pSQ->Atm[ii].mass*(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]);
			k += 1 ;
		}
	}

	k = 0 ;

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->vels[k] += pSQ->time_step*0.5*(pSQ->accels[k] - pSQ->xi_nose*pSQ->vels[k]);
			pSQ->vels[k+pSQ->n_atm] += pSQ->time_step*0.5*(pSQ->accels[k+pSQ->n_atm] - pSQ->xi_nose*pSQ->vels[k+pSQ->n_atm]);
			pSQ->vels[k+2*pSQ->n_atm] += pSQ->time_step*0.5*(pSQ->accels[k+2*pSQ->n_atm] - pSQ->xi_nose*pSQ->vels[k+2*pSQ->n_atm]);

			if(pSQ->MDstepCount==1) // non-mapped atom coords required for charge extrapolation
			{
				pSQ->Atm[ii].Rx[jj].main_atm_nm = pSQ->Atm[ii].Rx[jj].main_atm;
				pSQ->Atm[ii].Ry[jj].main_atm_nm = pSQ->Atm[ii].Ry[jj].main_atm;
				pSQ->Atm[ii].Rz[jj].main_atm_nm = pSQ->Atm[ii].Rz[jj].main_atm;
			}

			pSQ->Atm[ii].Rx[jj].main_atm += pSQ->time_step*pSQ->vels[k];
			pSQ->Atm[ii].Ry[jj].main_atm += pSQ->time_step*pSQ->vels[k+pSQ->n_atm];
			pSQ->Atm[ii].Rz[jj].main_atm += pSQ->time_step*pSQ->vels[k+2*pSQ->n_atm];

			pSQ->Atm[ii].Rx[jj].main_atm_nm += pSQ->time_step*pSQ->vels[k];
			pSQ->Atm[ii].Ry[jj].main_atm_nm += pSQ->time_step*pSQ->vels[k+pSQ->n_atm];
			pSQ->Atm[ii].Rz[jj].main_atm_nm += pSQ->time_step*pSQ->vels[k+2*pSQ->n_atm];

			k+=1;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NR(DS_SQ* pSQ)
{ 
	double *vel_temp, xin_nose,xio,delxi,*vonose,*hnose,*binose,dnose ;
	int k = 0 ;
	vel_temp = new double [3*pSQ->n_atm]();
	vonose = new double [3*pSQ->n_atm]();
	hnose = new double [3*pSQ->n_atm]();
	binose = new double [3*pSQ->n_atm]();
	xin_nose = pSQ->xi_nose ;

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			vel_temp[k] = pSQ->vels[k] ;
			vel_temp[k+pSQ->n_atm] = pSQ->vels[k+pSQ->n_atm];
			vel_temp[k+2*pSQ->n_atm] = pSQ->vels[k+2*pSQ->n_atm] ;
			k += 1 ;
		}
	}

	int ready = 0 ; // 0- false 1- true
	do
	{
		xio = xin_nose ;
		delxi = 0.0 ;
		k = 0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
			{
				vonose[k] = vel_temp[k] ;
				vonose[k+pSQ->n_atm] = vel_temp[k+pSQ->n_atm] ;
				vonose[k+2*pSQ->n_atm] = vel_temp[k+2*pSQ->n_atm] ;
				hnose[k] = -1.0*pSQ->time_step*0.5*(pSQ->accels[k] - xio*vonose[k]) - (pSQ->vels[k] - vonose[k]);
				hnose[k+pSQ->n_atm] = -1.0*pSQ->time_step*0.5*(pSQ->accels[k+pSQ->n_atm] - xio*vonose[k+pSQ->n_atm]) - (pSQ->vels[k+pSQ->n_atm] - vonose[k+pSQ->n_atm]);
				hnose[k+2*pSQ->n_atm] = -1.0*pSQ->time_step*0.5*(pSQ->accels[k+2*pSQ->n_atm] - xio*vonose[k+2*pSQ->n_atm]) - (pSQ->vels[k+2*pSQ->n_atm] - vonose[k+2*pSQ->n_atm]);
				binose[k] = vonose[k]*pSQ->time_step*pSQ->Atm[ii].mass/pSQ->Ms ;
				delxi += hnose[k]*binose[k] ; 
				binose[k+pSQ->n_atm] = vonose[k+pSQ->n_atm]*pSQ->time_step*pSQ->Atm[ii].mass/pSQ->Ms ;
				delxi += hnose[k+pSQ->n_atm]*binose[k+pSQ->n_atm] ;
				binose[k+2*pSQ->n_atm] = vonose[k+2*pSQ->n_atm]*pSQ->time_step*pSQ->Atm[ii].mass/pSQ->Ms ;
				delxi += hnose[k+2*pSQ->n_atm]*binose[k+2*pSQ->n_atm] ;
				k += 1 ;
			}
		}

		dnose = -1.0*(xio*pSQ->time_step/2.0 + 1.0) ;
		delxi += -1.0*dnose*((-1.0*pSQ->v2nose + pSQ->dof*pSQ->kB*pSQ->Treq)*pSQ->time_step*0.5/pSQ->Ms - (pSQ->xi_nose - xio)) ;
		delxi /= (-1.0*pSQ->time_step*pSQ->time_step*0.5*pSQ->v2nose/pSQ->Ms + dnose) ;
		pSQ->v2nose = 0.0 ;
		k = 0 ;                    
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				vel_temp[k] += (hnose[k] + pSQ->time_step*0.5*vonose[k]*delxi)/dnose ;
				vel_temp[k+pSQ->n_atm] += (hnose[k+pSQ->n_atm] + pSQ->time_step*0.5*vonose[k+pSQ->n_atm]*delxi)/dnose ;
				vel_temp[k+2*pSQ->n_atm] += (hnose[k+2*pSQ->n_atm] + pSQ->time_step*0.5*vonose[k+2*pSQ->n_atm]*delxi)/dnose ;
				pSQ->v2nose += pSQ->Atm[ii].mass*(vel_temp[k]*vel_temp[k]+vel_temp[k+pSQ->n_atm]*vel_temp[k+pSQ->n_atm]+vel_temp[k+2*pSQ->n_atm]*vel_temp[k+2*pSQ->n_atm]);
				k += 1 ;
			}
		}

		xin_nose = xio + delxi ;
		ready = 1 ;

		//Test for convergence
		int kk = -1,jj = 0 ;
		do
		{
			kk = kk + 1 ;
			if(kk >= pSQ->n_atm)
			{
				kk = 0 ;
				jj += 1 ;
			}
			if(kk < pSQ->n_atm && jj<3)
			{
				if(abs((vel_temp[kk + jj*pSQ->n_atm] - vonose[kk + jj*pSQ->n_atm])/vel_temp[kk + jj*pSQ->n_atm]) > pow(10.0,-7.0))
					ready = 0 ;
			}
			else
			{
				if(abs((xin_nose - xio)/xin_nose) > pow(10.0,-7.0))
					ready = 0 ; 
			}
		} while(kk < pSQ->n_atm && jj<3 && ready);
	} while (ready == 0);  

	pSQ->xi_nose = xin_nose ;
	k = 0 ;

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{
			pSQ->vels[k] = vel_temp[k] ;
			pSQ->vels[k+pSQ->n_atm] = vel_temp[k+pSQ->n_atm];
			pSQ->vels[k+2*pSQ->n_atm] = vel_temp[k+2*pSQ->n_atm] ;
			k += 1 ;
		}
	}

	delete [] vel_temp ;
	delete [] vonose ;
	delete [] binose ;
	delete [] hnose ;
}

//********************************* Atomic Relaxation ***************************************//
///////////////////////////////////////////////////////////////////////////////////////////////
//               NLCG: Nonlinear Conjugate Gradient method for ground state minimization     //
//             Reference: An Introduction to the Conjugate Gradient Method Without           //
//                         Agonizing Pain, Jonathan Richard Shewchuk                         //
///////////////////////////////////////////////////////////////////////////////////////////////
void NLCG(DS_SQ *pSQ)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int i=0,j,k=0,imax=100,jmax=6, n=30, jsum=0; // TODO: MAX_ITS = imax
	double deltaNew,deltad,deltaOld,deltaMid,tol1,tol2=1e-10,sigma0=0.5,alpha,etaPrev,eta,beta;
	double *r,*d,*s,*y,*F;  
	int Np=3*pSQ->n_atm;
	r = new double [Np]();
	d = new double [Np]();
	s = new double [Np]();
	y = new double [Np]();
	F = new double [Np]();
	int inCtr;

	tol1 = (1e-10)*3*pSQ->n_atm; // TODO: NLCG_TOL*3*Natm = tol1

	if(rank==0)
		printf("WARNING: NLCG max iterations and tol are hard-coded. (See md.cpp, lines 1144 and 1155.) \n");

	/*
	 * update electron density, energy and forces
	 */  
	SQDFT_forces(pSQ);

	deltaNew=0.0;
	int count=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
		{// r=forces
			r[count]=pSQ->forces[count];
			r[count+pSQ->n_atm]=pSQ->forces[count+pSQ->n_atm];
			r[count+2*pSQ->n_atm]=pSQ->forces[count+2*pSQ->n_atm];

			// s=r
			s[count]=r[count];
			s[count+pSQ->n_atm]=r[count+pSQ->n_atm];
			s[count+2*pSQ->n_atm]=r[count+2*pSQ->n_atm];

			// d=s
			d[count]=s[count];
			d[count+pSQ->n_atm]=s[count+pSQ->n_atm];
			d[count+2*pSQ->n_atm]=s[count+2*pSQ->n_atm];

			// deltaNew=dot(r,d)
			deltaNew+=r[count]*d[count]+r[count+pSQ->n_atm]*d[count+pSQ->n_atm]+r[count+2*pSQ->n_atm]*d[count+2*pSQ->n_atm];
			count+=1;
		}
	}
	inCtr=0;   

	while((i<imax) && (deltaNew >tol1))
	{      

		j=0;
		deltad=0.0;
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				deltad+=d[count]*d[count]+d[count+pSQ->n_atm]*d[count+pSQ->n_atm]+d[count+2*pSQ->n_atm]*d[count+2*pSQ->n_atm];
				count+=1;
			}
		}

		alpha = -sigma0;
		/*
		 * perturb atomic positions
		 */
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				y[count] = pSQ->Atm[ii].Rx[jj].main_atm + sigma0*d[count]; 
				y[count+pSQ->n_atm] = pSQ->Atm[ii].Ry[jj].main_atm + sigma0*d[count+pSQ->n_atm]; 
				y[count+2*pSQ->n_atm] = pSQ->Atm[ii].Rz[jj].main_atm + sigma0*d[count+2*pSQ->n_atm]; 

				// x_new=y
				pSQ->Atm[ii].Rx[jj].main_atm = y[count]; 
				pSQ->Atm[ii].Ry[jj].main_atm = y[count+pSQ->n_atm]; 
				pSQ->Atm[ii].Rz[jj].main_atm = y[count+2*pSQ->n_atm]; 

				count+=1;
			}
		}

		/*
		 * if an atom has moved out of the simulation domain, map it back.
		 * this is only applicable for periodic boundary conditions
		 */
		WrapAround_PBC(pSQ);  

		/*
		 * update electron density, energy and forces
		 */  
		SQDFT_forces(pSQ); 
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{// F=forces
				F[count]=pSQ->forces[count];
				F[count+pSQ->n_atm]=pSQ->forces[count+pSQ->n_atm];
				F[count+2*pSQ->n_atm]=pSQ->forces[count+2*pSQ->n_atm];

				count+=1;
			}
		}

		/*
		 * replace back the original atomic positions
		 */   
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				// x = y - sigma0*d
				pSQ->Atm[ii].Rx[jj].main_atm = y[count] - sigma0*d[count]; 
				pSQ->Atm[ii].Ry[jj].main_atm = y[count+pSQ->n_atm] - sigma0*d[count+pSQ->n_atm]; 
				pSQ->Atm[ii].Rz[jj].main_atm = y[count+2*pSQ->n_atm] - sigma0*d[count+2*pSQ->n_atm]; 
				count+=1;
			}
		}

		/*
		 * if an atom has moved out of the simulation domain, map it back.
		 * this is only applicable for periodic boundary conditions
		 */
		WrapAround_PBC(pSQ); 
		etaPrev=0.0;
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				// etaPrev=dot(F,d)
				etaPrev+=F[count]*d[count]+F[count+pSQ->n_atm]*d[count+pSQ->n_atm]+F[count+2*pSQ->n_atm]*d[count+2*pSQ->n_atm];
				count+=1;
			}
		}
		etaPrev = -etaPrev;

		/*
		 * line search
		 */
		do
		{     
			if(inCtr==0)
			{
				eta=0.0;
				count=0;
				for(int ii=0;ii<pSQ->n_typ;ii++)
				{
					for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
					{
						// eta=dot(r,d)
						eta+=r[count]*d[count]+r[count+pSQ->n_atm]*d[count+pSQ->n_atm]+r[count+2*pSQ->n_atm]*d[count+2*pSQ->n_atm];
						count+=1;
					}
				}
				eta = -eta;
			}
			else{
				/*
				 * update electron density, energy and forces
				 */
				SQDFT_forces(pSQ); 
				count=0;
				for(int ii=0;ii<pSQ->n_typ;ii++)
				{
					for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
					{// F=forces
						F[count]=pSQ->forces[count];
						F[count+pSQ->n_atm]=pSQ->forces[count+pSQ->n_atm];
						F[count+2*pSQ->n_atm]=pSQ->forces[count+2*pSQ->n_atm];

						count+=1;
					}
				}     
				eta=0.0;
				count=0;
				for(int ii=0;ii<pSQ->n_typ;ii++)
				{
					for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
					{
						// eta=dot(F,d)
						eta+=F[count]*d[count]+F[count+pSQ->n_atm]*d[count+pSQ->n_atm]+F[count+2*pSQ->n_atm]*d[count+2*pSQ->n_atm];
						count+=1;
					}
				}
				eta = -eta;
			}

			alpha = alpha*(eta/(etaPrev-eta));      
			/*
			 * perturb atomic positions
			 */
			count=0;
			for(int ii=0;ii<pSQ->n_typ;ii++)
			{
				for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
				{
					// x=x+alpha*d
					pSQ->Atm[ii].Rx[jj].main_atm += alpha*d[count]; 
					pSQ->Atm[ii].Ry[jj].main_atm += alpha*d[count+pSQ->n_atm]; 
					pSQ->Atm[ii].Rz[jj].main_atm += alpha*d[count+2*pSQ->n_atm]; 

					count+=1;
				}
			}  

			/*
			 * if an atom has moved out of the simulation domain, map it back.
			 * this is only applicable for periodic boundary conditions
			 */
			WrapAround_PBC(pSQ);   

			etaPrev = eta;
			j++; inCtr++;

		}while((j<jmax) && (alpha*alpha*deltad>tol2));
		jsum = jsum+j;

		/*
		 * update electron density, energy and forces
		 */
		SQDFT_forces(pSQ); 
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{// r=forces
				r[count]=pSQ->forces[count];
				r[count+pSQ->n_atm]=pSQ->forces[count+pSQ->n_atm];
				r[count+2*pSQ->n_atm]=pSQ->forces[count+2*pSQ->n_atm];

				count+=1;
			}
		}
		inCtr=0;

		deltaOld = deltaNew;

		deltaMid=0.0;
		deltaNew=0.0;
		count=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
			{
				// deltaMid=dot(r,s)
				deltaMid+=r[count]*s[count]+r[count+pSQ->n_atm]*s[count+pSQ->n_atm]+r[count+2*pSQ->n_atm]*s[count+2*pSQ->n_atm];

				// s=r
				s[count]=r[count];
				s[count+pSQ->n_atm]=r[count+pSQ->n_atm];
				s[count+2*pSQ->n_atm]=r[count+2*pSQ->n_atm];

				// deltaNew=dot(r,s)
				deltaNew+=r[count]*s[count]+r[count+pSQ->n_atm]*s[count+pSQ->n_atm]+r[count+2*pSQ->n_atm]*s[count+2*pSQ->n_atm];

				count+=1;
			}
		}

		if(rank==0)
		{
			//printf("%========================================%\n");
			printf("NLCG:outer iter no.=%d,deltaNew=%g \n",i,deltaNew);
			//printf("%========================================%\n");
		}
		beta = (deltaNew-deltaMid)/deltaOld;
		k++;

		if((k==n) || (beta<=0))
		{
			count=0;
			for(int ii=0;ii<pSQ->n_typ;ii++)
			{
				for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
				{
					// d=s
					d[count]=s[count];
					d[count+pSQ->n_atm]=s[count+pSQ->n_atm];
					d[count+2*pSQ->n_atm]=s[count+2*pSQ->n_atm];
					count+=1;
				}
			}
			k=0;
		}
		else
		{
			count=0;
			for(int ii=0;ii<pSQ->n_typ;ii++)
			{
				for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
				{
					// d=s
					d[count]=s[count];
					d[count+pSQ->n_atm]=s[count+pSQ->n_atm];
					d[count+2*pSQ->n_atm]=s[count+2*pSQ->n_atm];
					count+=1;
				}
			}

			count=0;
			for(int ii=0;ii<pSQ->n_typ;ii++)
			{
				for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)     
				{
					// d=s
					d[count]=s[count]+beta*d[count];
					d[count+pSQ->n_atm]=s[count+pSQ->n_atm]+beta*d[count+pSQ->n_atm];
					d[count+2*pSQ->n_atm]=s[count+2*pSQ->n_atm]+beta*d[count+2*pSQ->n_atm];
					count+=1;
				}
			}
		}
		i++;     
	}

	delete [] F;
	delete [] r;
	delete [] d;
	delete [] s;
	delete [] y;

	return;
}


