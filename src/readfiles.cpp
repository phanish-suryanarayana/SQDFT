/** \file readfiles.cpp
  \brief This file contains functions that read and store the given input files. 

  Input files include .input file with method parameters, .atoms file with atoms information and .psd pseudopotential files. 

*/

#include "headers.h"
#include "ds_sq.h"
#include "func_sq.h"
#include <unistd.h>

/** \def M_PI value of \f$ \pi\f$.
*/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/** \def TEMP_TOL a very small value used in some cases involving ceil and floor.
*/
#ifndef TEMP_TOL
#define TEMP_TOL 1e-12
#endif

// function to check inputs - processor domain number etc.
void CheckInputs(DS_SQ* pSQ, char ** argv)    
{
	int rank,count;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char temp; string inptfile;

	if(rank==0)
	{
		printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
		printf("       Spectral Quadrature Density Functional Theory (SQDFT) \n");
		printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
		cout << " "<<endl;
		//timestamp ( );

		char* c_time_str;
		time_t current_time=time(NULL);
		c_time_str=ctime(&current_time);  
		printf("SQDFT version-0 \n");
		printf("Starting time: %s \n",c_time_str);

	}

	double np_temp;
	if(rank==0)
		printf("Assuming CUBICAL DOMAIN. nprocx=nprocy=nprocz \n");
	np_temp = pow(pSQ->nproc,1.0/3.0);

	if(abs(pow(round(np_temp),3)-pSQ->nproc)>TEMP_TOL)
	{
		//if(rank==0)
		printf("Requested number of processors:%u, is not a perfect cube. Exiting. \n",pSQ->nproc);
		pSQ->nprocx=floor(np_temp);pSQ->nprocy=floor(np_temp);pSQ->nprocz=floor(np_temp);
		exit(0);
	}else
	{
		//if(rank==0)
		//printf("Requested number of processors:%u, is a perfect cube. \n",pSQ->nproc);
		pSQ->nprocx=round(np_temp);pSQ->nprocy=round(np_temp);pSQ->nprocz=round(np_temp);
	}
	if(rank==0)
		printf("Choosing %u processors in each direction. Total nproc= %u \n",pSQ->nprocx,pSQ->nprocx*pSQ->nprocy*pSQ->nprocz);
	if(pSQ->nprocx*pSQ->nprocy*pSQ->nprocz>pSQ->nproc)
	{//if(rank==0)
		printf("Error:Inconsistency with no. of processors.");
		exit(0);
	}

	count=1;
	temp = argv[1][5+count];
	while(temp!=0)
	{
		inptfile.append(1u,temp);
		count=count+1;
		temp = argv[1][5+count]; // first index [1] corresponds to 2nd argument of the mpirun command in the pbs file, which has the input file name
	}
	pSQ->input_file = inptfile;
	if(rank==0)
		cout << " "<<endl;
}

// function to read the input file
void Read_input(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	pSQ->mpi_time=0.0;
	double tcomm1,tcomm2;
	char str[80];
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	FILE *input_file;
	char inpt_filename[100]="./";
	const char * cchar = pSQ->input_file.c_str(); // convert string to const char
	strcat(inpt_filename,cchar);
	strcat(inpt_filename,".input");
	if(rank==0) // /////////////
		input_file = fopen(inpt_filename,"r");
	if(rank==0)
	{cout << " "<<endl;
		printf("Reading input file... \n");}

	if(rank==0) // /////////////
	{
		while (!feof(input_file))
		{
			fscanf(input_file,"%s",str);
			if(strcmp(str,"vel_dstr")==0) 
			{
				fscanf(input_file,"%u",&pSQ->vds);
			}else if(strcmp(str,"FD_order")==0) 
			{
				fscanf(input_file,"%u",&pSQ->FDn);
				pSQ->FDn = (pSQ->FDn)/2; // store half order                                              
			}else if(strcmp(str,"domain")==0) 
			{
				fscanf(input_file,"%lf", &pSQ->domain[0]);
				fscanf(input_file,"%lf", &pSQ->domain[1]);
				fscanf(input_file,"%lf", &pSQ->domain[2]);
			}else if(strcmp(str,"n_int")==0)  
			{
				fscanf(input_file,"%u", &pSQ->n_int[0]);
				fscanf(input_file,"%u", &pSQ->n_int[1]);
				fscanf(input_file,"%u", &pSQ->n_int[2]);
			}else if(strcmp(str,"atoms_file")==0)
			{
				fscanf(input_file,"%s",&pSQ->atoms_file[0]);
			}else if(strcmp(str,"ncell")==0)
			{
				fscanf(input_file,"%d",&pSQ->ncell);
			}else if(strcmp(str,"perturb")==0)
			{
				fscanf(input_file,"%lf",&pSQ->perturb);
			}else if(strcmp(str,"rand_seed")==0)
			{
				fscanf(input_file,"%d",&pSQ->rand_seed);
			}else if(strcmp(str,"T")==0)
			{
				fscanf(input_file,"%lf",&pSQ->T);
			}else if(strcmp(str,"npl")==0)
			{
				fscanf(input_file,"%u",&pSQ->npl);
			}else if(strcmp(str,"Rcut")==0)
			{
				fscanf(input_file,"%lf",&pSQ->Rcut);
			}else if(strcmp(str,"poisson_tol")==0)
			{
				fscanf(input_file,"%lf",&pSQ->poisson_tol);
			}else if(strcmp(str,"poisson_maxiter")==0)
			{
				fscanf(input_file,"%d",&pSQ->poisson_maxiter);
			}else if(strcmp(str,"lanczos_tol")==0)
			{
				fscanf(input_file,"%lf",&pSQ->lanczos_tol);
			}else if(strcmp(str,"fermi_tol")==0)
			{
				fscanf(input_file,"%lf",&pSQ->fermi_tol);
			}else if(strcmp(str,"scf_tol")==0)
			{
				fscanf(input_file,"%lf",&pSQ->scf_tol);
			}else if(strcmp(str,"scf_miniter")==0)
			{
				fscanf(input_file,"%d",&pSQ->scf_miniter);
			}else if(strcmp(str,"scf_maxiter")==0)
			{
				fscanf(input_file,"%d",&pSQ->scf_maxiter);
			}else if(strcmp(str,"beta_aaj")==0)
			{
				fscanf(input_file,"%lf",&pSQ->beta_aaj);
			}else if(strcmp(str,"beta_scf")==0)
			{
				fscanf(input_file,"%lf",&pSQ->beta_scf);
			}else if(strcmp(str,"m_aaj")==0)
			{
				fscanf(input_file,"%d",&pSQ->m_aaj);
			}else if(strcmp(str,"p_aaj")==0)
			{
				fscanf(input_file,"%d",&pSQ->p_aaj);
			}else if(strcmp(str,"p_scf")==0)
			{
				fscanf(input_file,"%d",&pSQ->p_scf);
			}else if(strcmp(str,"m_scf")==0)
			{
				fscanf(input_file,"%d",&pSQ->m_scf);
			}else if(strcmp(str,"non_blocking")==0)
			{
				fscanf(input_file,"%d",&pSQ->non_blocking);
			}else if(strcmp(str,"prnt_atoms")==0)
			{
				fscanf(input_file,"%d",&pSQ->prnt_atoms);
			}else if(strcmp(str,"MaxMDsteps")==0)
			{
				fscanf(input_file,"%d",&pSQ->MaxMDsteps);
			}else if(strcmp(str,"time_step")==0)
			{
				fscanf(input_file,"%lf",&pSQ->time_step);
			}else if(strcmp(str,"ChgExtrap")==0)
			{
				fscanf(input_file,"%d",&pSQ->ChgExtrap);
			}else if(strcmp(str,"restart_scf")==0)
			{
				fscanf(input_file,"%d",&pSQ->restart_scf);
			}else if(strcmp(str,"prnt_scf")==0)
			{
				fscanf(input_file,"%d",&pSQ->prnt_scf);
			}else if(strcmp(str,"restart_md")==0)
			{
				fscanf(input_file,"%d",&pSQ->restart_md);
			}else if(strcmp(str,"prnt_md")==0)
			{
				fscanf(input_file,"%d",&pSQ->prnt_md);
			}else if(strcmp(str,"RelaxAtoms")==0)
			{
				fscanf(input_file,"%d",&pSQ->RelaxAtoms);
			}else if(strcmp(str,"ensemble")==0)
			{
				fscanf(input_file,"%d",&pSQ->ensemble);
			}else if(strcmp(str,"name")==0)
			{
				fscanf(input_file,"%s",&pSQ->name[0]);
			}else if(strcmp(str,"fermi_alg")==0)
			{
				fscanf(input_file,"%d",&pSQ->fermi_alg);
			}else if(strcmp(str,"qmass")==0)
			{
				fscanf(input_file,"%lf",&pSQ->Ms);
			}else if(strcmp(str,"Correction")==0)
			{
				fscanf(input_file,"%d",&pSQ->Correction);
			}

		}
		fclose(input_file);
	}

	tcomm1 = MPI_Wtime();

	// ----------------- Bcast ints together and doubles together (two MPI_Bcast 's) --------------------
	int bcast_int[27]={pSQ->FDn,pSQ->n_int[0],pSQ->n_int[1],pSQ->n_int[2],pSQ->npl,pSQ->poisson_maxiter,pSQ->scf_maxiter, pSQ->m_aaj,pSQ->p_aaj,pSQ->p_scf,pSQ->non_blocking,pSQ->ncell,pSQ->m_scf,pSQ->prnt_atoms,pSQ->rand_seed,pSQ->MaxMDsteps,pSQ->ChgExtrap,pSQ->restart_scf,pSQ->prnt_scf,pSQ->restart_md,pSQ->prnt_md,pSQ->RelaxAtoms,pSQ->vds,pSQ->scf_miniter,pSQ->ensemble,pSQ->fermi_alg,pSQ->Correction};
	double bcast_double[14]={pSQ->domain[0],pSQ->domain[1],pSQ->domain[2],pSQ->T,pSQ->Rcut,pSQ->poisson_tol,pSQ->lanczos_tol,pSQ->fermi_tol,pSQ->scf_tol,pSQ->beta_aaj,pSQ->beta_scf,pSQ->perturb,pSQ->time_step,pSQ->Ms};
	MPI_Bcast(bcast_int,27,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(bcast_double,14,MPI_DOUBLE,0,MPI_COMM_WORLD);

	pSQ->FDn=bcast_int[0] ;
	pSQ->n_int[0]=bcast_int[1] ;
	pSQ->n_int[1]=bcast_int[2] ;
	pSQ->n_int[2]=bcast_int[3] ;
	pSQ->npl=bcast_int[4] ;
	pSQ->poisson_maxiter=bcast_int[5];
	pSQ->scf_maxiter=bcast_int[6] ;
	pSQ->m_aaj=bcast_int[7] ;
	pSQ->p_aaj=bcast_int[8] ;
	pSQ->p_scf=bcast_int[9] ;
	pSQ->non_blocking=bcast_int[10] ;
	pSQ->ncell=bcast_int[11] ;
	pSQ->m_scf=bcast_int[12];
	pSQ->prnt_atoms=bcast_int[13];
	pSQ->rand_seed=bcast_int[14];
	pSQ->MaxMDsteps=bcast_int[15];
	pSQ->ChgExtrap=bcast_int[16];
	pSQ->restart_scf=bcast_int[17];
	pSQ->prnt_scf=bcast_int[18];
	pSQ->restart_md=bcast_int[19];
	pSQ->prnt_md=bcast_int[20];
	pSQ->RelaxAtoms=bcast_int[21];
	pSQ->vds=bcast_int[22];
	pSQ->scf_miniter=bcast_int[23];
	pSQ->ensemble = bcast_int[24];
	pSQ->fermi_alg = bcast_int[25];
	pSQ->Correction = bcast_int[26];

	pSQ->domain[0]=bcast_double[0] ;
	pSQ->domain[1]=bcast_double[1] ;
	pSQ->domain[2]=bcast_double[2] ;
	pSQ->T=bcast_double[3] ;
	pSQ->Rcut=bcast_double[4] ;
	pSQ->poisson_tol=bcast_double[5] ;
	pSQ->lanczos_tol=bcast_double[6] ;
	pSQ->fermi_tol=bcast_double[7] ;
	pSQ->scf_tol=bcast_double[8] ;
	pSQ->beta_aaj=bcast_double[9] ;
	pSQ->beta_scf=bcast_double[10] ;
	pSQ->perturb=bcast_double[11] ;
	pSQ->time_step=bcast_double[12] ;
	pSQ->Ms=bcast_double[13] ;

	MPI_Bcast(pSQ->atoms_file,10,MPI_CHAR,0,MPI_COMM_WORLD);

	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;

	// Compute mesh size delta
	double delta_x = pSQ->domain[0]/pSQ->n_int[0];
	double delta_y = pSQ->domain[1]/pSQ->n_int[1];
	double delta_z = pSQ->domain[2]/pSQ->n_int[2];
	if((fabs(delta_x-delta_y) >=TEMP_TOL) || (fabs(delta_x-delta_z) >=TEMP_TOL) || (fabs(delta_y-delta_z) >=TEMP_TOL))
	{      if(rank==0)
		printf("Error: mesh spacing MUST be same in all directions. \n");
		exit(0);
	}else
	{ 
		pSQ->delta = delta_x;
	}

	pSQ->nloc = ceil(pSQ->Rcut/pSQ->delta);
	pSQ->Rcut = pSQ->nloc*pSQ->delta;

	pSQ->kB = 8.617343e-5; // Boltzmann constant in eV/K
	pSQ->Ceh= 27.211384523; // conversion from eV to Ha, 1 Ha=27.211384523 eV, so 1eV=1/Ceh Ha
	pSQ->kB = pSQ->kB/pSQ->Ceh ; 

	if(pSQ->T<10000.0 && rank==0)
		printf("WARNING: Need to use Brent's method for Fermi energy calculation for low smearing.\n"); // Using Newton-Raphson at low temperatures might fail since the Fermi-Dirac becomes steeper and derivative might be inaccurate. See sq.cpp, line 862.

	pSQ->bet = (1/((pSQ->kB)*pSQ->T)); // Inverse of smearing in 1/Ha
	pSQ->latconst = pSQ->domain[0]; // lattice constant of unit cell 
	pSQ->domain[0] = (double)pSQ->domain[0]*pSQ->ncell;
	pSQ->domain[1] = (double)pSQ->domain[1]*pSQ->ncell;
	pSQ->domain[2] = (double)pSQ->domain[2]*pSQ->ncell;
	pSQ->n_int[0] = pSQ->n_int[0]*pSQ->ncell;
	pSQ->n_int[1] = pSQ->n_int[1]*pSQ->ncell;
	pSQ->n_int[2] = pSQ->n_int[2]*pSQ->ncell;

	if(rank==0)
	{
		printf("  Atoms       : %s \n",pSQ->atoms_file);
		printf("  ncell       : %d \n",pSQ->ncell); 
		printf("  domain      : [%.8f %.8f %.8f] Bohr \n",pSQ->domain[0],pSQ->domain[1],pSQ->domain[2]);
		printf("  n_int       : [%d %d %d] \n",pSQ->n_int[0],pSQ->n_int[1],pSQ->n_int[2]);
		printf("  delta       : %.14f \n",pSQ->delta);
		printf("  unit cell   : %.4f Bohr\n",pSQ->latconst);
		printf("  max perturb : %.3f Bohr\n",pSQ->perturb);
		printf("  rand_seed   : %d \n",pSQ->rand_seed); 
		printf("  T (K)       : %.8f \n",pSQ->T);
		printf("  npl         : %u \n",pSQ->npl); 
		printf("  Rcut        : %.6f Bohr\n",pSQ->Rcut); 
		printf("  nloc        : %d \n",pSQ->nloc);
		printf("  FD order    : %u \n",2*pSQ->FDn);
		printf("  Poisson Tol : %g \n",pSQ->poisson_tol);
		printf("  Poisson MaxIter: %d \n",pSQ->poisson_maxiter);
		printf("  Lanczos Tol : %g \n",pSQ->lanczos_tol);
		printf("  Fermi Tol   : %g \n",pSQ->fermi_tol);
		printf("  SCF Tol     : %g \n",pSQ->scf_tol);
		printf("  SCF MinIter : %d \n",pSQ->scf_miniter);
		printf("  SCF MaxIter : %d \n",pSQ->scf_maxiter);
		printf("  beta_aaj    : %.2f \n",pSQ->beta_aaj);
		printf("  m_aaj       : %d \n",pSQ->m_aaj); 
		printf("  p_aaj       : %d \n",pSQ->p_aaj);
		printf("  beta_scf    : %.2f \n",pSQ->beta_scf);
		printf("  m_scf       : %d \n",pSQ->m_scf); 
		printf("  p_scf       : %d \n",pSQ->p_scf); 
		printf("  non_blocking: %d \n",pSQ->non_blocking); 
		printf("  prnt_atoms  : %d \n",pSQ->prnt_atoms); 
		printf("  time_step   : %.4f fs\n",pSQ->time_step); 
		printf("  MaxMDsteps  : %d \n",pSQ->MaxMDsteps);
		printf("  MD time     : %.4f fs\n",pSQ->time_step*pSQ->MaxMDsteps);
		printf("  ChgExtrap   : %d \n",pSQ->ChgExtrap);
		printf("  restart_scf : %d \n",pSQ->restart_scf);
		printf("  prnt_scf    : %d \n",pSQ->prnt_scf);
		printf("  restart_md  : %d \n",pSQ->restart_md);
		printf("  prnt_md     : %d \n",pSQ->prnt_md);
		printf("  RelaxAtoms  : %d \n",pSQ->RelaxAtoms);
		if (pSQ->vds == 1)
			printf("  Vel Dstr    : Maxwell-Boltzmann\n");
		else if (pSQ->vds ==2)
			printf("  Vel Dstr    : Uniform\n");
	}   

	// Compute Finite Difference coefficients for Laplacian (Del^2)
	pSQ->coeff_lap = new double[pSQ->FDn+1]; // need to de-allocate later
	pSQ->coeff_lap[0] = 0;
	for(int a=1; a<=pSQ->FDn; a++)
		pSQ->coeff_lap[0]+= -(2.0/(a*a));
	pSQ->coeff_lap[0]*=1/(pSQ->delta*pSQ->delta);

	for(int p=1; p<=pSQ->FDn; p++)
	{
		pSQ->coeff_lap[p] = (2*pow(-1,p+1)*fract(pSQ->FDn,p)/(p*p*pSQ->delta*pSQ->delta)); 

	}

	// Compute Finite Difference coefficients for Gradient (Del,only RHS stencil)
	pSQ->coeff_grad = new double[pSQ->FDn+1]; // need to de-allocate later
	pSQ->coeff_grad[0] = 0; // NOTE: THIS COEFF IS NOT NEEDED, but using it as zero for convenience
	for(int p=1; p<=pSQ->FDn; p++)
	{
		pSQ->coeff_grad[p] = (pow(-1,p+1)*fract(pSQ->FDn,p)/(p*pSQ->delta));
	}

}

// function required to compute the factorial terms of FD coeffs
double fract(int n,int k)
{
	int i;
	double Nr=1, Dr=1, val;

	for(i=n-k+1; i<=n; i++)
		Nr*=i;
	for(i=n+1; i<=n+k; i++)
		Dr*=i;
	val = Nr/Dr;

	return (val);
}

// function to read the atoms file
void Read_atoms(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double tcomm1,tcomm2;
	double *atom_pos;
	double randx,randy,randz;
	int natm_temp,pp,qq,rr,jj,kk,count,bcast_temp[3];
	char str[80];
	FILE *input_file=NULL;
	char inpt_filename[100]="./";
	strcat(inpt_filename,pSQ->atoms_file);
	strcat(inpt_filename,".atoms");
	if(rank==0) // ////////////////
		input_file = fopen(inpt_filename,"r");

	if(rank==0) // ////////////////
	{
		printf("Reading %s file... \n",inpt_filename);
		fscanf(input_file,"%u",&pSQ->n_atm); // store number of atoms (first line)  
		printf("%d atoms in unit cell. Replicating unit cell %d times in each direction... \n",pSQ->n_atm,pSQ->ncell);      
		pSQ->n_atm = pSQ->n_atm*pow(pSQ->ncell,3);
		printf("Total number of atoms = %d \n",pSQ->n_atm);
	}

	tcomm1 = MPI_Wtime();
	MPI_Bcast(&pSQ->n_atm,1,MPI_INT,0,MPI_COMM_WORLD);
	tcomm2 = MPI_Wtime();
	pSQ->mpi_time+=tcomm2-tcomm1;

	pSQ->dof = 3*(pSQ->n_atm - 1) ;
	double latconst=1.0;
	pSQ->frac_coord=0; // default=0

	int flag=0;
	while (flag==0) //(!feof(input_file))
	{
		if(rank==0) // ////////////////////
			fscanf(input_file,"%s",str);
		tcomm1 = MPI_Wtime();
		MPI_Bcast(str,80,MPI_CHAR,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;
		if(strcmp(str,"n_typ")==0)  
		{
			if(rank==0) // ///////////////////
			{
				fscanf(input_file,"%u",&pSQ->n_typ);
				printf("n_typ: %u \n",pSQ->n_typ);

			}
			tcomm1 = MPI_Wtime();
			MPI_Bcast(&pSQ->n_typ,1,MPI_INT,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;
			pSQ->Atm = new DS_Atm [pSQ->n_typ]; // need to de-allocate later   
		}

		if(strcmp(str,"frac_coord")==0) 
		{
			if(rank==0) // ///////////////////
			{
				fscanf(input_file,"%u",&pSQ->frac_coord);
				printf("frac_coord: %u \n",pSQ->frac_coord);
			}
			tcomm1 = MPI_Wtime();
			MPI_Bcast(&pSQ->frac_coord,1,MPI_INT,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;   
			if(pSQ->frac_coord==1)
			{latconst=pSQ->latconst;
				if(rank==0)
					printf("Input is in fractional coordinates! \n");
			}
		}

		if(strcmp(str,"Atoms")==0)  
		{ 
			for(int ii=0;ii<pSQ->n_typ;ii++)
			{
				if(rank==0) // ////////////////////
				{
					fscanf(input_file,"%s",pSQ->Atm[ii].symb);       
					fscanf(input_file,"%u",&pSQ->Atm[ii].Z);  
					fscanf(input_file,"%lf",&pSQ->Atm[ii].mass);      
					fscanf(input_file,"%u",&pSQ->Atm[ii].natm);       
					fscanf(input_file,"%u",&pSQ->Atm[ii].lloc);       
					fscanf(input_file,"%lf",&pSQ->Atm[ii].rb);

					fscanf(input_file,"%s",pSQ->Atm[ii].pseudo_path);
				}
				tcomm1 = MPI_Wtime();
				MPI_Bcast(&pSQ->Atm[ii].symb,10,MPI_CHAR,0,MPI_COMM_WORLD);
				MPI_Bcast(&pSQ->Atm[ii].mass,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Bcast(&pSQ->Atm[ii].rb,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Bcast(&pSQ->Atm[ii].pseudo_path,1000,MPI_CHAR,0,MPI_COMM_WORLD);

				bcast_temp[0]=pSQ->Atm[ii].Z;bcast_temp[1]=pSQ->Atm[ii].natm;bcast_temp[2]=pSQ->Atm[ii].lloc;
				MPI_Bcast(bcast_temp,3,MPI_INT,0,MPI_COMM_WORLD);
				tcomm2 = MPI_Wtime();
				pSQ->mpi_time+=tcomm2-tcomm1;
				pSQ->Atm[ii].Z=bcast_temp[0];
				pSQ->Atm[ii].natm=bcast_temp[1];
				pSQ->Atm[ii].lloc=bcast_temp[2];
				pSQ->Atm[ii].rb = ceil(pSQ->Atm[ii].rb/pSQ->delta)*pSQ->delta;
				if(rank==0)
				{printf("%s %u %.8f %u %u %.4f\n",pSQ->Atm[ii].symb,pSQ->Atm[ii].Z,pSQ->Atm[ii].mass,pSQ->Atm[ii].natm,pSQ->Atm[ii].lloc,pSQ->Atm[ii].rb);printf("%s \n",pSQ->Atm[ii].pseudo_path);}

				atom_pos = new double [3*pSQ->Atm[ii].natm];

				for(jj=0;jj<pSQ->Atm[ii].natm;jj++)
				{
					if(rank==0)
					{
						fscanf(input_file,"%lf",&atom_pos[jj]);     
						fscanf(input_file,"%lf",&atom_pos[jj+pSQ->Atm[ii].natm]);
						fscanf(input_file,"%lf",&atom_pos[jj+2*pSQ->Atm[ii].natm]);
					}
				}    
				tcomm1 = MPI_Wtime();
				MPI_Bcast(atom_pos,3*pSQ->Atm[ii].natm,MPI_DOUBLE,0,MPI_COMM_WORLD);  
				tcomm2 = MPI_Wtime();
				pSQ->mpi_time+=tcomm2-tcomm1;

				// ---------------- replicate unit cell ncell times in each direction -------------
				natm_temp = pSQ->Atm[ii].natm;
				pSQ->Atm[ii].natm = pSQ->Atm[ii].natm*pow(pSQ->ncell,3);

				pSQ->Atm[ii].Rx = new DS_coord[pSQ->Atm[ii].natm]; // need to de-allocate later
				pSQ->Atm[ii].Ry = new DS_coord[pSQ->Atm[ii].natm]; // need to de-allocate later
				pSQ->Atm[ii].Rz = new DS_coord[pSQ->Atm[ii].natm]; // need to de-allocate later

				srand(ii+pSQ->rand_seed); // random seed, NOTE: Cannot use "srand(time(NULL))" as we want the seed to be the same across all the processors. Instead can have an option to input a different integer (rand_seed) through .inpt file and do "srand(ii+rand_seed)". NOTE: Call srand() before the following loop, otherwise, different procs are generating different sequence

				if(pSQ->frac_coord==1 && rank==0)
				{printf("Note: Atom positions are in fractional coordinates. \n");}

				jj=0;
				for(pp=0;pp<pSQ->ncell;pp++)
				{
					for(qq=0;qq<pSQ->ncell;qq++)
					{
						for(rr=0;rr<pSQ->ncell;rr++)
						{
							count=0;
							for(kk=0;kk<natm_temp;kk++)
							{
								randx = (2*((double)(rand()) / (double)(RAND_MAX))-1); // random between -1 and +1
								randy = (2*((double)(rand()) / (double)(RAND_MAX))-1);
								randz = (2*((double)(rand()) / (double)(RAND_MAX))-1);

								pSQ->Atm[ii].Rx[jj].main_atm = (double)(pp)*pSQ->latconst + latconst*atom_pos[count] + randx*fabs(pSQ->perturb);
								pSQ->Atm[ii].Ry[jj].main_atm = (double)(qq)*pSQ->latconst + latconst*atom_pos[count+natm_temp] + randy*fabs(pSQ->perturb);;
								pSQ->Atm[ii].Rz[jj].main_atm = (double)(rr)*pSQ->latconst + latconst*atom_pos[count+2*natm_temp] + randz*fabs(pSQ->perturb);;
								count+=1;
								if(rank==0)
								{
									if(pSQ->n_atm<=8)
									{ printf("%.12f %.12f %.12f \n",pSQ->Atm[ii].Rx[jj].main_atm/latconst,pSQ->Atm[ii].Ry[jj].main_atm/latconst,pSQ->Atm[ii].Rz[jj].main_atm/latconst);
									}else
									{
										if(jj==0)
											printf("Total number of atoms is greater than 8. Not printing atom coordinates.\n");
									}
								}
								jj+=1;
							}

						}
					}
				}

				delete [] atom_pos;

			}
			break; // After the "Atoms" string, read all atom data and then break out of the while loop
		}

		if(rank==0 && feof(input_file))
			flag=1;   

		tcomm1 = MPI_Wtime();      
		MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;
	}

	if(rank==0) // ////////////////////
	{
		fclose(input_file);
		cout <<"Completed reading atoms!" << "\n";
		cout << " "<<endl;  
	}
}

// function to read the pseudopotential files
void Read_psd(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	double *bcast_vec;
	double bcast_temp[3];
	double tcomm1,tcomm2;
	pSQ->Psd = new DS_Psd [pSQ->n_typ]; // need to de-allocate later 
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		if(rank==0)
			printf("Reading pseudopotential file for %s... ",pSQ->Atm[ii].symb);

		FILE *input_file=NULL;
		char filename[1000];
		sprintf(filename,"%s",pSQ->Atm[ii].pseudo_path);
		if(rank==0) // ////////////////////
			input_file = fopen(filename,"r");

		int count=0; // number of elements in the radial grid 
		char str[60];
		//first check if pseudopotential file is correct
		if(rank==0) // ///////////////////
			fscanf(input_file,"%s",str);

		if(rank==0) // ///////////////////
			if(strcmp(str,pSQ->Atm[ii].symb)!=0)
			{
				if(rank==0)
					cout << "Error: Pseudopotential file does not match with input atom type: " << ii <<endl;
				exit(1);
			}

		if(rank==0) // ///////////////////
		{ 
			do
			{
				fgets(str,60,input_file);            
			}while(strcmp(str," Radial grid follows\n"));
			do
			{
				fscanf(input_file,"%s",str);
				count++;
			}while(strcmp(str,"Pseudopotential")!=0);

		} 

		tcomm1 = MPI_Wtime();     
		MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;

		pSQ->Psd[ii].size = count;  // store the size of the pseudopotential list for future
		if(rank==0) // ///////////////////
			fclose(input_file);

		// set rc_s,p,d,f to zero for each at type
		pSQ->Psd[ii].rc_s=0.0;
		pSQ->Psd[ii].rc_p=0.0;
		pSQ->Psd[ii].rc_d=0.0;
		pSQ->Psd[ii].rc_f=0.0;

		//allocate memory for arrays storing radial grid, pseudopotentials and pseudowavefunctions 
		bcast_vec = new double [11*pSQ->Psd[ii].size](); // need to de-allocate later 
		pSQ->Psd[ii].RadialGrid = new double [pSQ->Psd[ii].size](); // need to de-allocate later 
		if(pSQ->Psd[ii].RadialGrid == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].RadialGrid"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rVloc = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].rVloc == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Vloc"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rVs = new double [pSQ->Psd[ii].size](); // need to de-allocate later  
		if(pSQ->Psd[ii].rVs == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Vs"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rVp = new double [pSQ->Psd[ii].size](); // need to de-allocate later 
		if(pSQ->Psd[ii].rVp == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Vp"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rVd = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].rVd == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Vd"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rVf = new double [pSQ->Psd[ii].size](); // need to de-allocate later 
		if(pSQ->Psd[ii].rVf == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Vf"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].rUs = new double [pSQ->Psd[ii].size](); // need to de-allocate later 
		if(pSQ->Psd[ii].rUs == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Us"<< endl;
			exit(1);
		} 

		pSQ->Psd[ii].rUp = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].rUp == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Up"<< endl;
			exit(1);
		} 

		pSQ->Psd[ii].rUd = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].rUd == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Ud"<< endl;
			exit(1);
		} 

		pSQ->Psd[ii].rUf = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].rUf == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].Uf"<< endl;
			exit(1);
		}

		pSQ->Psd[ii].uu = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		if(pSQ->Psd[ii].uu == NULL)
		{
			if(rank==0)
				cout << "Memory allocation failed in pSQ->Psd[ii].uu"<< endl;
			exit(1);
		}
		////////////////////////////////////////////////////////////////////////
		// open file again and read the pseudopotentials and pseudo wave functions now 
		////////////////////////////////////////////////////////////////////////  

		if(rank==0) // //////////////////
			input_file = fopen(filename,"r");

		double value=0.0;
		int l;
		double rc;
		int ir;
		int lineCtr;
		double Pi=M_PI;

		if(rank==0) // /////////////////
			do
			{
				fgets(str,60,input_file);            
			}while(strcmp(str," Radial grid follows\n"));

		pSQ->Psd[ii].RadialGrid[0]=0; // include the point r=0
		for(ir=1;ir<pSQ->Psd[ii].size;ir++)
		{   
			if(rank==0)
				fscanf(input_file,"%lf",&value);        
			//MPI_Bcast(&value,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			pSQ->Psd[ii].RadialGrid[ir]=value;  
			// cout << pSQ->Psd[ii].RadialGrid[ir]<<endl;           
		}

		lineCtr =0;
		if(rank==0) // ///////////////////
			while(strcmp(str," Pseudopotential follows (l on next line)\n"))       
			{          
				fgets(str,60,input_file);
				lineCtr++;
			}
		//MPI_Bcast(&lineCtr,1,MPI_INT,0,MPI_COMM_WORLD);
		tcomm1 = MPI_Wtime();         
		MPI_Bcast(str,60,MPI_CHAR,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;

		while(strcmp(str," Pseudopotential follows (l on next line)\n")==0)       
		{          
			if(rank==0) // ///////////////////
				fscanf(input_file,"%u",&l);    
			tcomm1 = MPI_Wtime(); 
			MPI_Bcast(&l,1,MPI_INT,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;
			if(l==0) //s orbital
			{
				pSQ->Psd[ii].lmax=0;              
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // /////////////////
						fscanf(input_file,"%lf",&value);
					value = 0.5*value; // To convert Rydberg constant into Hartree energy
					pSQ->Psd[ii].rVs[ir]=value; 

					if(pSQ->Atm[ii].lloc==0)
						pSQ->Psd[ii].rVloc[ir] = value;
				}
				pSQ->Psd[ii].rVs[0]=pSQ->Psd[ii].rVs[1];
				if(pSQ->Atm[ii].lloc==0)
					pSQ->Psd[ii].rVloc[0] = pSQ->Psd[ii].rVs[0];
			}
			if(l==1) //p orbital
			{
				pSQ->Psd[ii].lmax=1;              
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					value = 0.5*value; // /(pSQ->Psd[ii].RadialGrid[ir]);
					pSQ->Psd[ii].rVp[ir]=value; 
					if(pSQ->Atm[ii].lloc==1)
						pSQ->Psd[ii].rVloc[ir] = value;
				}
				pSQ->Psd[ii].rVp[0]=pSQ->Psd[ii].rVp[1];
				if(pSQ->Atm[ii].lloc==1)
					pSQ->Psd[ii].rVloc[0] = pSQ->Psd[ii].rVp[0];
			}
			if(l==2) //d orbital
			{
				pSQ->Psd[ii].lmax=2;              
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					value = 0.5*value; // /(pSQ->Psd[ii].RadialGrid[ir]);
					pSQ->Psd[ii].rVd[ir]=value; 
					if(pSQ->Atm[ii].lloc==2)
						pSQ->Psd[ii].rVloc[ir] = value;
				}
				pSQ->Psd[ii].rVd[0]=pSQ->Psd[ii].rVd[1];
				if(pSQ->Atm[ii].lloc==2)
					pSQ->Psd[ii].rVloc[0] = pSQ->Psd[ii].rVd[0];
			}
			if(l==3) //f orbital
			{
				pSQ->Psd[ii].lmax=3;              
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					value = 0.5*value; // /(pSQ->Psd[ii].RadialGrid[ir]);
					pSQ->Psd[ii].rVf[ir]=value; 
					if(pSQ->Atm[ii].lloc==3)
						pSQ->Psd[ii].rVloc[ir] = value;
				}
				pSQ->Psd[ii].rVf[0]=pSQ->Psd[ii].rVf[1];
				if(pSQ->Atm[ii].lloc==3)
					pSQ->Psd[ii].rVloc[0] = pSQ->Psd[ii].rVf[0];
			}
			// detect string in file
			if(rank==0) // ////////////////////
			{
				for(ir=0;ir<lineCtr;ir++)
					fgets(str,60,input_file); 
			}
			tcomm1 = MPI_Wtime(); 
			MPI_Bcast(str,60,MPI_CHAR,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;
		}
		// read until valence charge block is found
		if(rank==0) // /////////////////////
			while(strcmp(str," Valence charge follows\n"))       
			{       
				fgets(str,60,input_file);                   
			}

		// valence charge read
		if(rank==0) // ///////////////////
			while(strcmp(str," Valence charge follows\n")==0)       
			{          
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					value = value/(4*Pi*pSQ->Psd[ii].RadialGrid[ir]*pSQ->Psd[ii].RadialGrid[ir]);
					pSQ->Psd[ii].uu[ir]=value;
				}
				pSQ->Psd[ii].uu[0]=pSQ->Psd[ii].uu[1];

				// detect string in file
				if(rank==0) // //////////////////////
					for(ir=0;ir<lineCtr;ir++)
						fgets(str,60,input_file); 
			}
		tcomm1 = MPI_Wtime(); 
		MPI_Bcast(str,60,MPI_CHAR,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;
		// pseudowavefunction read
		while(strcmp(str," Pseudo-wave-function follows (l, zelect, rc)\n")==0)       
		{          
			if(rank==0) // //////////////////////
			{
				fscanf(input_file,"%d",&l);              
				fscanf(input_file,"%lf",&value);
				fscanf(input_file,"%lf",&rc); // rc    
			}
			bcast_temp[0]=l;bcast_temp[1]=value;bcast_temp[2]=rc;
			tcomm1 = MPI_Wtime(); 
			MPI_Bcast(bcast_temp,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;
			l=bcast_temp[0];value=bcast_temp[1];rc=bcast_temp[2];

			// MPI_Barrier(MPI_COMM_WORLD);
			if(l==0) //s orbital
			{
				pSQ->Psd[ii].rc_s=rc;
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);

					pSQ->Psd[ii].rUs[ir]=value;                         
				}
				pSQ->Psd[ii].rUs[0]=pSQ->Psd[ii].rUs[1];
			}
			if(l==1) //p orbital
			{
				pSQ->Psd[ii].rc_p=rc;
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					pSQ->Psd[ii].rUp[ir]=value;                         
				}
				pSQ->Psd[ii].rUp[0]=pSQ->Psd[ii].rUp[1];
			}
			if(l==2) //d orbital
			{
				pSQ->Psd[ii].rc_d=rc;
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					pSQ->Psd[ii].rUd[ir]=value;                         
				}
				pSQ->Psd[ii].rUd[0]=pSQ->Psd[ii].rUd[1];
			}
			if(l==3) //f orbital
			{
				pSQ->Psd[ii].rc_f=rc;
				for(ir=1;ir<pSQ->Psd[ii].size;ir++)
				{
					if(rank==0) // //////////////////////
						fscanf(input_file,"%lf",&value);
					pSQ->Psd[ii].rUf[ir]=value;                         
				}
				pSQ->Psd[ii].rUf[0]=pSQ->Psd[ii].rUf[1];
			}

			// detect string in file
			if(rank==0) // //////////////////////
				for(ir=0;ir<lineCtr;ir++)
				{
					if(feof(input_file))
						break;
					fgets(str,60,input_file);
				}
			tcomm1 = MPI_Wtime(); 
			MPI_Bcast(str,60,MPI_CHAR,0,MPI_COMM_WORLD);
			tcomm2 = MPI_Wtime();
			pSQ->mpi_time+=tcomm2-tcomm1;
		}

		if(rank==0)
			fclose(input_file);

		// ---------------------- Bcast ------------------------- 
		lineCtr=0;
		for(ir=0;ir<=(pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].RadialGrid[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=1*pSQ->Psd[ii].size;ir<=(2*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rVloc[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=2*pSQ->Psd[ii].size;ir<=(3*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rVs[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=3*pSQ->Psd[ii].size;ir<=(4*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rVp[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=4*pSQ->Psd[ii].size;ir<=(5*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rVd[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=5*pSQ->Psd[ii].size;ir<=(6*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rVf[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=6*pSQ->Psd[ii].size;ir<=(7*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rUs[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=7*pSQ->Psd[ii].size;ir<=(8*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rUp[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=8*pSQ->Psd[ii].size;ir<=(9*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rUd[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=9*pSQ->Psd[ii].size;ir<=(10*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].rUf[lineCtr];lineCtr++;
		}
		lineCtr=0;
		for(ir=10*pSQ->Psd[ii].size;ir<=(11*pSQ->Psd[ii].size-1);ir++)
		{
			bcast_vec[ir]=pSQ->Psd[ii].uu[lineCtr];lineCtr++;
		}

		tcomm1 = MPI_Wtime(); 
		MPI_Bcast(bcast_vec,11*pSQ->Psd[ii].size,MPI_DOUBLE,0,MPI_COMM_WORLD);
		tcomm2 = MPI_Wtime();
		pSQ->mpi_time+=tcomm2-tcomm1;
		lineCtr=0;
		for(ir=0;ir<=(pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].RadialGrid[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=1*pSQ->Psd[ii].size;ir<=(2*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rVloc[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=2*pSQ->Psd[ii].size;ir<=(3*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rVs[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=3*pSQ->Psd[ii].size;ir<=(4*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rVp[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=4*pSQ->Psd[ii].size;ir<=(5*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rVd[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=5*pSQ->Psd[ii].size;ir<=(6*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rVf[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=6*pSQ->Psd[ii].size;ir<=(7*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rUs[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=7*pSQ->Psd[ii].size;ir<=(8*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rUp[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=8*pSQ->Psd[ii].size;ir<=(9*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rUd[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=9*pSQ->Psd[ii].size;ir<=(10*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].rUf[lineCtr]=bcast_vec[ir];lineCtr++;
		}
		lineCtr=0;
		for(ir=10*pSQ->Psd[ii].size;ir<=(11*pSQ->Psd[ii].size-1);ir++)
		{
			pSQ->Psd[ii].uu[lineCtr]=bcast_vec[ir];lineCtr++;
		}

		delete [] bcast_vec;

		// ---------------------------------------------------------

		// compute rz-cutoff radius at which Vloc become Z/r
		ir=0;
		while(fabs(pSQ->Psd[ii].rVloc[ir]+(double)pSQ->Atm[ii].Z)>TEMP_TOL) // use '+' since Z is +ve but rVloc becomes -ve Z
		{
			ir+=1;
		}
		pSQ->Psd[ii].rz = pSQ->Psd[ii].RadialGrid[ir-1];

		// Compute exact denominator for non-local and store interpolation information
		pSQ->Psd[ii].rc = max(max(pSQ->Psd[ii].rc_s,pSQ->Psd[ii].rc_p),max(pSQ->Psd[ii].rc_d,pSQ->Psd[ii].rc_f));      
		pSQ->Psd[ii].rc = ceil(pSQ->Psd[ii].rc/pSQ->delta)*pSQ->delta;
		pSQ->Psd[ii].rz = (ceil(pSQ->Psd[ii].rz/pSQ->delta))*pSQ->delta;

		if(rank==0)
			printf("max rc is: %f and rz=%f \n",pSQ->Psd[ii].rc,pSQ->Psd[ii].rz);

		double dr=1e-4;
		int nr=1+ceil((pSQ->Psd[ii].rc)/dr);
		double rr;
		double Urr;
		double Vrr;
		double Vlrr;
		double *DU,*DV,*DVloc;//,**DVlJ,**DUlJ;

		pSQ->Psd[ii].DVloc = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DVsJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DUsJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DVpJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DUpJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DVdJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DUdJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DVfJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		pSQ->Psd[ii].DUfJ = new double [pSQ->Psd[ii].size](); // need to de-allocate later

		DU = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		DV = new double [pSQ->Psd[ii].size](); // need to de-allocate later
		DVloc = new double [pSQ->Psd[ii].size](); // need to de-allocate later

		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc,pSQ->Psd[ii].DVloc,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc,DVloc,pSQ->Psd[ii].size); //derivatives for spline

		// Denom_s
		pSQ->Psd[ii].Denom_s=0.0;      
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUs,pSQ->Psd[ii].DUsJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVs,pSQ->Psd[ii].DVsJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUs,DU,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVs,DV,pSQ->Psd[ii].size); //derivatives for spline
		for(ir=0;ir<nr;ir++)
		{
			rr=ir*dr;
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUs, pSQ->Psd[ii].size, &rr,&Urr,1,DU);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVs, pSQ->Psd[ii].size, &rr,&Vrr,1,DV);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc, pSQ->Psd[ii].size, &rr,&Vlrr,1,DVloc);
			if(rr!=0)
				pSQ->Psd[ii].Denom_s += Urr*(Vrr-Vlrr)*Urr*dr/(rr*pSQ->delta*pSQ->delta*pSQ->delta);
		}

		// Denom_p
		pSQ->Psd[ii].Denom_p=0.0;      
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUp,pSQ->Psd[ii].DUpJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVp,pSQ->Psd[ii].DVpJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUp,DU,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVp,DV,pSQ->Psd[ii].size); //derivatives for spline
		for(ir=0;ir<nr;ir++)
		{
			rr=ir*dr;
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUp, pSQ->Psd[ii].size, &rr,&Urr,1,DU);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVp, pSQ->Psd[ii].size, &rr,&Vrr,1,DV);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc, pSQ->Psd[ii].size, &rr,&Vlrr,1,DVloc);
			if(rr!=0)
				pSQ->Psd[ii].Denom_p += Urr*(Vrr-Vlrr)*Urr*dr/(rr*pSQ->delta*pSQ->delta*pSQ->delta);
		}

		// Denom_d
		pSQ->Psd[ii].Denom_d=0.0;   
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUd,pSQ->Psd[ii].DUdJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVd,pSQ->Psd[ii].DVdJ,pSQ->Psd[ii].size); //derivatives for spline   
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUd,DU,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVd,DV,pSQ->Psd[ii].size); //derivatives for spline
		for(ir=0;ir<nr;ir++)
		{
			rr=ir*dr;
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUd, pSQ->Psd[ii].size, &rr,&Urr,1,DU);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVd, pSQ->Psd[ii].size, &rr,&Vrr,1,DV);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc, pSQ->Psd[ii].size, &rr,&Vlrr,1,DVloc);
			if(rr!=0)
				pSQ->Psd[ii].Denom_d += Urr*(Vrr-Vlrr)*Urr*dr/(rr*pSQ->delta*pSQ->delta*pSQ->delta);
		}

		// Denom_f
		pSQ->Psd[ii].Denom_f=0.0;   
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUf,pSQ->Psd[ii].DUfJ,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVf,pSQ->Psd[ii].DVfJ,pSQ->Psd[ii].size); //derivatives for spline   
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUf,DU,pSQ->Psd[ii].size); //derivatives for spline
		getYD_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVf,DV,pSQ->Psd[ii].size); //derivatives for spline
		for(ir=0;ir<nr;ir++)
		{
			rr=ir*dr;
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rUf, pSQ->Psd[ii].size, &rr,&Urr,1,DU);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVf, pSQ->Psd[ii].size, &rr,&Vrr,1,DV);
			ispline_gen(pSQ->Psd[ii].RadialGrid, pSQ->Psd[ii].rVloc, pSQ->Psd[ii].size, &rr,&Vlrr,1,DVloc);
			if(rr!=0)
				pSQ->Psd[ii].Denom_f += Urr*(Vrr-Vlrr)*Urr*dr/(rr*pSQ->delta*pSQ->delta*pSQ->delta);
		}

		if(rank==0)
			cout << "done! "<<endl;  

		delete [] DU;
		delete [] DV;
		delete [] DVloc;
	} // end of for
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// function to print atoms and forces
void PrintAtoms(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{

	FILE *aout_file=NULL;
	char aout_filename[100]="./";
	const char * cchar = pSQ->input_file.c_str(); // convert string to const char
	strcat(aout_filename,cchar);
	strcat(aout_filename,".aout"); 

	aout_file = fopen(aout_filename,"w");
	fprintf(aout_file,"%u\n",pSQ->n_atm); // write number of atoms (first line) 
	fprintf(aout_file,"n_typ %u\n",pSQ->n_typ);
	fprintf(aout_file,"frac_coord %u\n",pSQ->frac_coord);


	double latconst=1.0;
	if(pSQ->frac_coord==1)
	{latconst=pSQ->latconst;
		fprintf(aout_file,"Note: Atom positions below are in fractional coordinates. \n");
	}

	fprintf(aout_file,"\n%s\n","Atoms");

	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		fprintf(aout_file,"%s ",pSQ->Atm[ii].symb);
		fprintf(aout_file,"%u ",pSQ->Atm[ii].Z);
		fprintf(aout_file,"%lf ",pSQ->Atm[ii].mass);
		fprintf(aout_file,"%u ",pSQ->Atm[ii].natm);
		fprintf(aout_file,"%u ",pSQ->Atm[ii].lloc);
		fprintf(aout_file,"%lf\n",pSQ->Atm[ii].rb);
		fprintf(aout_file,"%s\n",pSQ->Atm[ii].pseudo_path);
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			fprintf(aout_file,"%.16f %.16f %.16f\n",pSQ->Atm[ii].Rx[jj].main_atm/latconst,pSQ->Atm[ii].Ry[jj].main_atm/latconst,pSQ->Atm[ii].Rz[jj].main_atm/latconst);
		}

	}

	fprintf(aout_file," \n");
	fprintf(aout_file,"Atomic Forces (Ha/Bohr): \n");
	fprintf(aout_file," \n");

	int k=0;
	for(int JJ_typ=0;JJ_typ<pSQ->n_typ;JJ_typ++) 
	{
		fprintf(aout_file,"%s\n",pSQ->Atm[JJ_typ].symb);
		for(int JJ=0;JJ<pSQ->Atm[JJ_typ].natm;JJ++) 
		{
			fprintf(aout_file,"%.14f %.14f %.14f\n",pSQ->forces[k],pSQ->forces[k+pSQ->n_atm],pSQ->forces[k+2*pSQ->n_atm]);
			k=k+1;
		}
	}

	fclose(aout_file);
	cout <<"Atoms and forces info written to .aout file!" << "\n";
	cout << " "<<endl;
}


// function to print MD velocities
void PrintVels(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{

	FILE *aout_file=NULL;
	char aout_filename[100]="./";
	const char * cchar = pSQ->input_file.c_str(); // convert string to const char
	strcat(aout_filename,cchar);
	strcat(aout_filename,".vout"); 

	if(pSQ->MDstepCount==1)
	{
		aout_file = fopen(aout_filename,"w");
	}else
	{
		aout_file = fopen(aout_filename,"a");
	}
	fprintf(aout_file,"%d\n",pSQ->MDstepCount); 
	double vel=0.0;
	int k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			vel=sqrt(pSQ->vels[k]*pSQ->vels[k]+pSQ->vels[k+pSQ->n_atm]*pSQ->vels[k+pSQ->n_atm]+pSQ->vels[k+2*pSQ->n_atm]*pSQ->vels[k+2*pSQ->n_atm]);
			fprintf(aout_file,"%.6f\n",vel);
			k+=1;
		}

	}

	fclose(aout_file);
	cout <<"Velocities info written to .vout file!" << "\n";
	cout << " "<<endl;
}


// function to print electron density in an SCF
void PrintSCF(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{  
	int rank,i,j,k,count;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  

	// Store local rho into an array on the main domain
	double *Rho;
	int Np=pSQ->n_int[0]*pSQ->n_int[1]*pSQ->n_int[2]; //=pSQ->nproc*(pSQ->np_x*pSQ->np_y*pSQ->np_z);
	Rho = new double [Np]();
	for(int P=0;P<pSQ->nproc;P++) // rank of processor
	{
		if(rank==P)
		{
			count=0;
			for(k=0;k<pSQ->np_z;k++)
			{
				for(j=0;j<pSQ->np_y;j++)
				{
					for(i=0;i<pSQ->np_x;i++)
					{
						Rho[P*(pSQ->np_x*pSQ->np_y*pSQ->np_z)+count]=pSQ->rho[k][j][i];
						count+=1;
					}
				}
			}  

		}
	}

	// All reduce array on main domain, from each processor
	MPI_Allreduce(MPI_IN_PLACE, Rho, Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if(rank==0)
	{
		FILE *aout_file=NULL;
		char aout_filename[100]="./";
		const char * cchar = pSQ->input_file.c_str(); // convert string to const char
		strcat(aout_filename,cchar);
		strcat(aout_filename,".restart"); 
		aout_file = fopen(aout_filename,"w");

		fprintf(aout_file,"%d\n",pSQ->nproc); // first line is no. of procs (since this dictates the order in which we store rho(x)
		for(i=0;i<Np;i++)
		{
			fprintf(aout_file,"%.14f\n",Rho[i]);
		}

		fclose(aout_file);
		cout <<"Electron density written to .restart file!" << "\n";

	}

	delete [] Rho;

}


// function to read electron density for an SCF
void RestartSCF(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{  
	int rank,i,j,k,count;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  

	double *Rho;
	int Np=pSQ->n_int[0]*pSQ->n_int[1]*pSQ->n_int[2]; //=pSQ->nproc*(pSQ->np_x*pSQ->np_y*pSQ->np_z);
	Rho = new double [Np]();

	if(rank==0)
	{
		FILE *aout_file=NULL;
		char aout_filename[100]="./";
		const char * cchar = pSQ->input_file.c_str(); // convert string to const char
		strcat(aout_filename,cchar);
		strcat(aout_filename,".restart"); 
		aout_file = fopen(aout_filename,"r");

		int np;
		fscanf(aout_file,"%d\n",&np);
		if(np!=pSQ->nproc){printf("Error: The restart file was written by different number of processors. Choose the same number as the first line on .restart file. \n\n");exit(0);}
		for(i=0;i<Np;i++)
		{          
			fscanf(aout_file,"%lf\n",&Rho[i]);
		}

		fclose(aout_file);
		cout <<"Electron density read from .restart file!" << "\n";

	}

	// All reduce array on main domain, from each processor
	MPI_Allreduce(MPI_IN_PLACE, Rho, Np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for(int P=0;P<pSQ->nproc;P++) // rank of processor
	{
		if(rank==P)
		{
			count=0;
			for(k=0;k<pSQ->np_z;k++)
			{
				for(j=0;j<pSQ->np_y;j++)
				{
					for(i=0;i<pSQ->np_x;i++)
					{
						pSQ->rho[k][j][i]=Rho[P*(pSQ->np_x*pSQ->np_y*pSQ->np_z)+count];
						count+=1;
					}
				}
			}  

		}
	}

	delete [] Rho;
}

// function to print MD into a restart file ".restartMD"
void PrintMD(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{  
	int rank,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  

	FILE *aout_file=NULL;
	char aout_filename[100]="./";
	const char * cchar = pSQ->input_file.c_str(); // convert string to const char
	strcat(aout_filename,cchar);
	strcat(aout_filename,".restartMD"); 

	aout_file = fopen(aout_filename,"w");

	fprintf(aout_file,"MDstepCount\n"); 
	fprintf(aout_file,"%d\n",pSQ->MDstepCount+pSQ->MDrestartCount); 

	fprintf(aout_file,"MD velocities\n"); 
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			fprintf(aout_file,"%.15f %.15f %.15f\n",pSQ->vels[k],pSQ->vels[k+pSQ->n_atm],pSQ->vels[k+2*pSQ->n_atm]);
			k+=1;
		}

	}

	fprintf(aout_file,"MD forces\n"); 
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			fprintf(aout_file,"%.15f %.15f %.15f\n",pSQ->forces[k],pSQ->forces[k+pSQ->n_atm],pSQ->forces[k+2*pSQ->n_atm]);
			k+=1;
		}

	}

	fprintf(aout_file,"Atom positions\n");
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			fprintf(aout_file,"%.15f %.15f %.15f\n",pSQ->Atm[ii].Rx[jj].main_atm,pSQ->Atm[ii].Ry[jj].main_atm,pSQ->Atm[ii].Rz[jj].main_atm);
		}

	}

	fprintf(aout_file,"MD stats\n");
	fprintf(aout_file,"%.15f %.15f\n",pSQ->mean_TE,pSQ->std_TE);
	fprintf(aout_file,"%.15f %.15f\n",pSQ->mean_PE,pSQ->std_PE);
	fprintf(aout_file,"%.15f %.15f\n",pSQ->mean_KE,pSQ->std_KE);
	fprintf(aout_file,"%.15f %.15f\n",pSQ->mean_T,pSQ->std_T);

	fclose(aout_file);
	cout <<"MD restart info written to .restartMD file!" << "\n";

}

// function to read MD restart file ".restartMD"
void RestartMD(DS_SQ* pSQ)    // pSQ is the pointer to the SQ object of the DS_SQ data structure
{  
	int rank,k;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  

	double *atom_pos;
	int Np=3*pSQ->n_atm;
	atom_pos = new double [Np]();

	if(rank==0)
	{
		FILE *aout_file=NULL;
		char aout_filename[100]="./";
		const char * cchar = pSQ->input_file.c_str(); // convert string to const char
		strcat(aout_filename,cchar);
		strcat(aout_filename,".restartMD"); 
		aout_file = fopen(aout_filename,"r");
		char str[100];
		do
		{
			fgets(str,60,aout_file);      
		}while(strcmp(str,"MDstepCount\n"));
		fscanf(aout_file,"%d",&pSQ->MDrestartCount);

		do
		{
			fgets(str,60,aout_file);            
		}while(strcmp(str,"MD velocities\n"));
		k=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
			{
				fscanf(aout_file,"%lf",&pSQ->vels[k]);
				fscanf(aout_file,"%lf",&pSQ->vels[k+pSQ->n_atm]);
				fscanf(aout_file,"%lf\n",&pSQ->vels[k+2*pSQ->n_atm]);
				k+=1;
			}
		}

		do
		{
			fgets(str,60,aout_file);            
		}while(strcmp(str,"MD forces\n"));
		k=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
			{
				fscanf(aout_file,"%lf",&pSQ->forces[k]);
				fscanf(aout_file,"%lf",&pSQ->forces[k+pSQ->n_atm]);
				fscanf(aout_file,"%lf\n",&pSQ->forces[k+2*pSQ->n_atm]);
				k+=1;
			}
		}

		do
		{
			fgets(str,60,aout_file);            
		}while(strcmp(str,"Atom positions\n"));
		k=0;
		for(int ii=0;ii<pSQ->n_typ;ii++)
		{
			for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
			{
				fscanf(aout_file,"%lf",&pSQ->Atm[ii].Rx[jj].main_atm);
				fscanf(aout_file,"%lf",&pSQ->Atm[ii].Ry[jj].main_atm);
				fscanf(aout_file,"%lf\n",&pSQ->Atm[ii].Rz[jj].main_atm);
				atom_pos[k]=pSQ->Atm[ii].Rx[jj].main_atm;
				atom_pos[k+pSQ->n_atm]=pSQ->Atm[ii].Ry[jj].main_atm;
				atom_pos[k+2*pSQ->n_atm]=pSQ->Atm[ii].Rz[jj].main_atm;
				k+=1;
			}

		}

		do
		{
			fgets(str,60,aout_file);            
		}while(strcmp(str,"MD stats\n"));
		fscanf(aout_file,"%lf",&pSQ->mean_TE);
		fscanf(aout_file,"%lf\n",&pSQ->std_TE);
		fscanf(aout_file,"%lf",&pSQ->mean_PE);
		fscanf(aout_file,"%lf\n",&pSQ->std_PE);
		fscanf(aout_file,"%lf",&pSQ->mean_KE);
		fscanf(aout_file,"%lf\n",&pSQ->std_KE);
		fscanf(aout_file,"%lf",&pSQ->mean_T);
		fscanf(aout_file,"%lf\n",&pSQ->std_T);

		fclose(aout_file);
		cout <<"MD info read from .restartMD file!" << "\n";
		//cout << " "<<endl;
	}


	MPI_Bcast(&pSQ->MDrestartCount,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(pSQ->vels,3*pSQ->n_atm,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(pSQ->forces,3*pSQ->n_atm,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(atom_pos,3*pSQ->n_atm,MPI_DOUBLE,0,MPI_COMM_WORLD);
	k=0;
	for(int ii=0;ii<pSQ->n_typ;ii++)
	{
		for(int jj=0;jj<pSQ->Atm[ii].natm;jj++)
		{
			pSQ->Atm[ii].Rx[jj].main_atm=atom_pos[k];
			pSQ->Atm[ii].Ry[jj].main_atm=atom_pos[k+pSQ->n_atm];
			pSQ->Atm[ii].Rz[jj].main_atm=atom_pos[k+2*pSQ->n_atm];
			k+=1;
		}
	}
	double bcast_double[8];
	bcast_double[0]=pSQ->mean_TE;
	bcast_double[1]=pSQ->std_TE;
	bcast_double[2]=pSQ->mean_PE;
	bcast_double[3]=pSQ->std_PE;
	bcast_double[4]=pSQ->mean_KE;
	bcast_double[5]=pSQ->std_KE;
	bcast_double[6]=pSQ->mean_T;
	bcast_double[7]=pSQ->std_T;
	MPI_Bcast(bcast_double,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
	pSQ->mean_TE=bcast_double[0];
	pSQ->std_TE=bcast_double[1];
	pSQ->mean_PE=bcast_double[2];
	pSQ->std_PE=bcast_double[3];
	pSQ->mean_KE=bcast_double[4];
	pSQ->std_KE=bcast_double[5];
	pSQ->mean_T=bcast_double[6];
	pSQ->std_T=bcast_double[7];

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

	delete [] atom_pos;

}



