/** \file main.cpp
  \brief The file containing the main() function.

*/

/**
  \mainpage Parallel Spectral Quadrature code for Density Functional Theory (SQDFT)

  This code performs DFT calculations using Spectral Quadrature method. It can do periodic infinite cell calculations on three dimensional cubical domain of atoms.

  \authors Phanisri Pradeep Pratapa, Georgia Institute of Technology
  \authors Abhiraj Sharma, Georgia Institute of Technology
  \authors Phanish Suryanarayana, Georgia Institute of Technology

*/

#include "headers.h" // standard header files
#include "ds_sq.h"   // header file with data structures
#include "func_sq.h" // header file with functions

#ifndef TEMP_TOL
#define TEMP_TOL 1e-12
#endif

/** \brief The main function that calls all the other functions for SQ code.

  \param argc first input argument
  \param argv second input argument
  */
int main(int argc, char ** argv)
{ 
	DS_SQ SQ={};  /// Declare an object called "SQ", details of data structure come from header file DS_SQ.h.

	int rank,psize;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&psize);
	SQ.nproc = psize; // no. of processors
	SQ.poiss_time=0.0; SQ.tpoiss_mpi=0.0;
	SQ.sq_time=0.0; SQ.tsq_mpi=0.0;
	SQ.engy_time=0.0; SQ.tengy_mpi=0.0;
	SQ.forcs_time=0.0; SQ.tforcs_mpi=0.0;

	double t_begin,t_end,t_read,t_barr,t_comm,t_init,t_free;
	double tcomm_inp,tcomm_int;
	t_begin = MPI_Wtime();

	CheckInputs(&SQ,argv);

	if(rank==0)
		printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");	

	/// Read from an input file, store the data in data structure  and write the data to the output file. 
	SQ.mpi_time=0.0;   
	Read_input(&SQ);   
	Read_atoms(&SQ);
	Read_psd(&SQ);
	tcomm_inp = SQ.mpi_time;  
	t_barr = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD); 
	t_read = MPI_Wtime();

	/// Initialize quantities.
	Processor_domain(&SQ);
	AllocateArrays(&SQ);

	SQ.mpi_time=0.0;   
	Comm_topologies(&SQ);  
	tcomm_int = SQ.mpi_time;  
	Laplacian_Comm_Indices(&SQ); // compute communication indices information for Laplacian
	SQ_Comm_Indices(&SQ);        // compute communication indices information for SQ 
	t_comm = MPI_Wtime();

	t_init = MPI_Wtime();
	if(rank==0){
		cout << " "<<endl;
		printf("Note: Times printed in parentheses () beside the full time taken indicate the MPI communications part of the calculation. \n");
		cout << " "<<endl;
		printf("Time spent in initialization = %.4f (%.4f) seconds. \n",t_init-t_begin,tcomm_inp+tcomm_int+(t_read-t_barr));
		printf("Break down: \n");
		printf("Reading input files time     = %.4f (%.4f) seconds. \n",t_barr-t_begin,tcomm_inp);
		printf("MPI_Barrier after reading    = %.4f seconds. \n",t_read-t_barr);
		printf("Init & create comm topos     = %.4f (%.4f) seconds. \n",t_comm-t_read,tcomm_int);
		cout << " "<<endl;}

	/// Initialize MD
	if(SQ.MaxMDsteps > 1)
		Initialize_MD(&SQ);
	double mean_TE=0.0,std_TE=0.0,mean_PE=0.0,std_PE=0.0,mean_KE=0.0,std_KE=0.0,mean_T=0.0,std_T=0.0,mean_TEext=0.0,std_TEext=0.0; // default
	SQ.MDrestartCount=0; // default
	if(SQ.restart_md==1)
	{
		RestartMD(&SQ); // vels, accels, atom positions, stats and MDrestartCount are updated
		mean_TE=SQ.mean_TE;std_TE=SQ.std_TE;
		mean_PE=SQ.mean_PE;std_PE=SQ.std_PE;
		mean_KE=SQ.mean_KE;std_KE=SQ.std_KE;
		mean_T=SQ.mean_T;std_T=SQ.std_T;
		if(SQ.MDrestartCount==SQ.MaxMDsteps)
		{if(rank==0){printf("Error: MaxMDsteps is already reached from the restart file. Nothing to be done.\n");}exit(0);}
	}

	SQ.MDstepCount=1;
	if(SQ.MaxMDsteps==1)
	{
		if(SQ.RelaxAtoms!=1)
			SQDFT_forces(&SQ);

		if(SQ.RelaxAtoms==1)
			NLCG(&SQ); // relaxation
	}

	while((SQ.MDstepCount+SQ.MDrestartCount)<=SQ.MaxMDsteps && SQ.MaxMDsteps > 1)
	{
		if(rank==0)
		{
			printf("MD step %d \n",SQ.MDstepCount+SQ.MDrestartCount);
		}
		t_init = MPI_Wtime();
		if(SQ.ensemble == 1)
			NVT_NH_mov8(&SQ) ; // Temperature Constant
		else if(SQ.ensemble == 2)
			NVE_MD(&SQ); // Energy constant
		mean_TE+=SQ.TE;
		std_TE+=SQ.TE*SQ.TE;
		mean_PE+=SQ.PE;
		std_PE+=SQ.PE*SQ.PE;
		mean_KE+=SQ.KE;
		std_KE+=SQ.KE*SQ.KE;
		mean_T+=SQ.T_MD;
		std_T+=SQ.T_MD*SQ.T_MD;
		mean_TEext+=SQ.TEext;
		std_TEext+=SQ.TEext*SQ.TEext;
		SQ.mean_TE=mean_TE;SQ.std_TE=std_TE;
		SQ.mean_PE=mean_PE;SQ.std_PE=std_PE;
		SQ.mean_KE=mean_KE;SQ.std_KE=std_KE;
		SQ.mean_T=mean_T;SQ.std_T=std_T;
		SQ.mean_TEext=mean_TEext;SQ.std_TEext=std_TEext;

		if(SQ.prnt_md==1 && rank==0)
			PrintMD(&SQ);

		t_end = MPI_Wtime();
		if(rank==0)
		{
			printf("TE (mean, std) (Ha/atom)     : %.14f %.14f  \n",mean_TE/(SQ.MDstepCount+SQ.MDrestartCount),sqrt(std_TE/(SQ.MDstepCount+SQ.MDrestartCount)-(mean_TE/(SQ.MDstepCount+SQ.MDrestartCount))*(mean_TE/(SQ.MDstepCount+SQ.MDrestartCount))));
			printf("PE (mean, std) (Ha/atom)     : %.14f %.14f  \n",mean_PE/(SQ.MDstepCount+SQ.MDrestartCount),sqrt(std_PE/(SQ.MDstepCount+SQ.MDrestartCount)-(mean_PE/(SQ.MDstepCount+SQ.MDrestartCount))*(mean_PE/(SQ.MDstepCount+SQ.MDrestartCount))));
			printf("KE (mean, std) (Ha/atom)     : %.14f %.14f  \n",mean_KE/(SQ.MDstepCount+SQ.MDrestartCount),sqrt(std_KE/(SQ.MDstepCount+SQ.MDrestartCount)-(mean_KE/(SQ.MDstepCount+SQ.MDrestartCount))*(mean_KE/(SQ.MDstepCount+SQ.MDrestartCount))));
			printf("T (mean, std) (K)            : %.14f %.14f  \n",mean_T/(SQ.MDstepCount+SQ.MDrestartCount),sqrt(fabs(std_T/(SQ.MDstepCount+SQ.MDrestartCount)-(mean_T/(SQ.MDstepCount+SQ.MDrestartCount))*(mean_T/(SQ.MDstepCount+SQ.MDrestartCount)))));
			printf("MD step time      = %.4f sec \n\n",t_end-t_init);

		}

		SQ.MDstepCount+=1;
	} // end while loop

	Deallocate_memory(&SQ);  ///< De-allocate memory.
	t_free = MPI_Wtime();

	if(rank==0)
	{
		printf("Total wall time   = %.4f \n",t_free-t_begin);
		printf("============================== \n");

		cout << " "<<endl;
		printf("References for SQDFT: \n");

		printf("   (1) Spectral Quadrature method for accurate O(N) electronic structure calculations of metals and insulators, P.P. Pratapa, P. Suryanarayana, and J.E. Pask, Comput. Phys. Commun. 200, 96--107 (2016). \n");
		printf("   (2) Anderson acceleration of the Jacobi iterative method: an efficient alternative to Krylov methods for large, sparse linear systems, P.P. Pratapa, P. Suryanarayana, J.E. Pask, J. Comput. Phys. 306, 43--54 (2016). \n");
		printf("   (3) Periodic Pulay method for robust and efficient convergence acceleration of self-consistent field iterations, A.S. Banerjee, P. Suryanarayana, J.E. Pask, Chem. Phys. Lett. 647, 31--35 (2016). \n");
		cout << " "<<endl;
		printf("End of run for %d atom system on %d processors. \n",SQ.n_atm,SQ.nproc);
		char* c_time_str;
		time_t current_time=time(NULL);
		c_time_str=ctime(&current_time);  
		printf("Ending time: %s",c_time_str);

	} 
	if(rank==0)
		printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

	MPI_Barrier(MPI_COMM_WORLD);   
	MPI_Finalize();
	return 0;
} 


