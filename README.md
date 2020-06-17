# SQDFT
  SQDFT is a code for performing high-temperature Born-Oppenheimer Quantum Molecular Dynamics 
  (QMD) in the framework of Kohn-Sham Density Functional Theory (DFT).

1. External Libraries:
    SQDFT uses the following external library:  
    MVAPICH2 2.1, <http://mvapich.cse.ohio-state.edu/news/>  
    The code has been tested using mvapich2/2.1 and intel/15.0 compilers, openmpi/1.8 and 
	gcc/4.9.0, and IBM mpicxx-4.7.2. 

2. Compilation:
    SQDFT can be compiled from within the "SQDFT" folder using the commands: 

     `make clean`  
     `make`

    A successful compilation will lead to the creation of the executable "SQDFT/lib/sqdft"

3. Input files:
    SQDFT requires the following files as input:

    - ".input"  - User options and parameters.
    - ".atoms"  - Atomic information.
      In addition, SQDFT requires Troullier-Martins pseudopotential files as generated by the 
      code "atom" (http://bohr.inesc-mn.pt/~jlm/pseudo.html). 

4. Execution:
    SQDFT can be executed in parallel using the `mpirun` command. Sample PBS script files are
    available in "SQDFT/tests" folder. The ".input" and ".atoms" input files should be present 
	in the root folder, and the corresponding pseudopotential files in the "./pseudopotential" 
	subfolder. Example input files can be found in the "SQDFT/tests" folder.
        
    Example: To run a simulation on 8 processors with input files as "filename1.input" and 
	"filename2.atoms", use the following command:
  
    `mpirun -np 8 ./lib/sqdft -name filename1`

5. User Guide:
    Please refer to the SQDFT User Guide ("SQDFTUserGuide.pdf" provided in "SQDFT/doc" subfolder) for
    more information on input and output files contents.
    
6. Tests: 
    The "SQDFT/tests" folder contains two tests:
    
    - Al 
    - LiH  
      These tests can be executed using the mpirun command as described above. The alternative is to
      use the PBS script files that have been provided. The energies and forces in the resulting 
	  output files should be compared with the reference output files: results_Al and results_LiH 
	  given in Al and LiH tests folders, respectively. Note that all simulations should be executed
  from the root folder. 
    
7. Brief Description of all files in SQDFT:
   
   - In "SQDFT/doc" folder
   
     **SQDFTUserGuide.pdf**      : The user guide for SQDFT
   
   - In "SQDFT/inc" folder
   
     ds_sq.h                   : Header file containing information of all data structures
     func_sq.h                : Header file declaring all the functions used in this code
     headers.h               : Header file declaring all the standard header files required for this code

   - In "SQDFT/lib" folder

     **sqdft**                   : Executable created after the successful compilation of the code

   - In "SQDFT/pseudopotentials" folder

     Subfolder containing pseudopotential files for Aluminium, Hydrogen and Lithium
   
   - In "SQDFT/src" folder  
   
     **anderson.cpp**            : Functions for performing Anderson extrapolation   
     **deallocate.cpp**           : Function for de-allocating the memory    
     **energy.cpp**                 : Functions for computing energy    
     **forces.cpp**                  : Functions for computing the atomic forces   
     **initialize.cpp**              : Functions for initializion    
     **main.cpp**                    : The main function    
     **md.cpp**                       : Functions required for Quantum Molecular Dynamics (QMD)    
     **nonlocal.cpp**              : Functions for calculation related to the nonlocal projectors    
     **poisson.cpp**               : Functions for solving the Poisson's equation    
     **readfiles.cpp**             : Functions for reading and storing the inputs from the given input files  
     **scf.cpp**                        : Functions required in the SCF iteration   
     **spline.cpp**                  : Functions for doing the 1d cubic spline interpolation   
     **sq.cpp**                        : Functions for performing Spectral Quadrature (SQ)    
   
   - In "SQDFT/tests" folder

     **Al**                                 : Input and output files for Al  
     **LiH**                              : input and output files for LiH  
     **sqdft.pbs**                   : PBS script file    
   
   - In "SQDFT" folder

     **README.md**              : This file  
     **makefile**                     : The makefile to compile SQDFT
