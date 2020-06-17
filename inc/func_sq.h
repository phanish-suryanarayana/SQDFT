/** \file func_sq.h
\brief Header file declaring all the functions used in this code.
 */

#ifndef FUNC_SQ_H // if FUNC_SQ_H already created then need not compile this header again
#define FUNC_SQ_H  // define FUNC_SQ_H if compiling this header for the first time
#include "headers.h"
#include "ds_sq.h"


//////////////////////////////////////////////////////////////////
/** \brief Function to read the input number of processors and check if it is appropriate for this code.

\param pSQ pointer to the DS_SQ global data structure
\param argv input argument from \ref main function that has input file name information
 */
void CheckInputs(DS_SQ* pSQ,char ** argv); // check input proc info  

/** \brief Function to read the input file and store the input information for the problem

\param pSQ pointer to the DS_SQ global data structure
 */
void Read_input(DS_SQ* pSQ); // function to read the input file

/** \brief Function to compute \f$ (n!n!)/((n-k)!(n+k)!)\f$, where \f$n,k\f$ are integers.

\param n first argument
\param k second argument
 */
double fract(int n,int k);   // function required to compute the factorial terms of FD coeffs

/** \brief Function to read and store the atom data input file

\param pSQ pointer to the DS_SQ global data structure
 */
void Read_atoms(DS_SQ* pSQ); // function to read the atoms file

/** \brief Function to read and store the pseudopotential files for different atom types

\param pSQ pointer to the DS_SQ global data structure
 */
void Read_psd(DS_SQ* pSQ);   // function to read pseudopotential files

/** \brief Function to compute the positions of the periodic replica atoms within a distance of \f$r_b\f$ from the boundary of processor domain

\param pSQ pointer to the DS_SQ global data structure
 */
void Replica_atoms(DS_SQ* pSQ);

/** \brief Function to compute the positions of the periodic replica atoms within a distance of \f$R_{cut}+r_c\f$ from the boundary of processor domain. This is required for non-local calculations.

\param pSQ pointer to the DS_SQ global data structure
 */
void Replica_atoms_Rcut(DS_SQ* pSQ); 

/** \brief Function to distribute the main domain among the processors and compute the start/end nodes of the processor domains.

\param pSQ pointer to the DS_SQ global data structure
 */
void Processor_domain(DS_SQ* pSQ);

/** \brief Function to compute the communication topologies that will be used by MPI collectives

\param pSQ pointer to the DS_SQ global data structure
 */
void Comm_topologies(DS_SQ* pSQ);

//////////////////////////////////////////////////////////////////
// spline interpolation

/** \brief Function that does one dimensional cubic spline interpolation

\param X1 reference radial grid
\param Y1 reference function value
\param len1 size of above reference one dimensional arrays
\param X2 radial points at which interpolation needs to be done
\param Y2 interpolated function value at the required radial points
\param len2 size of one dimensional arrays at which interpolation needs to be done
\param YD derivative information computed using getYD_gen function
 */
void ispline_gen(double *X1,double *Y1,int len1,double *X2,double *Y2,int len2,double *YD);
/** \brief Function that computes derivative information required for cubic spline interpolation by ispline_gen

\param X reference radial grid
\param Y reference function value
\param YD derivative information
\param len size of above reference one dimensional arrays
 */
void getYD_gen(double *X, double *Y, double *YD,int len);
/** \brief Function that is needed by getYD_gen to solve a tridiagonal matrix

\param A,B,C bands of tridiagonal matrix
\param D derivative information
\param len size of reference arrays used for interpolation
 */
void tridiag_gen(double *A,double *B,double *C,double *D,int len);

// auxiliary functions
/** \brief Function to compute 2-norm of a vector (across all the main domain)

\param Vec vector (one-dimensional array)
\param len size of the vector
\param ResVal scalar output, 2-norm of the vector
*/
void Vector2Norm(double* Vec, int len, double* ResVal); // function to compute the vector 2-norm

/** \brief Function to print electron density in SCF into a .restart file..

*/
void PrintSCF(DS_SQ* pSQ); 

/** \brief Function to read electron density in SCF from a .restart file..

*/
void RestartSCF(DS_SQ* pSQ); 

//////////////////////////////////////////////////////////////////
/** \brief Function to initialize the velocities and accelerations for Molecular Dynamics (MD).

*/
void Initialize_MD(DS_SQ* pSQ); // MD initialization

/** \brief Function to wrap around the atom positions across the main domain (periodic boundary condtions).

*/
void WrapAround_PBC(DS_SQ* pSQ); 

/** \brief Function to evaluate nuclear (ionic) charge density, guess electron density in the processor domain and energy correction.

This function also computes total nuclear charge, the self-energy and the potential associated with energy correction.
\param pSQ pointer to the DS_SQ global data structure
*/
void ChargeDensity(DS_SQ* pSQ); // b-calculation

/** \brief Function to compute correction in energy due to repulsion between overlapping nuclei.

*/
void OverlapCorrection(DS_SQ* pSQ); // energy corrections
/** \brief Function to compute correction in forces due to repulsion between overlapping nuclei.

*/
void OverlapCorrection_forces(DS_SQ* pSQ); // force corrections


// /** \brief Function to compute initial guess electron density in main domain for the first SCF iteration
//
//This function also computes total number of electrons.
// \param pSQ pointer to the DS_SQ global data structure
//*/
//void GuessElectronDensity(DS_SQ* pSQ); // rho-calculation

/** \brief Function to allocate memory for various arrays (b, rho, phi, Veff etc.)

\param pSQ pointer to the DS_SQ global data structure
*/
void AllocateArrays(DS_SQ* pSQ); 

/** \brief Function to perform NVE (microcanonical) Molecular Dynamics (MD).

*/
void NVE_MD(DS_SQ* pSQ); // MD

/** \brief Function to perform NVT (Canonical) Molecular Dynamics (MD) using Nose-Hoover thermostat.

*/

void NVT_NH_mov8(DS_SQ* pSQ);
void vv_1(DS_SQ* pSQ);
void NR(DS_SQ* pSQ);

/** \brief Function to update the atomic positions using Leapfrog algorithm in Molecular Dynamics (MD).

*/
void Leapfrog_part1(DS_SQ* pSQ); // MD

/** \brief Function to update the atomic velocities using Leapfrog algorithm in Molecular Dynamics (MD).

*/
void Leapfrog_part2(DS_SQ* pSQ); // MD

/** \brief Function to compute Total Energy (TE=PE+KE), Potential Energy (PE) and Kinetic Energy (KE) in Molecular Dynamics (MD).

*/
void MD_energies(DS_SQ* pSQ); // MD

/** \brief Function to perform charge extrapolation to provide better rho_guess for MD steps.

*/
void ChargeExtrapolation_MD(DS_SQ* pSQ); // MD


/** \brief Function to perform SCF and compute DFT forces in every MD step.

*/
void SQDFT_forces(DS_SQ* pSQ); // MD

/** \brief Function to print latest MD information into a .restartMD file..

*/
void PrintMD(DS_SQ* pSQ); 

/** \brief Function to read MD information from a .restartMD file..

*/
void RestartMD(DS_SQ* pSQ); 

/** \brief Function to perform atomic relaxation (or) energy minimization using Non-linear conjugate gradient method.

*/
void NLCG(DS_SQ* pSQ); // MD

/** \function to store atomic positions...

*/
void datacollect(DS_SQ* pSQ);


//////////////////////////////////////////////////////////////////
/** \brief Function to compute compute communication information for Laplacian by calling the \sa EdgeIndicesForPoisson() function.

\param pSQ pointer to the DS_SQ global data structure
*/
void Laplacian_Comm_Indices(DS_SQ* pSQ); 

//void EdgeIndicesForPoisson(DS_SQ* pSQ, int eout_s[3][6],int eout_e[3][6], int ein_s[3][6],int ein_e[3][6], int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind); // compute edge indices for residual calculation in Poisson solver
/** \brief Function to compute start/end indices of communication regions for residual calculation in Poisson solver

\param eout_s Two dimensional array of start node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t processor domain. This array corresponds to the outgoing array in MPI communication from the current processor.
\param eout_e Two dimensional array of end node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t processor domain. This array corresponds to the outgoing array in MPI communication from the current processor.
\param ein_s Two dimensional array of start node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t a domain which is larger and is finite difference stencil width away from the processor domain in all directions. This array corresponds to the incoming array in MPI communication to the current processor.
\param ein_e Two dimensional array of end node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t a domain which is larger and is finite difference stencil width away from the processor domain in all directions. This array corresponds to the incoming array in MPI communication to the current processor.
\param ereg_s, ereg_e are two dimensional arrays with start/end node indices of communication regions when number of nodes in the processor domain is larger than finite difference order. The regions comprise of 6 faces, 12 edges and 8 corners. This is useful in order to make the partial stencil updates after the Laplcian is partially applied for all finite difference nodes in the processor domain. These arrays include all the indices of nodes that will need partial stencil update.
\param stencil_sign is a two dimensional array that indicates whether the partial stencil update is done on the right (+1) or left (-1).
\param edge_ind is a two dimensional array of indices which indicate the starting point of the partial stencil update in each communication region.
\param displs_send,displs_recv are arrays of relative displacement indices of the first element of each outgoing/incoming array from each communication region defined w.r.t the first communication region's array. This is needed as an input for the vector form of the MPI communication collective.
\param ncounts_send,ncounts_recv are arrays of number of nodes to be communicated with each neighboring processor. This is needed as an input for the vector form of the MPI communication collective.
*/
void EdgeIndicesForPoisson(DS_SQ* pSQ, int **eout_s,int **eout_e, int **ein_s,int **ein_e, int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv); // compute edge indices for residual calculation in Poisson solver 

//void PoissonResidual(DS_SQ* pSQ,double ***phi_old,double ***phi_res,int iter,MPI_Comm comm_dist_graph_cart); // compute residual in Poisson solver
//void PoissonResidual(DS_SQ* pSQ,double ***phi_old,double ***phi_res,int iter,MPI_Comm comm_dist_graph_cart, int eout_s[3][6],int eout_e[3][6], int ein_s[3][6],int ein_e[3][6], int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind); // compute residual in Poisson solver
/** \brief Function to compute the residual at each iteration in the Poisson solver

\param pSQ pointer to the global data structure DS_SQ
\param phi_old is a three dimensional array of the current iterate in the Poisson solver. The array is defined on a domain which is larger than the processor domain by finite difference stencil width. The size of the array is (np_x+2FDn) x (np_y+2FDn) x (np_z+2FDn). 
\param phi_res is a three dimensional array of the residual of the current iteration of the Poisson solver. The array is defined on a domain which is larger than the processor domain by finite difference stencil width. The size of the array is (np_x+2FDn) x (np_y+2FDn) x (np_z+2FDn). 
\param iter is the current iteration number
\param comm_dist_graph_cart MPI graph topology of processors for communication involved in applying finite difference Laplacian 
\param eout_s, eout_e, ein_s, ein_e, ereg_s, ereg_e, stencil_sign, edge_ind, displs_send, displs_recv, ncounts_send, ncounts_recv are entities related to the communication regions. \sa EdgeIndicesForPoisson()
*/
void PoissonResidual(DS_SQ* pSQ,double ***phi_old,double ***phi_res,int iter,MPI_Comm comm_dist_graph_cart, int **eout_s,int **eout_e, int **ein_s,int **ein_e, int ereg_s[3][26],int ereg_e[3][26],int **stencil_sign,int **edge_ind,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv); // compute residual in Poisson solver 

/** \brief Function to compute average time taken by Poisson solver across all processors.

\param dt is the time taken by current processor.
\param nproc is the total number of processors in the domain.
*/
void AvgTime_Poisson(double dt, int nproc);


/** \brief Function to compute the Pythagoras sum (\f$ \sqrt{(a^2 + b^2)} \f$) of two scalars.

\param a,b are the two scalars. 
*/
double pythag(double a, double b); // computes (a^2 + b^2)^0.5

/** \brief Function to compute Singular Value Decomposition of an input matrix. This is needed for Pseudoinverse calculation. \f$A=UWV^T\f$.  

\param a is 2D array representing the input matrix of size m x n. \f$U\f$ replaces matrix \f$a\f$ on output.
\param m,n are the number of rows and columns of the matrix respectively.
\param w is an array of n elements representing the singular values.
\param v is a 2D array representing the output matrix \f$V\f$ of size n x n.
*/  
void SingularValueDecomp(double **a,int m,int n, double *w, double **v); // a is matrix (array) size m x n, A=UWV'. U replaces "a" on output. w is an array of singular values, size 1 x n. V is output as matrix v of size n x n. 

/** \brief Function to compute the pseudoinverse times vector of a square matrix. This is used in the AndersonExtrapolation function. \sa AndersonExtrapolation().

\param Ac is the matrix \f$A\f$ input in column major format.
\param b is the vector with which the pseudoinverse of matrix A should be multiplied.
\param x is the output vector.
\param m is the number of rows/columns in matrix A. 
*/  
void PseudoInverseTimesVec(double *Ac,double *b,double *x,int m); // returns x=pinv(A)*b, matrix A is m x m  but given in column major format Ac (so m^2 x 1) and vector b is m x 1

//void MPI_vec_sum(double *invec, double *inoutvec, int *len, MPI_Datatype *dtype);
/** \brief Function to compute the Anderson update vector for a fixed point iteration. This is used in both AAJ Poisson solver as well as SCF iteration.

\param DX is a two dimensional array of size anderson_history x (np_x*np_y*np_z). It stores the previous anderson_history+1 iterates of the Poisson solver
\param DF is a two dimensional array of size anderson_history x (np_x*np_y*np_z). It stores the previous anderson_history+1 residuals of the Poisson solver
\param phi_res_vec is a one dimensional array of size (np_x*np_y*np_z) x 1. It is the vector form of the three dimensional array \ref phi_res.
\param beta_mix is the relaxation parameter of Anderson mixing
\param anderson_history is an integer which is one less than the total number of iterates stored in the history arrays
\param N is the size of phi_res_vec array
\param am_vec is the Anderson update vector. It is of the same size of phi_res_vec
\param mpi_time is scalar to store time taken by mpi communication in this function
*/
void AndersonExtrapolation(double **DX, double **DF, double *phi_res_vec, double beta_mix, int anderson_history, int N, double *am_vec, double *mpi_time); // Anderson update
//////////////////////////////////////////////////////////////////
/** \brief Function to solve the Poisson's equation for electrostatic potential \f$ \phi\f$, using AAJ linear solver.

\param pSQ pointer to the global data structure DS_SQ
*/
void PoissonSolver_AAJ(DS_SQ* pSQ); // Alternating Anderson-Jacobi (AAJ) Iteration for Poisson's equation

/** \brief Function to compute the exchange correlation potential Vxc.

\param pSQ pointer to the global data structure DS_SQ
*/
void ExchangeCorrelationPotential(DS_SQ* pSQ); // Compute Vxc

/** \brief Function to compute the effective potential \f$V_{eff}=\phi+V_{xc}+V_{nloc}\f$.

\param pSQ pointer to the global data structure DS_SQ
*/
void EvaluateEffectivePotential(DS_SQ* pSQ); // Compute Vxc

/** \brief Function to compute communication information for SQ. \sa RcutRegionIndicesForSQ().

\param pSQ pointer to the global data structure DS_SQ
*/
void SQ_Comm_Indices(DS_SQ* pSQ); 


/** \brief Function to compute the start/end indices of truncation region (Rcut) required for communication in SQ method.

\param pSQ pointer to the global data structure DS_SQ
\param eout_s,eout_e Two dimensional array of start/end node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t processor domain. This array corresponds to the outgoing array in MPI communication from the current processor.
\param ein_s,ein_e Two dimensional array of start node indices with one of the dimension for the three directions and one for each of the communication region associated with each neighboring processor domain. So it is a 2d array of size 3 x nneigh, where nneigh is the total number of neighboring processor domains involved in communication. The indices stored in this array correspond to the nodes in the processor domain defined w.r.t a domain which is larger and is nloc width away from the processor domain in all directions. This array corresponds to the incoming array in MPI communication to the current processor.
\param displs_send,displs_recv are arrays of relative displacement indices of the first element of each outgoing/incoming array from each communication region defined w.r.t the first communication region's array. This is needed as an input for the vector form of the MPI communication collective.
\param ncounts_send,ncounts_recv are arrays of number of nodes to be communicated with each neighboring processor. This is needed as an input for the vector form of the MPI communication collective.
 */
void RcutRegionIndicesForSQ(DS_SQ* pSQ, int **eout_s,int **eout_e, int **ein_s,int **ein_e,int *displs_send,int *displs_recv,int *ncounts_send,int *ncounts_recv); // compute edge indices for SQ  MPI communication

/** \brief Function to compute sub Non-local times vector (matrix free).

\param pSQ pointer to the global data structure DS_SQ
\param vec is the array with which the matrix elements need to be multiplied.
\param i,j,k are the indices of the the finite difference node corresponding to the sub-Hamiltonian. The indices are w.r.t a domain which extends until Rcut distance from the processor domain.
\param Vv is the output array which is equal to the matrix times the input vector (array).
 */
void VnlocsubTimesVec(DS_SQ* pSQ,double ***vec,int i,int j,int k,double ***Vv); 

/** \brief Function to compute sub Non-local (J^th) times vector (matrix free) for force calculation.

\param pSQ pointer to the global data structure DS_SQ
\param vec is the array with which the matrix elements need to be multiplied.
\param i,j,k are the indices of the the finite difference node corresponding to the sub-Hamiltonian. The indices are w.r.t a domain which extends until Rcut distance from the processor domain.
\param Vv is the output array which is equal to the matrix times the input vector (array).
\param j_atm,j_typ are the parameters required for non-local force calculation which is done for j_atm.
 */
void VnlocsubTimesVec_J(DS_SQ* pSQ,double ***vec,int i,int j,int k,double ***Vv, int j_atm, int j_typ); 


/** \brief Function to compute Spherical Harmonics required in non-local.

\param x x-coordinate of the point at which Spherical Harmonic has to be evaluated.
\param y y-coordinate of the point at which Spherical Harmonic has to be evaluated.
\param z z-coordinate of the point at which Spherical Harmonic has to be evaluated.
\param r radial distance of the point at which Spherical Harmonic has to be evaluated.
\param l,m quantum numbers for which Spherical Harmonic has to be evaluated. l=0,1,2,3 and m = -l to +l.
 */
double SphericalHarmonics(double x,double y, double z, double r, int l, int m); 

/** \brief Function to compute non-local projectors for all atoms associated with the processor domain.

\param pSQ pointer to the global data structure DS_SQ.
*/
void NonlocalProjectors(DS_SQ* pSQ); 

/** \brief Function to compute sub-Hamiltonian times vector (matrix free).

\param pSQ pointer to the global data structure DS_SQ
\param vec is the array with which the matrix elements need to be multiplied.
\param i,j,k are the indices of the the finite difference node corresponding to the sub-Hamiltonian. The indices are w.r.t a domain which extends until Rcut distance from the processor domain.
\param Hv is the output array which is equal to the matrix times the input vector (array).
 */
void HsubTimesVec(DS_SQ* pSQ, double ***vec, int i, int j, int k, double ***Hv); // compute edge indices for SQ  MPI communication

/** \brief Function to compute square root of sum of squares of elements of a three dimensional array.

\param vec is the 3D array
\param nx is the number of elements of the array in each direction/dimension.
\param val is the output scalar value.
 */
void Vector2Norm_local(double ***vec,int nx,double *val);

/** \brief Function to compute sum of products of individual elements of two three dimensional arrays.

\param v1 is the first 3D array
\param v2 is the second 3D array
\param nx is the number of elements of the array in each direction/dimension. This is same for both the arrays.
\param val is the output scalar value.
 */
void VectorDotProduct_local(double ***v1,double ***v2,int nx,double *val);

/** \brief Function to compute eigenvalues of a real,symmetric tridiagonal matrix using QL algorithm with implicit shifts.

\param diag is the array of diagonal elements of the tridiagonal matrix.
\param subdiag is the array of off-diagonal elements of the tridiagonal matrix.
\param n is the number of elements in each row/column of the matrix.
\param lambda_min is the minimum eigenvalue of the tridiagonal matrix.
\param lambda_max is the maximum eigenvalue of the tridiagonal matrix.
 */

void TridiagEigenSolve(double *diag,double *subdiag,int n,double *lambda_min,double *lambda_max);

/** \brief Function to compute extremal eigenvalues of sub-Hamiltonia.

\param pSQ pointer to the global data structure DS_SQ
\param vec is a 3D array with 2*nloc+1 elements in each direction and is the initial guess for Lanczos iteration, where nloc is the the number of finite difference nodes in truncation region.
\param i,j,k are the indices of the the finite difference node corresponding to the sub-Hamiltonian. The indices are w.r.t a domain which extends until Rcut distance from the processor domain.
\param lambda_min is the minimum eigenvalue of the sub-Hamiltonian.
\param lambda_max is the maximum eigenvalue of the sub-Hamiltonian.
\param nd is the current nodes index in the processor. 
*/

void LanczosAlgorithm(DS_SQ* pSQ, double ***vec, int i, int j, int k, double *lambda_min, double *lambda_max,int nd);


/** \brief Function to evaluate Fermi-Dirac function g(t).

\param t is the argument at which to evaluate the Fermi-Dirac function.
\param lambda_f is a parameter in the Fermi-Dirac function.
\param bet is the smearing which is a parameter in the Fermi-Dirac function.
 */
double FermiDirac(double t,double lambda_f,double bet); 

/** \brief Function to evaluate t*g(t) where g(t) is the Fermi-Dirac function.

\param t is the argument at which to evaluate the Fermi-Dirac function.
\param lambda_f is a parameter in the Fermi-Dirac function.
\param bet is the smearing which is a parameter in the Fermi-Dirac function.
 */
double UbandFunc(double t,double lambda_f,double bet);

/** \brief Function to evaluate the entropy function g*log(g)+(1-g)*log(1-g) where g(t) is the Fermi-Dirac function.

\param t is the argument at which to evaluate the Fermi-Dirac function.
\param lambda_f is a parameter in the Fermi-Dirac function.
\param bet is the smearing which is a parameter in the Fermi-Dirac function.
 */
double EentFunc(double t,double lambda_f,double bet);


/*NOTE: (TODO) THIS FUNCTION CURRENTLY USES DISCRETE ORTHOGONALITY PROPERTY OF CHEBYSHEV POLYNOMIALS WHICH ONLY GIVES "approximate" COEFFICIENTS. SO IF THE FUNCTION IS NOT SMOOTH ENOUGH THE COEFFICIENTS WILL NOT BE ACCURATE. SO THIS FUNCTION HAS TO BE REPLACED WITH ITS CONTINUOUS COUNTER PART WHICH EVALUATES OSCILLATORY INTEGRALS AND USES FFT.
*/
/** \brief Function to compute Chebyshev expansion coefficients in Clenshaw-Curtis Spectral Quadrature method. 

\param npl is the degree of the Chebyshev fit
\param lambda_fp is the Fermi energy
\param bet is the inverse of smearing
\param Ci,di are arrays to store the Chebyshev coefficients.
\param func is the function handle to which the Chebyshev fit is computed
 */
void ChebyshevCoeff(int npl,double lambda_fp,double bet,double *Ci,double *di,double (*func)(double,double,double));
//void ChebyshevCoeff(int npl,double lambda_fp,double bet,double *Ci,double (*func)(double,double,double));


/** \brief Function to evaluate the constraint on number of electrons.

\param pSQ pointer to the global data structure DS_SQ.
\param lambda_f is a parameter in the Fermi-Dirac function. It would be the Fermi energy when the constraint is satisfied.
\param Ci,di are arrays to store the Chebyshev coefficients.
\param tferm_mpi is a scalar to store mpi communication time
 */
double NeConstraint(DS_SQ* pSQ,double lambda_f,double *Ci,double *di,double *tferm_mpi); 
//double NeConstraint(DS_SQ* pSQ,double lambda_f); 

/** \brief Function to find the root of NeConstraint() using Brent's algorithm.

\param pSQ pointer to the global data structure DS_SQ.
\param x1,x2 are initial guesses.
\param Ci,di are arrays to store the Chebyshev coefficients.
\param count is the SCF iteration number. Useful to provide better guesses from second iteration onwards.
\param tferm_mpi is a scalar to store mpi communication time
 */
double BrentsAlgorithm(DS_SQ* pSQ,double x1,double x2,double *Ci,double *di,int count,double *tferm_mpi); 




/** \brief Function to evaluate derivative of Fermi-Dirac function g(t) w.r.t lambda_f

\param t is the argument at which to evaluate the Fermi-Dirac function.
\param lambda_f is a parameter in the Fermi-Dirac function.
\param bet is the smearing which is a parameter in the Fermi-Dirac function.
 */
double FermiDirac_derv(double t,double lambda_f,double bet); 

/** \brief Function to evaluate the derivative of and the constraint on number of electrons.

\param pSQ pointer to the global data structure DS_SQ.
\param lambda_f is a parameter in the Fermi-Dirac function. It would be the Fermi energy when the constraint is satisfied.
\param Ci,di are arrays to store the Chebyshev coefficients.
\param fa,dfa are the trace evaluations of density matrix and the derivative terms.
\param tferm_mpi is a scalar to store mpi communication time
 */
double NeConstraint_derv(DS_SQ* pSQ,double lambda_f,double *Ci,double *di,double *fa,double *dfa,double *tferm_mpi); 

/** \brief Function to find the root of NeConstraint() using Newton-Raphson algorithm.

\param pSQ pointer to the global data structure DS_SQ.
\param Ci,di are arrays to store the Chebyshev coefficients.
\param tferm_mpi is a scalar to store mpi communication time
 */
double NewtonRaphson(DS_SQ* pSQ,double *Ci,double *di,double *tferm_mpi); 



/** \brief Function to compute Chebyshev expansion components in Clenshaw-Curtis Spectral Quadrature method.

\param pSQ pointer to the global data structure DS_SQ
\param count is the SCF iteration number
 */
void ClenshawCurtisSpectralQuadrature(DS_SQ* pSQ,int count); // SQ for electron density
  
/** \brief Function to compute the total energy of the system

\param pSQ pointer to the global data structure DS_SQ
*/
void EvaluateTotalEnergy(DS_SQ* pSQ); // compute energy


/** \brief Function to perform SCF iteration.

\param pSQ pointer to the global data structure DS_SQ
*/
void SCF_iteration(DS_SQ *pSQ); 

/** \brief Function to compute the local component of the forces and force corrections.

\param pSQ pointer to the global data structure DS_SQ
*/
void ForcesLocal(DS_SQ* pSQ); 

/** \brief Function to compute the non-local component of the forces.

\param pSQ pointer to the global data structure DS_SQ
*/
void ForcesNonLocal(DS_SQ* pSQ); 

/** \brief Function to compute total forces (including correction).

\param pSQ pointer to the global data structure DS_SQ
*/
void ForcesTotal(DS_SQ* pSQ); 

/** \brief Function to print atoms and forces information into .aout file.

\param pSQ pointer to the global data structure DS_SQ
*/
void PrintAtoms(DS_SQ* pSQ); 

/** \brief Function to print MD velocities information into .vout file.

\param pSQ pointer to the global data structure DS_SQ
*/
void PrintVels(DS_SQ* pSQ); 

//////////////////////////////////////////////////////////////////
/** \brief Function to de-allocate or free the memory that has been dynamically allocated during run-time.

\param pSQ pointer to the global data structure DS_SQ
*/
void Deallocate_memory(DS_SQ* pSQ); // function to de-allocate memory



#endif
