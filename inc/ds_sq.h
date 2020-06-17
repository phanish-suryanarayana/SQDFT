/** \file ds_sq.h
\brief Header file containing information of all data structures.

A single global data structure with all the variables is created for the SQ code. Separate data structures for atom data and pseudopotential data are also created which are parts of the global data structure.
*/

#ifndef DS_SQ_H  // if DS_SQ_H already created then need not compile this header again
#define DS_SQ_H  // define DS_SQ_H if compiling this header for the first time
#include "headers.h"

/** \brief Data structure for each atom associated with the processor domain.

This data structure is defined for each coordinate data structure DS_coord (namely \f$ R_x \f$, \f$ R_y \f$ and \f$ R_z \f$) which in turn is defined for each atom type in the main domain. It stores the cartesian coordinate of the atom and the start and end nodes of the intersection of main domain and a domain which is \f$r_b\f$ distance from the atom position.
*/
typedef struct   // Data structure for atom associated with processor domain
{
  double coord; ///< Cartesian coordinate
  int start;    ///< Start node index of intersection domain (index is w.r.t main domain)
  int end;      ///< End node index of intersection domain (index is w.r.t main domain)
}DS_ProcAtm;

/** \brief Data structure for each atom associated with the processor+Rcut domain.

This data structure is defined for each coordinate data structure DS_coord (namely \f$ R_x \f$, \f$ R_y \f$ and \f$ R_z \f$) which in turn is defined for each atom type in the main domain. It stores the cartesian coordinate of the atom and the nonlocal projectors. 
*/
typedef struct   // Data structure for atom associated with processor+Rcut domain
{
  double coord; ///< Cartesian coordinate
  double **UdVtm;
}DS_ProcAtmRcut;

/** \brief Coordinate Data structure defined for each atom type.

This data structure is defined for each atom type in the main domain. The replica atom information associated with the main atom is stored as members of this data structure.
 */

typedef struct   // data structure for each atom type
{
  double main_atm;     ///< Cartesian coordinate of main domain atom at current step (after Wraparound --> (t+dt))
  double main_atm_nm;  ///< Cartesian coordinate of main domain atom at current step (not Wrapped around --> (t+dt), for chg extrap)
  double main_atm_nm_old;  ///< Cartesian coordinate of main domain atom at current step (for PredicCorr) (not Wrapped around --> (t), for chg extrap)
  double main_atm_old; ///< Cartesian coordinate of main domain atom in the previous MD step (t)
  double main_atm_0dt; ///< Cartesian coordinate of main domain atom in the previous MD step (t) (not Wrapped around )
  double main_atm_1dt; ///< Cartesian coordinate of main domain atom in the 2nd previous MD step (t-1*dt) (not Wrapped around )
  double main_atm_2dt; ///< Cartesian coordinate of main domain atom in the 2nd previous MD step (t-2*dt) (not Wrapped around )
  DS_ProcAtm *ProcAtm;    ///< Coordinates & other info of atoms within rb distance from processor domain 
  DS_ProcAtmRcut *ProcAtmRcut; ///< Coordinates & other info of atoms associated with processor+Rcut domain
  int n_replica;      ///< Number of replica atoms (having rb-domain intersecting) associated with proc domain (includes main domain atom)
  int n_replica_Rcut;  ///< Number of replica atoms associated with proc+Rcut domain (includes main domain atom)
}DS_coord;

/** \brief Data structure defined for each atom type.

For each atom atom type, this data structure stores the element symbol, charge, coordinate data structure DS_coord etc.

 */

typedef struct   // data structure for each atom type
{
  char symb[10]; ///< Type of atom (symbol)
  int natm;      ///< Number of atoms of this type
  int Z;         ///< Valence charge of the nucleus  
  int lloc;      ///< Local component of the pseudopotential for this atom type
  double rb;      ///< Cut-off radius for b-calculation
  char pseudo_path[1000]; ///< Type of atom (symbol)
  double mass;    /// atomic mass in amu
  DS_coord *Rx, *Ry, *Rz; ///< Object for atom coords (both main & replicas)
}DS_Atm;

/** \brief Data structure that stores the pseudopotential information for each element.

The pseudopotential information stored in this data structure consists of the radial grid, different components of the pseudopotential, pseudowavefunctions, denominator values to be used in calculations involving non-local pseudopotential. The pseudopotentials and pseudowavefunctions are stored as the radial distance times the value at the point.

 */

typedef struct   // data structure for each pseudopotential type (corresponding to each atom type)
{
  double* rVloc; ///< Local component of the pseudopotential
  double* rVs;  ///< s-component of the pseudpotential
  double* rVp;  ///< p-component of the pseudpotential
  double* rVd;  ///< d-component of the pseudpotential
  double* rVf;  ///< f-component of the pseudpotential
  double* rUs;  ///< s-component of the pseudowavefunction
  double* rUp;  ///< p-component of the pseudowavefunction
  double* rUd;  ///< d-component of the pseudowavefunction
  double* rUf;  ///< f-component of the pseudowavefunction
  double* uu;  ///< Isolated atom electron density
  double* RadialGrid;  ///< Radial grid on which pseudopotential data is defined
  double rc_s;  ///< Pseudopotential cut-off radius for s-component
  double rc_p;  ///< Pseudopotential cut-off radius for p-component
  double rc_d;  ///< Pseudopotential cut-off radius for d-component
  double rc_f;  ///< Pseudopotential cut-off radius for f-component
  double rc;  ///< Maximum of cut-off radii across all components
  double rz;  ///< Cut-off radius at which Vloc becomes Z/r
  int lmax;  ///< Quantum number corresponding to maximum component present in the pseudopotential file. \f$l=0,1,2,3\f$ for s,p,d,f components respectively.
  int size;  ///< Size of the one-dimensional array RadialGrid (and corresponding pseudopotentials and pseudowavefunctions)
  double Denom_s; ///< Denominator value for s-component
  double Denom_p; ///< Denominator value for p-component
  double Denom_d; ///< Denominator value for d-component
  double Denom_f; ///< Denominator value for f-component
  double *DVloc,*DVsJ,*DUsJ,*DVpJ,*DUpJ,*DVdJ,*DUdJ,*DVfJ,*DUfJ; ///< Arrays to store the derivative information for spline fits
}DS_Psd;

/** \brief Data structure that stores the indices of the communication regions associated with Laplacian in Poisson solve.

*/

typedef struct   // data structure for each atom type
{
  int ereg_s[3][26],ereg_e[3][26]; // 6 faces, 12 eges, 8 corners so 26 regions disjointly make the total communication region of proc domain. 3 denotes x,y,z dirs
 int **eout_s,**eout_e,**ein_s,**ein_e,**stencil_sign,**edge_ind,*displs_send,*ncounts_send,*displs_recv,*ncounts_recv;
}DS_LapInd;

/** \brief Data structure that stores the indices of the communication regions associated with SQ.

*/
typedef struct   // data structure for each atom type
{
 int **eout_s,**eout_e,**ein_s,**ein_e,*displs_send,*ncounts_send,*displs_recv,*ncounts_recv;
}DS_SqInd;

/** \brief Data structure to store all the variables of SQDFT code.

This is a global data structure which contains all the variables and other data structures required for the SQDFT code.

 */

typedef struct   // data structure for SQ code
{
  string input_file;             ///< Input file name
  int nproc;                     ///< Total number of processors
  int nprocx,nprocy,nprocz;      ///< Number of processors in each direction of domain
  int pnode_s[3];                ///< Processor domain start nodes' indices w.r.t main (starts from 0)
  int pnode_e[3];                ///< Processor domain end nodes' indices w.r.t main (starts from 0)
  int np_x,np_y,np_z;            ///< Number of finite difference nodes in each direction of processor domain
  int fermi_alg;                 ///< Variable for choice between Newton_raphson and Brent's Algorithm
  int Correction;                ///< Variable for choice between Old and New Energy/force correction
  int vds;                       ///< Variable for choice of initial velocity distribution for MD
  int ensemble;                  ///< Variable for choice between NVE and NVT ensemble
  char name[20];                 ///< String to store the name of the output file during restart
  int FDn;                       ///< Half-the order of finite difference
  char atoms_file[10];           ///< Name to read .atoms file
  int frac_coord;                ///< Integer specifying whether atom positions are given as fractional coordinates in .atoms file
  double latconst;               ///< Size of single unit cell in each direction (lattice constant) (Bohr)
  double domain[3];              ///< Main domain size in each direction (Bohr)
  int n_int[3];                  ///< Number of finite difference intervals in each direction of main domain
  int ncell;                     ///< Number of unit cells in each direction (required for scaling tests to replicate unit cell from an atoms file)
  double perturb;                ///< All the atoms would be randomly perturbed by a max of perturb in each direction (Bohr)
  int rand_seed;                 ///< Natural number for random used for perturbing atoms
  double T;                      ///< Electronic temperature in Kelvin
  double bet;                    ///< Inverse of smearing (1/Ha)
  int npl;                       ///< Degree of Chebyshev polynomial to be used in SQ for Clenshaw-Curtis quadrature
  double Rcut;                   ///< Truncation (localization) radius of the density matrix
  int nloc;                      ///< Number of finite difference intervals in Rcut distance
  int n_atm;                     ///< Total number of atoms in the main domain
  int n_typ;                     ///< Total number of types of atoms in the main domain
  DS_Atm* Atm;                   ///< Object for each atom type 
  DS_Psd* Psd;                   ///< Object for each pseudopotential file type
  double delta;                  ///< Mesh size (Bohr)
  double* coeff_lap;             ///< Finite difference coefficients for Laplacian
  double* coeff_grad;            ///< Finite difference coefficients for gradient
  double ***b;                   ///< Nuclear (Ionic) charge density \f$b(x)\f$ in processor domain
  double Eself;                  ///< Total self energy from nuclei
  double Ncharge;                ///< Total charge of the nuclei which is equal to integral(b) over main domain
  double ***rho;                 ///< Electron density \f$\rho(x)\f$ in processor domain
  double ***rhs;                 ///< Right hand side of the Poisson's equation
  double Nelectron;              ///< Total number of electrons which is equal to integral(rho) over main domain
  double rc_tilda;               ///< Cut-off radius for reference charge denstiy used in overlap corrections
  double ***b_tilda;             ///< Reference charge density used to compute energy/force corrections
  double ***Vc;                  ///< Potential calculated in overlap corrections
  double Ecorr;                  ///< Energy correction due to overlapping charge densities
  DS_LapInd LapInd;              ///< Object for Laplacian communication information
  double poisson_tol;            ///< Convergence tolerance for Poisson solver
  int poisson_maxiter;           ///< Maximum number of iterations allowed in the Poisson solver
  int *neighs_lap;               ///< Array of neighboring processor ranks in the Laplacian stencil width from current processor
  double lanczos_tol;            ///< Convergence tolerance for Lanczos algorithm
  double fermi_tol;              ///< Convergence tolerance for Fermi energy calculation
  double scf_tol;                ///< Convergence tolerance for SCF
  int scf_miniter;               ///< Minimum number of SCF iterations
  int scf_maxiter;               ///< Maximum number of SCF iterations
  double ***phi_guess;           ///< Initial guess for Poisson solver
  double beta_aaj;               ///< Anderson mixing/or Weighted Jacobi relaxation parameter
  int m_aaj;                     ///< Number of iterates in Anderson mixing = m_aaj+1
  int p_aaj;                     ///< AAJ parameter. Anderson update done every p_aaj\f$^{th}\f$ iteration of Poisson solver
  double beta_scf;               ///< Periodic Pulay (for SCF) relaxation parameter
  int m_scf;                     ///< Number of iterates in Anderson mixing = m_scf+1 (for SCF)
  int p_scf;                     ///< AAJ parameter. Anderson update done every p_scf\f$^{th}\f$ SCF iteration
  int non_blocking;              ///< Option that indicates using non-blocking version of MPI command. 1=TRUE or 0=FALSE (for MPI collectives)
  double ***phi;                 ///< Electrostatic potential \f$\phi(x)\f$ in processor domain
  double Etot;                   ///< Total potential energy of atoms in main domain
  MPI_Comm comm_laplacian;       ///< Communicator topology for Laplacian
  MPI_Comm comm_sq;              ///< Communicator toplogy for Spectral Quadrature
  DS_SqInd SqInd;                ///< Object for SQ communication information
  double ***Vxc;                 ///< Exchange correlation potential
  double ***Veff;                ///< Effective potential (\f$V_{eff}=\phi+V_{xc}+V_{nloc}\f$)
  int *neighs_sq;                ///< Array of neighboring processor ranks within Rcut distance from current processor
  int nneigh_sq;                 ///< No. of neighboring processors in neighs_sq array
  double lambda_f;               ///< Fermi energy
  double **rho_pj;               ///< Components of Chebyshev expansion in SQ, for all nodes in the processor domain and for all degrees until npl
  double *chi;                   ///< Array of scaling parameters for sub-Hamiltonains for all nodes in the processor domain
  double *zee;                   ///< Array of scaling parameters for sub-Hamiltonains for all nodes in the processor domain
  double **Ci;                   ///< Chebyshev polynomial coefficients in the SQ method to find electron density
  double Ebs;                    ///< Band structure energy 
  double Eent;                   ///< Entropy energy	
  double *flocx,*flocy,*flocz;   ///< Local component of forces on all atoms
  double *fnlocx,*fnlocy,*fnlocz;///< Non-local component of forces on all atoms
  double *fcorrx,*fcorry,*fcorrz;///< Force correction due to overlapping charge densities
  double *fx,*fy,*fz,*forces;    ///< Total forces on all atoms
  double poiss_time,sq_time,engy_time,forcs_time; ///< Scalar variables to record time taken
  double tpoiss_mpi,tsq_mpi,tengy_mpi,tforcs_mpi; ///< Scalar variables to record mpi communication time taken
  double mpi_time;               ///< Dummy variable to compute time taken by mpi communication in each of the functions
  int prnt_atoms;                ///< Print atoms info into a .aout file
  double *vels,*vels_old;        ///< Velocities of all atoms (Bohr/fs)
  double *accels;                ///< Accelerations of all atoms (Bohr/(fs)^2)
  double time_step;              ///< Molecular Dynamics time step (in femto sec (fs))
  int MaxMDsteps;                ///< Maximum number of MD steps
  int MDstepCount;               ///< MD step count index
  double kB;                     ///< Boltzmann constant
  double Ceh;                    ///< Conversion factor from eV to Ha
  int dof;                       ///<  Degree of freedom of the system
  double TE;                     ///< Total Energy of the atoms (KE+PE)
  double KE;                     ///< Kinetic Energy of the atoms
  double PE;                     ///< Potential Energy of the atoms
  double ***rho_at;              ///< Superposition of isolated atom electron densities
  double ***drho_new,***drho,***drho_dt,***drho_2dt; ///< Arrays to be used for charge extrapolation in MD (drho=rho-rho_at)
  int ChgExtrap;                 ///< Using charge extrapolation for rho_guess in each MD step
  int restart_scf;                ///< Reads electron density guess from a .restart file
  int prnt_scf;                   ///< Prints electron density at latest scf iteration into .restart file
  double Treq;                   ///< Required temperature (K) for NVT MD (same as T in input file)
  double Ms,v2nose,xi_nose,snose;///< coupling parameter for NVT MD using NH
  double TEext;                  ///< Total energy of the extended system (including heat bath) for NVT MD using NH
  double T_MD;                   ///< Inonic temperature (can be different from electronic temperature)
  int MDrestartCount;            ///< MD step count in the latest restart file
  double mean_TE,std_TE,mean_PE,std_PE,mean_KE,std_KE,mean_T,std_T,mean_TEext,std_TEext; ///< MD statistics
  int restart_md;                ///< Reads MD info from a .restartMD file
  int prnt_md;                   ///< Prints latest MD info into .restartMD file
  int RelaxAtoms;                ///< Perform atomic relaxation (energy minimization)
  }DS_SQ;
#endif
