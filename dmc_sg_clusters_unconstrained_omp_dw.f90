!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS IN INPUT FILE/GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nparticles: Number of particles (atoms or rigid molecules) in the system
! N0: Initial number of random walkers
! dtau: Time step
! tmax: Maximum projection time for single DMC run
! Kaver: Controls projection times at which the energies are written out and the configurations are stored for restarting the simulation
! Vref: Input Vref value if the simulation is not started from the global minimum; also the instantaneous estimator of the DMC ground state energy
! Ndw: Number of times walker weights are computed for descendant weighting
! dw_tau_max: Projection time at which the observables are computed with DW and the DW parameters are reset
! dw_tau_incr: Increment that controls the number of DMC steps that elapse before the weights for descendant weighting are computed again
! quench_equilib: Number of DMC steps after which the instantaneous isomer multiplicities do not write out
! hist_pairs_min,hist_pairs_max,nb_pairs: Determines the minimum value, maximum value, and number of bins for the pair correlation function
! hist_CM_max,nb_CM: Determines the maximum value, and number of bins for the density profile function
! hist_Q4_min,hist_Q4_max,nb_Q4: Determines the minimum value, maximum value, and number of bins for the Q4 distribution
! hist_Q6_min,hist_Q6_max,nb_Q6: Determines the minimum value, maximum value, and number of bins for the Q6 distribution
! rc_sphere: Radius of the constraining sphere
! coord_config0: String that gives the name of the initial starting configuration to be read
! atom_type: Character that determines the single type of atom comprising the Lennard Jones cluster
! pseudo_atom_mass: Mass of the ficitious/fales atom used in the DMC simulation used to determine the quantum delocalization lambda
! start_configs: String that controls how the configurations are initialized
! restart: String that controls if the simulation is restarted by changing the random walker population or by changing the atom mass
! structure_parms: Logical variable that determines if the pair correlation function, the density profile, and the Q4 and Q6 distributions are computed
! quenching: Logical variable that determines if quenching and isomer assignments are performed
! write_configs: Logical variable that determines if all of the configurations are written out to the disk at the end of each DMC run
! Nt_incr: Increment that controls the number of configurations that get quenched
! NMMC: Number of Monte Carlo steps used in the Metropolis algorithm
! incr_MMC: Increment that determines when configurations from the Metropolis algorithm are written out to an xyz file
! rc_MMC: Radius of the constraining sphere for the Metropolis Monte Carlo runs
! Temp: Temperature used during the Metropolis Monte Carlo runs
! NDMC: Number of DMC simulations that are run with different target random walker populations
! Vref_parm: Determines if the initial Vref is taken directly from the input file or is initialized from the global minimum configuration
! N0_target_max: The maximum target random walker population for the last DMC run
! Nisomers_max: Maximum number of configurations/isomers that can be stored after quenching
! Nisomers_thresh: Maximum number of configurations/isomers that are written out to the disk
! Nisomers_configs: Number of configurations actually written out with xyz coordinates
! step_quench: DMC steps at which the configurations are quenched
! r_energy_thresh, r_Q4_thresh, r_Q6_thresh: Threshold values for assigning the quenched energy, Q4, and Q6 values to isomers stored in the library
! grad_tol: Tolerance for the congugate gradient subroutine that determines the number of gradient evaluations
! grad_maxfac: Factor that determines the maximum number of gradient evaluations
! energy_thresh: Energy threshold that determines if a configuration is unphysical
! num_threads: Total number of threads that can be used in the parallel regions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS IN INPUT FILE/GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigma_lj: Core radius of the LJ potential fit to the SG potential for the LJ cluster
! epsilon_lj: Well depth of the LJ potential fit to the SG potential for the LJ cluster
! Dim: Dimensionality of the configuration space
! Npairs: Number of particle pairs in the system
! N0_part: Integer fraction that determines the number of random walkers that need to be added to the population if N0 is large
! Nmax: Maximum number of random walkers that can be present in the system
! Nt: Instantaneous number of random walkers
! Nt1: Number of random walkers stored prior to descendant weighting (DW)
! N_iter_max: Maximum number of DMC steps calculated from tmaxand dtau
! m: Current DMC simulation number
! t: Current iteration number in loop 1,N_iter_max (outer loop of main program)
! Nt_quench: Total number of random walkers that need to be quenched
! dw_step_max: Number of steps at which the observables are computed and the configurations for DW are reset 
! dw_step0: Specifies the iteration number at which the DW weights start being accumulated
! dw_step: Iteration number at which one set of weights is accumulated and stored
! dw_incr: Determines the number of time steps that elapse until the next set of weights are computed for descendant weighting
! Epartaver: Partial average of the DMC energy over Kaver steps
! seed(num_threads): Array that stores the random number generator seeds for each individual thread
! sigma: Width of the Gaussian distribution for the DMC moves
! x(Dim,Nmax): Array that stores the instantaneous position in configuration space for each configuration/random walker
! dw_count: Counts the number of non-zero walkers that were used to compute the observeables for each value of tau dw
! dw_walkers(Nmax): Array of walker indicies for descendant weighting
! dw_weight(Ndw,Nmax): Array of weights for each walker from descendant weighting used to compute observables
! dw_Nt(Ndw): Array that stores the current walker numbers for use in calculating observables
! dw_x(Dim,Nmax): Array of confiugrations for descendant weighting used to compute observables
! d_pairs: Distance increment for the pair correlation function
! d_CM: Distance increment for the density profile function
! d_Q4: Distance increment for the Q4 cluster distribution
! d_Q6: Distance increment for the Q6 cluster distribution
! pairs_hist(0:nb_pairs,Ndw): Pair correlation function of the full parahydrogen cluster
! CM_hist(0:nb_CM,Ndw): Density profile of the full parahydrogen cluster
! Q4_hist(0:nb_Q4,Ndw): Q4 distribution of the full parahydrogen cluster
! Q6_hist(0:nb_Q6,Ndw): Q6 distribution of the full parahydrogen cluster
! file_num: Number that specifies which DMC simulation the current file belongs to
! file: String that gives the name of the file for the current DMC simulation
! Nisomers: Running tally of the number of isomers stored in the library
! DMC_num: String for the DMC simulation number that is written on the isomer library directory name
! count_config(Nisomers_max): Populations/isomer fractions of each of the isomers stored in the library for each value of tau dw
! lib_energy(Nisomers_max): Stores the classical energies of each isomer in the library
! lib_Q4(Nisomers_max): Stores the Q4 values of each isomer in the library
! lib_Q6(Nisomers_max): Stores the Q6 values of each isomer in the library
! lib_energy_HA(Nisomers_max): Stores the harmonic approximation ground state energies for each isomer
! lib_logW(Nisomers_max): Stores the log of the product of the frequncies from the HA for each isomer
! lib_order(Nisomers_max): Stores the point group order for each isomer
! lib_pgsym(Nisomers_max): Stores the point group symbol for each isomer 
! lib_count_config(Nisomers_max): Stores the running averaged populations of each isomer in the library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE FOR METROPOLIS MONTE CARLO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! U(Nparticles,Nparticles): Array that stores the energies for all of the pairs of particles
! U_move(Nparticles): Array that stores the energies for single pairs of the particles for a particle that has undergone a trial move
! coord_MMC: String for the name of the xyz file from Metropolis Monte Carlo
! max_displacement: Maximum displacement for a particle in the cluster
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE FOR METROPOLIS MONTE CARLO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN DIFFUSION SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel: Logical variable that is true if the program is running in parallel and false if not
! i,k: Loop indicies
! n: Index for loop over random walkers
! thread: Current thread being used in the parallel region; index of the seed array for the random number generator
! jsr: Instantaneous value of the seed for the random number generator
! Nreplicate: Number of random walkers designated to be replicated during branching: population=2 or higher
! Nkill: Number of random walkers designated to be killed during branching: populations=0
! Nreplicate1: Number of random walkers to designated to be replicated when the number of walkers to be killed is larger than the number to be replicated
! population(Nmax): Array that determines which random walkers are replicated or killed
! replicate(Nmax): Array of indicies for random walkers that are replicated during branching
! kill(Nmax): Array of indicies for random walkers that are killed during branching
! u: Uniform random number on the interval (0,1)
! rho1,rho2: Two uniform random numbers on the interval (0,1) used to compute a Gaussian random number
! Vaver: Average potential energy computed from population of random walkers
! Pb: Expoential function used to determine which random walkers are replicated or killed
! energy: Potential energy of a configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN DIFFUSION SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES in HISTOGRAMS SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel: Logical variable that says whether or not the subroutine is running in parallel
! i,j,k,l,n: Loop indicies
! info: 0 if Q4 and Q6 values are computed correctly; 1 if Q4 and Q6 values for the configuration should be ignored
! ignore_config: Number of configurations whose Q4 and Q6 values were ignored
! pair_dist: Single pair distance for a given configuration
! Q4: Q4 structure parameter for the full cluster/given configuration
! Q6: Q6 structure parameter for the full cluster/given configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES in HISTOGRAMS SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN QUENCHING SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel: Logical variable that says whether or not the subroutine is running in parallel
! i,k,n: Loop indicies
! info, iter, nfunc, ngrad: Parameters output by the conjugate gradient subroutine
! nquench, ngrad_sum, count_unphysical: Stores the number of successful quenches, number of gradient evaluations, and number of unphysical configurations
! info_order: Parameter output by the subroutine that computes Q4 and Q6 values
! x_quench(Nt_quench): Array that stores the quenched configurations with non-zero weights
! energy_opt(Nt_quench): Stores the classical energies of quenched configurations wiht non-zero weights
! Q4(Nt_quench): Stores the Q4 values of the quenched configurations with non-zero weights
! Q6(Nt_quench): Stores the Q6 values of the quenched configurations with non-zero weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN QUENCHING SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! i,j,n: Loop indicies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARALLEL REGIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random number generator found on John Burkardt's website.
! URL: https://people.sc.fsu.edu/~jburkardt/f_src/ziggurat_openmp/ziggurat_openmp.f90
! Parallel Region 1: Computation of energies of each random walker and determination of which random walkers will be replicated or killed in the diffusion subroutine
! Parallel Region 2: Replacement of random walkers to be killed with those to be replicated in the diffusion subroutine
! Parallel Region 3: If the number of random walkers to be killed is larger than those to be replicated, then replace the killed ones by those on the top of the replicate array in the diffusion subroutine
! Parallel Region 4: Compute the pair correlation functions, density profiles, Q4 and Q6 distributions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARALLEL REGIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OTHER IMPORTANT NOTES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code has no restart capability as the random walker populations used for the SG clusters are too large to write out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OTHER IMPORTANT NOTES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Module containing all of the global variables and subroutines for the DMC true ground state code
module dmc_module
  Implicit None
  Double Precision, Parameter :: bohr=0.52917721092
  Double Precision, Parameter :: autokcalmol=627.51 
  Double Precision, Parameter :: autokJmol=2625.5002
  Double Precision, Parameter :: autoK=315775.13
  Double Precision, Parameter :: melectron=1822.88839
  Double Precision, Parameter :: deg=180/dacos(-1d0)
  Double Precision, Parameter :: pi=dacos(-1d0)
  Double Precision, Parameter :: days=86400,hours=3600
  Double Precision, Parameter :: sigma_lj=3.02179 ! angstroms
  Double Precision, Parameter :: epsilon_lj=1.05768620853707D-4 ! atomic units  
  Double Precision, Parameter :: H2mass=2*1.00782503223*melectron
  Double Precision, Parameter :: D2mass=2*2.01410177812*melectron
  Double Precision, Parameter :: Hemass=4.00260325413*melectron
  Double Precision, Parameter :: Nemass=19.9924401762*melectron
  Double Precision, Parameter :: Armass=39.9623831237*melectron
  Double Precision, Parameter :: Krmass=83.9114977282*melectron
  Double Precision, Parameter :: Xemass=129.903509349*melectron
  Double Precision, Parameter :: cmtoau=4.5563352527D-6
  Double Precision, Parameter :: freq_cutoff=2*cmtoau
  
  ! Global variables for main program and/or DMC code
  Integer :: Nparticles,Dim,Npairs,Nmax,N0,Nt,Nt1,Nt_incr,Nt_quench,N0_target_max
  Integer :: N_iter_max,Kaver,m,t,num_threads,NDMC,N0_part,quench_equilib,step_quench,Nisomers_configs
  Double Precision :: dtau,tmax,sigma,Vref,Epartaver,rc_sphere,lambda,pseudo_atom_mass
  Integer, Allocatable :: seed(:)
  Real, Allocatable :: x(:,:)
  Character(len=2) :: atom_type
  Character(len=50) :: coord_config0,file,file_num,DMC_num,start_configs,restart
  Logical :: structure_parms,quenching,write_configs,Vref_parm

  ! Global variables for Metropolis Monte Carlo
  Integer :: NMMC,incr_MMC
  Double Precision :: rc_MMC,Temp
  Character(len=50) :: coord_MMC
  Double Precision :: max_displacement=0.1
  Double Precision, Allocatable :: U(:,:),U_move(:)

  ! Global variables for descendant weighting
  Integer :: Ndw,dw_step_max,dw_step0,dw_step,dw_incr,dw_count
  Double Precision :: dw_tau_max,dw_tau0,dw_tau_incr
  Integer, Allocatable :: dw_walkers(:),dw_weight(:,:),dw_Nt(:)
  Real, Allocatable :: dw_x(:,:)
  
  ! Global variables for the pair correlation function, structure parameters (Q4 and Q6), and potential energy
  Integer :: nb_pairs,nb_CM,nb_Q4,nb_Q6
  Double Precision :: hist_pairs_min,hist_Q4_min,hist_Q6_min
  Double Precision :: hist_pairs_max,hist_CM_max,hist_Q4_max,hist_Q6_max
  Double Precision :: d_pairs,d_CM,d_Q4,d_Q6
  Integer(kind=8), Allocatable :: pairs_hist(:,:),CM_hist(:,:),Q4_hist(:,:),Q6_hist(:,:)

  ! Global variables for quenching subroutine
  Integer :: Nisomers,Nisomers_max,Nisomers_thresh,grad_maxfac
  Double Precision :: grad_tol,energy_thresh,r_energy_thresh,r_Q4_thresh,r_Q6_thresh
  Character(len=4), Allocatable :: lib_pgsym(:)
  Integer, Allocatable :: count_config(:),lib_count_config(:),lib_order(:)
  Double Precision, Allocatable :: lib_energy(:),lib_Q4(:),lib_Q6(:),lib_energy_HA(:),lib_logW(:)

contains
  
  ! Silvera-Goldman pairwise additive potential for para_H2
  function pot_SG_opt(x,Dim)
    Implicit None 
    
    ! Constants for the repulsive exponential term at short range in angstroms
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41    
    
    ! Variables used to compute the potential energy
    Integer, Intent(In) :: Dim
    Real, Intent(In) :: x(Dim)
    Integer :: i,j
    Double Precision :: pot_SG_opt,r,r2 ! angstroms for r and r2; pot_SG_opt in atomic units
    
    pot_SG_opt=0d0
    do i=2,Nparticles
       do j=1,i-1
          r2=sum((x(i*3-2:i*3)-x(j*3-2:j*3))**2)
          r=dsqrt(r2)
          r2=1/r2         
          if (r<cutoff) then ! damp the dispersion interactions at short range
             pot_SG_opt=pot_SG_opt+dexp(alpha-(beta+gamma*r)*r)-dexp(-((cutoff/r)-1)**2)*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          else ! use the undamped dispersion interactions at long range
             pot_SG_opt=pot_SG_opt+dexp(alpha-(beta+gamma*r)*r)-((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          endif
       enddo
    enddo
  end function pot_SG_opt

  ! Silvera-Goldman pairwise additive potential for para_H2 and conjugate gradient subroutine
  subroutine pot_SG_opt_cg(energy,x,Dim)
    Implicit None 
    
    ! Constants for the repulsive exponential term at short range in angstroms
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41    
    
    ! Variables used to compute the potential energy
    Integer, Intent(In) :: Dim
    Double Precision, Intent(In) :: x(Dim)
    Integer :: i,j
    Double Precision :: r,r2 ! angstroms for r and r2
    Double Precision, Intent(Out) :: energy ! energy in atomic units
    
    energy=0d0
    do i=2,Nparticles
       do j=1,i-1
          r2=sum((x(i*3-2:i*3)-x(j*3-2:j*3))**2) 
          r=dsqrt(r2)
          r2=1/r2         
          if (r<cutoff) then ! damp the dispersion interactions at short range
             energy=energy+dexp(alpha-(beta+gamma*r)*r)-dexp(-((cutoff/r)-1)**2)*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          else ! use the undamped dispersion interactions at long range
             energy=energy+dexp(alpha-(beta+gamma*r)*r)-((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          endif
       enddo
    enddo
  end subroutine pot_SG_opt_cg

  ! Gradients for Silvera-Goldman potential
  subroutine grad_SG_opt_cg(grad,x,Dim)
    Implicit None
    
    ! Constants for the repulsive exponential term at short range in angstroms
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41    
    
    ! Variables used to compute the gradient
    Integer, Intent(In) :: Dim
    Double Precision, Intent(In) :: x(Dim)
    Integer :: i,j
    Double Precision :: r,r2,xij(3) ! angstroms for r, r2, and xij
    Double Precision, Intent(Out) :: grad(Dim) ! atomic units per angstrom
    
    grad=0d0
    do i=2,Nparticles
       do j=1,i-1
          xij(:)=x(i*3-2:i*3)-x(j*3-2:j*3)
          r2=sum(xij(:)**2)
          r=dsqrt(r2)
          r2=1/r2         
          if (r<cutoff) then ! damp the dispersion interactions at short range
             xij(:)=(-(2*gamma*r+beta)*dexp(alpha-(beta+gamma*r)*r) &
                  +dexp(-((cutoff/r)-1)**2)*((((6*C6/r)+((8*C8/r)+((10*C10/r)-9*C9)*r2)*r2)*r2**3) &
                  +2*cutoff*r2*(1-(cutoff/r))*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3))) &
                  *xij(:)/r
          else ! use the undamped dispersion interactions at long range to compute the gradient function
             xij(:)=(-(2*gamma*r+beta)*dexp(alpha-(beta+gamma*r)*r) &
                  +(((6*C6/r)+((8*C8/r)+((10*C10/r)-9*C9)*r2)*r2)*r2**3)) &
                  *xij(:)/r
          endif
          grad(i*3-2:i*3)=grad(i*3-2:i*3)+xij(:)
          grad(j*3-2:j*3)=grad(j*3-2:j*3)-xij(:)
       enddo
    enddo
  end subroutine grad_SG_opt_cg

  ! Silvera-Goldman potential: single pair array
  subroutine U_array(x,Dim)
    Implicit None 
    
    ! Constants for the repulsive exponential term at short range in angstroms
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41    
    
    ! Variables used to compute the potential energy
    Integer, Intent(In) :: Dim
    Real, Intent(In) :: x(Dim)
    Integer :: i,j
    Double Precision :: r,r2 ! angstroms for r and r2
    
    U=0d0
    U_move=0d0
    do i=2,Nparticles
       do j=1,i-1
          r2=sum((x(i*3-2:i*3)-x(j*3-2:j*3))**2)
          r=dsqrt(r2)
          r2=1/r2         
          if (r<cutoff) then ! damp the dispersion interactions at short range
             U(i,j)=dexp(alpha-(beta+gamma*r)*r)-dexp(-((cutoff/r)-1)**2)*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          else ! use the undamped dispersion interactions at long range
             U(i,j)=dexp(alpha-(beta+gamma*r)*r)-((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
          endif
          U(j,i)=U(i,j)
       enddo
    enddo
  end subroutine U_array

  ! Lennard-Jones potential: single pair
  function U_pair(x,Dim,k,i)
    Implicit None

    ! Constants for the repulsive exponential term at short range in angstroms
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41    
    
    ! Variables used to compute the potential energy
    Integer, Intent(In) :: i,k,Dim 
    Real, Intent(In) :: x(Dim)
    Integer :: j
    Double Precision :: U_pair,r,r2 ! angstroms for r2, atomic units for U_pair
    
    r2=sum((x(k*3-2:k*3)-x(i*3-2:i*3))**2)
    r=dsqrt(r2)
    r2=1/r2         
    if (r<cutoff) then ! damp the dispersion interactions at short range
       U_pair=dexp(alpha-(beta+gamma*r)*r)-dexp(-((cutoff/r)-1)**2)*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
    else ! use the undamped dispersion interactions at long range
       U_pair=dexp(alpha-(beta+gamma*r)*r)-((C6+(C8+C10*r2-(C9/r))*r2)*r2**3)
    endif
  end function U_pair

  subroutine CM(x,x_CM,flag)
    Implicit None
    Logical, Intent(In) :: flag
    Real, Intent(Inout) :: x(Dim)
    Double Precision, Intent(Out) :: x_CM(3)
    Integer :: i

    x_CM=0d0
    do i=1,Nparticles
       x_CM(:)=x_CM(:)+x(3*i-2:3*i)
    enddo
    x_CM=x_CM/Nparticles
    if (flag) then
       do i=1,Nparticles
          x(3*i-2:3*i)=x(3*i-2:3*i)-x_CM(:)
       enddo
    endif
  end subroutine CM

  subroutine Cluster_R(x,R)
    Implicit None
    Real, Intent(Inout) :: x(Dim)
    Double Precision, Intent(Out) :: R
    Double Precision :: x_CM(3),dist
    Integer :: j

    call CM(x,x_CM,.false.)
    R=0d0
    do j=1,Nparticles
       dist=sum((x(j*3-2:j*3)-x_CM(:))**2)
       dist=dsqrt(dist)
       if (R<dist) R=dist
    enddo
  end subroutine Cluster_R

  subroutine Metropolis_MC(NMMC,x,beta,Rc)
    Implicit None
    Integer, Intent(In) :: NMMC
    Double Precision, Intent(In) :: beta,Rc ! beta in 4*epsilon/kT, Rc in angstroms
    Real, Intent(Inout) :: x(Dim) ! coordinates in angstroms  
    Integer :: i,j,k,Nrestart,accept,all
    Real :: x0(3)
    Double Precision :: Delta_E,x_CM(3),R ! energies in reduced units, distances in angstroms
    Double Precision :: r1,t,s(3) ! harvest for random numbers
    
    ! Initialize all counters and the energies for all single pairs in the cluster
    Nrestart=min(1000,NMMC)
    accept=0
    all=0
    
    do j=1,NMMC
       call random_number(r1) ! used to decide which particle undergos displacement
       i=floor(r1*Nparticles)+1 ! displace particle i
       if (i>Nparticles) i=Nparticles
       call random_number(t) ! used to accept the trial move
       call random_number(s) ! used for trial displacement
       s=2*s-1

       x0(1:3)=x(3*i-2:3*i) ! hold the coordinates of the i-th atom
       all=all+1
       x(3*i-2:3*i)=x(3*i-2:3*i)+max_displacement*s
       call Cluster_R(x,R)
       if (R<Rc) then ! Compute the energy difference involving only the ith particle
          Delta_E=0d0
          do k=1,Nparticles
             if (k.ne.i) then
                U_move(k)=U_pair(x,Dim,k,i)
                Delta_E=Delta_E+U(k,i)-U_move(k)
             endif
          enddo
          
          ! Compute the energies of the configurations for only the particle that underwent the displacement
          if (dexp(beta*Delta_E)>=t) then 
             ! accept the new configuration
             do k=1,Nparticles
                if (k.ne.i) then
                   U(k,i)=U_move(k)
                   U(i,k)=U(k,i)
                endif
             enddo
             accept=accept+1
             goto 113
          endif
       endif
       ! reject the new configuration
       x(3*i-2:3*i)=x0(1:3)
113    continue
          
       if (mod(j,Nrestart)==0) then
          call CM(x,x_CM,.true.) ! Translate the center of mass
          ! Adjust the displacement parameters every 1000 moves for ~50% acceptance rates
          if (dble(accept)/all<0.5) then 
             max_displacement=max_displacement*0.9
          else
             max_displacement=max_displacement*1.1
          endif
          accept=0
          all=0
       endif
    enddo
  end subroutine Metropolis_MC
    
  function grad_pair(Q1,Q2)
    Implicit None
    Double Precision, Parameter :: alpha=1.713
    Double Precision, Parameter :: beta=1.5671/bohr 
    Double Precision, Parameter :: gamma=0.00993/bohr**2    
    
    ! Constants for the multipole expansion of the dispersion interaction in angstroms
    Double Precision, Parameter :: C6=12.14*bohr**6     
    Double Precision, Parameter :: C8=215.2*bohr**8
    Double Precision, Parameter :: C9=143.1*bohr**9
    Double Precision, Parameter :: C10=4813.9*bohr**10  
    
    ! Cutoff for the dipsersion interaction damping function in angstroms
    Double Precision, Parameter :: cutoff=1.28*3.41
    Double Precision, Intent(In) :: Q1(3),Q2(3)
    Double Precision :: Q12(3),grad_pair(3),r,r2 ! angstroms for r, r2, and Q12
    
    Q12=Q1-Q2
    r2=sum(Q12(:)**2)
    r=dsqrt(r2)
    r2=1/r2         
    if (r<cutoff) then ! damp the dispersion interactions at short range
       grad_pair(:)=(-(2*gamma*r+beta)*dexp(alpha-(beta+gamma*r)*r) &
            +dexp(-((cutoff/r)-1)**2)*((((6*C6/r)+((8*C8/r)+((10*C10/r)-9*C9)*r2)*r2)*r2**3) &
            +2*cutoff*r2*(1-(cutoff/r))*((C6+(C8+C10*r2-(C9/r))*r2)*r2**3))) &
            *Q12(:)/r
    else ! use the undamped dispersion interactions at long range to compute the gradient function
       grad_pair(:)=(-(2*gamma*r+beta)*dexp(alpha-(beta+gamma*r)*r) &
            +(((6*C6/r)+((8*C8/r)+((10*C10/r)-9*C9)*r2)*r2)*r2**3)) &
            *Q12(:)/r
    endif
  end function grad_pair

  subroutine Hessian_pair(x,H)
    ! Compute the full mass scaled Hessian by a symmetrized finite difference for a pairwise monoatomic potential with pair gradient function grad_pair(r1,r2)
    Implicit None
    Double Precision, Parameter :: s=1d-5 
    Double Precision, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: H(Dim,Dim)
    Double Precision :: r1(3),r2(3)
    Integer :: i,j,k
    
    H=0d0
    do i=1,Nparticles-1
       r1(:)=x(3*i-2:3*i)
       r2(:)=r1(:)
       do k=1,3          
          r1(k)=x(k+(i-1)*3)+s/2
          r2(k)=x(k+(i-1)*3)-s/2
          do j=i+1,Nparticles
             H(3*(i-1)+k,3*j-2:3*j)=(grad_pair(r2,x(3*j-2:3*j))-grad_pair(r1,x(3*j-2:3*j)))/s
          enddo
          r1(k)=x(k+(i-1)*3)
          r2(k)=x(k+(i-1)*3)
       enddo
    enddo
    do i=1,Nparticles-1
       do j=i+1,Nparticles
          H(3*j-2:3*j,3*i-2:3*i)=H(3*i-2:3*i,3*j-2:3*j)
          H(3*i-2:3*i,3*i-2:3*i)=H(3*i-2:3*i,3*i-2:3*i)-H(3*i-2:3*i,3*j-2:3*j)
          H(3*j-2:3*j,3*j-2:3*j)=H(3*j-2:3*j,3*j-2:3*j)-H(3*i-2:3*i,3*j-2:3*j)
       enddo
    enddo
    if (atom_type=='Uk') then
       H=H/pseudo_atom_mass
    else if (atom_type=='H') then
       H=H/atom_mass('H2')
    else if (atom_type=='D') then
       H=H/atom_mass('D2')
    else if (atom_type=='He') then
       H=H/atom_mass('He')
    else if (atom_type=='Ne') then 
       H=H/atom_mass('Ne')
    else if (atom_type=='Ar') then
       H=H/atom_mass('Ar')
    else if (atom_type=='Kr') then 
       H=H/atom_mass('Kr')
    else if (atom_type=='Xe') then
       H=H/atom_mass('Xe')
    endif
  end subroutine Hessian_pair

  subroutine harm_approx_gs(x,energy,logW,IERR)
    Implicit None
    Double Precision, Intent(In) :: x(Dim)
    Integer, Intent(Out) :: IERR
    Double Precision, Intent(Out) :: energy,logW
    Integer :: i
    Double Precision :: H(Dim,Dim),omega(Dim),FV1(Dim),FV2(Dim),dummy
    
    call Hessian_pair(x,H)
    ! Diagonalize the mass-scaled Hessian
    call RS(Dim,Dim,H,omega,0,dummy,FV1,FV2,IERR)
    do i=Dim,1,-1
       omega(i)=sign(dsqrt(dabs(omega(i))),omega(i))*bohr ! Frequencies
    enddo
    call pot_SG_opt_cg(energy,x,Dim)
    logW=0d0
    IERR=0    
    do i=1,Dim
       if (omega(i)>freq_cutoff) then 
          energy=energy+omega(i)/2
          logW=logW+log(omega(i))
       else
          IERR=IERR+1
       endif
    enddo
  end subroutine harm_approx_gs

  function atom_mass(atom)
    Implicit None
    Double Precision :: atom_mass
    Character(len=2), Intent(In) :: atom
    
    if (atom=='H2') then
       atom_mass=H2mass
    else if (atom=='D2') then
       atom_mass=D2mass
    else if (atom=='He') then
       atom_mass=Hemass
    else if (atom=='Ne') then
       atom_mass=Nemass
    else if (atom=='Ar') then
       atom_mass=Armass
    else if (atom=='Kr') then
       atom_mass=Krmass
    else if (atom=='Xe') then
       atom_mass=Xemass
    else
       write(16,*) 'atom ', atom, ' is not recognized'
       stop 
    endif
    return
  end function atom_mass

  subroutine read_config(xyz_file,x)
    Implicit None
    Character(len=30), Intent(In) :: xyz_file
    Character(len=2) :: dummy
    Integer :: i
    Real, Intent(Out) :: x(Dim)
    
    open(3,file=xyz_file)
    read(3,*) 
    read(3,*)
    x=0d0
    do i=1,Nparticles
       read(3,*) dummy,x(3*i-2:3*i)
    enddo
    close(3)
  end subroutine read_config
  
  ! Integer generator to produce the random number seed for each omp thread
  subroutine shr3(jsr,rand_int)
    Implicit None
    Integer :: jsr_input
    Integer, Intent(Inout) :: jsr 
    Integer, Intent(Out) :: rand_int
    
    jsr_input=jsr
    
    jsr=ieor(jsr,ishft(jsr,13))
    jsr=ieor(jsr,ishft(jsr,-17))
    jsr=ieor(jsr,ishft(jsr,5))
    
    rand_int=jsr_input+jsr
  end subroutine shr3
  
  ! Random number generator for uniform distribution in the range (0,1): Does not include 0 and 1 explicitly
  subroutine r4_uni(jsr,rand_num)
    Implicit None
    Integer :: jsr_input
    Integer, Intent(Inout) :: jsr
    Double Precision, Intent(Out) :: rand_num
    
    jsr_input=jsr
    
    jsr=ieor(jsr,ishft(jsr,13))
    jsr=ieor(jsr,ishft(jsr,-17))
    jsr=ieor(jsr,ishft(jsr,5))
    
    ! This random number generator does not include 0
    rand_num=0.5E+00+0.2328306E-09*real(jsr_input+jsr,kind=8)
    
    ! This random number generator is commented out because it does include 0. This will make the diffusion step blow up
    ! rand_num=0.5E+00+real(jsr_input+jsr,kind=8) & 
    ! /real (65536,kind=8)/real(65536,kind=8)
  end subroutine r4_uni

  subroutine particle_pair_distances(x,pair_dist)
    Implicit None
    Integer :: i,j,l
    Real, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: pair_dist(Npairs)
    
    ! Given a configuration of particles, compute all of the pair distances
    pair_dist=0d0
    l=0
    do i=2,Nparticles
       do j=1,i-1
          l=l+1
          pair_dist(l)=sum((x(i*3-2:i*3)-x(j*3-2:j*3))**2)
       enddo
    enddo
    pair_dist(:)=dsqrt(pair_dist(:))
  end subroutine particle_pair_distances

  function cluster_size(x)
    Implicit None
    Integer :: j
    Real :: x(Dim)
    Double Precision :: CM(3),dist
    Logical :: cluster_size
    
    ! Determine the distances of each particles in the cluster to the center of mass
    CM=0d0
    do j=1,Nparticles
       CM(:)=CM(:)+x(3*j-2:3*j)
    enddo
    CM=CM/Nparticles
    cluster_size=.true.
    do j=1,Nparticles
       dist=sum((x(3*j-2:3*j)-CM(:))**2)
       if (dist>rc_sphere**2) then
          cluster_size=.false.
          return
       endif
    enddo

    ! Translate the cluster such that the center of mass is at the origin
    do j=1,Nparticles
       x(3*j-2:3*j)=x(3*j-2:3*j)-CM(:)
    enddo
  end function cluster_size

  subroutine diffusion()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,k,n,thread,jsr,Nreplicate,Nkill,Nreplicate1
    Integer :: population(Nmax),replicate(Nmax),kill(Nmax)
    Double Precision :: u,rho1,rho2,Vaver,Pb,energy

    ! If this is a new simulation at the very first iteration, then compute a separate random number generator seed for each omp thread
    if (m==0.and.t==1) then
       jsr=123456789
       do i=0,num_threads-1 
          call shr3(jsr,seed(i))
       enddo
    endif

    ! Initialize Vaver, the population array, and the radius of the clusters
    Vaver=0d0
    population(1:Nt)=1

    ! Parallel Region: Diffusive displacements/random walk
    !$omp parallel &
    !$omp shared (Dim,dtau,Nt,x,Vref,sigma,population,seed) &                 
    !$omp private (n,k,thread,jsr,energy,Pb,rho1,rho2,u)
    
    ! On the first iteration only, determine if the code is running in parallel and write out the number of threads that will be used 
    if (t==1) then
       num_threads=omp_get_num_threads() 
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif
    
    ! Loop over the random walkers
    !$omp do reduction (+:Vaver)
    do n=1,Nt
       ! Get the number of the current thread in use and the instantaneous value of the seed
       thread=omp_get_thread_num()
       jsr=seed(thread)
       do k=1,Dim
          ! Compute two uniform random numbers to calculate a Gaussian random number
          call r4_uni(jsr,rho1)
          if (rho1==0d0) call r4_uni(jsr,rho1) ! Ensures that 0 will not go into the formula to compute the Gaussian random number
          call r4_uni(jsr,rho2)
          x(k,n)=x(k,n)+dsqrt(-2*dlog(rho1))*dcos(2*pi*rho2)*sigma
       enddo

       ! If a particle has a distance to the center of mass that exceeds the radius of the constraining sphere, kill the random walker
       if (.not.cluster_size(x(:,n))) population(n)=0
       
       ! Compute Pb and compare to a random number: assign walkers to be killed, retained, or replicated
       if (population(n)==1) then
          energy=pot_SG_opt(x(:,n),Dim) ! x(:,n) in angstroms, pot_SG_opt in atomic units
          Pb=dexp(-(energy-Vref)*dtau)
          population(n)=floor(Pb)
          call r4_uni(jsr,u)
          if (Pb-population(n)>u) population(n)=population(n)+1 
          if (population(n).ne.0) Vaver=Vaver+energy*population(n)
       endif
       seed(thread)=jsr
    enddo
    !$omp enddo
    !$omp end parallel

    ! Designate which random walkers will be replicated and which ones will be killed 
    ! replicate(i) (i=1,Nreplicate1) contains indices of the walkers to be replicated
    ! kill(i) (i=1,Nreplicate1) contains indices of walkers to be replaced
    Nreplicate=0
    Nkill=0
    do i=1,Nt
       do k=2,population(i) ! Determine the indicies of random walkers that will be replicated
          Nreplicate=Nreplicate+1
          replicate(Nreplicate)=i
       enddo
       if (population(i)==0) then ! Determine the indicies of random walkers that will be killed
          Nkill=Nkill+1
          kill(Nkill)=i
       endif
    enddo
    if (Nreplicate>Nkill) then ! replicate the remaining walkers to Nt+1,Nt+2.... 
       do i=1,Nreplicate-Nkill
          kill(Nkill+i)=Nt+i
       enddo
    else if (Nkill>Nreplicate) then ! replace the remaining walkers to be killed by the walkers at the top (Nt-1,Nt-2,...)
       Nreplicate1=Nreplicate
       do i=Nt,Nt-(Nkill-Nreplicate)+1,-1
          if (population(i)>0) then  ! if (population(i)==0) then skip (no need to move it as it is already at the top)
             Nreplicate1=Nreplicate1+1
             replicate(Nreplicate1)=i
          endif
       enddo
    endif
    Nt=Nt+Nreplicate-Nkill
    if (Nt>Nmax) then
       write(16,*) 'Nt=',Nt,'>',Nmax
       write(16,*) 'Projection Time=',t*dtau
       stop 
    endif
    
    ! Replace the killed configurations with the ones that were replicated
    ! Replace the indicies of the killed walkers with those of the replicated walkers
    ! Parallel Region
    !$omp parallel &
    !$omp shared (Nreplicate,x,replicate,kill,dw_walkers) &                  
    !$omp private (i)
    !$omp do
    do i=1,Nreplicate
       x(:,kill(i))=x(:,replicate(i))
       dw_walkers(kill(i))=dw_walkers(replicate(i))
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Replace the killed configurations with the ones that were replicated
    ! Replace the indicies of the killed walkers with those of the replicated walkers
    if (Nkill>Nreplicate) then ! replace the remaining walkers to be killed by the walkers at the top (Nt-1,Nt-2,...)
       !$omp parallel &
       !$omp shared (Nreplicate,x,replicate,kill,dw_walkers) &                  
       !$omp private (i)
       !$omp do
       do i=Nreplicate+1,Nreplicate1
          x(:,kill(i))=x(:,replicate(i))
          dw_walkers(kill(i))=dw_walkers(replicate(i))
       enddo
       !$omp enddo
       !$omp end parallel
    endif

    ! Compute Vref and Epartaver for the new value of Nt
    Vref=Vaver/Nt-((Nt-N0)/(N0*dtau))
    Epartaver=Epartaver+Vref 
     
    ! Write out the projection time, energy, and number of walkers
    if (mod(t,Kaver)==0) then
       write(1,*) t*dtau,(Epartaver*autoK)/(Nparticles*Kaver),Nt
       Epartaver=0d0
       call flush(1)
    endif
  end subroutine diffusion

  subroutine dw_reset(j)
    Implicit None
    Integer, Intent(Out) :: j
    Integer :: n
    
    ! Prior to computing the observables, reset the indicies and store the configurations and current walker number for descendant weighting
    dw_x=x
    Nt1=Nt
    
    ! Reset the weight index counter, the DMC step at which the weights are computed, and the weights
    j=1
    dw_step=dw_step0
    dw_weight=0
    dw_Nt=0
    
    ! Initialize/reset the indicies for descendant weighting
    do n=1,Nt
       dw_walkers(n)=n
    enddo
    
    ! Compute the tau_dw=0 results for observables when Ndw is larger than 1
    if (Ndw>1.and.dw_step0==0) call descendant_weighting(j)
  end subroutine dw_reset

  subroutine descendant_weighting(j)
    Implicit None
    Integer, Intent(Inout) :: j
    Integer :: i,n
    
    ! Only reset if the value of Ndw is 1 and the tau_dw=0 results should be computed
    if (Ndw==1.and.dw_step0==0) call dw_reset(j)
    
    ! Compute the weights of each walker for the j-th time
    do n=1,Nt
       dw_weight(j,dw_walkers(n))=dw_weight(j,dw_walkers(n))+1
    enddo
    ! Store the current values of Nt to compute the observables later
    dw_Nt(j)=Nt 
    
    ! Determine the number of walkers with non-zero weight
    dw_count=0
    do n=1,Nt1
       if (dw_weight(j,n).ne.0) then 
          dw_count=dw_count+1
       endif
    enddo
    write(30,*) t*dtau,j,dw_count
    
    ! Update the number of time steps at which the next set of weights is computed, and set the new value of the weight index counter
    if (j<Ndw) then 
       dw_step=dw_step+dw_incr  
       j=j+1
       ! Compute the observables once all of the weights are determined  
    else if (j==Ndw) then
       call flush(30)
       if (structure_parms) call histograms()
       if (Ndw==1.and.dw_step0==0) return ! Go back to the main program without resetting the DW indicies, etc.
       call dw_reset(j)
    endif
  end subroutine descendant_weighting

  subroutine histograms()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,j,k,l,n,info,ignore_config(Ndw)
    Double Precision :: pair_dist,CM_dist,Q4,Q6
    
    pairs_hist=0
    CM_hist=0
    Q4_hist=0
    Q6_hist=0
    ignore_config=0

    ! Parallel Region: Compute the pair correlation function, density profile, and the Q4 and Q6 distributions
    !$omp parallel &
    !$omp shared (Nt1,Nparticles,Ndw,dw_x,dw_weight) &
    !$omp shared (hist_pairs_min,hist_Q4_min,hist_Q6_min) &
    !$omp shared (nb_pairs,nb_CM,nb_Q4,nb_Q6) &
    !$omp shared (d_pairs,d_CM,d_Q4,d_Q6) &
    !$omp private (i,j,k,l,n,pair_dist,CM_dist,Q4,Q6,info)

    ! For the first time only, determine if the code is running in parallel and write out the number of threads that will be used 
    if (t==dw_step_max) then
       num_threads=omp_get_num_threads() 
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif

    !$omp do reduction (+:pairs_hist,CM_hist,Q4_hist,Q6_hist,ignore_config)
    do n=1,Nt1
       if (dw_weight(1,n).ne.0) then
          ! Compute the pair distances and determine the correlation function
          do i=2,Nparticles
             do j=1,i-1
                pair_dist=sum((dw_x(3*i-2:3*i,n)-dw_x(3*j-2:3*j,n))**2) 
                l=(dsqrt(pair_dist)-hist_pairs_min)/d_pairs
                if (l>nb_pairs) l=nb_pairs
                if (l<0) l=0
                do k=1,Ndw
                   pairs_hist(l,k)=pairs_hist(l,k)+dw_weight(k,n)
                enddo
             enddo
          enddo

          ! Compute the distances of each particle to the center of mass of each configuration and append to the density profile histogram: center of mass is already at the origin
          do j=1,Nparticles
             CM_dist=sum((dw_x(3*j-2:3*j,n))**2)
             l=dsqrt(CM_dist)/d_CM
             if (l>nb_CM) l=nb_CM
             do k=1,Ndw
                CM_hist(l,k)=CM_hist(l,k)+dw_weight(k,n)
             enddo
          enddo
          
          ! Determine the Q4 and Q6 cluster distributions for the cluster
          call order_par(Nparticles,dble(dw_x(:,n)),Q4,Q6,info)
          
          if (info==0) then
             l=(Q4-hist_Q4_min)/d_Q4
             if (l>nb_Q4) l=nb_Q4
             if (l<0) l=0
             do k=1,Ndw
                Q4_hist(l,k)=Q4_hist(l,k)+dw_weight(k,n)
             enddo
             
             l=(Q6-hist_Q6_min)/d_Q6
             if (l>nb_Q6) l=nb_Q6
             if (l<0) l=0
             do k=1,Ndw
                Q6_hist(l,k)=Q6_hist(l,k)+dw_weight(k,n)
             enddo

          else if (info==1) then ! ignore the configuration
             do k=1,Ndw
                ignore_config(k)=ignore_config(k)+dw_weight(k,n)
             enddo
          endif
       endif
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Write out the pair correlation function, density profile, Q4 and Q6 distributions
    do k=1,Ndw
       do i=0,nb_pairs
          write(8,*) hist_pairs_min+(i+0.5)*d_pairs,pairs_hist(i,k)/(dble(Npairs)*dble(dw_Nt(k))*d_pairs)
       enddo
       do i=0,nb_CM
          write(9,*) (i+0.5)*d_CM,CM_hist(i,k)/(4*pi*dble(dw_Nt(k))*dble((i+0.5)**2)*d_CM**3)
       enddo
       do i=0,nb_Q4 
          write(10,*) hist_Q4_min+(i+0.5)*d_Q4,Q4_hist(i,k)/((dw_Nt(k)-ignore_config(k))*d_Q4)
       enddo
       do i=0,nb_Q6 
          write(11,*) hist_Q6_min+(i+0.5)*d_Q6,Q6_hist(i,k)/((dw_Nt(k)-ignore_config(k))*d_Q6)
       enddo
       write(16,*) t*dtau,k,ignore_config(k)
    enddo
    call flush(8)
    call flush(9)
    call flush(10)
    call flush(11)
    call flush(16)
  end subroutine histograms

  subroutine quench()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: k,n,info,iter,nfunc,ngrad,nquench,count_unphysical,info_order,order,IERR
    Integer(kind=8) :: ngrad_sum
    Character(len=4) :: point_group
    Double Precision :: x_quench(Dim,Nt_quench),energy_opt(Nt_quench),Q4(Nt_quench),Q6(Nt_quench)
    Double Precision :: energy_HA,logW
    
    ! Initialize the necessary variables and arrays
    energy_opt=0d0
    Q4=0d0
    Q6=0d0
    count_config=0
    nquench=0
    ngrad_sum=0
    count_unphysical=0
    
    ! Open the isomer/configuration file to append with new isomers
    open(18,file='lib_isomers'//trim(adjustl(DMC_num))//'.dat',form='unformatted',position='append') 
    
    ! Parallel Region
    !$omp parallel &
    !$omp shared (Dim,Nparticles,Nt_quench,Nt_incr,grad_tol,grad_maxfac,x,energy_thresh) &
    !$omp shared (x_quench,energy_opt,Q4,Q6) &
    !$omp private (n,info,iter,nfunc,ngrad,info_order)

    ! Determine if the quenching subroutine is running in parallel for only the first time it is accessed
    if (t==step_quench) then
       num_threads=omp_get_num_threads()
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif

    ! Quench the configurations in parallel and compute their corresponding Q4 and Q6 distributions
    !$omp do reduction (+:nquench,ngrad_sum,count_unphysical)
    do n=1,Nt_quench
       ! Assign the current config/walker to a different array for quenching
       x_quench(:,n)=dble(x(:,n*Nt_incr))
       
       ! Quench the configuration
       call cg_descent(grad_tol,grad_maxfac*Dim,x_quench(:,n),Dim,pot_SG_opt_cg,grad_SG_opt_cg,&
            info,energy_opt(n),iter,nfunc,ngrad,energy_thresh)
       
       ! If the configuration is physical then calculate its Q4 and Q6
       if (info==0.or.info==2) then
          call order_par(Nparticles,x_quench(:,n),Q4(n),Q6(n),info_order)
          nquench=nquench+1
          ngrad_sum=ngrad_sum+ngrad
       else if (info==1.or.info.ge.3) then ! Assign the configs to "unphysical" and don't bother with the assignment
          count_unphysical=count_unphysical+1
       endif
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Assign the current quenched configurations to one of the isomers already contained in the library based on the classical energy, Q4, and Q6 values
    ! If no match is found, add the configuration to the library and write out the corresponding configuration to an unformatted file
    do n=1,Nt_quench
       call harm_approx_gs(x_quench(:,n),energy_HA,logW,IERR)
       if (IERR==6) then ! the config is a minimum and try to assign the quenched configuration to one of the isomers already in the library
          do k=1,Nisomers
             if (dabs(energy_opt(n)-lib_energy(k))<r_energy_thresh & 
                  .and.dabs(Q4(n)-lib_Q4(k))<r_Q4_thresh.and.dabs(Q6(n)-lib_Q6(k))<r_Q6_thresh) then ! Assign the current config to one of the isomers already in the library
                lib_count_config(k)=lib_count_config(k)+1
                count_config(k)=count_config(k)+1
                goto 100 ! Continue and do not try to add another configuration to the library
             endif
          enddo
          
          Nisomers=Nisomers+1
          if (Nisomers>Nisomers_max) then ! write out the final results to the library and stop
             ! Open the library to replace the already existing isomer fraction information in the library with new data
             open(20,file='lib_info'//trim(adjustl(DMC_num))//'.dat',status='replace') 
             write(20,*) '# Column 1: Isomer Index'
             write(20,*) '# Column 2: Classical Energy (K per Particle)'
             write(20,*) '# Column 3: Quenched Q4'
             write(20,*) '# Column 4: Quenched Q6'
             write(20,*) '# Column 5: HA Ground State Energy (K per Particle)'
             write(20,*) '# Column 6: Normalized Log of the Sqrt Frequency Product'
             write(20,*) '# Column 7: Instantaneous Isomer Multiplicity'
             write(20,*) '# Column 8: Point Group Symbol'
             write(20,*) '# Column 9: Point Group Order'

             ! Write out the projection time, isomer index, quenched energies, Q4, Q6, harmonic approximation parameters, isomer multiplicity, and the point group symbol and order
             do k=1,Nisomers_max
                if (t.le.quench_equilib.and.count_config(k)>0) write(17,*) k,(lib_energy(k)*autoK)/Nparticles,&
                     lib_Q4(k),lib_Q6(k),(lib_energy_HA(k)*autoK)/Nparticles,(lib_logW(k)-lib_logW(1))/2,count_config(k),&
                     lib_pgsym(k),lib_order(k)
                write(20,*) k,(lib_energy(k)*autoK)/Nparticles,&
                     lib_Q4(k),lib_Q6(k),(lib_energy_HA(k)*autoK)/Nparticles,(lib_logW(k)-lib_logW(1))/2,lib_count_config(k),& 
                     lib_pgsym(k),lib_order(k)
             enddo
             close(18)
             close(20)
             
             ! Write the current number of isomers, successful quenching events, average number of gradient evaluations, and the number of unphysical configurations
             write(19,*) t*dtau,nquench,dble(ngrad_sum)/nquench,count_unphysical
             call flush(17)
             call flush(19)
             
             write(16,*) 'Stop: Nisomers=',Nisomers,'>',Nisomers_max ! Nconfigs exceeds Nisomers_max
             stop
          else   
             ! Add the current isomer's parameters to the library and add the coordinats of the isomer to an unformatted file
             call pgsym(Nparticles,x_quench(:,n),order,point_group)
             lib_energy(Nisomers)=energy_opt(n)
             lib_Q4(Nisomers)=Q4(n)
             lib_Q6(Nisomers)=Q6(n)
             lib_energy_HA(Nisomers)=energy_HA
             lib_logW(Nisomers)=logW
             lib_pgsym(Nisomers)=point_group
             lib_order(Nisomers)=order
             lib_count_config(Nisomers)=lib_count_config(Nisomers)+1
             count_config(Nisomers)=count_config(Nisomers)+1
             ! Write out only if the number of isomers is below a threshold value
             if (Nisomers.le.Nisomers_thresh) write(18) real(x_quench(:,n)) 
          endif
       endif
100    continue
    enddo
    
    ! Open the library to replace the already existing isomer fraction information in the library with new data
    open(20,file='lib_info'//trim(adjustl(DMC_num))//'.dat',status='replace') 
    write(20,*) '# Column 1: Isomer Index'
    write(20,*) '# Column 2: Classical Energy (K per Particle)'
    write(20,*) '# Column 3: Quenched Q4'
    write(20,*) '# Column 4: Quenched Q6'
    write(20,*) '# Column 5: HA Ground State Energy (K per Particle)'
    write(20,*) '# Column 6: Normalized Log of the Sqrt Frequency Product'
    write(20,*) '# Column 7: Instantaneous Isomer Multiplicity'
    write(20,*) '# Column 8: Point Group Symbol'
    write(20,*) '# Column 9: Point Group Order'
    
    ! Write out the projection time, isomer index, quenched energies, Q4, Q6, harmonic approximation parameters, isomer multiplicity, and the point group symbol and order
    do k=1,Nisomers
       if (t.le.quench_equilib.and.count_config(k)>0) write(17,*) k,(lib_energy(k)*autoK)/Nparticles,&
            lib_Q4(k),lib_Q6(k),(lib_energy_HA(k)*autoK)/Nparticles,(lib_logW(k)-lib_logW(1))/2,count_config(k),&
            lib_pgsym(k),lib_order(k)
       write(20,*) k,(lib_energy(k)*autoK)/Nparticles,&
            lib_Q4(k),lib_Q6(k),(lib_energy_HA(k)*autoK)/Nparticles,(lib_logW(k)-lib_logW(1))/2,lib_count_config(k),& 
            lib_pgsym(k),lib_order(k)
    enddo
    close(18)
    close(20)
    
    ! Write out the current number of isomers, successful quenching events, average number of gradient evaluations, and the number of unphysical configurations
    write(19,*) t*dtau,nquench,dble(ngrad_sum)/nquench,count_unphysical
    call flush(17)
    call flush(19)
  end subroutine quench

  subroutine restart_DMC()
    Implicit None
    Integer :: i,j,n

    ! Check the parameters for restarting the simulation
    open(50,file='input_restart.dat')
    read(50,*) ! Header
    do j=1,m
       read(50,*)
    enddo
    read(50,*) N0,dtau,tmax,Kaver,lambda,step_quench,Ndw,dw_tau_max,dw_tau_incr,NDMC
    close(50)
    
    ! Initialize the new parameters
    Epartaver=0d0
    N_iter_max=tmax/dtau
    if (structure_parms) then
       dw_step_max=dw_tau_max/dtau
       dw_incr=dw_tau_incr/dtau
       if (Ndw==1.and.dw_incr==0) then
          dw_step0=0 ! Compute the observables using the tau_dw=0 result only
       else
          dw_step0=dw_step_max-(Ndw-1)*dw_incr ! Compute DW results for observables with tau_dw spaced out in regular intervals
       endif
    endif
    
    if (restart=='random_walkers'.and.N0.le.N0_target_max) then ! change the target random walker population for the next DMC run
       ! Add Nt random walkers N0_part-1 times to the pre-existing random walker population
       N0_part=N0/Nt ! N0>Nt as N0 is the new target random walker population for the next round: integer division
       do i=1,N0_part-1
          do n=1,Nt
             x(:,n+i*Nt)=x(:,n)
          enddo
       enddo
       
       ! Add the remaining random walkers to the population to reach the N0 target value for the next DMC run
       do n=1,mod(N0,Nt)
          x(:,n+N0_part*Nt)=x(:,n)
       enddo
       
       ! Reset/adjust the maximum projection time and the target random walker population for the subsequent DMC run
       Nt=N0
    else if (restart=='lambda') then ! change quantum delocalization parameter for the next DMC run
       pseudo_atom_mass=melectron/lambda**2 ! atomic units
       sigma=dsqrt(dtau/pseudo_atom_mass)*bohr ! reset the width of the Gaussian distribution using the new mass
       if (atom_type.ne.'Uk') atom_type='Uk'
    endif
           
    ! Reset the libraries for quenching and isomer assignments
    if (quenching) then
       Nisomers=0
       lib_energy=0d0
       lib_Q4=0d0
       lib_Q6=0d0
       lib_energy_HA=0d0
       lib_logW=0d0
       lib_order=0
       lib_pgsym=''
       lib_count_config=0
    endif
    close(1)
    if (structure_parms) then
       close(30)
       close(8)
       close(9)
       close(10)
       close(11)
    endif
    if (quenching) then
       close(17)
       close(19)
    endif
        
    ! Initialize the files for the next DMC run
    write(file_num,'(i7)') m+1
    file='energy'//trim(adjustl(file_num))//'.dat'
    open(1,file=file) 
    write(1,*) '# Nparticles=',Nparticles,'N0=',N0,'dtau=',dtau
    write(1,*) '# Projection_time=',tmax,'Kaver=',Kaver
    if (atom_type=='H') then
       write(1,*) '# H2mass=',atom_mass('H2')/melectron
       lambda=dsqrt(melectron/atom_mass('H2'))
    else if (atom_type=='D') then
       write(1,*) '# D2mass=',atom_mass('D2')/melectron
       lambda=dsqrt(melectron/atom_mass('D2'))
    else if (atom_type=='He') then
       write(1,*) '# Hemass=',atom_mass('He')/melectron 
       lambda=dsqrt(melectron/atom_mass('He'))
    else if (atom_type=='Ne') then
       write(1,*) '# Nemass=',atom_mass('Ne')/melectron
       lambda=dsqrt(melectron/atom_mass('Ne'))
    else if (atom_type=='Ar') then
       write(1,*) '# Armass=',atom_mass('Ar')/melectron
       lambda=dsqrt(melectron/atom_mass('Ar'))
    else if (atom_type=='Kr') then
       write(1,*) '# Krmass=',atom_mass('Kr')/melectron
       lambda=dsqrt(melectron/atom_mass('Kr'))
    else if (atom_type=='Xe') then
       write(1,*) '# Xemass=',atom_mass('Xe')/melectron
       lambda=dsqrt(melectron/atom_mass('Xe'))
    else if (atom_type=='Uk') then
       write(1,*) '# Pseudo Atom Mass=',pseudo_atom_mass/melectron
    endif
    write(1,*) '# Incr Opt=',Kaver
    write(1,*) '# Number of Threads=',num_threads
    write(1,*) '# Potential=Silvera Goldman'
    write(1,*) '# Quantum Parameter=',lambda 
    write(1,*) '# Initial Vref per Particle (K)=',(Vref*autoK)/Nparticles
    
    if (structure_parms) then
       file='dw_walker_count'//trim(adjustl(file_num))//'.dat'
       open(30,file=file)
       write(30,*) '# Column 1: Projection Time'
       write(30,*) '# Column 2: Weight Index'
       write(30,*) '# Column 3: Number of Walkers with Non-zero Weight'

       file='pairs_dist'//trim(adjustl(file_num))//'.dat'
       open(8,file=file)  
       write(8,*) '# Histogram Particle Pairs=',hist_pairs_min,hist_pairs_max
       write(8,*) '# Histogram Particle Pairs=',d_pairs
       write(8,*) '# Column 1=Particle Pair Distances'
       write(8,*) '# Column 2=Normalized Bin Count'

       file='CM_dist'//trim(adjustl(file_num))//'.dat'
       open(9,file=file) 
       write(9,*) '# Histogram CM Distances=',hist_CM_max
       write(9,*) '# Histogram CM Distances=',d_CM
       write(9,*) '# Column 1=CM Distances'
       write(9,*) '# Column 2=Normalized Bin Count'
     
       file='Q4_dist'//trim(adjustl(file_num))//'.dat'
       open(10,file=file)
       write(10,*) '# Histogram Q4 Cluster Distance=',hist_Q4_min,hist_Q4_max
       write(10,*) '# Histogram Q4 Cluster Distance=',d_Q4
       write(10,*) '# Column 1=Q4 Cluster Distances'
       write(10,*) '# Column 2=Normalized Bin Count'
       
       file='Q6_dist'//trim(adjustl(file_num))//'.dat'
       open(11,file=file)
       write(11,*) '# Histogram Q6 Cluster Distance=',hist_Q6_min,hist_Q6_max
       write(11,*) '# Histogram Q6 Cluster Distance=',d_Q6
       write(11,*) '# Column 1=Q6 Cluster Distances'
       write(11,*) '# Column 2=Normalized Bin Count'
    endif

    if (quenching) then
       file='populations'//trim(adjustl(file_num))//'.dat'
       open(17,file=file)
       write(17,*) '# Column 1: Isomer Index'
       write(17,*) '# Column 2: Classical Energy (K per Particle)'
       write(17,*) '# Column 3: Quenched Q4'
       write(17,*) '# Column 4: Quenched Q6'
       write(17,*) '# Column 5: HA Ground State Energy (K per Particle)'
       write(17,*) '# Column 6: Normalized Log of the Sqrt Frequency Product'
       write(17,*) '# Column 7: Instantaneous Isomer Multiplicity'
       write(17,*) '# Column 8: Point Group Symbol'
       write(17,*) '# Column 9: Point Group Order'
       
       file='grad_evals'//trim(adjustl(file_num))//'.dat'
       open(19,file=file)
       write(19,*) '# Column 1: Projection Time'
       write(19,*) '# Column 2: Number of Configs Quenched Successfully'
       write(19,*) '# Column 3: Average Number of Gradient Evaluations'
       write(19,*) '# Column 4: Number of Unphysical Quenched Configs'
    endif
    call flush(1)
    if (structure_parms) then
       call flush(30)
       call flush(8)
       call flush(9)
       call flush(10)
       call flush(11)
    endif
    if (quenching) then
       call flush(17)
       call flush(19)
    endif
  end subroutine restart_DMC
end module dmc_module

Program dmc_true_gs
  
  use dmc_module
  use omp_lib
  
  Implicit None
  Integer :: i,j,n
  Character(len=2) :: dummy

  ! All distances are in angstroms
  open(2,file='input.dat')
  read(2,*) Nparticles
  read(2,*) N0 ! Target random walker population for preliminary equilibration or target population for single DMC run
  read(2,*) dtau
  read(2,*) tmax ! Maximum projection time for preliminary equilibration or maximum projection time for single DMC run
  read(2,*) Kaver
  read(2,*) Vref ! Kelvin per particle
  read(2,*) Ndw,dw_tau_max,dw_tau_incr,quench_equilib
  read(2,*) hist_pairs_min,hist_pairs_max,nb_pairs
  read(2,*) hist_CM_max,nb_CM
  read(2,*) hist_Q4_min,hist_Q4_max,nb_Q4
  read(2,*) hist_Q6_min,hist_Q6_max,nb_Q6
  read(2,*) rc_sphere
  read(2,*) coord_config0
  read(2,*) atom_type,lambda ! quantum delocalization parameter in amu^(-1/2)
  read(2,*) start_configs,restart ! Strings that determine where the initial configs come from, and how the DMC simulation should be restarted, respectively
  read(2,*) structure_parms,quenching,write_configs,Nt_incr
  read(2,*) NMMC,incr_MMC,rc_MMC,Temp ! rc_MMC in reduced units of sigma and Temp in reduced units of epsilon
  read(2,*) NDMC,Vref_parm,N0_target_max
  read(2,*) Nisomers_max,Nisomers_thresh,Nisomers_configs,step_quench
  read(2,*) r_energy_thresh,r_Q4_thresh,r_Q6_thresh ! r_energy_thresh in Kelvin per particle
  read(2,*) grad_tol,grad_maxfac,energy_thresh ! Kelvin per particle
  read(2,*) num_threads
  close(2)

  ! Initialize the parameters for dimensionality, the maximum numberof iterations, the number of particle-particle pairs,the well depth of the LJ fit for the LJ potential in atomic units, the mass of the pseudo atom in atomic units, the mass decrementation in atomic units, and the quantum parameter lambda
  Dim=3*Nparticles
  N_iter_max=tmax/dtau
  Npairs=(Nparticles*(Nparticles-1))/2
  rc_sphere=rc_sphere*sigma_lj
  pseudo_atom_mass=melectron/lambda**2 ! atomic units
  if (restart=='random_walkers') then ! Set the largest target random walker population: set Nmax to just larger than the largest N0 target value
     Nmax=N0_target_max*1.1
  else ! Nmax is set to just larger than N0
     Nmax=N0*1.1
  endif
  
  ! Define the distance parameter for the histograms and allocate memory to each array
  if (structure_parms) then
     d_pairs=(hist_pairs_max-hist_pairs_min)/(nb_pairs+1)
     d_CM=hist_CM_max/(nb_CM+1)
     d_Q4=(hist_Q4_max-hist_Q4_min)/(nb_Q4+1)
     d_Q6=(hist_Q6_max-hist_Q6_min)/(nb_Q6+1)

     ! Initialize the first value of dw_step and the increment for computing the weights
     dw_step_max=dw_tau_max/dtau
     dw_incr=dw_tau_incr/dtau
     if (Ndw==1.and.dw_incr==0) then
        dw_step0=0 ! Compute the observables using the tau_dw=0 result only
     else
        dw_step0=dw_step_max-(Ndw-1)*dw_incr ! Compute DW results for observables with tau_dw spaced out in regular intervals
     endif
  endif
  
  allocate(seed(0:num_threads-1),x(Dim,Nmax),dw_walkers(Nmax))
  if (structure_parms) then
      allocate(dw_weight(Ndw,Nmax),dw_x(Dim,Nmax),dw_Nt(Ndw))
     allocate(pairs_hist(0:nb_pairs,Ndw),CM_hist(0:nb_CM,Ndw),Q4_hist(0:nb_Q4,Ndw),Q6_hist(0:nb_Q6,Ndw))
  endif
  if (quenching) then
     allocate(count_config(Nisomers_max),lib_count_config(Nisomers_max))
     allocate(lib_energy(Nisomers_max),lib_Q4(Nisomers_max),lib_Q6(Nisomers_max))
     allocate(lib_energy_HA(Nisomers_max),lib_logW(Nisomers_max),lib_pgsym(Nisomers_max),lib_order(Nisomers_max))
  endif

  ! Initialize the libraries for quenching and isomer assignments
  if (quenching) then 
     energy_thresh=(energy_thresh*Nparticles)/autoK
     r_energy_thresh=(r_energy_thresh*Nparticles)/autoK
     Nisomers=0
     lib_energy=0d0
     lib_Q4=0d0
     lib_Q6=0d0
     lib_energy_HA=0d0
     lib_logW=0d0
     lib_order=0
     lib_pgsym=''
     lib_count_config=0
  endif
 
  if (NDMC>0) then
     write(file_num,'(i7)') 0
     file='energy'//trim(adjustl(file_num))//'.dat'
     open(1,file=file) 
  else
     open(1,file='energy.dat')
  endif
  write(1,*) '# Nparticles=',Nparticles,'N0=',N0,'dtau=',dtau
  write(1,*) '# Projection_time=',tmax,'Kaver=',Kaver
  if (atom_type=='H') then
     write(1,*) '# H2mass=',atom_mass('H2')/melectron
     lambda=dsqrt(melectron/atom_mass('H2'))
  else if (atom_type=='D') then
     write(1,*) '# D2mass=',atom_mass('D2')/melectron
     lambda=dsqrt(melectron/atom_mass('D2'))
  else if (atom_type=='He') then
     write(1,*) '# Hemass=',atom_mass('He')/melectron 
     lambda=dsqrt(melectron/atom_mass('He'))
  else if (atom_type=='Ne') then
     write(1,*) '# Nemass=',atom_mass('Ne')/melectron
     lambda=dsqrt(melectron/atom_mass('Ne'))
  else if (atom_type=='Ar') then
     write(1,*) '# Armass=',atom_mass('Ar')/melectron
     lambda=dsqrt(melectron/atom_mass('Ar'))
  else if (atom_type=='Kr') then
     write(1,*) '# Krmass=',atom_mass('Kr')/melectron
     lambda=dsqrt(melectron/atom_mass('Kr'))
  else if (atom_type=='Xe') then
     write(1,*) '# Xemass=',atom_mass('Xe')/melectron
     lambda=dsqrt(melectron/atom_mass('Xe'))
  else if (atom_type=='Uk') then
     write(1,*) '# Pseudo Atom Mass=',pseudo_atom_mass/melectron
  endif
  write(1,*) '# Incr Opt=',Kaver
  write(1,*) '# Number of Threads=',num_threads
  write(1,*) '# Potential=Silvera Goldman'
  write(1,*) '# Quantum Parameter=',lambda

  if (structure_parms) then
     if (NDMC>0) then
        file='dw_walker_count'//trim(adjustl(file_num))//'.dat'
        open(30,file=file)
     else
        open(30,file='dw_walker_count.dat')
     endif
     write(30,*) '# Column 1: Projection Time'
     write(30,*) '# Column 2: Weight Index'
     write(30,*) '# Column 3: Number of Walkers with Non-zero Weight'
     
     if (NDMC>0) then
        file='pairs_dist'//trim(adjustl(file_num))//'.dat'
        open(8,file=file) 
     else
        open(8,file='pairs_dist.dat') 
     endif
     write(8,*) '# Histogram H2-H2 Pairs=',hist_pairs_min,hist_pairs_max
     write(8,*) '# Histogram H2-H2 Pairs=',d_pairs
     write(8,*) '# Column 1=H2-H2 Pair Distances'
     write(8,*) '# Column 2=Normalized Bin Count'
     
     if (NDMC>0) then
        file='CM_dist'//trim(adjustl(file_num))//'.dat'
        open(9,file=file) 
     else
        open(9,file='CM_dist.dat')
     endif
     write(9,*) '# Histogram CM Distances=',hist_CM_max
     write(9,*) '# Histogram CM Distances=',d_CM
     write(9,*) '# Column 1=CM Distances'
     write(9,*) '# Column 2=Normalized Bin Count'
     
     if (NDMC>0) then
        file='Q4_dist'//trim(adjustl(file_num))//'.dat'
        open(10,file=file)
     else
        open(10,file='Q4_dist.dat')
     endif
     write(10,*) '# Histogram Q4 Cluster Distance=',hist_Q4_min,hist_Q4_max
     write(10,*) '# Histogram Q4 Cluster Distance=',d_Q4
     write(10,*) '# Column 1=Q4 Cluster Distances'
     write(10,*) '# Column 2=Normalized Bin Count'
     
     if (NDMC>0) then
        file='Q6_dist'//trim(adjustl(file_num))//'.dat'
        open(11,file=file)
     else
        open(11,file='Q6_dist.dat')
     endif
     write(11,*) '# Histogram Q6 Cluster Distance=',hist_Q6_min,hist_Q6_max
     write(11,*) '# Histogram Q6 Cluster Distance=',d_Q6
     write(11,*) '# Column 1=Q6 Cluster Distances'
     write(11,*) '# Column 2=Normalized Bin Count'
  endif

  if (quenching) then
     if (NDMC>0) then
        file='populations'//trim(adjustl(file_num))//'.dat'
        open(17,file=file)
     else
        open(17,file='populations.dat')
     endif
     write(17,*) '# Column 1: Isomer Index'
     write(17,*) '# Column 2: Classical Energy (K per Particle)'
     write(17,*) '# Column 3: Quenched Q4'
     write(17,*) '# Column 4: Quenched Q6'
     write(17,*) '# Column 5: HA Ground State Energy (K per Particle)'
     write(17,*) '# Column 6: Normalized Log of the Sqrt Frequency Product'
     write(17,*) '# Column 7: Instantaneous Isomer Multiplicity'
     write(17,*) '# Column 8: Point Group Symbol'
     write(17,*) '# Column 9: Point Group Order'

     if (NDMC>0) then
        file='grad_evals'//trim(adjustl(file_num))//'.dat'
        open(19,file=file)
     else
        open(19,file='grad_evals.dat')
     endif
     write(19,*) '# Column 1: Projection Time'
     write(19,*) '# Column 2: Number of Configs Quenched Successfully'
     write(19,*) '# Column 3: Average Number of Gradient Evaluations'
     write(19,*) '# Column 4: Number of Unphysical Quenched Configs'
  endif
  
  open(16,file='parallel_test.dat')
  write(16,*) '# Column 1: Parallel On or Off'
  write(16,*) '# Column 2: Number of Threads in Use'
  write(16,*) '# Projection Time and Number of Ignored Configs'
  call flush(16)

  if (structure_parms) then
     call flush(30)
     call flush(8)
     call flush(9)
     call flush(10)
     call flush(11)
  endif
  if (quenching) then
     call flush(17)
     call flush(19)
  endif
  
  if (start_configs=='MMC_configs') then ! use Metropolis Monte Carlo to generate randomized configurations
     allocate(U(Nparticles,Nparticles),U_move(Nparticles))
     call read_config(coord_config0,x(:,1))
     call U_array(x(:,1),Dim)
     rc_MMC=rc_MMC*sigma_lj
     Temp=Temp*epsilon_lj
     do n=2,N0
        x(:,n)=x(:,n-1)
        call Metropolis_MC(NMMC,x(:,n),1/Temp,rc_MMC)
        if (mod(n,incr_MMC)==0) then
           write(coord_MMC,'(i7)') n/incr_MMC
           coord_MMC='MMC_configs/'//trim(adjustl(coord_MMC))//'.xyz'
           open(31,file=coord_MMC)
           write(31,*) Nparticles
           write(31,*) (pot_SG_opt(x(:,n),Dim)*autoK)/Nparticles,'|K'
           if (atom_type=='Uk') then
              do i=1,Nparticles
                 write(31,*) 'H',x(3*i-2:3*i,n)
              enddo
           else
              do i=1,Nparticles
                 write(31,*) atom_type,x(3*i-2:3*i,n)
              enddo
           endif
           close(31)
        endif
     enddo
     deallocate(U,U_move)
  else if (start_configs=='single_config') then ! place all walkers in the basin of attraction corresponding to the initial configuration
     call read_config(coord_config0,x(:,1))
     do n=2,N0
        x(:,n)=x(:,1)
     enddo
  else if (start_configs=='file_configs_unquenched') then ! all walkers are initialized from a file containing unquenched DMC configs
     open(40,file='coord_configs.dat',form='unformatted')
     do n=1,N0
        do i=1,Nparticles
           read(40) dummy,x(3*i-2:3*i,n)
        enddo
        read(40)
     enddo
     close(40)
  else if (start_configs=='file_configs_quenched') then ! all walkers are initialized from a file containing quenched DMC configs
     open(40,file='coord_configs.dat',form='unformatted')
     N0_part=N0/Nisomers_configs ! integer division
     do i=1,Nisomers_configs
        read(40) x(:,1+(i-1)*N0_part)
        
        ! Place an equal share of random walkers at the position of the reference configuration in config space
        do n=2,N0_part
           x(:,n+(i-1)*N0_part)=x(:,1+(i-1)*N0_part)
        enddo
     enddo
     close(40)
     
     ! If N0 and Nlib_configs do not divide evenly, place the remaining walkers at the position of each refenerce config until no walkers remain
     do i=1,mod(N0,Nisomers_configs)            
        x(:,i+N0_part*Nisomers_configs)=x(:,1+(i-1)*N0_part)
     enddo
  endif
  
  ! Calculate the width (sigma) of the Gaussian distribution for the DMC diffusion moves
  ! Multiply by bohr to convert to sigma to angstroms
  if (atom_type=='H') then
     sigma=dsqrt(dtau/atom_mass('H2'))*bohr
  else if (atom_type=='D') then
     sigma=dsqrt(dtau/atom_mass('D2'))*bohr
  else if (atom_type=='He') then
     sigma=dsqrt(dtau/atom_mass('He'))*bohr
  else if (atom_type=='Ne') then
     sigma=dsqrt(dtau/atom_mass('Ne'))*bohr
  else if (atom_type=='Ar') then
     sigma=dsqrt(dtau/atom_mass('Ar'))*bohr
  else if (atom_type=='Kr') then
     sigma=dsqrt(dtau/atom_mass('Kr'))*bohr
  else if (atom_type=='Xe') then
     sigma=dsqrt(dtau/atom_mass('Xe'))*bohr
  else if (atom_type=='Uk') then
     sigma=dsqrt(dtau/pseudo_atom_mass)*bohr
  else
     write(16,*) 'unexpected atom type=',atom_type
     stop
  endif
  write(16,*) 'Atom_type=',atom_type

  ! Initialize the constants in the order_par.f subroutine
  if (structure_parms.or.quenching) call SHINIT()

  ! Initialize Vref, Nt, and Epartaver
  if (Vref_parm) then
     write(1,*) '# Initial Vref per Particle (K)=',Vref 
     Vref=(Vref*Nparticles)/autoK
  else
     Vref=pot_SG_opt(x(:,1),Dim)
     write(1,*) '# Initial Vref per Particle (K)=',(Vref*autoK)/Nparticles
  endif
  call flush(1)
  Nt=N0
  Epartaver=0d0

  ! Initialize the indicies for descendant weighting
  do n=1,Nt
     dw_walkers(n)=n
  enddo
  
  ! Set the number of threads for all of the parallel regions embedded within the subroutines
  ! If not set, the maximum number of threads on the machine is used by default
  call omp_set_num_threads(num_threads)
  
  ! Outer loop: Number of DMC simulations that will be run: 0th run is for preliminary equilibration
  ! Inner loop: Number of time steps for each DMC run
  do m=0,NDMC
     write(DMC_num,'(i7)') m
     if (structure_parms) call dw_reset(j)
     do t=1,N_iter_max
        ! Perform the diffusion steps and branching
        call diffusion()

        ! Compute the weight of each walker at different DMC steps specified by dw_step, and compute the observables by the dw method
        if (structure_parms) then
           if (mod(t,dw_step_max)==dw_step.or.mod(t,dw_step_max)==0) call descendant_weighting(j)
        endif

        ! Quench a fraction of the configs without dw
        if (quenching.and.mod(t,step_quench)==0) then
           Nt_quench=Nt/Nt_incr
           call quench()
        endif
        
        ! Write out the configurations at the last iteration of each DMC run if desired 
        if (write_configs.and.m==NDMC.and.t==N_iter_max) then
           write(file_num,'(i7)') m
           file='coord_configs'//trim(adjustl(file_num))//'.dat'
           open(14,file=file,form='unformatted')
           do n=1,Nt
              if (mod(n,Nt_incr)==0) write(14) x(:,n) ! write out a given fraction of the configurations at random
           enddo
           close(14)
        endif
        
        ! After N_max_iter DMC steps have completed, increase the random walker population by a specified amount and reset the time step loop for the next round
        if (m<NDMC.and.t==N_iter_max) call restart_DMC()
     enddo
  enddo
end Program dmc_true_gs

