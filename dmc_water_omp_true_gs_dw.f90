!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS IN INPUT FILE/GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Natoms: Number of atoms in system
! N0: Initial number of random walkers
! dtau: Time step
! tmax: Maximum projection time
! Kaver: Controls projection times at which the energies are written out and the configurations are stored for restarting the simulation
! Ndw: Number of times walker weights are computed for descendant weighting
! dw_tau_max: Parameter that controls the number of DMC steps that elapse before the weights for descendant weighting are reset
! dw_tau_incr: Increment that controls the number of DMC steps that elapse before the weights for descendant weighting are computed again
! hist_dist_min,histangle_min,hist_pairs_OO_min: Lower bound for histogram data on x axis
! hist_dist_max,hist_angle_max,hist_pairs_OO_max: Upper bound for histogram data on x axis
! nb_dist,nb_angle,nb_pairs_OO: Number of bins used to generate histgrams
! histograms: Logical parameter that controls if histograms are computed
! quenching: Logical parameter that controls if the quenching procedure is performed
! inertia: Logical parameter that determines if the moments of inertia are computed if no quenching is done
! write_configs: Logical parameter that controls if the configurations from DMC at a given time step are written out to a file
! step_configs: Time step at which the configurations are written out
! Nt_incr: Increment of random walkers at which the configurations are selected to write out
! isotope: Default: Assigns 'H' masses to all 'H' atoms; 'D_only' if xyz file has 'H' labels but simulation is performed with 'D' masses; 'H_and_D' every other 'H' changed to 'D'; 'one_D' if only one atom is changed to 'D' mass; 'one_H' if only one atom stays 'H' while all others are changed to 'D'
! potential: 'tip4p', 'ttm3', or 'mbpol' characters that specify which potential energy surface (PES) is used in the simulation
! Vref_parm: Specifies the basin of attraction from which Vref is intially obtained
! r_energy_thresh: Energy threshold for assignment of configurations to library reference configurations
! r_pairs_thresh: O-O pair distance threshold for assignment of configurations to library reference configurations
! r_cos_thresh: Cosine of bond angles threshold for assignment of configurations to library reference configurations 
! grad_tol: Controls the accuracy of the quenching procedure
! grad_maxfac: Controls the total number of gradient evaulations
! energy_thresh: Lower bound on energy of a quenched structure below which the configuration is ignored
! restart_parm: Logical variable that determines if configurations are written out for restarting the program
! num_threads: Total number of threads that can be used in the parallel regions
! file_type: Character strings "new" or "restart" determine if the simulation is new or needs to be restarted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS IN INPUT FILE/GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dim: Dimensionality of the configuration space
! NH2O: Number of oxygen atoms in system
! Npairs: Number of pairs of oxygen atoms
! Nmax: Maximum number of random walkers that can be present in the system
! N0_part: Number of random walkers placed in each basin of attraction
! Nt: Instantaneous number of random walkers
! Nt1: Stored number of random walkers for use in computing observables from descendant weighting
! Nisomers: Number of reference configurations in library
! N_iter_max: Maximum number of DMC steps calculated from tmax and dtau
! t: Current iteration number in loop 1,N_iter_max (outer loop of main program)
! t1: Main loop index at which to start or restart the DMC simulation
! dw_step: Specifies the current time step at which the weights are computed for descendant weighting
! dw_count: Counts the number of non-zero walkers that were used to compute the observeables
! Vref: Instantaneous DMC energy: estimator of the true quantum mechanical ground state energy
! Epartaver: Partial average of the DMC energy over Kaver steps
! sum_mass: Sum of masses over each atom in the system
! seed(num_threads): Array that stores the random number generator seeds for each individual thread
! mass(Natoms): Array that stores the mass of each atom
! sigma(Dim): Array of widths from Gaussian distribution for DMC moves
! x(Dim,Nmax): Array that stores the instantaneous position in configuration space for each configuration/random walker
! file_num: Number that designates each file stored in the library_opt_configs directory
! coord_config: String that gives the name of the reference configuration to be read
! atom_type(Natoms): Character array of chemical symbols for atoms appearing in an xyz coordinate file
! dw_walkers(Nmax): Array of walker indicies for descendant weighting
! dw_weight(Ndw,Nmax): Array of weights for each walker from descendant weighting used to compute observables
! dw_Nt(Ndw): Array that stores the current walker numbers for use in calculating observables
! dw_x(Dim,Nmax): Array of confiugrations for descendant weighting used to compute observables
! d_dist,d_angle,d_pairs_OO: Length of each histogram bin
! OH_hist(0:nb_dist,Ndw),HOH_hist(0:nb_angle,Ndw),pairs_OO_hist(0:nb_pairs_OO,Ndw): Arrays containing data for each histogram
! count_config(Nisomers+2,Ndw): Array of populations for the reference configurations: includes the population of the configs assigned to "others" and "unphysical"
! min_energy(Nisomers): Array that stores the classical energies of each reference configuration
! lib_config_pd(Npairs,Nisomers): Array to store the pair distances of each reference configuration
! lib_config_cos(Npairs,Nisomers): Array to store the value of the cosines of each reference configuration
! MI_eigenvals_aver(3,Nisomers,Ndw): Running sum of the eigenvalues for the moment of inertia tensor; used to compute the average moments of inertia for each configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GLOBAL VARIABLES IN DMC MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
! energy: Potential energy of a configuration from one of the active PESs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN DIFFUSION SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN HISTOGRAM SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel: Logical variable that is true if the program is running in parallel and false if not
! i,j,k,l: Loop indicies
! n: Index for loop over random walkers
! OH_dist(2): Array of the two OH bond lengths for each water molecule
! HH_dist: Distance between two hydrogen atoms for each water molecule
! HOH_angle: Angle between the atoms in each water moelcule
! pairs_OO_dist: Pair distances between oxygen atoms computed for the pair distance histogram
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN HISTOGRAM SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN QUENCH SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel: Logical variable that is true if the program is running in parallel and false if not
! i,j,k: Indices for loops during quenching
! n: Index for loop over random walkers
! info: Number output by the cg_descent program: 0 or 2 if the configurations quenched correctly; other numbers mean an error occurred during quenching
! iter: Number of iterations performed by the cg_descent program
! nfunc: Number of potential energy function evaluations performed by the cg_descent program
! ngrad: Number of gradient evaluations performed by the cg_descent program
! ngrad_sum: Keeps track of number of gradient evaluations during crude quenching
! nquench: Counter that determines the number of successful quenching events for a given set of configurations
! energy_opt: Energy of a quenched configuration
! x_quench(Dim): Array that stores a configuration designated for quenching
! pair_dist(Npairs): Array of pair distances for an instantaneous configuration
! cosines(Npairs): Value of the cosine for the OH bond angles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN QUENCH SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! i,j,n: Loop indicies
! MI_eigenvals(3): Array the contains the eigenvalues of the moment of inertia tensor for the reference isomers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLES IN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARALLEL REGIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random number generator found on John Burkardt's website.
! URL: https://people.sc.fsu.edu/~jburkardt/f_src/ziggurat_openmp/ziggurat_openmp.f90
! Parallel Region 1: Computation of energies of each random walker and determination of which random walkers will be replicated or killed in the diffusion subroutine
! Parallel Region 2: Replacement of random walkers to be killed with those to be replicated in the diffusion subroutine
! Parallel Region 3: If the number of random walkers to be killed is larger than those to be replicated, then replace the killed ones by those on the top of the replicate array in the diffusion subroutine
! Parallel Region 4: Compute the histograms for the bond lengths, bond angles, and pair correlation functions in the histograms subroutine
! Parallel Region 5: Assign instantaneous configurations to one of the configurations in the library after quenching in the quench subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARALLEL REGIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Module containing all of the global variables and subroutines for the DMC true ground state code
module dmc_module
  Implicit None
  Double Precision, Parameter :: bohr=0.52917721092
  Double Precision, Parameter :: autokcalmol=627.51 
  Double Precision, Parameter :: autokJmol=2625.5002
  Double Precision, Parameter :: autoMHz=6579683920.711
  Double Precision, Parameter :: melectron=1822.88839
  Double Precision, Parameter :: deg=180/dacos(-1d0)
  Double Precision, Parameter :: pi=dacos(-1d0)
  Double Precision, Parameter :: days=86400,hours=3600    
  Double Precision, Parameter :: Hmass=1.00782503223*melectron
  Double Precision, Parameter :: Dmass=2.01410177812*melectron
  Double Precision, Parameter :: Omass=15.99491461957*melectron
  
  ! Global variables for main program and/or DMC code
  Integer :: Natoms,Dim,NH2O,Npairs,Nmax,N0,N0_part,Nt,Nt1,Nisomers
  Integer :: N_iter_max,Kaver,t,t1,Vref_parm,num_threads
  Integer :: Ndw,dw_step0,dw_step_max,dw_step,dw_incr,dw_count,step_configs,Nt_incr
  Double Precision :: dtau,tmax,dw_tau_max,dw_tau_incr,Vref,Epartaver,sum_mass
  Integer, Allocatable :: seed(:)
  Double Precision, Allocatable :: mass(:),sigma(:),x(:,:)
  Character(len=5) :: potential
  Character(len=7) :: isotope
  Character(len=7) :: file_type
  Character(len=50) :: file_num,coord_config
  Character(len=2), Allocatable :: atom_type(:)
  Logical :: histograms,quenching,inertia,restart_parm,write_configs
  
  ! Global variables for descendant weighting
  Integer, Allocatable :: dw_walkers(:),dw_weight(:,:),dw_Nt(:)
  Double Precision, Allocatable :: dw_x(:,:)

  ! Global variables for histograms
  Integer :: nb_dist,nb_angle,nb_pairs_OO
  Double Precision :: hist_dist_min,hist_angle_min,hist_pairs_OO_min
  Double Precision :: hist_angle_max,hist_dist_max,hist_pairs_OO_max
  Double Precision :: d_dist,d_angle,d_pairs_OO
  Integer, Allocatable :: OH_hist(:,:),HOH_hist(:,:),pairs_OO_hist(:,:) 
  
  ! Global variables for quenching
  Integer :: grad_maxfac
  Double Precision :: grad_tol,r_energy_thresh,r_cos_thresh,r_pairs_thresh,energy_thresh
  Integer, Allocatable :: count_config(:,:)
  Double Precision, Allocatable :: min_energy(:),lib_config_pd(:,:),lib_config_cos(:,:)
  Double Precision, Allocatable :: MI_eigenvals_aver(:,:,:)
 
contains
  
  subroutine pot_energy(energy,x,Dim)
    use iso_c_binding
    Implicit None
    Integer, Intent(In) :: Dim 
    Double Precision, Intent(In) :: x(Dim) ! Bring in the coordinates (angstroms)
    Double Precision, Intent(Out) :: energy
    
    ! Call the Fortran 90/C++ interface for the TIP4P potential with energy only
    interface
       subroutine qtip4pf_energy(nH2O, q, energy) bind(C)
         use iso_c_binding
         integer(c_size_t), value :: nH2O
         real(c_double), intent(in) :: q(*)
         real(c_double), intent(out) :: energy
       end subroutine qtip4pf_energy
    end interface
    
    if (potential=='tip4p') then
       call qtip4pf_energy(INT8(Dim/9),x,energy) ! x in angstrom, energy in kcal/mol
    else if (potential=='ttm3') then 
       call calcttm3f(INT8(Dim/9),energy,x) ! x in angstrom, energy in kcal/mol
    else if (potential=='mbpol') then
       call calcpot(INT8(Dim/9),energy,x) ! x in angstrom, energy in kcal/mol
    else
       stop 'unrecognized potential'
    endif
  end subroutine pot_energy

  subroutine gradient(grad,x,Dim)
    use iso_c_binding
    Implicit None
    Integer, Intent(In) :: Dim                 
    Double Precision, Intent(In) :: x(Dim) ! Bring in the coordinates (angstroms)
    Double Precision, Intent(Out) :: grad(Dim)
    Double Precision :: energy
    
    ! Call the Fortran 90/C++ interface for the TIP4P potential with energy and gradients
    interface
       subroutine qtip4pf_gradient(nH2O, q, energy, grad) bind(C)
         use iso_c_binding
         Integer(c_size_t), Value :: nH2O
         Real(c_double), Intent(in) :: q(*)
         Real(c_double), Intent(out) :: energy, grad(*)
       end subroutine qtip4pf_gradient
    end interface
    
    if (potential=='tip4p') then
       call qtip4pf_gradient(INT8(Dim/9),x,energy,grad) ! x in angstrom, energy in kcal/mol, grad in kcal/mol/angstrom
    else if (potential=='ttm3') then 
       call calcttm3fg(INT8(Dim/9),energy,x,grad) ! x in angstrom, energy in kcal/mol, grad in kcal/mol/angstrom
    else if (potential=='mbpol') then
       call calcpotg(INT8(Dim/9),energy,x,grad) ! x in angstrom, energy in kcal/mol, grad in kcal/mol/angstrom
    else
       stop 'unrecognized potential'
    endif
  end subroutine gradient

  function atom_mass(atom)
    Implicit None
    Double Precision :: atom_mass
    Character(len=2), Intent(In) :: atom
    
    if (atom=='O') then
       atom_mass=Omass
    else if (atom=='H') then
       atom_mass=Hmass
    else if (atom=='D') then
       atom_mass=Dmass
    else
       write(16,*) 'atom ', atom, ' is not recognized'
       stop 
    endif
    return
  end function atom_mass

  subroutine read_config(xyz_file,x,energy)
    Implicit None
    Character(len=30), Intent(In) :: xyz_file 
    Integer :: i
    Double Precision, Intent(Out) :: x(Dim),energy

    open(3,file=xyz_file)
    read(3,*) 
    read(3,*) energy
    x=0d0
    do i=1,Natoms
       read(3,*) atom_type(i),x(3*i-2:3*i)
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

  function distance(N,a,b)
    Implicit None
    Integer :: N
    Double Precision :: a(N),b(N),distance
    distance=sum((a-b)**2)
    distance=sqrt(distance)
  end function distance

  subroutine OO_distances(x,pair_dist)
    Implicit None
    Integer :: i,j,l
    Double Precision, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: pair_dist(Npairs)
    
    ! Given a configuration of NH2O water molecules, compute O-O pair distances.
    pair_dist=0d0
    l=0
    do i=2,NH2O
       do j=1,i-1
          l=l+1
          pair_dist(l)=sum((x(9*i-8:9*i-6)-x(9*j-8:9*j-6))**2) 
       enddo
    enddo
    pair_dist=dsqrt(pair_dist)
    if (Npairs>1) call hpsort(Npairs,pair_dist)
  end subroutine OO_distances

  subroutine OH_bonds(x,cosines)
    Implicit None
    Integer :: i,j,l
    Double Precision, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: cosines(Npairs)
    Double Precision :: vec1(3),vec2(3),vec(3),vecn(3,NH2O),len,len1,len2
  
    ! Determine the value of the cosine for each pair
    l=0
    do i=1,NH2O
       len1=1000
       len2=2000
       do j=2,Natoms
          if (mod(j,3).ne.1) then
             vec=x(3*j-2:3*j)-x(9*i-8:9*i-6)
             len=dot_product(vec,vec)
             if (len<len1) then
                len2=len1
                len1=len
                vec2=vec1
                vec1=vec
             else if (len<len2) then
                len2=len
                vec2=vec
             endif
          endif
       enddo
       vec=vec1+vec2
       vecn(:,i)=vec/sqrt(dot_product(vec,vec))
       do j=1,i-1
          l=l+1
          cosines(l)=dot_product(vecn(:,i),vecn(:,j))
       enddo
    enddo
    if (Npairs>1) call hpsort(Npairs,cosines)
  end subroutine OH_bonds

  subroutine CM(Natoms,x)
    Implicit None
    Integer, Intent(In) :: Natoms
    Double Precision, Intent(Inout) :: x(3,Natoms)
    Integer :: i
    Double Precision :: vec(3)
    vec=0
    do i=1,Natoms
       vec(:)=vec(:)+mass(i)*x(:,i)
    enddo
    vec=vec/sum_mass
    do i=1,Natoms
       x(:,i)=x(:,i)-vec(:)
    enddo
  end subroutine CM

  subroutine moment_of_inertia(x,MI_eigenvals)   
    Implicit None
    Integer :: i,k,IERR
    Double Precision, Intent(In) :: x(3,Natoms)
    Double Precision, Intent(Out) :: MI_eigenvals(3)
    Double Precision :: IT(3,3),MI_eigenvects(3,3),FV1(3),FV2(3)
    
    ! Compute the elements of the inertia tensor
    IT=0d0
    do i=1,Natoms
       ! Diagonal elements
       IT(1,1)=IT(1,1)+mass(i)*(x(2,i)**2+x(3,i)**2) ! I_xx
       IT(2,2)=IT(2,2)+mass(i)*(x(3,i)**2+x(1,i)**2) ! I_yy
       IT(3,3)=IT(3,3)+mass(i)*(x(1,i)**2+x(2,i)**2) ! I_zz
       
       ! Off-diagonal elements, Symmetric matrix
       IT(1,2)=IT(1,2)-mass(i)*x(1,i)*x(2,i) ! I_xy
       IT(2,3)=IT(2,3)-mass(i)*x(2,i)*x(3,i) ! I_yz
       IT(3,1)=IT(3,1)-mass(i)*x(3,i)*x(1,i) ! I_zx
    enddo
    
    IT(2,1)=IT(1,2)
    IT(3,2)=IT(2,3)
    IT(1,3)=IT(3,1)
    
    ! Diagonalize the inertia tensor and find the eigenvalues
    MI_eigenvals=0d0
    call RS(3,3,IT,MI_eigenvals,0,MI_eigenvects,FV1,FV2,IERR)
    MI_eigenvals=MI_eigenvals/melectron
  end subroutine moment_of_inertia

  ! Sorts arrays in ascending order
  subroutine hpsort(N,RA)
    Implicit None
    Integer, Intent(In) :: N
    Double Precision, Intent(Inout) :: RA(N)
    Integer :: I,IR,J,L
    Double Precision :: RRA
    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if (L > 1) then
       L=L-1
       RRA=RA(L)
    else
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       if (IR.eq.1) then
          RA(1)=RRA
          return
       endif
    endif
    I=L
    J=L+L
20  if (J.le.IR) then
       if (J < IR) then
          if (RA(J) < RA(J+1)) J=J+1
       endif
       if (RRA < RA(J)) then
          RA(I)=RA(J)
          I=J; J=J+J
       else
          J=IR+1
       endif
       goto 20
    endif
    RA(I)=RRA
    goto 10
  end subroutine hpsort

  subroutine diffusion()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,k,n,thread,jsr,Nreplicate,Nkill,Nreplicate1
    Integer :: population(Nmax),replicate(Nmax),kill(Nmax)
    Double Precision :: u,rho1,rho2,Vaver,Pb,energy

    ! If this is a new simulation at the very first iteration, then compute a separate random number generator seed for each omp thread
    if (t==1) then
       jsr=123456789
       do i=0,num_threads-1 
          call shr3(jsr,seed(i))
       enddo
    endif

    ! Initialize Vaver, the population array, and the counter for configs with energies below the global minimum
    Vaver=0d0
    population(1:Nt)=1
    
    ! Parallel Region: Diffusive displacements/random walk
    !$omp parallel &
    !$omp shared (Natoms,Dim,dtau,Nt,x,Vref,sigma,population,seed,energy_thresh) &                 
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
          x(k,n)=x(k,n)+dsqrt(-2*dlog(rho1))*dcos(2*pi*rho2)*sigma(k)
       enddo

       ! Translate the configuration such that the center of mass is at the origin
       call CM(Natoms,x(:,n))

       ! Compute the potential energy
       call pot_energy(energy,x(:,n),Dim) ! x(:,n) in angstrom, energy in kcal/mol
       
       ! If the potential energy of the configuration falls just below the global minimum then remove this walker from the population
       if (energy<energy_thresh) population(n)=0

       ! Compute Pb and compare to a random number: assign walkers to be killed, retained, or replicated
       if (population(n)==1) then ! proceed with computing Pb and Vaver
          energy=energy/autokcalmol
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
       write(1,*) t*dtau,Epartaver*autokcalmol/Kaver,Nt
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
       if (histograms) call histogram()
       if (quenching) call quench()
       if (inertia) call MOI() 
       if (Ndw==1.and.dw_step0==0) return ! Go back to the main program without resetting the DW indicies, etc.
       call dw_reset(j)
    endif
  end subroutine descendant_weighting

  subroutine histogram()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,j,k,l,m,n,j1,j2
    Double Precision :: OH_dist_temp,OH_dist(2),HH_dist,HOH_angle,pair_OO_dist

    OH_hist=0
    HOH_hist=0
    pairs_OO_hist=0

    ! Parallel Region
    !$omp parallel &
    !$omp shared (Nt1,NH2O,Ndw,dw_x,hist_dist_min,hist_angle_min,hist_pairs_OO_min) &
    !$omp shared (d_dist,d_angle,d_pairs_OO,nb_dist,nb_angle,nb_pairs_OO,dw_weight) &
    !$omp private (i,j,k,l,m,n,j1,j2,OH_dist_temp,OH_dist,HH_dist,HOH_angle,pair_OO_dist)
    
    ! Determine if the the histogram subroutine is running in parallel for only the first time it is accessed
    if (t==dw_step_max) then
       num_threads=omp_get_num_threads()   
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif
    
    ! Loop over the configurations with 3 corresponding weights and compute the bond distances, bond angles, and pair distances only once 
    !$omp do reduction (+:OH_hist,HOH_hist,pairs_OO_hist)
    do n=1,Nt1
       if (dw_weight(1,n).ne.0) then
          ! Compute the OH bond lengths and HOH bond angles for each configuration and simultaneously update the histograms
          do i=1,NH2O
             OH_dist=100 ! Initialize the OH_dist array: OH dist is an array of 2 bond lengths per water molecule
             do j=1,Natoms
                if (mod(j,3).ne.1) then
                   OH_dist_temp=sum((dw_x(3*(j-1)+1:3*(j-1)+3,n)-dw_x(9*i-8:9*i-6,n))**2) ! Distance of all the hydrogens to an oxygen
                   if (OH_dist(1)>OH_dist_temp) then ! find the two closest hydrogens
                      OH_dist(2)=OH_dist(1)
                      j2=j1
                      OH_dist(1)=OH_dist_temp
                      j1=j
                   else if (OH_dist(2)>OH_dist_temp) then
                      OH_dist(2)=OH_dist_temp
                      j2=j
                   endif
                endif
             enddo
             do m=1,2 ! Update the OH histogram
                l=(dsqrt(OH_dist(m))-hist_dist_min)/d_dist
                if (l>nb_dist) l=nb_dist
                if (l<0) l=0
                do k=1,Ndw
                   OH_hist(l,k)=OH_hist(l,k)+dw_weight(k,n)
                enddo
             enddo
             
             ! Compute the HH distance used to determine the HOH bond angle
             HH_dist=sum((dw_x(3*j1-2:3*j1,n)-dw_x(3*j2-2:3*j2,n))**2)
             
             ! Compute the HOH bond angle using the law of cosines
             HOH_angle=dacos((OH_dist(2)+OH_dist(1)-HH_dist)/(2*dsqrt(OH_dist(2)*OH_dist(1))))*deg
             l=(HOH_angle-hist_angle_min)/d_angle
             if (l>nb_angle) l=nb_angle
             if (l<0) l=0
             do k=1,Ndw
                HOH_hist(l,k)=HOH_hist(l,k)+dw_weight(k,n)
             enddo
          enddo
          
          ! Compute the pair distances and append to the histogram
          do i=2,NH2O
             do j=1,i-1 
                pair_OO_dist=sum((dw_x(9*i-8:9*i-6,n)-dw_x(9*j-8:9*j-6,n))**2)
                l=(dsqrt(pair_OO_dist)-hist_pairs_OO_min)/d_pairs_OO
                if (l>nb_pairs_OO) l=nb_pairs_OO
                if (l<0) l=0
                do k=1,Ndw
                   pairs_OO_hist(l,k)=pairs_OO_hist(l,k)+dw_weight(k,n)
                enddo
             enddo
          enddo
       endif
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Write out the histogram data to the specified files
    do j=1,Ndw
       do i=0,nb_dist 
          write(7,*) hist_dist_min+(i+0.5)*d_dist,OH_hist(i,j)/(2*NH2O*dw_Nt(j)*d_dist)
       enddo
       do i=0,nb_angle 
          write(8,*) hist_angle_min+(i+0.5)*d_angle,HOH_hist(i,j)/(NH2O*dw_Nt(j)*d_angle)
       enddo
       do i=0,nb_pairs_OO
          write(9,*) hist_pairs_OO_min+(i+0.5)*d_pairs_OO,pairs_OO_hist(i,j)/(Npairs*dw_Nt(j)*d_pairs_OO)
       enddo
    enddo
    call flush(7)
    call flush(8)
    call flush(9)
  end subroutine histogram

  subroutine quench()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,j,k,n,info,iter,nfunc,ngrad,ngrad_sum,nquench
    Double Precision :: energy_opt,x_quench(Dim),pair_dist(Npairs),cosines(Npairs),MI_eigenvals(3)
    
    ! Initialize the necessary variables and arrays
    nquench=0
    ngrad_sum=0
    count_config=0
    
    ! Parallel Region
    !$omp parallel &
    !$omp shared (grad_tol,grad_maxfac,dw_x,Dim,min_energy,Nt1,Npairs,Nisomers,Ndw,energy_thresh) &
    !$omp shared (dw_weight,lib_config_pd,lib_config_cos,r_energy_thresh,r_cos_thresh,r_pairs_thresh) &
    !$omp private (i,j,k,n,info,iter,nfunc,ngrad,energy_opt,x_quench,pair_dist,cosines)
    
    ! Determine if the quenching subroutine is running in parallel for only the first time it is accessed
    if (t==dw_step_max) then
       num_threads=omp_get_num_threads()
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif

    ! Loop over the configurations with 3 corresponding weights and quench them only once   
    !$omp do reduction (+:count_config,nquench,ngrad_sum)
    do n=1,Nt1
       if (dw_weight(1,n).ne.0) then
          ! Assign the current config/walker to a different array for quenching
          x_quench(:)=dw_x(:,n) 
          
          ! Quench the configuration
          call cg_descent(grad_tol,grad_maxfac*Dim,x_quench,Dim,pot_energy,gradient,&
               info,energy_opt,iter,nfunc,ngrad,energy_thresh)
          
          ! If the configuration is physical then calculate its pair distances, bond cosines, and moments of inertia
          if (info==0.or.info==2) then 
             ngrad_sum=ngrad_sum+ngrad
             nquench=nquench+1
             
             ! Compute the pair distances and cosines of the quenched config
             ! Compute the moments of inertia of the unquenched config
             call OO_distances(x_quench,pair_dist)
             call OH_bonds(x_quench,cosines)               
             
             do k=1,Nisomers
                if (dabs(energy_opt-min_energy(k))<r_energy_thresh &   
                     .and.distance(Npairs,pair_dist,lib_config_pd(:,k))<r_pairs_thresh &
                     .and.distance(Npairs,cosines,lib_config_cos(:,k))<r_cos_thresh) then
                   ! Assign config "n" with its coresponding weight to one of the isomers in the library
                   do j=1,Ndw
                      count_config(k,j)=count_config(k,j)+dw_weight(j,n) 
                   enddo
                   goto 100 ! Once the config is assigned, move on to check for the next non-zero weight
                endif
             enddo
             
             ! If none of the criteria are satisified for any of the configs, then assign the config to "others"
             do j=1,Ndw
                count_config(Nisomers+1,j)=count_config(Nisomers+1,j)+dw_weight(j,n)
             enddo
             ! Otherwise, if the cg_descent program encounters internal errors or gives an unphysical config don't bother with the assignment process.
          else if (info==1.or.info.ge.3) then
             do j=1,Ndw
                count_config(Nisomers+2,j)=count_config(Nisomers+2,j)+dw_weight(j,n)
             enddo
          endif
       endif
100    continue
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Compute the populations/fractions of each isomer at time t
    do j=1,Ndw
       do k=1,Nisomers+2
          write(4,*) j,k,count_config(k,j)/dble(dw_Nt(j))
       enddo
    enddo
    
    ! Write out the number of gradient evaluations and the number of times configs were assigned to "others or unphysical"
    write(11,*) t*dtau,nquench,dble(ngrad_sum)/nquench
    call flush(4)
    call flush(11)
  end subroutine quench

  subroutine MOI()
    use omp_lib
    Implicit None
    Logical :: parallel
    Integer :: i,j,n
    Double Precision :: MI_eigenvals(3)
    
    ! Initialize the average moment of inertia
    MI_eigenvals_aver=0d0
    
    ! Parallel Region
    !$omp parallel &
    !$omp shared (Nt1,Ndw,dw_weight,dw_x) &
    !$omp private (j,n,MI_eigenvals)
    
    ! Determine if the quenching subroutine is running in parallel for only the first time it is accessed
    if (t==dw_step_max) then
       num_threads=omp_get_num_threads()
       parallel=omp_in_parallel()
       write(16,*) parallel,num_threads
       call flush(16)
    endif
    
    !$omp do reduction (+:MI_eigenvals_aver)
    do n=1,Nt1
       if (dw_weight(1,n).ne.0) then
          call moment_of_inertia(dw_x(:,n),MI_eigenvals)
          do j=1,Ndw
             MI_eigenvals_aver(:,1,j)=MI_eigenvals_aver(:,1,j)+MI_eigenvals(:)*dw_weight(j,n)
          enddo
       endif
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Compute the average moments of inertia and rotational constants
    do j=1,Ndw
       do i=1,3
          write(17,*) j,1,MI_eigenvals_aver(i,1,j)/dw_Nt(j)
          write(19,*) j,1,(dw_Nt(j)*autoMHz*bohr**2)/(2*MI_eigenvals_aver(i,1,j)*melectron)
       enddo
    enddo
    call flush(17)
    call flush(19)
  end subroutine MOI
end module dmc_module

Program dmc_true_gs
  
  use dmc_module
  use omp_lib
  
  Implicit None
  Integer :: i,j,k,n
  Double Precision :: MI_eigenvals(3)

  ! All distances are in angstroms
  open(2,file='input.dat')
  read(2,*) Natoms
  read(2,*) N0
  read(2,*) dtau
  read(2,*) tmax
  read(2,*) Kaver
  read(2,*) Ndw,dw_tau_max,dw_tau_incr
  read(2,*) hist_dist_min,hist_dist_max,nb_dist 
  read(2,*) hist_angle_min,hist_angle_max,nb_angle
  read(2,*) hist_pairs_OO_min,hist_pairs_OO_max,nb_pairs_OO
  read(2,*) histograms,quenching,inertia
  read(2,*) write_configs,step_configs,Nt_incr
  read(2,*) isotope
  read(2,*) potential
  read(2,*) Vref_parm 
  read(2,*) r_energy_thresh,r_pairs_thresh,r_cos_thresh
  read(2,*) grad_tol,grad_maxfac
  read(2,*) energy_thresh ! Threshold energy in kcal/mol taken to be just below the PES global minimum
  read(2,*) restart_parm
  read(2,*) num_threads
  read(2,*) file_type ! String that determines if the simulation is new or needs to be restarted
  close(2)

  ! Read in the number of reference configurations from the readme file
  open(12,file='library_opt_configs/readme.dat')
  read(12,*)
  read(12,*)
  read(12,*) Nisomers

  ! Initialize the parameters for dimensionality, the maximum number of walkers, and the maximum number of iterations
  Dim=3*Natoms
  Nmax=N0*1.1
  N_iter_max=tmax/dtau
  
  ! Compute the number of water molecules and number of OO pairs
  NH2O=Natoms/3
  Npairs=(NH2O*(NH2O-1))/2

  ! Define the distance parameter for the histograms and allocate memory to each array
  if (histograms) then
     d_dist=(hist_dist_max-hist_dist_min)/(nb_dist+1)
     d_angle=(hist_angle_max-hist_angle_min)/(nb_angle+1)
     d_pairs_OO=(hist_pairs_OO_max-hist_pairs_OO_min)/(nb_pairs_OO+1)
  endif

  ! Initialize the first value of dw_step and the increment for computing the weights
  dw_step_max=dw_tau_max/dtau
  dw_incr=dw_tau_incr/dtau
  if (Ndw==1.and.dw_incr==0) then
     dw_step0=0 ! Compute the observables using the tau_dw=0 result only
  else
     dw_step0=dw_step_max-(Ndw-1)*dw_incr ! Compute DW results for observables with tau_dw spaced out in regular intervals
  endif

  allocate(atom_type(Natoms),mass(Natoms),seed(0:num_threads-1),x(Dim,Nmax),sigma(Dim))
  allocate(dw_walkers(Nmax),min_energy(Nisomers))
  if (histograms.or.quenching.or.inertia) then 
     allocate(dw_weight(Ndw,Nmax),dw_x(Dim,Nmax),dw_Nt(Ndw),MI_eigenvals_aver(3,Nisomers,Ndw))
  endif
  if (histograms) then
     allocate(OH_hist(0:nb_dist,Ndw),HOH_hist(0:nb_angle,Ndw),pairs_OO_hist(0:nb_pairs_OO,Ndw))
  endif
  if (quenching) allocate(lib_config_pd(Npairs,Nisomers),lib_config_cos(Npairs,Nisomers),count_config(Nisomers+2,Ndw))
  
  ! If the file is new, then write in all of the file headers
  ! Otherwise, if the file is restart, then open up all of the old files and start appending from the nearest Kaver step where the simulation left off
  if (file_type=='new') then
     open(1,file='energy.dat') 
     write(1,*) '# Natoms=',Natoms,'N0=',N0,'dtau=',dtau
     write(1,*) '# Projection_time=',tmax,'Kaver=',Kaver
     if (isotope=='D_only') then 
        write(1,*) '# Dmass=',Dmass/melectron,'Omass=',Omass/melectron
     else if (isotope=='H_and_D'.or.isotope=='one_D'.or.isotope=='one_H') then 
        write(1,*) '# Hmass=',Hmass/melectron,'Dmass=',Dmass/melectron 
        write(1,*) '# Omass=',Omass/melectron
     else
        write(1,*) '# Hmass=',Hmass/melectron,'Omass=',Omass/melectron
     endif
     write(1,*) '# Incr Opt=',Kaver
     write(1,*) '# Number of Reference Configs=',Nisomers
     write(1,*) '# Number of Threads=',num_threads
     if (potential=='tip4p') write(1,*) '# Potential=TIP4P'
     if (potential=='ttm3') write(1,*) '# Potential=TTM3'
     if (potential=='mbpol') write(1,*) '# Potential=MBPOL'
     
     open(16,file='parallel_test.dat')
     write(16,*) '# Column 1: Parallel On or Off'
     write(16,*) '# Column 2: Number of Threads in Use'
     
     if (restart_parm) then 
        open(21,file='completed_storage.dat')
        write(21,*) '# Column 1: Projection Time'
        write(21,*) '# Column 2: Storage of Data Completed'
        call flush(21)
     endif
     
     open(31,file='moments_of_inertia_isomers.dat')
     write(31,*) '# Column 1: Isomer Number'
     write(31,*) '# Columns 2-4: Moments of Inertia (amu*A**2)'
     write(31,*) '# Columns 2-4: Rotational Constants (MHz)'
     call flush(31)
     
     if (histograms.or.quenching.or.inertia) then
        open(30,file='dw_walker_count.dat')
        write(30,*) '# Column 1: Projection Time'
        write(30,*) '# Column 2: Number of Walkers with Non-zero Weight'
        call flush(30)
     endif
        
     if (inertia) then
        open(17,file='moments_of_inertia.dat')
        write(17,*) '# Column 1=DW Weight Index'
        write(17,*) '# Column 2=Isomer Index'
        write(17,*) '# Column 3=Average Moment of Inertia (amu*A**2)'
        
        open(19,file='rotational_constants.dat')
        write(19,*) '# Column 1=DW Weight Index'
        write(19,*) '# Column 2=Configuration Index'
        write(19,*) '# Column 3=Average Rotational Constants (MHz)'
        
        call flush(17)
        call flush(19)
     endif

     if (histograms) then
        open(7,file='histogram_OH_dist.dat') 
        write(7,*) '# Histogram OH Distance=',hist_dist_min,hist_dist_max
        write(7,*) '# Histogram OH Distance=',d_dist
        write(7,*) '# Column 1=OH Distances'
        write(7,*) '# Column 2=Normalized Bin Count'
  
        open(8,file='histogram_HOH_angle.dat') 
        write(8,*) '# Histogram HOH Angles=',hist_angle_min,hist_angle_max
        write(8,*) '# Histogram HOH Angles=',d_angle
        write(8,*) '# Column 1=HOH Angles'
        write(8,*) '# Column 2=Normalized Bin Count'
    
        open(9,file='histogram_pairs_OO_dist.dat') 
        write(9,*) '# Histogram OO Pairs=',hist_pairs_OO_min,hist_pairs_OO_max
        write(9,*) '# Histogram OO Pairs=',d_pairs_OO
        write(9,*) '# Column 1=OO Pair Distances'
        write(9,*) '# Column 2=Normalized Bin Count'

        call flush(7)
        call flush(8)    
        call flush(9)
     endif

     if (quenching) then
        open(4,file='populations.dat')
        write(4,*) '# Column 1=DW Weight Index'
        write(4,*) '# Column 2=Isomer Number'
        write(4,*) '# Column 3=Isomer Fraction'
        
        open(11,file='grad_evals.dat')
        write(11,*) '# Column 1=Projection Time'
        write(11,*) '# Column 2=Number of Successful Gradient Evals'
        write(11,*) '# Column 3=Average Number of Gradiant Evals'
        
        call flush(4)
        call flush(11)
     endif
     
  else if (file_type=='restart') then ! open up the files as "old" and start appending from the nearest energy partial average where the simulation terminated
     open(1,file='energy.dat',status='old',position='append')
     open(16,file='parallel_test.dat',status='old',position='append')
     if (restart_parm) open(21,file='completed_storage.dat',status='old',position='append') 
     open(18,file='restart_configs.dat',form='unformatted',status='old')
     open(31,file='moments_of_inertia_isomers.dat',status='old',position='append')
     
     if (histograms.or.quenching.or.inertia) then 
        open(30,file='dw_walker_count.dat',status='old',position='append')
     endif

     if (inertia) then
        open(17,file='moments_of_inertia.dat',status='old',position='append')
        open(19,file='rotational_constants.dat',status='old',position='append')
     endif

     if (histograms) then
        open(7,file='histogram_OH_dist.dat',status='old',position='append')
        open(8,file='histogram_HOH_angle.dat',status='old',position='append')
        open(9,file='histogram_pairs_OO_dist.dat',status='old',position='append')
     endif

     if (quenching) then
        open(4,file='populations.dat',status='old',position='append')
        open(11,file='grad_evals.dat',status='old',position='append')
     endif
  endif
  
  ! Read in the reference configurations based on their file numbering pattern and place an equal number of random walkers at the position of each reference config
  N0_part=N0/Nisomers   
  do i=1,Nisomers
     write(file_num,'(i7)') i
     coord_config='library_opt_configs/'//trim(adjustl(file_num))//'.xyz'
     
     ! Read in the reference configs in the library and write out the minimum energies of the reference configs to the readme file 
     call read_config(coord_config,x(:,1+(i-1)*N0_part),min_energy(i))
     write(12,*) trim(adjustl(file_num))//'.xyz',min_energy(i)
     
     ! Store the pair distances, cosines, and moments of inertia of the reference configs
     if (quenching) then 
        call OO_distances(x(:,1+(i-1)*N0_part),lib_config_pd(:,i))
        call OH_bonds(x(:,1+(i-1)*N0_part),lib_config_cos(:,i))
     endif
     
     ! Place an equal share of random walkers at the position of the reference configuration in config space
     do n=2,N0_part
        x(:,n+(i-1)*N0_part)=x(:,1+(i-1)*N0_part)
     enddo
  enddo
  close(12)
  
  ! If N0 and N_configs do not divide evenly, place the remaining walkers at the position of each refenerce config until no walkers remain
  do i=1,mod(N0,Nisomers)            
     x(:,i+N0_part*Nisomers)=x(:,1+(i-1)*N0_part)
  enddo

  ! Store the masses of each atom type in the mass array
  ! Calculate the width (sigma) of the Gaussian distribution for each type of atom and the total mass of the cluster
  ! Multiply by width by bohr to convert sigma to angstroms
  do i=1,Natoms
     if (isotope=='D_only'.and.atom_type(i)=='H') atom_type(i)='D' ! change all H to D
     mass(i)=atom_mass(atom_type(i))
     sigma(i*3-2:i*3)=dsqrt(dtau/mass(i))*bohr
  enddo
  sum_mass=sum(mass)
  write(16,*) 'Sum_mass=',sum_mass/melectron
  write(16,*) 'Atom_type=',atom_type
  call flush(16)

  ! Translate the configs so the center of mass is at the origin
  do n=1,N0
     call CM(Natoms,x(:,n))
  enddo
  
  ! Compute the moments of inertia and rotational constants for each isomer in the library
  do i=1,Nisomers
     call moment_of_inertia(x(:,1+(i-1)*N0_part),MI_eigenvals)
     write(31,*) i,(MI_eigenvals(j),j=1,3) 
     write(31,*) i,((autoMHz*bohr**2)/(2*MI_eigenvals(j)*melectron),j=1,3)
     write(31,*)
  enddo
  close(31)
     
  ! Initialize Vref, all of the histogram arrays, quenching arrays, and moment of inertia arrays if the simulation is "new"
  if (file_type=='new') then
     ! If file type is new, then get the initial Vref from one of the basins of attraction
     ! Compute the initial value of Vref: Can start Vref/random walk from different basins of attraction  
     call pot_energy(Vref,x(:,1+(Vref_parm-1)*N0_part),Dim)
     write(1,*) '# Initial Vref (kcal/mol)=',Vref
     write(1,*) '# Energy Thresh (kcal/mol)=',energy_thresh
     call flush(1)
     Vref=Vref/autokcalmol
     
     ! For a "new" simulation, then the initial Nt gets assigned to N0 and t1 is 1
     t1=1
     Nt=N0
     
     ! Otherwise, "restart" the simulation where it left off 
     ! Read in the stored values to restart the simulation
  else if (file_type=='restart') then
     x=0d0
     read(18) Nt
     read(18) Vref
     read(18) t1
     read(18) (seed(i),i=0,num_threads-1)
     read(18) (x(:,n),n=1,Nt)
     close(18)
     t1=t1+1
  endif
  Epartaver=0d0

  ! Initialize/reset the indicies for descendant weighting
  do n=1,Nt
     dw_walkers(n)=n
  enddo
  
  ! Set the number of threads for all of the parallel regions embedded within the subroutines
  ! If not set, the maximum number of threads on the machine is used by default
  call omp_set_num_threads(num_threads)
  
  ! Main/outer time step loop
  if (histograms.or.quenching.or.inertia) call dw_reset(j)
  do t=t1,N_iter_max
     ! Perform the diffusion steps and branching
     call diffusion()

     if (mod(t,Kaver)==0) then
        ! Overwrite the unformatted "restart_configs.dat" file
        ! Store and write out Vref, Nt, projection time, seed, and all of the configurations every Kaver steps
        if (restart_parm) then
           open(18,file='restart_configs.dat',form='unformatted',status='replace')   
           write(18) Nt
           write(18) Vref
           write(18) t
           write(18) (seed(i),i=0,num_threads-1)
           write(18) (x(:,n),n=1,Nt) ! Store the configs for the next iteration: Restarts can occur every Kaver steps when Epartaver is reset to 0
           close(18)
           write(21,*) t*dtau,'Storage=Completed'
           call flush(21)
        endif
     endif
        
     ! Compute the observables using the DW method
     if (histograms.or.quenching.or.inertia) then
        if (mod(t,dw_step_max)==dw_step.or.mod(t,dw_step_max)==0) call descendant_weighting(j)
     endif
   
     ! Write out all of the configurations at a specified time to be used in some capacity later
     if (write_configs.and.t==step_configs) then
        open(14,file='coord_configs.dat',form='unformatted')
        do n=1,Nt
           if (mod(n,Nt_incr)==0) write(14) real(x(:,n))
        enddo
        close(14)
     endif
  enddo
end Program dmc_true_gs

