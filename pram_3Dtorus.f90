module pram
!
! Input file for general simulation structure
! Each model parameter is written in each subroutines of mdgrmodel.f90
!
  implicit none

  character*65 :: jname='3D stationary torus' !- problem name
  character*15 :: date='May. 04. 2015' !- Date

  integer, parameter :: imax=128, jmax=64+6, kmax=128 !- Grid number  
  integer, parameter :: nv=9 !- Number of variables
  integer, parameter :: nkmax=20 !- Wave number for perturbation
  integer, parameter :: iprocs=2, jprocs=1, kprocs=1 !- CPU number in i-,j-, and  k- direction

!  integer, parameter :: ncpu=2 ! Number of CPUs (=npe) for MPI
!  integer, parameter :: kkm= (kmax-1)/ncpu+16+1 ! modified grid number for MPI 
!  integer, parameter :: kkm=64

  integer, parameter :: icpu=2000000 !- CPU time
  real(8), parameter :: dt0=0.01d0 !- not used
  integer, parameter :: itmax=2000000 !- Max number of iteration
  integer, parameter :: nshot=1 !- Number of Output
  integer, parameter :: icres=0 !- Researt parameter
  real(8), parameter :: dtmin=1.d-9 !- Minimum dt 
  real(8), parameter :: tmax=1.0d0 !- Maximum simulation time
  real(8), parameter :: cfl=0.3d0 !- CFL number
  integer, parameter :: isskip=1 !- Skip for researt data output
  integer, parameter :: iter=1000 !- Maxmum iteration number for Newton-Raphson

  real(8), parameter :: pi=3.14159265358979323d0 !- Pi
  integer, parameter :: iseed=2012 !- seed number for random number generation

  real(8), parameter :: pmin=1.0d-8 !- limiter for minimum pressure
  real(8), parameter :: pmax=1.0d6 !- limiter for maximum pressure
  real(8), parameter :: dmin=1.0d-5 !- limiter for minimum density
  
  real(8), parameter :: c0=1.0d0 !- light speed
  real(8), parameter :: gam=4.d0/3.d0 !- adiabatic index
  real(8), parameter :: gam0=4.d0/3.d0
  real(8), parameter :: kpol=1.d-3 !- polytropic constant
!
  real(8), parameter :: rbh=1.d0 !- black hole radius rg=GM/c^2
  real(8), parameter :: akm=0.5d0 !- black hole spin

  !
  real(8), parameter :: vlim=1.d0 !- limiter for velocity/c
  real(8), parameter :: adampa=0.d0  !- parameter for damping 
  real(8), parameter :: adamp=0.d0 !- parameter for damping
  real(8), parameter :: rdamp=0.d0 !- parameter for damping

  real(8), parameter :: xmin=2.1d0
  real(8), parameter :: xmax=40.0d0 !- min & max size of i-direc
  real(8), parameter :: ymin=0.0d0-3.*(2.*pi)/float(jmax-7)
  real(8), parameter :: ymax=2.*pi+3.*(2.*pi)/float(jmax-7)  !- min & max size of j-direc  
  real(8), parameter :: zmin=0.01d0
  real(8), parameter :: zmax=3.131592d0 !- min & max size of k-direc

  integer, parameter :: ix1=1, ix2=1, ix3=1 !- grid spacing 
  real(8), parameter :: R0=0.5d0, hslope=0.3d0 !- parameter for grid spacing
  real(8), parameter :: tilang=0.d0*(pi/180.d0) ! tilted angle

  integer, parameter :: model=9 !- simulation model
  integer, parameter :: metric=203 !- metric
  integer, parameter :: ieos=0 !- EoS
  integer, parameter :: icha=1 !- calculation method for wave speed

  integer, parameter :: iwvec=2 !- inversion procedure scheme
  integer, parameter :: irkt=3 !- time-advance scheme
  integer, parameter :: irec=1 !- reconstruction scheme
  integer, parameter :: ihll=1 !- Approximate Riemann solver scheme
  integer, parameter :: ict=1  !- Constrained Transport scheme

  integer, parameter :: iflag(7)=(/0, 0, 0, 0, 0, 0, 0/)
  integer :: iboux1in(9)=(/5, 5, 5, 5, 5, 5, 5, 5, 5/) !- boundary cond.x1-in
  integer :: iboux1ot(9)=(/5, 5, 5, 5, 5, 5, 5, 5, 5/) !- boundary cond.x1out
  integer :: iboux2in(9)=(/1, 1, 1, 1, 1, 1, 1, 1, 1/) !- boundary cond.x2-in
  integer :: iboux2ot(9)=(/1, 1, 1, 1, 1, 1, 1, 1, 1/) !- boundary cond.x2-out
  integer :: iboux3in(9)=(/5, 5, 5, 5, 5, 5, 5, 5, 5/) !- boundary cond.x3-in
  integer :: iboux3ot(9)=(/5, 5, 5, 5, 5, 5, 5, 5, 5/) !- boundary cond.x3-out 
!
!!!!!!! Parameter explanation !!!!!!!
!
! Researt Option
!
!     icres = 0  : restart off
!     icres = 1  : restart on
!
! Simulation Model
!
!    model = 1 : 1D Test
!
! Mertic (coordinates) 
!
!     metric =   100 :  Minkowski spacetime (x,y,z)
!     metric =   103 :  Schwarzschild spacetime in BL coord (r,\theta,\phi)
!     metric =   203 :  Kerr spacetime in BL coord (r, \theta, \phi)
!
! Equation of State
!
!     ieos=0 : Gamma-law constant equation of state    
!     ieos=1 : TM variable equation of state (cal enthalpy)
!     ieos=2 : TM variable equation of state (cal internal energy)
!
! Reconstruction scheme
!
!     irec = 1 : MC slope-limiter reconstruction (2nd order)
!     irec = 2 : Minmod slope-limiter reconstruction (2nd order)
!     irec = 3 : MUSCL method (3rd order), not slope-limiter
!     irec = 4 : Convex ENO reconstruction method (3rd order)
!     irec = 5 : rPPM reconstruction method (4th order)
!     irec = 6 : Monotonicity-Preserving 5th reconstruction (5th order)
!     irec = 7 : Weighted ENO 5th reconstruction (5th order)
!     irec = 8 : MP-WENO reconstruction
!     irec = 9 : WENO-Z 5th reconstruction
!     irec = 10: Mapped WENO 5th reconstruction
!     irec = 11: Limited Reconstruction (3rd order)
!
! Approximate Riemann solver scheme
!
!     ihll = 1 : HLLE single-state approximate Riemann solver
!     ihll = 2 : HLLC two-state approximate Riemann solver
!                (Mignone & Bodo 200?)
!     ihll = 3 : HLLC two-state approximate Riemann solver
!                (Honkkila & Janhunen 2007)
!     ihll = 4 : HLLD multi-state approximate Riemann solver
!                (underconstruction)
!     ihll = 5 : HLLD multi-state approximate Riemann solver
!                (Mignone et al. 2010)
!
! Constrained Transport scheme
!
!     ict = 0 : No contrained transport
!     ict = 1 : Flux CT sheme (Toth 2000)
!     ict = 2 : Modified Flux CT scheme (Gardiner & Stone 2005)
!     ict = 3 : Upwind Flux CT scheme (Gardiner & Stone 2005)
!
! RK time advance scheme
!
!     irkt = 2 : 2nd order Runge-Kutta time advance method
!     irkt = 3 : 3rd order Runge-Kutta time advance method
!
! Inversion procedure scheme
!
!     iwvec=1 : new 2 variable newton-raphson (from Noble et al. 2005)
!     iwvec=2 : 1 variable newton-raphson with variable EOS
!               (from Mignone & McKinney 2007) 
!
! Boundary Condiiton
!
!               1   :    periodic boundary
!               2   :    fixed boundary
!               3   :    Neumann boundary
!               4   :    free boundary
!               5   :    u_0 = u_1       
!               7   :    radiative boundary without eigenvalue
!              12   :    anti-symmetric boundary condition (u_0=-u_1)
!              13   :    regid wall boundary condition (u_0=0.0)
!              15   :    special reflecting boundary 
!                        for jet propagation (z-direction only)
!              17   :    special radiative boundary without eigenvalue
!                        for jet propagation (z-direction only)
!              19   :    special boundary of u_0=u_1
!                        for jet propagation (z-direction only)
!
end module pram

