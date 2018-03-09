module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 logical :: yin

logical :: brk

 character(len=50) :: prefix, scale      ! prefix for output filenames
 character(len=04) :: radlaw

 integer :: ndim
 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file marker


 real :: time, dt, timem, timep, svel 
 real :: gam, pi, gamm
 real, parameter :: courant = 0.30        ! timestep fraction of courant limit
 real, parameter :: xwig = 0.05          ! fraction of a zone to wiggle grid for dissipation
 real, parameter :: smallp = 1.0e-45     ! Set small values to prevent divide by zero
 real, parameter :: smallr = 1.0e-45
 real, parameter :: small  = 1.0e-45

! physical constants
 real, parameter :: GRV  = 6.67e-08	! gravitational constant
 real, parameter :: Rsun = 6.95e+10	! solar radius
 real, parameter :: AU   = 1.50e+13 ! 1 AU
 real, parameter :: Msun = 1.99e+33 ! solar mass
 real, parameter :: Lsun = 3.89e+33 ! solar luminosity
 real, parameter :: mp   = 1.67e-24	! proton mass
 real, parameter :: kB   = 1.67e-16	! Boltzmann constant
 real, parameter :: GM   = GRV*Msun

! rotation variables
 real :: sep, omega, rcm, opd
 real :: GMP, TmpP, dmP, LumP, Cs2P, ucP, rcP, rhoP, IIP, radP
 real :: GMS, TmpS, dmS, LumS, Cs2S, ucS, rcS, rhoS, IIS, radS

 !real, dimension(6) :: uin

 real :: sth, cth, sph, cph
 integer :: step
 real :: x1min, x1max, x2min, x2max, Z1, Z2, Z3
 real :: rad1, rad2, rad3, kap1, kap2, kap3
 real :: p_flr, T_flr, dmscal, x0min, minrad

 ! step 2/3 parameters
 integer                                 :: amax
 real, dimension(:),         allocatable :: alpha
 real, dimension(:,:,:),     allocatable :: mybound
 real, dimension(:,:,:,:),   allocatable :: myray
 real, dimension(:,:,:,:,:), allocatable :: myring
real :: assdw
integer :: nssdw

 real :: uinflo, dinflo, vinflo, winflo, pinflo, einflo 
 real :: uotflo, dotflo, votflo, wotflo, potflo, eotflo
      
end module 


module sweepsize
!=======================================================================
! Dimension of 1D sweeps.  maxsweep must be as long as the longest of the 
! 3D arrays PLUS the ghost zones:  maxsweep = max(imax,jmax,kmax) + 12
!----------------------------------------------------------------------

 integer, parameter :: maxsweep=1036 

end module sweepsize


module sweeps      
!=======================================================================
!  data structures used in 1D sweeps, dimensioned maxsweep (in sweep_size)
!-----------------------------------------------------------------------

use sweepsize

   character(len=1) :: sweep                                    ! direction of sweep: x,y,z
   integer :: nmin, nmax, ngeom                                 ! number of first and last real zone  
   real, dimension(maxsweep) :: r, p, e, q, u, v, w, c          ! fluid variables
   real, dimension(maxsweep) :: xa, xa0, dx, dx0, dvol          ! coordinate values
   real, dimension(maxsweep) :: f, flat                         ! flattening parameter
   real, dimension(maxsweep,5) :: para                          ! parabolic interpolation coefficients
   real :: radius, theta, stheta, ctheta, sphi, cphi

end module sweeps


module noise
!=======================================================================
! Perlin Noise Generation
!----------------------------------------------------------------------

   integer, parameter                            :: seed    = 183214	! pseudorandomization seed
   integer, parameter                            :: octaves = 2			! controls 'smoothness' of noise
																		! (don't set higher than 4)
   ! length scales
   integer, parameter                            :: x_lscale = 2.0*3.14
   integer, parameter                            :: y_lscale = 1.0*3.14
   integer, parameter                            :: z_lscale = 1.0

   ! tile size
   integer, parameter                            :: tile    = 512
   integer, parameter                            :: tiledim = 512
   integer, parameter                            :: repeats = 1			! generally don't repeat

   ! noise values
   real,    parameter                            :: persistence = 0.8
   real                                          :: tilesize, maxamplitude

   ! pseudo-random permutation
   integer, dimension(2*tiledim)                 :: pp
   integer, dimension(  tiledim)                 :: perm
   

end module noise


