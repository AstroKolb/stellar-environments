subroutine ppmlr

! Using the 1D arrays of rho, u, and P, perform a 1D lagrangian hydrodynamics evolution:
!   - set boundary conditions by filling ghost zones
!   - obtain parabolic interpolations of rho, u, P
!   - compute input states from these interpolations for Riemann problem
!   - call the Riemann solver to find the time averages umid, pmid
!   - evolve the finite difference equations to get updated values of rho, u, and E
! Then perform a conservative remap back to the Eulerian grid
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL :: xwag, alternate
REAL, DIMENSION(maxsweep) :: dr, du, dp, r6, u6, p6, rl, ul, pl
REAL, DIMENSION(maxsweep) :: rrgh, urgh, prgh, rlft, ulft, plft, umid, pmid
REAL, DIMENSION(maxsweep) :: xaf, dxf

!-----------------------------------------------------------------------

! Compute parabolic coefficients and volume elements
call paraset( nmin-4, nmax+5, para, dx, xa )

! Calculate flattening coefficients for smoothing near shocks
call flatten

! Interpolate parabolae for fluid variables 
call parabola(nmin-4, nmax+4, para, p, dp, p6, pl, flat)
call parabola(nmin-4, nmax+4, para, r, dr, r6, rl, flat)
call parabola(nmin-4, nmax+4, para, u, du, u6, ul, flat)

! Integrate parabolae to get input states for Riemann problem
call states( pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh )
   
! Call the Riemann solver to obtain the zone face averages, umid and pmid
call riemann( nmin-3, nmax+4, gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid )

! do lagrangian update using umid and pmid
call evolve( umid, pmid )

!#########################################################################
! EXTRA DISSIPATION TO REDUCE CARBUNCLE NOISE

alternate = 2*(ncycle - 2*(ncycle/2)) - 1
xwag = sum(flat)
if ((xwig*xwag /= 0.0).and.(sweep /= 'x')) then ! wiggle grid, remap, then remap back to Eulerian grid

 ! save Eulerian coordinates for second remap
 xaf = xa0
 dxf = dx0

 ! wiggle grid where there is a shock, including edges
 do n = nmin, nmax+1
  if (max(flat(n-1),flat(n)) > 0.0) xa0(n) = xa0(n) + xwig*dx0(n)*alternate
  dx0(n-1) = xa0(n) - xa0(n-1)
 enddo
 dx0(nmax+1) = xa0(nmax+2) - xa0(nmax+1)

 call remap

 ! put wiggled grid into xa(n), but keep ghost cells untouched
 do n = nmin, nmax+1
  xa(n)   = xa0(n)
  dx(n-1) = xa(n) - xa(n-1)
 enddo
 dx(nmax+1) = xa(nmax+2) - xa(nmax+1)

 call volume(nmin, nmax, ngeom, radius, xa, dx, dvol)

 ! put Eulerian grid back into xa0
 xa0 = xaf
 dx0 = dxf

endif
!#########################################################################

! remap onto original Eulerian grid
call remap

return
end

