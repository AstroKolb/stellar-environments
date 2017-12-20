      subroutine flatten

! Flaten looks for signs of strong shocks and sets the variable flat
! between 0.0 (smooth flow) and 1.0 (strong shock).
! Simplified method of C&W: eqn. A.1 and A.2.
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

! LOCALS
integer :: n
real, dimension(maxsweep) :: steep
real :: delp1, delp2, shock, temp1, temp2, mach2, old_flat

real, parameter :: omega1 = 0.5
real, parameter :: omega2 = 10.0
real, parameter :: epsilon = 1.0

!--------------------------------------------------------------------------
! Look for presence of a shock using pressure gradient and sign of
! velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
! Compute steepness parameter based on steepness of pressure jump IF 
! there is a shock.

do n = nmin-4, nmax+4

  delp1 = p(n+1) - p(n-1)
  delp2 = p(n+2) - p(n-2)
  if(abs(delp2) .lt. small) delp2 = small
  shock = abs(delp1)/min(p(n+1),p(n-1))-epsilon
  shock = max(0.0,shock)
  if(shock .gt. 0.0) shock = 1.0
  if(u(n-1).lt.u(n+1)) shock = 0.0
  temp1 = ( delp1 / delp2 - omega1 ) * omega2
  steep(n) = shock * max( 0., temp1 )

  !Do no apply any flattening if Mach > 100
  mach2 = (u(n)**2+v(n)**2+w(n)**2)/(gam*p(n)/r(n))
  temp2 = max(0.0, 1.0e+04 - mach2)
  temp2 = min(temp2, 1.0)
  steep(n) = steep(n) * temp2

enddo

! Set phony boundary conditions for the steepness parameter

steep(nmin-5) = steep(nmin-4)
steep(nmax+5) = steep(nmax+4)

! Set flatening coefficient based on the steepness in nieghboring zones

do n = nmin-4, nmax+4
  temp2   = max( steep(n-1), steep(n), steep(n+1) )
  flat(n) = max( 0.0, min( 0.5, temp2 ) )
enddo

! flat(n) should be set to old_flat if no shock in this direction

do n = nmin-3, nmax+3
  old_flat = f(n) - int(f(n))
  if (flat(n).gt.0.0) then
    flat(n) = max(flat(n),old_flat)
    f(n)    = max(flat(n) + (2.0*ndim-3.0), 0.0)   ! f(n)=0 for 1D, 1.f for 2D, 3.f for 3D
   else
    f(n)    = max(0.0, f(n) - 1.0)
    flat(n) = old_flat
  endif
enddo


return
end


