subroutine forces(xf,grav,fict)
  
! Calculate both real (ie, gravity and whatnot) 
! and fictitious (ie, coriolis and centrifugal) forces.
!---------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
integer :: n
real :: xf0, sinxf0
real :: rs1n2, rad      ! step 1
real :: rad2p, rad2s, r2s, sal, cal    ! step 2
real, dimension(maxsweep) :: grav, fict, xf

!----------------------------------------------------------------

if (step == 1) then
   if (sweep == 'x') then
      if (yin) then
         do n = nmin-4, nmax+5
            rad = xf(n)
            rs1n2 = (rad**2 + sep**2 - 2.0*rad*sep*sth*cph)**0.5

            grav(n) = - GMP / xf(n)**2 * c(n)                           !&		! primary
                     ! - GMS / rs1n2**3 * (rad-sep*sth*cph)   			   !&		! secondary
                     ! + LumP * exp(-c(n)) / (4.0*pi*3.0e+10*xf(n)**2)   		! radiation
            
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                         &		! grid
                      + omega**2*xf(n)*sth**2                           &		! centrifugal
                      - omega**2*rcm*sth*cph                            &		! origin shift
                      + 2.0*omega*sth*w(n)										   ! coriolis

            !grav(n) = 0.0
            !fict(n) = (w(n)**2+v(n)**2)/xf(n)
         enddo
      else
         do n = nmin-4, nmax+5
            rad = xf(n)
            rs1n2 = (rad**2 + sep**2 + 2.0*rad*sep*sth*cph)**0.5

            grav(n) = - GMP / xf(n)**2 * c(n)                           !&		! primary
                     ! - GMS / rs1n2**3 * (rad+sep*sth*cph)            	!&		! secondary
                     ! + LumP * exp(-c(n)) / (4.0*pi*3.0e+10*xf(n)**2)   		! radiation
            
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                         &		! grid
                      + omega**2*xf(n)*(cth**2*sph**2+cph**2)           &		! centrifugal
                      + omega**2*rcm*sth*cph                            &		! origin
                      - 2.0*omega*(sph*cth*w(n)-cph*v(n))						   ! coriolis

            !grav(n) = 0.0
            !fict(n) = (w(n)**2+v(n)**2)/xf(n)
         enddo
      endif
   else if (sweep == 'y') then
      if (yin) then
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius
            rs1n2 = (rad**2 + sep**2 - 2.0*rad*sep*sth*cph)**0.5

            !grav(n) = + GMS / rs1n2**3 * sep*cth*cph
            
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius     &
                      + omega**2*radius*sth*cth                         &
                      - omega**2*rcm*cth*cph                            &
                      + 2.0*omega*cth*v(n)

            grav(n) = 0.0
            !fict(n) = v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      else
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius
            rs1n2 = (rad**2 + sep**2 + 2.0*rad*sep*sth*cph)**0.5

            !grav(n) = - GMS / rs1n2**3 * sep*cth*cph
            
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius     &
                      - omega**2*radius*sph**2*cth*sth                  &
                      + omega**2*rcm*cth*cph                            &
                      - 2.0*omega*(cph*w(n)-sph*sth*v(n))

            grav(n) = 0.0
            !fict(n) = v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      endif
   else
      if (yin) then
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth
            rs1n2 = (rad**2 + sep**2 - 2.0*rad*sep*sth*cph)**0.5

            !grav(n) = - GMS / rs1n2**3 * sep*sph
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth     &
                      + omega**2*rcm*sph                                &
                      - 2.0*omega*(sth*v(n)+cth*w(n))

            grav(n) = 0.0
            !fict(n) = u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth
         enddo
      else
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth
            rs1n2 = (rad**2 + sep**2 + 2.0*rad*sep*sth*cph)**0.5

            !grav(n) = + GMS / rs1n2**3 * sep*sph
            
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth     &
                      - omega**2*radius*cph*sph                         &
                      - omega**2*rcm*sph                                &
                      - 2.0*omega*(sth*sph*w(n)-cth*sph*v(n)) 

            grav(n) = 0.0
            !fict(n) = u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth
         enddo
      endif
   endif
else if (step == 2) then
   r2s = sep - rcm
   sal = sin(omega*time)
   cal = cos(omega*time)

   if (sweep == 'x') then
      if (yin) then
         do n = nmin-4, nmax+5
            rad = xf(n)
            rad2p = (rad**2 + rcm**2 + 2.0*rad*rcm*sth*(cph*cal+sph*sal))**0.5
            rad2s = (rad**2 + r2s**2 - 2.0*rad*r2s*sth*(cph*cal+sph*sal))**0.5

            grav(n) = - GMP / rad2p**3 * (rad+rcm*sth*(cph*cal+sph*sal))    &
                      - GMS / rad2s**3 * (rad-r2s*sth*(cph*cal+sph*sal))    !&
                     ! + LumP * exp(-c(n)) / (4.0*pi*3.0e+10*xf(n)**2)   		! radiation
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                        
         enddo
      else
         do n = nmin-4, nmax+5
            rad = xf(n)
            rad2p = (rad**2 + rcm**2 - 2.0*rad*rcm*(sth*cph*cal-cth*sal))**0.5
            rad2s = (rad**2 + r2s**2 + 2.0*rad*r2s*(sth*cph*cal-cth*sal))**0.5

            grav(n) = - GMP / rad2p**3 * (rad-rcm*(sth*cph*cal-cth*sal))    &
                      - GMS / rad2s**3 * (rad+r2s*(sth*cph*cal-cth*sal))    !&
                     ! + LumP * exp(-c(n)) / (4.0*pi*3.0e+10*xf(n)**2)   		! radiation
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                         
         enddo
      endif
   else if (sweep == 'y') then
      if (yin) then
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius
            rad2p = (rad**2 + rcm**2 + 2.0*rad*rcm*sth*(cph*cal+sph*sal))**0.5
            rad2s = (rad**2 + r2s**2 - 2.0*rad*r2s*sth*(cph*cal+sph*sal))**0.5

            grav(n) = - GMP / rad2p**3 * rcm*cth*(cph*cal+sph*sal)      &
                      + GMS / rad2s**3 * r2s*cth*(cph*cal+sph*sal)
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      else
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius
            rad2p = (rad**2 + rcm**2 - 2.0*rad*rcm*(sth*cph*cal-cth*sal))**0.5
            rad2s = (rad**2 + r2s**2 + 2.0*rad*r2s*(sth*cph*cal-cth*sal))**0.5

            grav(n) = + GMP / rad2p**3 * rcm*(cth*cph*cal+sth*sal)      &
                      - GMS / rad2s**3 * r2s*(cth*cph*cal+sth*sal)
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      endif
   else
      if (yin) then
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth
            rad2p = (rad**2 + rcm**2 + 2.0*rad*rcm*sth*(cph*cal+sph*sal))**0.5
            rad2s = (rad**2 + r2s**2 - 2.0*rad*r2s*sth*(cph*cal+sph*sal))**0.5

            grav(n) = + GMP / rad2p**3 * rcm*(sph*cal-cph*sal)          &
                      - GMS / rad2s**3 * r2s*(sph*cal-cph*sal)
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth
         enddo
      else
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth
            rad2p = (rad**2 + rcm**2 - 2.0*rad*rcm*(sth*cph*cal-cth*sal))**0.5
            rad2s = (rad**2 + r2s**2 + 2.0*rad*r2s*(sth*cph*cal-cth*sal))**0.5

            grav(n) = - GMP / rad2p**3 * rcm*sph*cal         &
                      + GMS / rad2s**3 * r2s*sph*cal
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth 
         enddo
      endif
   endif
else
   if (sweep == 'x') then
      if (yin) then
         do n = nmin-4, nmax+5
            rad = xf(n)

            grav(n) = 0.0
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                        
         enddo
      else
         do n = nmin-4, nmax+5
            rad = xf(n)

            grav(n) = 0.0
            fict(n) = + (w(n)**2+v(n)**2)/xf(n)                         
         enddo
      endif
   else if (sweep == 'y') then
      if (yin) then
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius

            grav(n) = 0.0
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      else
         do n = nmin-4, nmax+5
            sth = sin(xf(n))
            cth = cos(xf(n))
            rad = radius

            grav(n) = 0.0
            fict(n) = + v(n)**2*cth/(radius*sth) - u(n)*w(n)/radius
         enddo
      endif
   else
      if (yin) then
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth

            grav(n) = 0.0
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth
         enddo
      else
         do n = nmin-4, nmax+5
            sph = sin(xf(n))
            cph = cos(xf(n))
            rad = radius / sth

            grav(n) = 0.0
            fict(n) = - u(n)*w(n)/radius*cth - u(n)*v(n)/radius*sth 
         enddo
      endif
   endif
endif

!----------------------------------------------------------------
return
end
