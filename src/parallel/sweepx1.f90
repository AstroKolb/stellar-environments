subroutine sweepx1(xexpand)

! This subroutine performs 1D hydro sweeps in the X direction, looping over j by k rows.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

! LOCALS
integer :: i, j, k, n, myj, myk
real :: rad, xexpand, costhetan, dis

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nmin   = 7
nmax   = imax + 6
radius = 1.0

do i = 1, imax
  n = i + 6
  xa0(n) = zxa(i) * (1.0 + xexpand)
  dx0(n) = zdx(i) * (1.0 + xexpand)
enddo

call boundaryR

! Now Loop over each row...
do k = 1, ks
 myk  = k+krow*ks
 cphi = cos(zzc(myk))
 sphi = sin(zzc(myk))

 sph  = sin(zzc(myk))
 cph  = cos(zzc(myk))

 do j = 1, js
   myj    = j+jcol*js
   ctheta = cos(zyc(myj))
   stheta = sin(zyc(myj))

   sth    = sin(zyc(myj))
   cth    = cos(zyc(myj))

   if (yin) then
     costhetan = ctheta
   else
     costhetan = stheta*sphi
   endif

   c = 0.0

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do i = 1,imax
     n = i + 6
     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zux(i,j,k)
     v  (n) = zuy(i,j,k)
     w  (n) = zuz(i,j,k)
     f  (n) = zfl(i,j,k)
     c  (n) = 1.0-kap(i)

     xa (n) = zxa(i)
     dx (n) = zdx(i)
     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5d0*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Set boundary conditions, compute volume elements, evolve flow, then remap

   do n = 1, 6
    dx (nmin-n)= dx (nmin)
    dx0(nmin-n)= dx0(nmin)
    xa(nmin-n) = xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)

    rad = xa(nmin-n)+0.5d0*dx(nmin-n)

    if (step == 1) then
      u (nmin-n) = uin(n,1,1)
      v (nmin-n) = 0.0
      w (nmin-n) = 0.0

      r (nmin-n) = rin(n,1,1)
      p (nmin-n) = pin(n,1,1)
      c (nmin-n) = c(nmin)
    else if (step == 2) then
      u (nmin-n) = uin(n,j,k)
      v (nmin-n) = vin(n,j,k)
      w (nmin-n) = win(n,j,k)

      r (nmin-n) = rin(n,j,k)
      p (nmin-n) = pin(n,j,k)
      c (nmin-n) = cin(n,j,k)
    else
      u (nmin-n) = rad / time
      v (nmin-n) = 0.0
      w (nmin-n) = 0.0

      r (nmin-n) = assdw / u(nmin-n)**nssdw
      p (nmin-n) = pinflo * r(nmin-n) * u(nmin-n)**2
      c (nmin-n) = 1.0
    endif


    dx (nmax+n)= dx (nmax)
    dx0(nmax+n)= dx0(nmax)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)

    rad = xa(nmax+n) + 0.5d0*dx(nmax+n)

    if (step == 1) then
      u (nmax+n) = u(nmax)
      v (nmax+n) = v(nmax)
      w (nmax+n) = w(nmax)
      
      r (nmax+n) = r(nmax)
      p (nmax+n) = p(nmax)
      c (nmax+n) = c(nmax)
    else if (step == 2) then
      u (nmax+n) = u(nmax)
      v (nmax+n) = v(nmax)
      w (nmax+n) = w(nmax)

      r (nmax+n) = r(nmax)
      p (nmax+n) = p(nmax)
      c (nmax+n) = 0.0
    else
      u (nmax+n) = 0.0!uin(n,j,k)
      v (nmax+n) = 0.0!vin(n,j,k)
      w (nmax+n) = 0.0!win(n,j,k)

      r (nmax+n) = rin(n,j,k)
      p (nmax+n) = pin(n,j,k)
      c (nmax+n) = 0.0
    endif

    f (nmin-n) = 0.0
    f (nmax+n) = 0.0

    e (nmin-n) = p(nmin-n)/(gamm*r(nmin-n))+0.5d0*(u(nmin-n)**2+v(nmin-n)**2+w(nmin-n)**2)
    e (nmax+n) = p(nmax+n)/(gamm*r(nmax+n))+0.5d0*(u(nmax+n)**2+v(nmax+n)**2+w(nmax+n)**2)
   enddo

   if (step == 0) then
    do i = 1, 6
      n = i + 6
      u(n) = uin(n,1,1)
      r(n) = rin(n,1,1)
      p(n) = pin(n,1,1)
    enddo
   endif

   if (ongrid(myj,myk) == 0) call ppmlr 

   ! Put updated values back into 3D arrays, dropping ghost zones
   do i = 1, imax
     n = i + 6
     send1(1,k,j,i) = r(n)
     send1(2,k,j,i) = p(n)
     send1(3,k,j,i) = u(n)
     send1(4,k,j,i) = v(n)
     send1(5,k,j,i) = w(n)
     send1(6,k,j,i) = f(n)
     send1(7,k,j,i) = c(n)
   enddo

 enddo
enddo

do i = 1, imax
  n = i + 6
  zxa(i) = xa0(n)
  zdx(i) = dx0(n)
  zxc(i) = xa0(n) + 0.5*dx0(n)
enddo
 
return
end

