subroutine sweepz

! This subroutine performs sweeps in the Y direction, looping over i by k rows.   
! Data is read in from RECV after call to MPI_ALLTOALL.
!    DATA is blocked in i.
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweeps

IMPLICIT NONE

include 'mpif.h'

! LOCALS
integer :: i, j, k, n, nn, m, a2abuffer_size, mpierr, myj

!-----------------------------------------------------------------------

do i = 1, isy
 do m = 1, npey
  n = i + isy*(m-1)
  do j = 1, js
   do k = 1, ks
     send3(1,j,k,n) = recv2(1,k,i,j,m)
     send3(2,j,k,n) = recv2(2,k,i,j,m)
     send3(3,j,k,n) = recv2(3,k,i,j,m)
     send3(4,j,k,n) = recv2(4,k,i,j,m)
     send3(5,j,k,n) = recv2(5,k,i,j,m)
     send3(6,j,k,n) = recv2(6,k,i,j,m)
     send3(7,j,k,n) = recv2(7,k,i,j,m)
   enddo
  enddo
 enddo
enddo

a2abuffer_size = ks * isz * js * 7
call MPI_ALLTOALL(send3,a2abuffer_size,VH1_DATATYPE,recv3,a2abuffer_size,VH1_DATATYPE,MPI_COMM_COL,mpierr)  
   
call boundaryZ

!-----------------------------------------------------------------------

sweep  = 'z'
ngeom  = ngeomz

! Now Loop over each column...

do j = 1, js
  myj = j + jcol*js
  nmin = knmin(myj)
  nmax = knmax(myj)
  theta  = zyc(myj)
  stheta = sin(theta)

  sth = sin(theta)
  cth = cos(theta)

  do i = 1, isz
   radius = zxc(i+mypez*isz)
   radius = radius * stheta

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do m = 1, npez
    do k = 1, ks
     n = k + ks*(m-1) + 6
     r(n) = recv3(1,j,k,i,m)
     p(n) = recv3(2,j,k,i,m)
     u(n) = recv3(5,j,k,i,m)
     v(n) = recv3(3,j,k,i,m)
     w(n) = recv3(4,j,k,i,m)
     f(n) = recv3(6,j,k,i,m)
     c(n) = recv3(7,j,k,i,m)
    enddo
   enddo

   do k = 1, kmax
     n = k + 6
     xa (n) = zza(k)
     dx (n) = zdz(k)
     xa0(n) = zza(k)
     dx0(n) = zdz(k)
     e  (n) = p(n)/(r(n)*gamm)+0.5d0*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Set boundary conditions, compute volume elements, evolve flow, then remap
   do n = 1, 6
    nn = n + 6
    
    dx (nmin-n) = dx (nmax+1-n)
    xa (nmin-n) = xa (nmin-n+1) - dx (nmin-n)
    dx0(nmin-n) = dx0(nmax+1-1)
    xa0(nmin-n) = xa0(nmin-n+1) - dx0(nmin-n)
    dx (nmax+n) = dx (nmin+n-1)
    xa (nmax+n) = xa (nmax+n-1) + dx (nmax+n-1)
    dx0(nmax+n) = dx0(nmin+n-1)
    xa0(nmax+n) = xa0(nmax+n-1) + dx0(nmax+n-1)

    r (nmin-n) = zbr(n,j,i)
    p (nmin-n) = zbp(n,j,i)
    u (nmin-n) = zbu(n,j,i)
    v (nmin-n) = zbv(n,j,i)
    w (nmin-n) = zbw(n,j,i)
    f (nmin-n) = zbf(n,j,i)
    c (nmin-n) = zbc(n,j,i)

    r (nmax+n) = zbr(nn,j,i)
    p (nmax+n) = zbp(nn,j,i)
    u (nmax+n) = zbu(nn,j,i)
    v (nmax+n) = zbv(nn,j,i)
    w (nmax+n) = zbw(nn,j,i)
    f (nmax+n) = zbf(nn,j,i)
    c (nmax+n) = zbc(nn,j,i)

    e (nmin-n) = p(nmin-n)/(gamm*r(nmin-n))+0.5d0*(u(nmin-n)**2+v(nmin-n)**2+w(nmin-n)**2)
    e (nmax+n) = p(nmax+n)/(gamm*r(nmax+n))+0.5d0*(u(nmax+n)**2+v(nmax+n)**2+w(nmax+n)**2)
   enddo

   call ppmlr 

   do m = 1, npez
    do k = 1, ks
     n = k + ks*(m-1) + 6
     recv3(1,j,k,i,m) = r(n)
     recv3(2,j,k,i,m) = p(n)
     recv3(5,j,k,i,m) = u(n)
     recv3(3,j,k,i,m) = v(n)
     recv3(4,j,k,i,m) = w(n)
     recv3(6,j,k,i,m) = f(n)
     recv3(7,j,k,i,m) = c(n)
    enddo
   enddo

  enddo
enddo

!######################################################################################
! Repeat a second hydro step

ncycle = ncycle + 1
call boundaryZ

do j = 1, js
  myj = j + jcol*js
  nmin = knmin(myj)
  nmax = knmax(myj)
  theta  = zyc(myj)
  stheta = sin(theta)

  sth = sin(theta)
  cth = cos(theta)

  do i = 1, isz
   radius = zxc(i+mypez*isz)
   radius = radius * stheta

   do m = 1, npez
    do k = 1, ks
     n = k + ks*(m-1) + 6
     r(n) = recv3(1,j,k,i,m)
     p(n) = recv3(2,j,k,i,m)
     u(n) = recv3(5,j,k,i,m)
     v(n) = recv3(3,j,k,i,m)
     w(n) = recv3(4,j,k,i,m)
     f(n) = recv3(6,j,k,i,m)
     c(n) = recv3(7,j,k,i,m)
    enddo
   enddo

   do k = 1, kmax
     n = k + 6
     xa (n) = zza(k)
     dx (n) = zdz(k)
     xa0(n) = zza(k)
     dx0(n) = zdz(k)
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Set boundary conditions, compute volume elements, evolve flow, then remap
   do n = 1, 6
    nn = n + 6
    
    dx (nmin-n) = dx (nmax+1-n)
    xa (nmin-n) = xa (nmin-n+1) - dx (nmin-n)
    dx0(nmin-n) = dx0(nmax+1-1)
    xa0(nmin-n) = xa0(nmin-n+1) - dx0(nmin-n)
    dx (nmax+n) = dx (nmin+n-1)
    xa (nmax+n) = xa (nmax+n-1) + dx (nmax+n-1)
    dx0(nmax+n) = dx0(nmin+n-1)
    xa0(nmax+n) = xa0(nmax+n-1) + dx0(nmax+n-1)

    r (nmin-n) = zbr(n,j,i)
    p (nmin-n) = zbp(n,j,i)
    u (nmin-n) = zbu(n,j,i)
    v (nmin-n) = zbv(n,j,i)
    w (nmin-n) = zbw(n,j,i)
    f (nmin-n) = zbf(n,j,i)
    c (nmin-n) = zbc(n,j,i)

    r (nmax+n) = zbr(nn,j,i)
    p (nmax+n) = zbp(nn,j,i)
    u (nmax+n) = zbu(nn,j,i)
    v (nmax+n) = zbv(nn,j,i)
    w (nmax+n) = zbw(nn,j,i)
    f (nmax+n) = zbf(nn,j,i)
    c (nmax+n) = zbc(nn,j,i)

    e (nmin-n) = p(nmin-n)/(gamm*r(nmin-n))+0.5d0*(u(nmin-n)**2+v(nmin-n)**2+w(nmin-n)**2)
    e (nmax+n) = p(nmax+n)/(gamm*r(nmax+n))+0.5d0*(u(nmax+n)**2+v(nmax+n)**2+w(nmax+n)**2)
   enddo

   call ppmlr 

   ! Put updated values into 2D arrays, dropping ghost zones

   do k = 1, kmax
     n = k + 6
     send4(1,j,i,k) = r(n)
     send4(2,j,i,k) = p(n)
     send4(3,j,i,k) = v(n)
     send4(4,j,i,k) = w(n)
     send4(5,j,i,k) = u(n)
     send4(6,j,i,k) = f(n)
     send4(7,j,i,k) = c(n)
   enddo

  enddo
enddo

call MPI_ALLTOALL(send4,a2abuffer_size,VH1_DATATYPE,recv4,a2abuffer_size,VH1_DATATYPE,MPI_COMM_COL,mpierr)  

! Load data from sweepz recv buffer into sweepy send buffer for next sweepy
do i = 1, isz
 do m = 1, npez
  n = i + isz*(m-1)
  do k = 1, ks
   do j = 1, js
      send1(1,k,j,n) = recv4(1,j,i,k,m)
      send1(2,k,j,n) = recv4(2,j,i,k,m)
      send1(3,k,j,n) = recv4(3,j,i,k,m)
      send1(4,k,j,n) = recv4(4,j,i,k,m)
      send1(5,k,j,n) = recv4(5,j,i,k,m)
      send1(6,k,j,n) = recv4(6,j,i,k,m)
      send1(7,k,j,n) = recv4(7,j,i,k,m)
   enddo
  enddo
 enddo
enddo

return
end

