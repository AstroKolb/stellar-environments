subroutine sweepy

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
INTEGER :: i, j, k, m, n, nn, a2abuffer_size, mpierr, myk

!-----------------------------------------------------------------------
! Transpose data

a2abuffer_size = ks * isy * js * 7 
call MPI_ALLTOALL(send1,a2abuffer_size,VH1_DATATYPE,recv1,a2abuffer_size,VH1_DATATYPE, MPI_COMM_ROW, mpierr)  
    
! load boundary data 
call boundaryY

!-----------------------------------------------------------------------

sweep  = 'y'
ngeom  = ngeomy

! Now Loop over each column...

do k = 1, ks
  myk = k + krow*ks
  nmin = jnmin(myk)
  nmax = jnmax(myk)
  sphi = sin(zzc(myk))

  sph  = sin(zzc(myk))
  cph  = cos(zzc(myk))

  do i = 1, isy
   radius = zxc(i+mypey*isy)

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do j = 1, js
    do m = 1, npey
     n = j + js*(m-1) + 6
     r (n) = recv1(1,k,j,i,m)
     p (n) = recv1(2,k,j,i,m)
     u (n) = recv1(4,k,j,i,m)
     v (n) = recv1(5,k,j,i,m)
     w (n) = recv1(3,k,j,i,m)
     f (n) = recv1(6,k,j,i,m)
     c (n) = recv1(7,k,j,i,m)
    enddo
   enddo

   do j = 1, jmax
     n = j + 6
     xa (n) = zya(j)
     dx (n) = zdy(j)
     xa0(n) = zya(j)
     dx0(n) = zdy(j)
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

    r (nmin-n) = ybr(n,k,i)
    p (nmin-n) = ybp(n,k,i)
    u (nmin-n) = ybu(n,k,i)
    v (nmin-n) = ybv(n,k,i)
    w (nmin-n) = ybw(n,k,i)
    f (nmin-n) = ybf(n,k,i)
    c (nmin-n) = ybc(n,k,i)

    r (nmax+n) = ybr(nn,k,i)
    p (nmax+n) = ybp(nn,k,i)
    u (nmax+n) = ybu(nn,k,i)
    v (nmax+n) = ybv(nn,k,i)
    w (nmax+n) = ybw(nn,k,i)
    f (nmax+n) = ybf(nn,k,i)
    c (nmax+n) = ybc(nn,k,i)

    e (nmin-n) = p(nmin-n)/(gamm*r(nmin-n))+0.5d0*(u(nmin-n)**2+v(nmin-n)**2+w(nmin-n)**2)
    e (nmax+n) = p(nmax+n)/(gamm*r(nmax+n))+0.5d0*(u(nmax+n)**2+v(nmax+n)**2+w(nmax+n)**2)
   enddo

   ! 1D Hydrodynamic update using PPMLR
   call ppmlr 

   ! Put updated values into 2D arrays, dropping ghost zones

   do j = 1, jmax
     n = j + 6
     send2(1,k,i,j) = r(n)
     send2(2,k,i,j) = p(n)
     send2(3,k,i,j) = w(n)
     send2(4,k,i,j) = u(n)
     send2(5,k,i,j) = v(n)
     send2(6,k,i,j) = f(n)
     send2(7,k,i,j) = c(n)
   enddo

 enddo
enddo

call MPI_ALLTOALL(send2,a2abuffer_size,VH1_DATATYPE,recv2,a2abuffer_size,VH1_DATATYPE, MPI_COMM_ROW, mpierr)  

return
end

