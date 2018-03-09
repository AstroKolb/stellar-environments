subroutine init

! Ejecta-Driven SNR with exponential density profile
! 07oct11 blondin
!=======================================================================

! GLOBALS
use global
use zone
use noise
use netcdf

include 'mpif.h'

! LOCALS
integer :: i, j, k, n, mpierr
real :: ridt, xvel, yvel, zvel, width, widthz, widthy
real :: zoom, rad, vthe, vphi, theta, phi, dshock, pshock, ushock
real :: density, pressure, radius, perturb, port, bump, costheta

real :: cx, cy, cz, tmpu, RHS
real :: gen_noise

real :: dely, delz

real :: xmin,xmax,ymin,ymax,zmin,zmax
real :: c1, c2, rodt


integer :: nabove, nbelow
real :: ialpha, Esn, Mej, d0, age, vt, capA, Rsh
real :: d_plateau, r_plateau, rscale, dscale, pscale, uscale
real :: rad_max, den_max, pre_max

real, dimension(2048) :: rss, dss, pss, uss, css



character(LEN=1) :: tmp3
character(LEN=4) :: tmp1, tmp2
character(LEN=50) :: filename
real, dimension(imax,jmax,kmax/pez) :: prin_buff

real, dimension(imax,6,js,ks,npey) :: recv_buff
real, dimension(imax,6,js,ks) :: send_buff
integer :: jj, nv, gathbuffer_size

!--------------------------------------------------------------------------------
! Set up grid geometry: This Yin-Yang code is strictly 3D spherical

ndim   = 3
ngeomx = 2 
ngeomy = 4
ngeomz = 5

pi   = 2.0 * asin(1.0)

dely =     pi/2.0/(jmax-2)/2.0
delz = 3.0*pi/2.0/(kmax-2)/2.0

ymin = +1.0*pi/4.0 - dely
ymax = +3.0*pi/4.0 + dely
zmin = -3.0*pi/4.0 - delz
zmax = +3.0*pi/4.0 + delz

! set time and cycle counters
time   = 0.0
timep  = 0.
timem  = 0.
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 1000



!=======================================================================
! Log parameters of problem in history file

if (mype == 0) then
   write(8,*) 'Stellar Environment Code'
   write(8,"(' Step = ', i1)") step
   write(8,"(' Grid dimensions: ', i4, ' x ', i4, ' x ', i4)") imax, jmax, kmax
   write(8,*)
   write(8,"(' Adiabatic index = ', f8.6)") gam
   write(8,"(' Mass ratio = ', f4.2, ' x ', f4.2)") GMP, GMS
   write(8,"(' Binary Separation = ', e6.3E2, a)") sep, scale
   write(8,*)
   write(8,"(' Inner grid dimensions: ', e6.3e2, a, ' to ', e6.3e2, a)") x1min, scale, x1max, scale
   write(8,"(' Outer grid dimensions: ', e6.3e2, a, ' to ', e6.3e2, a)") x1max, scale, x2max, scale
endif

!======================================================================
! Set up parameters for rotation problem

if (step == 3) gam = 5.0/3.0
gamm = gam - 1.0

GMP = GMP * GM
GMS = GMS * GM

sep = sep * Rsun
x0min = x0min * Rsun
x1min = x1min * Rsun
x2min = x1max
x1max = x1max * sep
x2max = x2max * sep + rcm

radP  = radP  * Rsun
rad1  = rad1  * radP
rad2  = rad2  * radP
rad3  = rad3  * radP


! binary properties
rcm   = sep*GMS/(GMP+GMS)
omega = ((GMP+GMS)/sep**3)**0.5
opd   = 2.0*pi/omega
P_flr = gam * kB * T_flr / mp


Z1 = (x1max/x1min)**(1.0/float(imax)) - 1.0
Z2 = (x2max/x2min)**(1.0/float(imax)) - 1.0

do while (x2min*(1.0+Z2/2.0) > x1max*(1.0+Z1/2.0)/(1.0+Z1)-rcm)
   x2min = x2min * 0.999
   Z2 = (x2max/x2min)**(1.0/float(imax)) - 1.0
enddo


! calculate some other stuff
Cs2P = gam * kB * TmpP / mp
ucP  = Cs2P**0.5
rcP  = GMP / (2.0*Cs2P)
rhoP = (dmP*Msun/3.15e7)/(4.0*pi*rcP**2*ucP)
IIP  = rhoP*ucP*rcP**2


!=======================================================================
! Set up grid coordinates and parabolic coefficients

if (step == 1) then
   xmin = x1min
   xmax = x1max
   zoom = Z1
else if (step == 2) then
   xmin = x2min
   xmax = x2max
   zoom = Z2
else
   xmin = x1min * 0.40! * 4
   xmax = x1min * 1.75! * 4
   zoom = (xmax/xmin)**(1.0/float(imax)) - 1.0
endif


call grid(imax,xmin,xmax,zxa,zxc,zdx,zoom)
call grid(jmax,ymin,ymax,zya,zyc,zdy,0.00)
call grid(kmax,zmin,zmax,zza,zzc,zdz,0.00)

!=======================================================================
! initialize grid:

!call initialize_noise     ! unused, rand() seems to cause an issue on stampede
call boundaryRsetup
call boundaryR

if (step == 3) then
   nabove = 2
   do while (alpha(nabove) < 0.5 * sep)
      nabove = nabove + 1
   enddo

   !write(*,*) myray(1,nabove,1,1)

   ! add remnant parameters
   nssdw = 12
   ialpha = 1.048

   Esn = 1.0e+51
   Mej = 8.0e+33
   d0  = 1.0e-15   ! print this..
   age = 2.5e-2 * 3.15e+07

   vt   = sqrt((10.*nssdw-50.)/(3.*nssdw-9.)*Esn/Mej)
   capA = (5.*nssdw-25.)/(2.*pi*nssdw)*Esn/vt**5.0
   Rsh  = ialpha * (capA * vt**nssdw/d0)**(1.0/nssdw)*age**((nssdw-3.0)/nssdw)

   d_plateau = capA / age**3
   r_plateau = vt * age
endif

!----------------------------------------------------------------------!
! read in 1D data for SSDW; scaling to above parameters
! input data is scaled to a radius of 1.0, age of 0.5, ambient density of 1.0

if (step == 1) then
   !open(27, file='src/assets/wind16.dat')
   open(27, file='src/assets/mywind.dat')
   do n = 1, 2048
      read(27,*) rss(n), dss(n), pss(n), uss(n), css(n)
   enddo
   close(27)
endif

if (mype==0) write(*,*) 'asdf', rss(1:5)/zxc(1)

if (step == 3) then
   open(27,file='x060c')
   do n = 1, 2048 
      read(27,*) rss(n), dss(n), pss(n), uss(n)
   enddo
   close(27)

   n = 2048
   do while (dss(n) < 1.1)
      n = n - 1
   enddo
   dscale = d0
   rscale = Rsh / rss(n)
   uscale = (rss(1)*Rsh/rss(n)) / age / uss(1)  ! force uss(1) = rss(1)/age
   pscale = dscale * uscale**2
   time   = age

   dinflo = d_plateau
   pinflo = pss(1) / dss(1) / uss(1)**2

   do n = 1, 2048
   rss(n) = rss(n) * rscale
   dss(n) = dss(n) * dscale
   pss(n) = pss(n) * pscale
   uss(n) = uss(n) * uscale
   enddo

   rad_max = rss(2048)
   den_max = dss(2048)
   pre_max = pss(2048)
   assdw   = dss(1) * uss(1)**nssdw
   if (mype == 0) write(*,*) 'blastwave initialized at R=', rad_max/sep, 'sep'
endif

!==============================================================================

do i = 1, imax
   if (zxc(i) < rad1) then
      kap(i) = kap1
   else if (zxc(i) < rad2) then
      kap(i) = (kap2-kap1)/(rad2-rad1)*(zxc(i)-rad1) + kap1
   else if (zxc(i) < rad3) then
      kap(i) = (kap3-kap2)/(rad3-rad2)*(zxc(i)-rad2) + kap2
   else
      kap(i) = kap3
   endif
enddo


!  set min rad
i = 1
do while (zxc(i) < sep)
   i = i + 1
enddo
minrad = zdx(i)/2.0


if (step == 1) then

   ! parker wind
   if (.false.) then
      do i = 1, imax
         RHS = 4.0*log(zxc(i)/rcP) + 4.0*rcP/zxc(i) - 3.0

         tmpu = ucP
         do while (tmpu**2/ucP**2 - log(tmpu**2/ucP**2) < RHS)
            if (zxc(i) < rcP) tmpu = 0.9999*tmpu
            if (zxc(i) > rcP) tmpu = 1.0001*tmpu
         enddo
         do k = 1, ks
            do j = 1, js
               zux(i,j,k) = tmpu
               zuy(i,j,k) = 0.0
               zuz(i,j,k) = 0.0

               zro(i,j,k) = IIP / (zux(i,j,k) * zxc(i)**2)
               zpr(i,j,k) = zro(i,j,k)*kB*TmpP/mp
               zfl(i,j,k) = 0.0
               zcl(i,j,k) = 0.0
            enddo
         enddo
      enddo


      ! determine inflow boundary conditions
      do i = 1, 6
         rad = zxc(1) - i*zdx(1)

         RHS = 4.0 * (log(rad/rcP) + rcP/rad - 1.0)

         tmpu = ucP
         do while (tmpu**2/ucP**2 - log(tmpu**2/ucP**2) - 1.0 < RHS)
            tmpu = 0.9999 * tmpu
         enddo

         uin(i,1,1) = tmpu
      enddo
   endif


   ! write from file
   if (.true.) then
      do i = 1, imax
         nabove = 2
         do while (rss(nabove) < zxc(i) .and. nabove < 1024)
            nabove = nabove + 1
         enddo
         nbelow = nabove - 1

         c1 = (rss(nabove)-zxc(i))/(rss(nabove)-rss(nbelow))
         c2 = (zxc(i)-rss(nbelow))/(rss(nabove)-rss(nbelow))

         if (nabove < 1024) then
            do k = 1, ks
               do j = 1, js
                  zux(i,j,k) = c1*uss(nbelow) + c2*uss(nabove)
                  zuy(i,j,k) = 0.0
                  zuz(i,j,k) = 0.0

                  zro(i,j,k) = c1*dss(nbelow) + c2*dss(nabove)
                  zpr(i,j,k) = c1*pss(nbelow) + c2*pss(nabove)
                  zcl(i,j,k) = 0.0
               enddo
            enddo
         else
            do k = 1, ks
               do j = 1, js
                  zux(i,j,k) = zux(i-1,j,k)
                  zuy(i,j,k) = zuy(i-1,j,k)
                  zuz(i,j,k) = zuz(i-1,j,k)

                  zro(i,j,k) = zro(i-1,j,k)
                  zpr(i,j,k) = zpr(i-1,j,k)
                  zcl(i,j,k) = 0.0
               enddo
            enddo
         endif
      enddo

      do i = 1, 6
         rad = zxc(1) - i*zdx(1)
         nabove = 2
         do while (rss(nabove) < rad .and. nabove < 2048)
            nabove = nabove + 1
         enddo
         nbelow = nabove - 1

         c1 = (rss(nabove)-rad)/(rss(nabove)-rss(nbelow))
         c2 = (rad-rss(nbelow))/(rss(nabove)-rss(nbelow))
if (mype==0) write(*,*) zdx(1)/zxc(1), rss(1)/rad, rad
         uin(i,1,1) = c1*uss(nbelow) + c2*uss(nabove)
         rin(i,1,1) = c1*dss(nbelow) + c2*dss(nabove)
         pin(i,1,1) = c1*pss(nbelow) + c2*pss(nabove)
      enddo
      if (mype == 0) write(*,*) 'uin:', uin(:,1,1)
      if (mype == 0) write(*,*) 'rin:', rin(:,1,1)
      if (mype == 0) write(*,*) 'pin:', pin(:,1,1)
   endif

else if (step == 2) then ! ------------------------------------------
   do k = 1, ks
      do j = 1, js
         do i = 1, imax
            zux(i,j,k) = uin(1,j,k)! * 10.0
            zuy(i,j,k) = vin(1,j,k)
            zuz(i,j,k) = win(1,j,k)

            zro(i,j,k) = rin(1,j,k)
            zpr(i,j,k) = pin(1,j,k)
            zcl(i,j,k) = cin(1,j,k)
         enddo
      enddo
   enddo
else ! --------------------------------------------------------------

   ! view ..
   if (.false.) then
      do k = 1, ks
         do j = 1, js
            do i = 1, imax
               zro(i,j,k) = myray(1,i,j,k)
               zpr(i,j,k) = myray(2,i,j,k)
               zcl(i,j,k) = myray(3,i,j,k)

               zux(i,j,k) = myray(4,i,j,k)
               zuy(i,j,k) = myray(5,i,j,k)
               zuz(i,j,k) = myray(6,i,j,k)
            enddo
         enddo
      enddo
   endif

   ! add blastwave
   if (.true.) then
      do k = 1, ks
         do j = 1, js
            do i = 1, imax
               if (zxc(i) < rss(1)) then ! plateau
                  zux(i,j,k) = zxc(i) / time
                  zuy(i,j,k) = 0.0
                  zuz(i,j,k) = 0.0

                  zro(i,j,k) = assdw / zux(i,j,k)**nssdw
                  zpr(i,j,k) = pinflo * zro(i,j,k) * zux(i,j,k)**2
                  zcl(i,j,k) = 0.0
               else if (zxc(i) > rad_max) then   ! ism
                  nabove = 2
                  do while (alpha(nabove) < zxc(i))
                     nabove = nabove + 1
                  enddo
                  nbelow = nabove - 1
                  
                  c1 = (alpha(nabove)-zxc(i)) / (alpha(nabove)-alpha(nbelow))
                  c2 = (zxc(i)-alpha(nbelow)) / (alpha(nabove)-alpha(nbelow))
                  !zro(i,j,k) = d0*(x1min/zxc(i))**2
                  !zpr(i,j,k) = 1e-10
                  zro(i,j,k) = myray(1,nbelow,j,k)*c1 + myray(1,nabove,j,k)*c2
                  zpr(i,j,k) = myray(2,nbelow,j,k)*c1 + myray(2,nabove,j,k)*c2
                  zcl(i,j,k) = 0.0
                  
                  zux(i,j,k) = 0.!myray(4,nbelow,j,k)*c1 + myray(4,nabove,j,k)*c2
                  zuy(i,j,k) = 0.!myray(5,nbelow,j,k)*c1 + myray(5,nabove,j,k)*c2
                  zuz(i,j,k) = 0.!myray(6,nbelow,j,k)*c1 + myray(6,nabove,j,k)*c2
               else  ! snr
                  nabove = 2
                  do while (rss(nabove) < zxc(i))
                     nabove = nabove + 1
                  enddo
                  nbelow = nabove - 1

                  c1 = (rss(nabove)-zxc(i)) / (rss(nabove)-rss(nbelow))
                  c2 = (zxc(i)-rss(nbelow)) / (rss(nabove)-rss(nbelow))

                  zro(i,j,k) = dss(nbelow)*c1 + dss(nabove)*c2
                  zpr(i,j,k) = pss(nbelow)*c1 + pss(nabove)*c2
                  zcl(i,j,k) = 0.0

                  zux(i,j,k) = uss(nbelow)*c1 + uss(nabove)*c2
                  zuy(i,j,k) = 0.0
                  zuz(i,j,k) = 0.0
               endif
            enddo
         enddo
      enddo
   endif
endif


!=======================================================================
! Log parameters of problem in history file

if (zxc(imax) > sep) then
  i = 1
  do while (zxc(i) < sep)
    i = i + 1
  enddo
else
  i = imax
endif
tmpu = zux(i,1,1)

! finish writing parameters to grid
if (mype == 0) then
  write (8,*)
  write (8,"('Orbital period:  ', f7.2, ' years')") opd/3.15e7
  write (8,"('Vorbit to Vwind: ', f5.3          )") omega*sep / tmpu
  write (8,*) 
endif

! print variables
if (mype == 0) then
  write(*,*) 'separation (Rsun):',      sep/Rsun
  write(*,*) 'angular velocity (1/s):', omega
  write(*,*) 'orbital period (yr):',    opd/3.15e7
  write(*,*) 'sound speed (cm/s):',     ucP
  write(*,*) 'V_orbit to V_wind:',      omega*sep / tmpu
  write(*,*) 'R_cm:',                   rcm
endif


!########################################################################
! Compute Courant-limited timestep

ridt = 0.
do k = 1, ks
 do j = 1, js
  do i = 1, imax
     widthy = zdy(j+jcol*js)
     widthz = zdz(k+krow*ks)
     widthy = widthy*zxc(i)
     widthz = widthz*zxc(i)*sin(zyc(j+jcol*js))
     width  = min(zdx(i),widthy,widthz)
     svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
     xvel = abs(zux(i,j,k)) / zdx(i)
     yvel = abs(zuy(i,j,k)) / widthy
     zvel = abs(zuz(i,j,k)) / widthz
     ridt = max(xvel,yvel,zvel,svel,ridt)
  enddo
 enddo
enddo

call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
dt = courant / rodt

return
end

!#######################################################################


subroutine grid( nzones, xmin, xmax, xa, xc, dx, zoom )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax, zoom

!=======================================================================

if (zoom==0.0d0) then
  dxfac = (xmax - xmin) / float(nzones)
  do n = 1, nzones
    xa(n) = xmin + (n-1)*dxfac
    dx(n) = dxfac
    xc(n) = xa(n) + 0.5*dx(n)
  enddo
else
  xa(1) = xmin
  dx(1) = zoom*xa(1)
  xc(1) = xa(1) + 0.5d0*dx(1)
  do n  = 2, nzones
    xa(n) = xa(n-1)+dx(n-1)
    dx(n) = zoom*xa(n)
    xc(n) = xa(n)+0.5d0*dx(n)
  enddo
endif

return
end
