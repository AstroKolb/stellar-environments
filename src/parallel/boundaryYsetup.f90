subroutine boundaryYsetup
! Set up arrays for passing boundary data and interpolating between Yin and Yang grids.
! Written by Michael Owen (2007), John Blondin (2010)
!--------------------------------------------------------------------------------------

use global
use zone
use yinyang

IMPLICIT NONE

include 'mpif.h'

! LOCALS
INTEGER :: i,j,k,n, m, md, myk, jjmin, jjmax, mpierr, ns, ms, ktarget, il, klocal, mdmax, nvar
INTEGER, DIMENSION(20000) :: jbndry, kbndry, nbndry, js0, ks0, targetpe, jyt, kyt, nyt
REAL :: cphi, ctheta, x, y, z, cthetan, cphin, cosphi, sinphi, costhe, sinthe

!---------------------------------------------------------------------------------------
! Loop over my boundary zones and find indicies and coefficients for interpolation
!
!     neighbor 1 -> (m,n)
!     neighbor 2 -> (m,ns)
!     neighbor 3 -> (ms,n)
!     neighbor 4 -> (ms,ns)
!     interpolated value = c3*( c1*N1 + c2*N2 ) + c4*( c1*N3 + c2*N4 )
!     theta velocity = -c5*ut - c6*vt
!     phi velocity   = +c6*ut - c5*vt
!
! Output arrays dimension(j=12,k=ks) :: cy1, cy2, cy3, cy4, cy5, cy6
!---------------------------------------------------------------------------------------
nvar = 7  ! number of variables used in boundary conditions

do k = 1, ks
 myk    = krow*ks + k
 cphi   = zzc(myk)
 cosphi = cos(cphi)
 sinphi = sin(cphi)
 jjmin  = jnmin(myk) - 6
 jjmax  = jnmax(myk) - 6
  
 do j = 1, 12  ! loop over the 12 ghost cells in the j direction (6 left, 6 right)

  ! find the theta value for the ghost cell and convert to (x,y,z)
  if (j<7) then
   ctheta = zyc(jjmin) - real(j)*zdy(jjmin)
  else
   ctheta = zyc(jjmax) + real(j-6)*zdy(jjmax)
  endif
  
  sinthe = sin(ctheta)
  costhe = cos(ctheta)
  x = -sinthe*cosphi
  z =  sinthe*sinphi
  y =  costhe

  ! find the theta and phi of this point on the 'other' grid
  cthetan = acos(z)
  cphin = atan2(y,x)
    
  ! find the indicies of nearest zones on the 'other' grid
  m = int((cthetan - zya(1))/zdy(1)) + 1
  n = int((cphin   - zza(1))/zdz(1)) + 1

  cy2(j,k) = (cphin-zzc(n))/zdz(kmax)
  if (cy2(j,k) < 0.0) then
   cy2(j,k) = -cy2(j,k)
  endif

  cy1(j,k) = 1.0 - cy2(j,k)

  cy4(j,k) = (cthetan-zyc(m))/zdy(jmax)
  if(cy4(j,k) < 0.0) then
   cy4(j,k) = -cy4(j,k)
  endif

  cy3(j,k) = 1.0 - cy4(j,k)
  
  cy5(j,k) = sinphi*sin(cphin)
  cy6(j,k) = cos(cphin)/sinthe
  
 enddo
enddo
 
   
!---------------------------------------------------------------------------------------
! Loop over full 2D slice and see if I have any zones needed for interpolation
! Output:   md               : total number of blocks I have to send out
!           jytarget(md)     : index on target PE of ghost zone (1-6,7-12)
!           kytarget(md)     : index on target PE of ghost zone (1-ks)
!           nytarget(md)     : index on target PE of ghost zone neighbors (1-4) used in interpolation
!           jysource(md)     : j index (1-jmax) of zone on my grid that I need to send
!           kysource(md)     : k index (1-ks) of zone on my grid that I need to send
!           sycount (npez)   : array of block counts I am sending to other PEs
!           rycount (npez)   : array of block counts I am receiving from other PEs
!           syoffset(npez)   : displacement vector used in ALLTOALLV
!           ryoffset(npez)   : displacement vector used in ALLTOALLV
!---------------------------------------------------------------------------------------

md = 0
do k = 1, kmax
 ktarget = (k-1)/ks + 1
 klocal  = k - (ktarget-1)*ks
 cphi   = zzc(k)
 cosphi = cos(cphi)
 sinphi = sin(cphi)
 jjmin = jnmin(k) - 6
 jjmax = jnmax(k) - 6
  
 do j = 1, 12  ! loop over the 12 ghost cell in the j directions

  ! find the theta value for the ghost cell and convert to (x,y,z)
  if (j<7) then
   ctheta = zyc(jjmin) - real(j)*zdy(jjmin)
  else
   ctheta = zyc(jjmax) + real(j-6)*zdy(jjmax)
  endif
  
  sinthe = sin(ctheta)
  costhe = cos(ctheta)
  x = -sinthe*cosphi
  z =  sinthe*sinphi
  y =  costhe

  ! find the theta and phi of this point on the 'other' grid
  cthetan = acos(z)
  cphin = atan2(y,x)
    
  ! find the k indicies of nearest zone on the 'other' grid
  n = int((cphin - zza(1))/zdz(1)) + 1
  ns = n + 1
  if ((cphin-zzc(n)) < 0.0) ns = n - 1

  ! are either n or ns on my grid?
  if ((n-1)/ks == krow) then
    m = int((cthetan - zya(1))/zdy(1)) + 1
    ms = m + 1
    if((cthetan-zyc(m)) < 0.0) ms = m - 1
    
    md = md + 1
    jbndry(md) = j
    kbndry(md) = klocal
    nbndry(md) = 1  ! neighbor = 1 if m,n; 2 if m,ns; 3 if ms,n; 4 if ms,ns
    js0(md) = m 
    ks0(md) = n - ks*krow  ! need local 'k'
    targetpe(md) = ktarget
    
    md = md + 1
    jbndry(md) = j
    kbndry(md) = klocal
    nbndry(md) = 3  ! neighbor = 1 if m,n; 2 if m,ns; 3 if ms,n; 4 if ms,ns
    js0(md) = ms 
    ks0(md) = n - ks*krow  ! need local 'k'
    targetpe(md) = ktarget
  endif

  if ((ns-1)/ks == krow) then
    m = int((cthetan - zya(1))/zdy(1)) + 1
    ms = m + 1
    if((cthetan-zyc(m)) < 0.0) ms = m - 1
    
    md = md + 1
    jbndry(md) = j
    kbndry(md) = klocal
    nbndry(md) = 2  ! neighbor = 1 if m,n; 2 if m,ns; 3 if ms,n; 4 if ms,ns
    js0(md) = m 
    ks0(md) = ns - ks*krow  ! need local 'k'
    targetpe(md) = ktarget
    
    md = md + 1
    jbndry(md) = j
    kbndry(md) = klocal
    nbndry(md) = 4  ! neighbor = 1 if m,n; 2 if m,ns; 3 if ms,n; 4 if ms,ns
    js0(md) = ms 
    ks0(md) = ns - ks*krow  ! need local 'k'
    targetpe(md) = ktarget
  endif
    
 enddo
enddo
    
! used to make sure source arrays are dimensioned large enough
!call MPI_REDUCE( md, mdmax, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, mpierr )
!if (mype==0) write(8,*) 'maximum mdy = ', mdmax

!---------------------------------------------------------------------------------------
! repack boundary data into send buffer sorted by target PE
! md is total number of blocks this PE is sending out

mdy = md
ALLOCATE (sybuff(nvar,isy,mdy))
ALLOCATE (sybufg(nvar,isy,mdy))
syc_tot = nvar*isy*mdy

sycount = 0  
il = 0  
do n = 1, npez
 do m = 1, mdy
 
  if (targetpe(m) == n) then ! this block is going to this PE
    il = il + 1
    jyt(il) = jbndry(m)
    kyt(il) = kbndry(m)
    nyt(il) = nbndry(m)
    sycount (n) = sycount(n) + 1
    jysource(il) = js0(m)
    kysource(il) = ks0(m)
  endif  ! should end with il = md
  
 enddo
enddo

CALL MPI_ALLTOALL(sycount,1,MPI_INTEGER,rycount,1,MPI_INTEGER,MPI_COMM_COL,mpierr)

if ( sum(rycount) /= 48*ks ) print *, 'rycount ', sum(rycount), ks, mype

syoffset = 0
ryoffset = 0
do n = 2, npez
 syoffset(n) = syoffset(n-1) + sycount(n-1)
 ryoffset(n) = ryoffset(n-1) + rycount(n-1)
enddo

call MPI_ALLTOALLV(jyt,sycount,syoffset,MPI_INTEGER,jytarget,rycount,ryoffset,MPI_INTEGER,MPI_COMM_COL,mpierr)
call MPI_ALLTOALLV(kyt,sycount,syoffset,MPI_INTEGER,kytarget,rycount,ryoffset,MPI_INTEGER,MPI_COMM_COL,mpierr)
call MPI_ALLTOALLV(nyt,sycount,syoffset,MPI_INTEGER,nytarget,rycount,ryoffset,MPI_INTEGER,MPI_COMM_COL,mpierr)

! now redo indexing for variable buffer
sycount = sycount * nvar *isy
rycount = rycount * nvar *isy
syoffset = 0
ryoffset = 0
do n = 2, npez
 syoffset(n) = syoffset(n-1) + sycount(n-1)
 ryoffset(n) = ryoffset(n-1) + rycount(n-1)
enddo

return
end

