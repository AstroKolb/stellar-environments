subroutine boundaryZsetup

use global
use zone
use yinyang

IMPLICIT NONE

include 'mpif.h'

! LOCALS
integer :: i,j,k,n, m, md, il, l, ll, nn, mm, myk, myj, ns, ms
integer :: nvar, jtarget, nv, kkmin, kkmax, mpierr, jlocal, mdmax
INTEGER, DIMENSION(20000) :: jbndry, kbndry, nbndry, js0, ks0, targetpe, jzt, kzt, nzt
REAL :: cphi, ctheta, x, y, z, cthetan, cphin, cosphi, sinphi, costhe, sinthe

!---------------------------------------------------------------------------------------
! Loop over my boundary zones and find indicies and coefficients for interpolation
!
!     neighbor 1 -> (m,n)
!     neighbor 2 -> (ms,n)
!     neighbor 3 -> (m,ns)
!     neighbor 4 -> (ms,ns)
!     interpolated value = c3*( c1*N1 + c2*N2 ) + c4*( c1*N3 + c2*N4 )
!     theta velocity = -c5*wt - c6*ut
!     phi velocity   = +c6*wt - c5*ut
!
! Output arrays dimension(j=js,k=12) :: cz1, cz2, cz3, cz4, cz5, cz6
!---------------------------------------------------------------------------------------
nvar = 7  ! number of variables used in boundary conditions
 
do j = 1, js
 myj = jcol*js + j
 ctheta = zyc(myj)
 sinthe = sin(ctheta)
 costhe = cos(ctheta)
 kkmin  = knmin(myj) - 6
 kkmax  = knmax(myj) - 6
  
 do k = 1, 12

  if (k<7) then
   cphi = zzc(kkmin) - real(k)*zdz(kkmin)
  else
   cphi = zzc(kkmax) + real(k-6)*zdz(kkmax)
  endif
  sinphi = sin(cphi)
  cosphi = cos(cphi)

  x = -sinthe*cosphi
  z =  sinthe*sinphi
  y =  costhe

  cthetan = acos(z)
  cphin = atan2(y,x)

  m = int((cthetan - zya(1))/zdy(1)) + 1
  n = int((cphin   - zza(1))/zdz(1)) + 1

  cz2(j,k) = (cphin-zzc(n))/zdz(kmax)
  if (cz2(j,k) < 0.0) then
    cz2(j,k) = -cz2(j,k)
  endif

  cz1(j,k) = 1.0 - cz2(j,k)

  cz4(j,k) = (cthetan-zyc(m))/zdy(jmax)
  if(cz4(j,k) < 0.0) then
    cz4(j,k) = -cz4(j,k)
  endif

  cz3(j,k) = 1.0 - cz4(j,k)
  
  cz5(j,k) = sinphi*sin(cphin)
  cz6(j,k) = cos(cphin)/sinthe
 enddo
enddo

!---------------------------------------------------------------------------------------
! Loop over full 2D slice and see if I have any zones needed for interpolation
! Output:   md               : total number of blocks I have to send out
!           sbuff(md)%targ_j : index on target PE of ghost zone (1-6,7-12)
!           sbuff(md)%targ_k : index on target PE of ghost zone (1-ks)
!           sbuff(md)%targ_n : index on target PE of ghost zone neighbors (1-4) used in interpolation
!           jsource(md)      : j index of zone on my grid that I need to send
!           ksource(md)      : k index of zone on my grid that I need to send
!           scount (npez)    : array of block counts I am sending to other PEs
!           rcount (npez)    : array of block counts I am receiving from other PEs
!           soffset(npez)    : displacement vector used in ALLTOALLV
!           roffset(npez)    : displacement vector used in ALLTOALLV
!---------------------------------------------------------------------------------------

md = 0
do j = 1, jmax
 jtarget = (j-1)/js + 1
 jlocal  = j - (jtarget-1)*js
 ctheta = zyc(j)
 sinthe = sin(ctheta)
 costhe = cos(ctheta)
 kkmin  = knmin(j) - 6
 kkmax  = knmax(j) - 6
  
 do k = 1, 12

  if (k<7) then
   cphi = zzc(kkmin) - real(k)*zdz(kkmin)
  else
   cphi = zzc(kkmax) + real(k-6)*zdz(kkmax)
  endif
  sinphi = sin(cphi)
  cosphi = cos(cphi)

  x = -sinthe*cosphi
  z =  sinthe*sinphi
  y =  costhe

  cthetan = acos(z)
  cphin = atan2(y,x)

  m = int((cthetan - zya(1))/zdy(1)) + 1
  n = int((cphin   - zza(1))/zdz(1)) + 1

  ms = m + 1
  if((cthetan-zyc(m)) < 0.0) ms = m - 1
  
  ! are either n or ns on my grid?
  if ((m-1)/js == jcol) then
    ns = n + 1
    if ((cphin-zzc(n)) < 0.0) ns = n - 1
    
    md = md + 1
    jbndry(md) = jlocal
    kbndry(md) = k
    nbndry(md) = 1  
    js0(md) = m - js*jcol
    ks0(md) = n
    targetpe(md) = jtarget
    
    md = md + 1
    jbndry(md) = jlocal
    kbndry(md) = k
    nbndry(md) = 2  
    js0(md) = m - js*jcol
    ks0(md) = ns
    targetpe(md) = jtarget
  endif

  if ((ms-1)/js == jcol) then
    ns = n + 1
    if ((cphin-zzc(n)) < 0.0) ns = n - 1
    
    md = md + 1
    jbndry(md) = jlocal
    kbndry(md) = k
    nbndry(md) = 3  
    js0(md) = ms - js*jcol
    ks0(md) = n
    targetpe(md) = jtarget
    
    md = md + 1
    jbndry(md) = jlocal
    kbndry(md) = k
    nbndry(md) = 4  
    js0(md) = ms - js*jcol
    ks0(md) = ns
    targetpe(md) = jtarget
  endif
  
 enddo
enddo
 
!used to make sure source arrays are dimensioned large enough   
!call MPI_REDUCE( md, mdmax, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, mpierr )
!if (mype==0) write(8,*) 'maximum mdz = ', mdmax

!---------------------------------------------------------------------------------------
! repack boundary data into send buffer sorted by target PE
! md is total number of blocks this PE is sending out

mdz = md
ALLOCATE (szbuff(nvar,isz,mdz))
ALLOCATE (szbufg(nvar,isz,mdz))
szc_tot = nvar*isz*mdz

szcount = 0  
il = 0  
do n = 1, npey
 do m = 1, mdz
 
  if (targetpe(m) == n) then ! this block is going to this PE=n
    il = il + 1
    jzt(il) = jbndry(m)
    kzt(il) = kbndry(m)
    nzt(il) = nbndry(m)
    szcount (n)  = szcount(n) + 1
    jzsource(il) = js0(m)
    kzsource(il) = ks0(m)
  endif  ! should end with il = md
  
 enddo
enddo

CALL MPI_ALLTOALL(szcount,1,MPI_INTEGER,rzcount,1,MPI_INTEGER,MPI_COMM_ROW,mpierr)

if ( sum(rzcount) /= 48*js ) print *, 'rzcount ', sum(rzcount), 48*js, mype

szoffset = 0
rzoffset = 0
do n = 2, npey
 szoffset(n) = szoffset(n-1) + szcount(n-1)
 rzoffset(n) = rzoffset(n-1) + rzcount(n-1)
enddo

call MPI_ALLTOALLV(jzt,szcount,szoffset,MPI_INTEGER,jztarget,rzcount,rzoffset,MPI_INTEGER,MPI_COMM_ROW,mpierr)
call MPI_ALLTOALLV(kzt,szcount,szoffset,MPI_INTEGER,kztarget,rzcount,rzoffset,MPI_INTEGER,MPI_COMM_ROW,mpierr)
call MPI_ALLTOALLV(nzt,szcount,szoffset,MPI_INTEGER,nztarget,rzcount,rzoffset,MPI_INTEGER,MPI_COMM_ROW,mpierr)

! now redo indexing for variable buffer
szcount = szcount * nvar * isz
rzcount = rzcount * nvar * isz
szoffset = 0
rzoffset = 0
do n = 2, npey
 szoffset(n) = szoffset(n-1) + szcount(n-1)
 rzoffset(n) = rzoffset(n-1) + rzcount(n-1)
enddo

return
end
    

