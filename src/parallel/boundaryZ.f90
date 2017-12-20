subroutine boundaryZ
! Arrays of ghost-cell values are loaded for the Yin (Yang) grid
! using linear interpolations of variables from the Yang (Yin) grid.
! Written by Michael Owen (2007), John Blondin (2010)
!--------------------------------------------------------------------------------------

use global
use zone
use yinyang

IMPLICIT NONE

include 'mpif.h'

! LOCALS
INTEGER :: i,j,k,n, m, l, ll, nreceive, mpierr, nv, yy_pe, mpitag, yystat(MPI_STATUS_SIZE)
REAL ::  wt, ut
REAL, DIMENSION(7,isz,4,js,12) :: bdata

!---------------------------------------------------------------------------------------

if (mdz > 0) then
 do n = 1, mdz  
  m = (kzsource(n)-1)/ks + 1
  k = kzsource(n) - ks*(m-1)
  do i = 1, isz
   do nv = 1, 7
     szbuff(nv,i,n) = recv3(nv,jzsource(n),k,i,m)
   enddo
  enddo
 enddo
endif

! swap buffers between Yin and Yang grids
mpitag = mypeyy  ! swapping pes should share the same value of mypeyy
if (yin) then
  yy_pe = mype + npe/2
  call MPI_SEND(szbuff, szc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, mpierr)
  call MPI_RECV(szbufg, szc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, yystat, mpierr)
else
  yy_pe = mypeyy
  call MPI_RECV(szbufg, szc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, yystat, mpierr)
  call MPI_SEND(szbuff, szc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, mpierr)
endif

! scatter data to target pes
call MPI_ALLTOALLV(szbufg,szcount,szoffset,VH1_DATATYPE,rzbuff,rzcount,rzoffset,VH1_DATATYPE,MPI_COMM_ROW,mpierr)

!---------------------------------------------------------------------------------------
! Loop over receive buffer and pack interpolation zone data into bdata array

nreceive = 12 * 4 * js  ! should equal sum(rcount) ?
do n = 1, nreceive
 do i = 1, isz
  do nv = 1, 7
    bdata(nv,i,nztarget(n),jztarget(n),kztarget(n)) = rzbuff(nv,i,n)
  enddo
 enddo
enddo

!---------------------------------------------------------------------------------------
! load up bc arrays for use in sweepy using interpolation data passed in bdata

do j = 1, js
 do k = 1, 12  ! 1,6 for upper boundary ; 7,12 for lower boundary
  do i = 1, isz
    
    zbr(k,j,i) = cz3(j,k)*(cz1(j,k)*bdata(1,i,1,j,k) + cz2(j,k)*bdata(1,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(1,i,3,j,k) + cz2(j,k)*bdata(1,i,4,j,k))
    zbp(k,j,i) = cz3(j,k)*(cz1(j,k)*bdata(2,i,1,j,k) + cz2(j,k)*bdata(2,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(2,i,3,j,k) + cz2(j,k)*bdata(2,i,4,j,k))
    zbv(k,j,i) = cz3(j,k)*(cz1(j,k)*bdata(3,i,1,j,k) + cz2(j,k)*bdata(3,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(3,i,3,j,k) + cz2(j,k)*bdata(3,i,4,j,k))
    wt         = cz3(j,k)*(cz1(j,k)*bdata(4,i,1,j,k) + cz2(j,k)*bdata(4,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(4,i,3,j,k) + cz2(j,k)*bdata(4,i,4,j,k))
    ut         = cz3(j,k)*(cz1(j,k)*bdata(5,i,1,j,k) + cz2(j,k)*bdata(5,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(5,i,3,j,k) + cz2(j,k)*bdata(5,i,4,j,k))
    zbf(k,j,i) = cz3(j,k)*(cz1(j,k)*bdata(6,i,1,j,k) + cz2(j,k)*bdata(6,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(6,i,3,j,k) + cz2(j,k)*bdata(6,i,4,j,k))
    zbc(k,j,i) = cz3(j,k)*(cz1(j,k)*bdata(7,i,1,j,k) + cz2(j,k)*bdata(7,i,2,j,k)) &
               + cz4(j,k)*(cz1(j,k)*bdata(7,i,3,j,k) + cz2(j,k)*bdata(7,i,4,j,k))
            
    zbw(k,j,i) = -cz5(j,k)*wt - cz6(j,k)*ut
    zbu(k,j,i) = +cz6(j,k)*wt - cz5(j,k)*ut
    
  enddo
 enddo
enddo

return
end


