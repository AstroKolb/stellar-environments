subroutine boundaryY
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
REAL ::  ut, vt
REAL, DIMENSION(7,isy,4,12,ks) :: bdata

!---------------------------------------------------------------------------------------
! load send buffers with data for boundary conditions
if (mdy > 0) then
 do n = 1, mdy  
  m = (jysource(n)-1)/js + 1
  j = jysource(n) - js*(m-1)
  do i = 1, isy
   do nv = 1, 7
     sybuff(nv,i,n) = recv1(nv,kysource(n),j,i,m)
   enddo
  enddo
 enddo
endif

! swap buffers between Yin and Yang grids
mpitag = mypeyy
if (yin) then
  yy_pe = mype + npe/2
  call MPI_SEND(sybuff, syc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, mpierr)
  call MPI_RECV(sybufg, syc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, yystat, mpierr)
else
  yy_pe = mypeyy
  call MPI_RECV(sybufg, syc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, yystat, mpierr)
  call MPI_SEND(sybuff, syc_tot, VH1_DATATYPE, yy_pe, mpitag, MPI_COMM_WORLD, mpierr)
endif

! scatter data to target pes
call MPI_ALLTOALLV(sybufg,sycount,syoffset,VH1_DATATYPE,rybuff,rycount,ryoffset,VH1_DATATYPE,MPI_COMM_COL,mpierr)

!---------------------------------------------------------------------------------------
! Loop over receive buffer and pack interpolation zone data into bdata array

nreceive = 12 * 4 * ks  ! should equal sum(rcount) ?
do n = 1, nreceive
 do i = 1, isy
  do nv = 1, 7
    bdata(nv,i,nytarget(n),jytarget(n),kytarget(n)) = rybuff(nv,i,n)
  enddo
 enddo
enddo

!---------------------------------------------------------------------------------------
! load up bc arrays for use in sweepy using interpolation data passed in bdata

do k = 1, ks
 do j = 1, 12  ! 1,6 for upper boundary ; 7,12 for lower boundary
  do i = 1, isy
    
    ybr(j,k,i) = cy3(j,k)*(cy1(j,k)*bdata(1,i,1,j,k) + cy2(j,k)*bdata(1,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(1,i,3,j,k) + cy2(j,k)*bdata(1,i,4,j,k))
    ybp(j,k,i) = cy3(j,k)*(cy1(j,k)*bdata(2,i,1,j,k) + cy2(j,k)*bdata(2,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(2,i,3,j,k) + cy2(j,k)*bdata(2,i,4,j,k))
    ybw(j,k,i) = cy3(j,k)*(cy1(j,k)*bdata(3,i,1,j,k) + cy2(j,k)*bdata(3,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(3,i,3,j,k) + cy2(j,k)*bdata(3,i,4,j,k))
    ut         = cy3(j,k)*(cy1(j,k)*bdata(4,i,1,j,k) + cy2(j,k)*bdata(4,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(4,i,3,j,k) + cy2(j,k)*bdata(4,i,4,j,k))
    vt         = cy3(j,k)*(cy1(j,k)*bdata(5,i,1,j,k) + cy2(j,k)*bdata(5,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(5,i,3,j,k) + cy2(j,k)*bdata(5,i,4,j,k))
    ybf(j,k,i) = cy3(j,k)*(cy1(j,k)*bdata(6,i,1,j,k) + cy2(j,k)*bdata(6,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(6,i,3,j,k) + cy2(j,k)*bdata(6,i,4,j,k))
    ybc(j,k,i) = cy3(j,k)*(cy1(j,k)*bdata(7,i,1,j,k) + cy2(j,k)*bdata(7,i,2,j,k)) &
               + cy4(j,k)*(cy1(j,k)*bdata(7,i,3,j,k) + cy2(j,k)*bdata(7,i,4,j,k))
            
    ybu(j,k,i) = -cy5(j,k)*ut - cy6(j,k)*vt
    ybv(j,k,i) = +cy6(j,k)*ut - cy5(j,k)*vt
    
  enddo
 enddo
enddo

return
end
