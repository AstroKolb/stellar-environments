subroutine prin

! collect all rows of data onto jcol=0, and then 
! write out pez netcdf data files containing all variables (plus time as an attribute).
!---------------------------------------------------------------------------

use NETCDF
use global
use zone

include 'mpif.h'

! LOCALS

INTEGER, PARAMETER :: nvars = 6

character(LEN=1) :: char, tmp3
character(LEN=4) :: tmp1, tmp2
character(LEN=50) :: filename
CHARACTER(LEN=16), DIMENSION(nvars) :: varname
 
INTEGER :: i, j, k
INTEGER :: status, ncid, gathbuffer_size, mpierr, m, nvar, jsk, nv
INTEGER :: xDimID, yDimID, zDimID
INTEGER, DIMENSION(nvars) :: varID
INTEGER :: XScale_varID, YScale_varID, ZScale_varID

REAL, DIMENSION(kmax/pez) :: zshort
REAL, DIMENSION(imax,jmax,kmax/pez) :: prin_buff
REAL, DIMENSION(imax,nvars,jmax/pey,kmax/pez) :: send_buff
REAL, DIMENSION(imax,nvars,jmax/pey,kmax/pez,pey) :: recv_buff

REAL :: xvl, yvl, zvl, xct, yct, the

!------------------------------------------------------------------------------
! everybody loads up a send buffer; data is gathered on jcol=0 procs.
nvar = nvars   ! # of 3D arrays to write into netcdf file

varname(1) = 'Density'
varname(2) = 'Pressure'
varname(3) = 'Color'
varname(4) = 'Xvelocity'
varname(5) = 'Yvelocity'
varname(6) = 'Zvelocity'

do k = 1, ks
   sph = sin(zzc(mypez*ks+k))
   cph = cos(zzc(mypez*ks+k))
do j = 1, js
   sth = sin(zyc(mypey*js+j))
   cth = cos(zyc(mypey*js+j))
do i = 1, imax
   send_buff(i,1,j,k) = zro(i,j,k)
   send_buff(i,2,j,k) = zpr(i,j,k)
   send_buff(i,3,j,k) = zcl(i,j,k)
   send_buff(i,4,j,k) = zux(i,j,k)
   send_buff(i,5,j,k) = zuy(i,j,k)
   send_buff(i,6,j,k) = zuz(i,j,k)

!   ! find angle of point
!   if (yin) then
!      xct = zxc(i)*sth*cph
!      yct = zxc(i)*sth*sph
!   else
!      xct =-zxc(i)*sth*cph
!      yct = zxc(i)*cth
!   endif

!   ! convert velocity to cartesian
!   xvl = zux(i,j,k)*sth*cph + zuy(i,j,k)*cth*cph - zuz(i,j,k)*sph
!   yvl = zux(i,j,k)*sth*sph + zuy(i,j,k)*cth*sph + zuz(i,j,k)*cph
!   zvl = zux(i,j,k)*cth     - zuy(i,j,k)*sth

!   ! write velocity less rotation
!   if (yin) then
!      send_buff(i,4,j,k) = xvl - yct*omega             !zxc(i)*omega*sin(the)
!      send_buff(i,5,j,k) = yvl + xct*omega - rcm*omega !zxc(i)*omega*cos(the)
!      send_buff(i,6,j,k) = zvl
!   else
!      send_buff(i,4,j,k) =-xvl - yct*omega             !zxc(i)*omega*sin(the)
!      send_buff(i,5,j,k) = zvl + xct*omega - rcm*omega !zxc(i)*omega*cos(the)
!      send_buff(i,6,j,k) = yvl
!   endif

enddo
enddo
enddo

    
gathbuffer_size = imax * js * ks * nvar 
call MPI_GATHER(send_buff, gathbuffer_size, VH1_DATATYPE, recv_buff, gathbuffer_size, VH1_DATATYPE, 0, MPI_COMM_ROW, mpierr)

!--------------------------------------------------------------------------------
! only the jcol=0 procs create nc files, unload receive buffer and write data to disk

if (jcol == 0) then

! Create filename from integers nfile and krow and prefix such that filename
! looks like prefx_1000.0000 where 1000 is the value of nfile and 0000 is krow

write(tmp1,"(i4)") nfile
write(tmp2,"(i4)") krow
if (yin) then
  tmp3 = 'n'
else
  tmp3 = 'g'
endif
do i = 1, 4
  if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
enddo
if (step == 1) then
  filename = 'output/' // trim(prefix) // '.s1_' // tmp1 // '.' // tmp3 // tmp2
else if (step == 2) then
  filename = 'output/' // trim(prefix) // '.s2_' // tmp1 // '.' // tmp3 // tmp2
else
  filename = 'output/' // trim(prefix) // '.s3_' // tmp1 // '.' // tmp3 // tmp2
endif

nfile = nfile + 1

! create the netCDF file
status = nf90_create(filename, nf90_Clobber, ncid)
if (status /= nf90_NoErr ) print *, NF90_STRERROR(status)

! Initialize Dimensions
status = nf90_def_dim(ncid, "x", imax, xDimID)      
status = nf90_def_dim(ncid, "y", jmax, yDimID)
status = nf90_def_dim(ncid, "z", ks  , zDimID)

! define the coordinate scales as 1D variables
status = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
status = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
status = nf90_def_var(ncid, "z", nf90_float,(/ zDimID /), ZScale_varID)

! define the simulation variables
do nv = 1, nvar
 status = nf90_def_var(ncid, varname(nv), nf90_float,(/ xDimID, yDimID, zDimID /), varID(nv))
enddo

! perhaps add variable attributes like units, min/max...

! define some global attributes
status = nf90_put_att(ncid, nf90_global, "time", time)
status = nf90_put_att(ncid, nf90_global, "gamma", gam)
 
! take dataset out of definition mode
status = NF90_ENDDEF(ncid)

do nv = 1, nvar

  do m = 1, npey
   do k = 1, ks
    do j = 1, js
     jsk = (m-1)*js + j
     do i = 1, imax
       prin_buff(i,jsk,k) = recv_buff(i,nv,j,k,m)
     enddo
    enddo
   enddo
  enddo
 status = nf90_put_var(ncid, varID(nv), prin_buff)

enddo

status = nf90_put_var(ncid, XScale_varID, zxc)
status = nf90_put_var(ncid, YScale_varID, zyc)
do k = 1, ks
 zshort(k) = zzc(krow*ks+k)
enddo
status = nf90_put_var(ncid, ZScale_varID, zshort)

! always a necessary statment to flush output 
status = nf90_close(ncid)

endif


if (mype==0) write(8,6000) filename, time, ncycle
6000 format('Wrote ',a14,' to disk at time =',1pe12.5,' (ncycle =', i6,')')

return
end



