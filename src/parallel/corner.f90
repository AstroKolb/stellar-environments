subroutine corner

! This subroutine performs cuts the corners off the grids
   
!-----------------------------------------------------------------------

! GLOBALS
use NETCDF
use global
use zone

IMPLICIT NONE

! LOCALS
CHARACTER(LEN=50) :: filename
 
INTEGER :: status, ncid
INTEGER :: xDimID, yDimID, zDimID, XScale_varID, YScale_varID, ZScale_varID,ghost_varID 
INTEGER :: i,j,k,l,m,n,imult,itheta,nphi,ntheta,kmin
REAL :: cphi, ctheta, cthetan, cphin,x,y,z,dphi,dtheta,phi0,theta0,thetap,phip,phi,theta

!-----------------------------------------------------------------------
! create an array, ongrid, that =1 if on the computational grid, =0 if off

ongrid = 0

dphi   = zdz(1)
dtheta = zdy(1)
phi0   = zza(1)
theta0 = zya(1)

imult = 1
do itheta = 1,3,2
  theta = real(itheta)*pi*0.25
  imult = -imult

  phi = dphi/2.0 - pi/2.0
  do while (phi <= dphi/2.0 + pi/2.0)
    thetap = acos(sin(theta)*sin(phi))
    phip   = asin(cos(theta)/sin(thetap))

    nphi   = 1 + int((phip-phi0+imult*pi)/dphi)
    ntheta = 1 + int((thetap-theta0)/dtheta)

    do m = -1,1,1
      do n = -1,1,1
        ongrid(min(jmax,max(1,ntheta+m)),min(kmax,max(1,nphi-n))) = 1
      enddo
    enddo

    phi = phi + dphi 
  enddo
enddo

! use ongrid(j,k) to set nmin, nmax for sweepz: knmin(j), knmax(j)
do j = 1, jmax/2
  k = 1
  do while((ongrid(j,k)==0).and.(k<(kmax/6)))
    k = k + 1
  enddo
  knmin(j) = k + 6
  knmax(j) = kmax - k + 7
  knmin(jmax+1-j) = knmin(j)
  knmax(jmax+1-j) = knmax(j)
  if(j==1) kmin = k-1
enddo

! repeat for nmin, nmax for sweepy: jnmin(k), jnmax(k), but only do the corners
! set default values
jnmin = 7
jnmax = jmax+6

! cut off corners
do k = 1, kmin
  j = 1
  do while((ongrid(j,k)==0).and.(j<jmax))
    j = j + 1
  enddo
  jnmin(k) = j + 6
  jnmax(k) = jmax - j + 7
  jnmin(kmax-k+1) = jnmin(k)
  jnmax(kmax-k+1) = jnmax(k)
enddo

! now repopulate ongrid with 0 if on the grid, 1 if in the corners
ongrid = 1
do j = 1, jmax
  do k = knmin(j)-6, knmax(j)-6
    ongrid(j,k) = 0
  enddo
enddo

if (mype==0) then
  ! open up a netcdf file and output array marking ghost cells
  if (step == 1) then
    filename = 'output/' // trim(prefix) // '.s1ghost.nc'
  else if (step == 2) then
    filename = 'output/' // trim(prefix) // '.s2ghost.nc'
  else
    filename = 'output/' // trim(prefix) // '.s3ghost.nc'
  endif

  status = nf90_create(filename, nf90_Clobber, ncid)
    if (status /= nf90_NoErr ) print *, NF90_STRERROR(status)

  status = nf90_def_dim(ncid, "x", imax, xDimID)
  status = nf90_def_dim(ncid, "y", jmax, yDimID)
  status = nf90_def_dim(ncid, "z", kmax, zDimID)
  status = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
  status = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
  status = nf90_def_var(ncid, "z", nf90_float,(/ zDimID /), ZScale_varID)
  status = nf90_def_var(ncid, "Ghosts",   nf90_int, (/ yDimID, zDimID /), ghost_varID)
  status = NF90_ENDDEF(ncid)

  status = nf90_put_var(ncid, XScale_varID, zxc)
  status = nf90_put_var(ncid, YScale_varID, zyc)
  status = nf90_put_var(ncid, ZScale_varID, zzc)
  status = nf90_put_var(ncid, ghost_varID,  ongrid)

  status = nf90_close(ncid)
endif

return
 
end

