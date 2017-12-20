program convert

! read in a slabs of netCDF data from a sliced Yin-Yang grid and write out in EnSight Gold format
!
! gfortran -o merge_slabs merge_slabs.f90 -I/usr/local/include -L/usr/local/lib -lnetcdf
!---------------------------------------------------------------------------

use NETCDF  ! netcdf modules

! LOCALS

LOGICAL :: skipgeo

CHARACTER(LEN=1) :: tmp1
CHARACTER(LEN=4) :: tmp4, tmp2
CHARACTER(LEN=80) :: varname
CHARACTER(LEN=50) :: prefix
CHARACTER(LEN=50) :: filename
CHARACTER(LEN=50) :: ghostfile
CHARACTER(LEN=50) :: geofile 
CHARACTER(LEN=50) :: casefile
CHARACTER(LEN=50) :: varfile

CHARACTER(LEN=80) :: label1a, label1b, label1c, label1d, label1e, label1f, label1g, label1h, label1i
CHARACTER(LEN=80) :: label2, label3, label4, coord1, coord2, coord3

INTEGER :: ipart, nsteps
INTEGER :: ncid, ncstat, ncid_slabs, rank

INTEGER :: n, i, j, k, ndim, nvar, nv
INTEGER :: imax, jmax, kmax
INTEGER, DIMENSION(4) :: start
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ghostcells
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: gs

REAL, DIMENSION(:,:,:), ALLOCATABLE :: var, slab, slab2
REAL, DIMENSION(:), ALLOCATABLE :: zxc, zyc, zzc
REAL, DIMENSION(0:500) :: time

!------------------------------------------------------------------------------
! set up labels for EnSight files
label1a = 'Fortran Binary'
label1b = 'cartesian geometry file from VH-1'
label1c = 'junk line 2'
label1d = 'node id off'
label1e = 'element id off'
label1f = 'part '
label1g = 'cartesian coordinates'
label1h = 'block curvilinear with_ghost'
label1i = 'ghost_flags'
ipart  = 1
label2 = 'part'
label3 = 'block'

start = 1

!------------------------------------------------------------------------------
write(6,*) 'Convert a series of netCDF data files to Ensight Gold format'
write(6,fmt="('Input PREFIX of netCDF data file: ')")
read (5,*) prefix
write(6,fmt="('Input first data file number (from zero): ')")
read (5,*) nstart
nfile = nstart 


! OPEN DATA FILE with ghost zones and geometry
ghostfile = trim(prefix) // 'ghost.nc'
print *, 'opening ', ghostfile

ncstat = nf90_open(ghostfile, nf90_nowrite, ncid)

 ncstat = nf90_inquire_dimension(ncid, 1, coord1, imax)
 ncstat = nf90_inquire_dimension(ncid, 2, coord2, jmax)
 ncstat = nf90_inquire_dimension(ncid, 3, coord3, kmax)

 ! allocate coordinate arrays
 ALLOCATE( zxc(imax) )
 ALLOCATE( zyc(jmax) )
 ALLOCATE( zzc(kmax) )
 ! allocate data array for variable
 ALLOCATE(ghostcells(jmax,kmax) )
 ALLOCATE(var(imax,jmax,kmax) )
 ALLOCATE(gs(imax-1,jmax-1,kmax-1) )

 print *, 'Allocated data arrays'
 
 ! read coordinates arrays
 ncstat = nf90_get_var(ncid, 1, zxc)
     if(ncstat.ne.0) write(*,*) NF90_STRERROR(ncstat)
 ncstat = nf90_get_var(ncid, 2, zyc)
     if(ncstat.ne.0) write(*,*) NF90_STRERROR(ncstat)
 ncstat = nf90_get_var(ncid, 3, zzc)
     if(ncstat.ne.0) write(*,*) NF90_STRERROR(ncstat)
 ncstat = nf90_get_var(ncid, 4, ghostcells)
     if(ncstat.ne.0) write(*,*) NF90_STRERROR(ncstat)
ncstat = nf90_close(ncid)

geofile  = 'E' // trim(prefix) // '.geo'
inquire (file=geofile, exist=skipgeo)

if (.not.skipgeo) then  ! write out geometry files

 open(unit=25,file=geofile,form='unformatted')
 geofile  = 'F' // trim(prefix) // '.geo'
 open(unit=26,file=geofile,form='unformatted')

 ! first shuffle right edge back one zone
 do j = 1, jmax
  do k = kmax/2, kmax-1
   ghostcells(j,k) = ghostcells(j,k+1)
  enddo
 enddo
 ! first shuffle bottom edge back one zone
 do k = 1, kmax
  do j = jmax/2, jmax-1
   ghostcells(j,k) = ghostcells(j+1,k)
  enddo
 enddo
 
 ! now copy ghostcells into a smaller 3D array (imax-1)*(jmax-1)*(kmax-1)
 do k = 1, kmax-1
  do j = 1, jmax-1
   do i = 1, imax-1
    gs(i,j,k) = ghostcells(j,k)
   enddo
  enddo
 enddo 


  write(25) label1a
  write(25) label1b
  write(25) label1c
  write(25) label1d
  write(25) label1e

  write(26) label1a
  write(26) label1b
  write(26) label1c
  write(26) label1d
  write(26) label1e

  ! Coordinates of Yin Grid
  write(25) label1f
  write(25) ipart
  write(25) label1g
  write(25) label1h
  write(25) imax, jmax, kmax

  write(26) label1f
  write(26) ipart
  write(26) label1g
  write(26) label1h
  write(26) imax, jmax, kmax

  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k) = zxc(i)*sin(zyc(j))*cos(zzc(k))
    enddo
   enddo
  enddo
  write(25) var

  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k) = -var(i,j,k)
    enddo
   enddo
  enddo
  write(26) var

  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k) = zxc(i)*sin(zyc(j))*sin(zzc(k))
    enddo
   enddo
  enddo
  write(25) var

  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k) = zxc(i)*cos(zyc(j))
    enddo
   enddo
  enddo
  write(26) var

  write(25) var
  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k) = zxc(i)*sin(zyc(j))*sin(zzc(k))
    enddo
   enddo
  enddo
  write(26) var

  write(25) label1i
  write(25) gs

  write(26) label1i
  write(26) gs

 close(25)
 close(26)

endif

DEALLOCATE (gs)

! open up first slab file to get dimensions and time
write(tmp4,"(i4)") nfile + 1000
n=0
write(tmp2,"(i4)") n
do i = 1, 4
  if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
enddo
filename = trim(prefix) // '_' // tmp4 // '.n' // tmp2  
ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
ncstat = nf90_inquire(ncid_slab, rank, nvar)
ncstat = nf90_inquire_dimension(ncid_slab, 1, coord1, is)
ncstat = nf90_inquire_dimension(ncid_slab, 2, coord2, js)
ncstat = nf90_inquire_dimension(ncid_slab, 3, coord3, ks)
npe = kmax / ks

ALLOCATE(slab(imax,jmax,ks) )
ALLOCATE(slab2(imax,jmax,ks) )

!#########################################################################
! Loop over all netCDF files (until open statement produces an error)
do while (ncstat /= 2)
 print *, 'processing file ', filename

 l = 1  

!%%%%%%%%%%%%%%%%%%%%%%   DENSITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nv = 4
var = 0.0
!already opened n=0...
  ncstat = nf90_get_att(ncid_slab, nf90_global, 'time', time(nfile))
  ncstat = nf90_get_var(ncid_slab, nv, slab)
  ncstat = nf90_get_var(ncid_slab,  5, slab2)   !  also read in pressure
  ncstat = nf90_close(ncid_slab)
  do k = 1, ks
   do j = 1, jmax
    do i = 1, imax
!      if (slab2(i,j,k)>0.02) var(i,j,k) = slab(i,j,k)  ! only copy density if pressure > 0.02 (shocked gas)
      var(i,j,k) = slab(i,j,k)
    enddo
   enddo
  enddo     
 
do n = 1, npe-1
  write(tmp2,"(i4)") n
  do i = 1, 4
    if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
  enddo
  filename = trim(prefix) // '_' // tmp4 // '.n' // tmp2  
  ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
  ncstat = nf90_get_var(ncid_slab, nv, slab)
  ncstat = nf90_get_var(ncid_slab,  5, slab2)
  ncstat = nf90_close(ncid_slab)
  koff = n*ks
  do k = 1, ks
   do j = 1, jmax
    do i = 1, imax
!      if (slab2(i,j,k)>0.02) var(i,j,k+koff) = slab(i,j,k)
      var(i,j,k+koff) = slab(i,j,k)
    enddo
   enddo
  enddo     
enddo


  varname = 'density'
  varfile = 'E' // trim(prefix) // tmp4 // '.den'
  open(unit=15,file=varfile,form='unformatted')
   write(15) varname
   write(15) label2
   write(15) l
   write(15) label3
   write(15) var
  close(15)

var = 0.0
do n = 0, npe-1
   write(tmp2,"(i4)") n
   do i = 1, 4
     if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
   enddo
   filename = trim(prefix) // '_' // tmp4 // '.g' // tmp2  
   ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
   ncstat = nf90_get_var(ncid_slab, nv, slab)
   ncstat = nf90_get_var(ncid_slab,  5, slab2)
   ncstat = nf90_close(ncid_slab)
  koff = n*ks
  do k = 1, ks
   do j = 1, jmax
    do i = 1, imax
!      if (slab2(i,j,k)>0.02) var(i,j,k+koff) = slab(i,j,k)
      var(i,j,k+koff) = slab(i,j,k)
    enddo
   enddo
  enddo     
enddo

  varfile = 'F' // trim(prefix) // tmp4 // '.den'
  open(unit=15,file=varfile,form='unformatted')
   write(15) varname
   write(15) label2
   write(15) l
   write(15) label3
   write(15) var
  close(15)

!%%%%%%%%%%%%%%%%%%%%%%   PRESSURE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nv = 5
do n = 0, npe-1
   write(tmp2,"(i4)") n
   do i = 1, 4
     if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
   enddo
   filename = trim(prefix) // '_' // tmp4 // '.n' // tmp2  
   ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
   ncstat = nf90_get_var(ncid_slab, nv, slab)
   ncstat = nf90_close(ncid_slab)
  koff = n*ks
  do k = 1, ks
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k+koff) = slab(i,j,k)
    enddo
   enddo
  enddo     
enddo
  varname = 'pressure'
  varfile = 'E' // trim(prefix) // tmp4 // '.pre'
  open(unit=15,file=varfile,form='unformatted')
   write(15) varname
   write(15) label2
   write(15) l
   write(15) label3
   write(15) var
  close(15)

do n = 0, npe-1
   write(tmp2,"(i4)") n
   do i = 1, 4
     if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
   enddo
   filename = trim(prefix) // '_' // tmp4 // '.g' // tmp2  
   ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
   ncstat = nf90_get_var(ncid_slab, nv, slab)
   ncstat = nf90_close(ncid_slab)
  koff = n*ks
  do k = 1, ks
   do j = 1, jmax
    do i = 1, imax
      var(i,j,k+koff) = slab(i,j,k)
    enddo
   enddo
  enddo     
enddo
  varfile = 'F' // trim(prefix) // tmp4 // '.pre'
  open(unit=15,file=varfile,form='unformatted')
   write(15) varname
   write(15) label2
   write(15) l
   write(15) label3
   write(15) var
  close(15)

 ! OPEN next netCDF FILE
 nfile = nfile + 1
 write(tmp4,"(i4)") nfile + 1000
 n=0
 write(tmp2,"(i4)") n
 do i = 1, 4
  if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
 enddo
 filename = trim(prefix) // '_' // tmp4 // '.n' // tmp2  
 ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
print *, filename, ncstat

enddo

nfile = nfile - 1
print *, 'did not find ', filename

!####################################################################
! Write out two case files

casefile = 'E' // trim(prefix) // '.case'
geofile  = 'E' // trim(prefix) // '.geo'
 open(unit=14,file=casefile)

  write(14,*) 'FORMAT'
  write(14,*) 'type: ensight gold'
  write(14,*) 'GEOMETRY'
  write(14,*) 'model: ', geofile
  write(14,*) 'VARIABLE'
  
  if (nfile == nstart) then  
  
   write(tmp4,"(i4)") nstart + 1000
   write(14,*) 'scalar per node:    density    ', 'E' // trim(prefix) // tmp4 // '.den'
   write(14,*) 'scalar per node:    pressure   ', 'E' // trim(prefix) // tmp4 // '.pre'
!   write(14,*) 'vector per node:    velocity   ', 'E' // trim(prefix) // tmp4 // '.vel'

  else

   write(14,*) 'scalar per node:    density    ', 'E' // trim(prefix) // '****.den'
   write(14,*) 'scalar per node:    pressure   ', 'E' // trim(prefix) // '****.pre'
!   write(14,*) 'vector per node:    velocity   ', 'E' // trim(prefix) // '****.vel'

   write(14,*) 'TIME'
   write(14,*) 'time set:    1'
   write(14,*) 'number of steps:    ', nfile -nstart+1
   write(14,*) 'filename start number:  ', 1000+nstart
   write(14,*) 'filename increment:    1'

   write(14,fmt="('time values: ')",advance="no") 
   framecount = 0
   do n = nstart, nfile
    write(14,fmt="(1pe12.4)",advance="no") time(n)
    framecount = framecount + 1
    if (framecount > 5) then
     write(14,*)
     framecount = 0
    endif
   enddo

  endif

 close(14)
 
geofile  = 'F' // trim(prefix) // '.geo'
casefile = 'F' // trim(prefix) // '.case'
 open(unit=14,file=casefile)

  write(14,*) 'FORMAT'
  write(14,*) 'type: ensight gold'
  write(14,*) 'GEOMETRY'
  write(14,*) 'model: ', geofile
  write(14,*) 'VARIABLE'
  
  if (nfile == nstart) then  
  
   write(tmp4,"(i4)") nstart + 1000
   write(14,*) 'scalar per node:    density    ', 'F' // trim(prefix) // tmp4 // '.den'
   write(14,*) 'scalar per node:    pressure   ', 'F' // trim(prefix) // tmp4 // '.pre'
!   write(14,*) 'vector per node:    velocity   ', 'F' // trim(prefix) // tmp4 // '.vel'

  else

   write(14,*) 'scalar per node:    density    ', 'F' // trim(prefix) // '****.den'
   write(14,*) 'scalar per node:    pressure   ', 'F' // trim(prefix) // '****.pre'
!   write(14,*) 'vector per node:    velocity   ', 'F' // trim(prefix) // '****.vel'

   write(14,*) 'TIME'
   write(14,*) 'time set:    1'
   write(14,*) 'number of steps:    ', nfile -nstart+1
   write(14,*) 'filename start number:  ', 1000+nstart
   write(14,*) 'filename increment:    1'
   write(14,fmt="('time values: ')",advance="no") 

   framecount = 0
   do n = nstart, nfile
    write(14,fmt="(1pe12.4)",advance="no") time(n)
    framecount = framecount + 1
    if (framecount > 5) then
     write(14,*)
     framecount = 0
    endif
   enddo

  endif

 close(14)
 
end


