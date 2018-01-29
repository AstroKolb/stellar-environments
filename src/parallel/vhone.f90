program vhone

!--------------------------------------------------------------------------
!
!          VV              VV   HH       HH         111
!           VV            VV    HH       HH        1111
!            VV          VV     HH       HH          11
!             VV        VV      HH       HH          11
!              VV      VV       HHHHHHHHHHH   ==     11
!               VV    VV        HH       HH          11
!                VV  VV         HH       HH          11
!                 VVVV          HH       HH          11
!                  VV           HH       HH        111111
!
! 
!                        VIRGINIA HYDRODYNAMICS #1
!
!--------------------------------------------------------------------------
!  The Virginia Numerical Bull Session ideal hydrodynamics PPMLR
!  
!  Version 2.0 July 1997
!    Many small modifications, mainly brought BCs into only one call to sweepbc 
!  Version 1.0 October 1991 
!    This version was removed from the ``loop'' in the continued development
!    of VH-1 on Sep 1 1990, but has followed a similar evolution until this
!    release date.
! 
! Output channels:
!                 4 = aaaaadaa,   binary dump files for restarting
!                 8 = aaaaahst,   history file
!                     aaaaaXY.nc,   netCDF movie of XY slice
!                     aaaaa_1000.nc, netCDF output of full dataset
!-------------------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweepsize

include 'mpif.h'

! LOCALS
character(len=50) :: hstfile
character(len=52) :: indat_chrs, icons_chrs
character(len=8)  :: dmpfile
character(len=4)  :: tmp4 
character(len=2)  :: rstrt 
character(len=8)  :: todayis
character(len=MPI_MAX_PROCESSOR_NAME) :: myname

integer :: i, j, k, nfile_start
integer :: ncycend, nprin, ndump, nmovie, tmpcyc, mpierr, yy
integer :: a2abuffer_size,start_cycl, namelength
integer, dimension(5) :: indat_ints

real :: olddt
real :: endtime,tprin,tmovie
real :: rshockmax, Rs_max, xexpand
real, dimension(3) :: indat_rels
real, dimension(10) :: icons_rels
real :: start_time, end_time, run_time, zones



namelist / hinput / step, rstrt, prefix, ncycend, ndump, nprin, nmovie, endtime, tprin, tmovie
namelist / hicons / scale, GMP, GMS, sep, x1min, x1max, x2max, gam, LumP, dmP, TmpP
!-----------------------------------------------------------------------------------------
! Initialize MPI
call MPI_INIT(mpierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, npe, mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, mype, mpierr)

! Set VH1_DATATYPE for MPI calls
if (SIZEOF(endtime)==8) then
  VH1_DATATYPE = MPI_DOUBLE_PRECISION
elseif (SIZEOF(endtime)==4) then
  VH1_DATATYPE = MPI_REAL
else
  if(mype==0) write(*,*) 'Mode not recognized. MPI commands will fail.'
  stop
endif

! Break up processes into two equal halves for Yin and Yang grids
yy = mype / (npe/2)
call MPI_COMM_SPLIT(MPI_COMM_WORLD, yy, mype, MPI_COMM_YY, mpierr)
  call MPI_COMM_RANK(MPI_COMM_YY, mypeyy, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_YY, npeyy, mpierr)

yin = .false.; if (yy==0) yin = .true.

! set dimensions of simulation grid and subgrids to array dimensions

js   = jmax / pey
ks   = kmax / pez
isy  = imax / pey
isz  = imax / pez

jcol = mod(mypeyy,pey)   ! processor rank in the Y direction
krow = mypeyy / pey      ! processor rank in the Z direction

! Create communicators for each value of krow (for XY transpose)
call MPI_COMM_SPLIT(MPI_COMM_YY, krow, mypeyy, MPI_COMM_ROW, mpierr)
  call MPI_COMM_RANK(MPI_COMM_ROW, mypey, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_ROW, npey, mpierr)

! Create communicators for each value of jcol (for XZ transpose)
call MPI_COMM_SPLIT(MPI_COMM_YY, jcol, mypeyy, MPI_COMM_COL, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_COL, npez, mpierr)
  call MPI_COMM_RANK(MPI_COMM_COL, mypez, mpierr)

! should have npey = pey, npez = pez, mypey = jcol, mypez = krow

!-----------------------------------------------------------------------------------------
! mype=0 does some checking on core counts and array dimensions
if (mype == 0) then

 ! Check that code is running with desired number of PE's
 if ((pey*pez*2) /= npe) then
  write(*,*) 'sorry, I was compiled for ',2*pey*pez,' PEs.'
  stop
 endif

 ! Check that arrays are large enough for desired number of physical zones

 if (max(imax,jmax,kmax)+12 > maxsweep) then
  write(*,*) 'maxsweep too small'
  stop
 endif

 ! Check that imax and jmax are multiples of npey
 if (isy*npey /= imax) then
  write(*,*) 'imax is not an integer multiple of pey', imax, pey
  stop
 endif
 if (js*npey /= jmax) then
  write(*,*) 'jmax is not an integer multiple of pey', jmax, pey
  stop
 endif

 ! Check that imax and kmax are multiples of npez
 if (isz*npez /= imax) then
  write(*,*) 'imax is not an integer multiple of pez', imax, pez
  stop
 endif
 if (ks*npez /= kmax) then
  write(*,*) 'kmax is not an integer multiple of pez', kmax, pez
  stop
 endif

endif

!-----------------------------------------------------------------------------------------
! Begin by reading input deck, and defining file names
!   only jcol=0 tasks read input, broadcasts to everybody else in their domain

if (jcol==0) then
 open (unit=15,file='indat',status='old',form='formatted')
  read (15,nml=hinput)
 close(15)
 indat_ints(1) = step
 indat_ints(2) = ncycend
 indat_ints(3) = ndump
 indat_ints(4) = nmovie
 indat_ints(5) = nprin
 indat_rels(1) = endtime
 indat_rels(2) = tprin
 indat_rels(3) = tmovie
 indat_chrs(1:2) = rstrt
 indat_chrs(3:52) = prefix
endif

call MPI_BCAST( indat_ints, 5, MPI_INTEGER, 0, MPI_COMM_ROW, mpierr )
 step    = indat_ints(1)
 ncycend = indat_ints(2)
 ndump   = indat_ints(3)
 nmovie  = indat_ints(4)
 nprin   = indat_ints(5) 
call MPI_BCAST( indat_rels, 3, VH1_DATATYPE, 0, MPI_COMM_ROW, mpierr )
 endtime = indat_rels(1)
 tmovie  = indat_rels(2)
 tprin   = indat_rels(3)
call MPI_BCAST( indat_chrs, 52, MPI_CHARACTER, 0, MPI_COMM_ROW, mpierr )
 rstrt  = indat_chrs(1:2)
 prefix = indat_chrs(3:52)


if (.false.) then
 ! load in icons
if (jcol==0) then
  open (unit=17, file='icons', status='old', form='formatted')
    read (17, nml=hicons)
  close(17)
  icons_chrs     = scale
  icons_rels(1)  = GMP
  icons_rels(2)  = GMS
  icons_rels(3)  = sep
  icons_rels(4)  = x1min
  icons_rels(5)  = x1max
  icons_rels(6)  = x2max
  icons_rels(7)  = gam
  icons_rels(8)  = LumP
  icons_rels(9)  = dmP
  icons_rels(10) = TmpP
endif

call MPI_BCAST( icons_rels, 10, VH1_DATATYPE, 0, MPI_COMM_ROW, mpierr )
  GMP   = icons_rels(1)
  GMS   = icons_rels(2)
  sep   = icons_rels(3)
  x1min = icons_rels(4)
  x1max = icons_rels(5)
  x2max = icons_rels(6)
  gam   = icons_rels(7)
  LumP  = icons_rels(8)
  dmP   = icons_rels(9)
  TmpP  = icons_rels(10)
call MPI_BCAST( icons_chrs, 52, MPI_CHARACTER, 0, MPI_COMM_ROW, mpierr )
  scale = icons_chrs
endif

call read_icons

!-----------------------------------------------------------------------------------------
write(tmp4,"(i4)") mype
do i = 1, 4
  if ((tmp4(i:i)) == ' ') tmp4(i:i) = '0'
enddo
dmpfile = 'daa.' // tmp4

if ((rstrt == 'NO').or.(rstrt == 'no')) then ! Initialize variables for new problem; 

  if (mype == 0) then
   call date_and_time(todayis)
   hstfile = 'output/' // trim(prefix) // '.hst'
   open (unit= 8,file=hstfile,form='formatted')
   write (8,*) 'History File for YY-VH1 simulation run on ', todayis(5:6), ' / ', todayis(7:8), ' / ', todayis(1:4)
   call MPI_GET_PROCESSOR_NAME(myname, namelength, mpierr)
   write (8,"(' machine id: ',a25)") myname
   write (8,"(' Running on ',i4,' processors, split ',i3, ' X ', i3,' X 2')") npe,npey,npez
   if(VH1_DATATYPE==MPI_DOUBLE_PRECISION) then
     write(8,*) 'Running in double precision mode'
   else
     write(8,*) 'Running in single precision mode'
   endif
   write (8,*)
  endif

  call init
  nfile_start = nfile
  call prin
  !call images

else ! Restart from old dumpfile...
  
  dmpfile = 'd' // rstrt(1:2) // dmpfile(4:8)
  if (mype == 0) then
    hstfile = 'output/' // trim(prefix) // '.hst'
    open (unit= 8,file=hstfile,form='formatted',position='append')
  endif
  call undump(dmpfile)
  if (timep >= tprin)  timep = tprin-2*dt    ! adjust so that dt is not set to zero on first step
  if (timem >= tmovie) timem = tmovie-2*dt
  nfile_start = nfile

endif

call corner
call boundaryYsetup
call boundaryZsetup

write(*,*) 'made it..'
if(mype == 0) then
  !start_time = MPI_WTIME()
  start_cycl = ncycle
  write(8,*) 'Starting on cycle number ',ncycle
endif


if (step == 1) then
  ncycend = 1000000
  nprin   = 1000000
  endtime = 2.00 * opd
  tprin   = 0.20 * opd
else if (step == 2) then
  ncycend = 1000000
  nprin   = 1000000
  endtime = 10.0 * opd
  tprin   = 00.5 * opd
else


endif

!############################################################################
!                         MAIN COMPUTATIONAL LOOP

do while (ncycle < ncycend)

  if (step < 3) then
    if (mype == 0) write(*,*) 'STEP = ', ncycle, dt/opd, time/opd
  else
    if (mype == 0) write(*,*) 'STEP = ', ncycle, zxc(i) / x2max
  endif

  ncycle = ncycle + 1
  ncycp  = ncycp  + 2
  ncycd  = ncycd  + 2
  ncycm  = ncycm  + 2
  olddt  = dt
  svel   = 0.

  if ( time + 2*dt .gt. endtime ) then ! set dt to land on endtime
    if(mype.eq.0) write(8,*) 'cutting to the end...', ncycle, ncycend
    dt = 0.5*(endtime - time)
    ncycend = ncycle-1
    ncycp   = nprin
    ncycd   = ndump
  else if ( timep+2*dt .gt. tprin ) then ! set dt to land on tprin
    dt = 0.5*(tprin - timep)
    ncycp = nprin
  else if ( timem+2*dt .gt. tmovie ) then ! set dt to land on tmovie
    dt = 0.5*(tmovie - timem)
    ncycm = nmovie
  endif

! =============================================================
! Find maximum radius of shock and terminate if it gets close to the edge
  rshockmax = 0.0
  do k = 1, ks
    do j = 1, js
      i = imax
      do while (zux(i,j,k) < 1e3)
        i = i - 1
      enddo
      rshockmax = max(rshockmax,zxa(i))
    enddo
  enddo
  call MPI_ALLREDUCE( rshockmax, Rs_max, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
  if (Rs_max > zxa(imax-6) .and. step == 3) then
    xexpand = 0.7*zdx(imax)/zxa(imax)
  else
    xexpand = 0.0
  endif
! =============================================================

! Alternate sweeps to approximate 2nd order operator splitting
  call sweepx1(xexpand)
  call sweepy
  call sweepz
  call sweepy
  call sweepx2(xexpand)

  time  = time  + 2.0*dt
  timep = timep + 2.0*dt
  timem = timem + 2.0*dt
  dt = olddt

  call dtcon   ! Check constraints on the timestep

  !if ((ncycle - 100*(ncycle/100)) == 0) 

  ! output movie images/datasets/dumpfiles based on simulation time or cycle #

  if ( ncycm.ge.nmovie .or. timem.ge.tmovie ) then
    !if (krow == npez/2) call imagesxy
    !if (jcol == npey/2) call imagesxz
    !if (krow == 0) call imagesxy
    !if (jcol == 0) call imagesxz
    !call images
    timem = 0.0
    ncycm = 0
  endif

  if ( ncycp >= nprin .or. timep >= tprin ) then
    call prin
    timep = 0.
    ncycp = 0
  endif

  if ( ncycd >= ndump ) then
    ncycd = 0
    call dump(dmpfile)
  endif

enddo
!                           END OF MAIN LOOP
!#########################################################################
      
if(mype == 0) then
  !end_time = MPI_WTIME()
  zones = (imax*js*ks)
  zones = zones*(ncycle - start_cycl)
  !run_time = end_time - start_time
  write(8,*) 'successful stop at cycle number', ncycle
  !if (run_time>4000.0) then
  !  write(8,*) 'elapsed time = ',run_time/3600., ' hours'
  !else
  !  write(8,*) 'elapsed time = ',run_time/60., ' minutes'
  !endif
  !write(8,"('speed = ',f5.1,' kz/s/pe')") 1.0e-3*zones/run_time
  close( 8 )

  ! open up a file to store data needed for post-processing
  open(9,file='postprocess')
  write(9,*) trim(prefix)
  write(9,*) nfile_start
  write(9,*) nfile - nfile_start
  write(9,*) npey, npez
  close(9)
endif
      
call MPI_FINALIZE(mpierr)

stop      
end 
