subroutine boundaryRsetup
! =========================================================
! set up arrays for rotating inner boundary
! grid dimensions are zxc(6) / zyc(pi/4, 3pi/4) / zzc(0:2pi), for yin+yang
! rotation is assumed around the z-axis
! =========================================================

  use zone
  use global

implicit none

  include 'mpif.h'

  integer :: i,j,k,jj,kk,ii
  integer :: n,nv,y,pej,pek
  integer :: scattersize, mpierr

  integer, parameter :: nvars = 6

  real, dimension(:,:,:,:,:,:), allocatable :: rings2, scatter_buffs2
  real, dimension(:,:,:,:,:),   allocatable :: rings3, scatter_buffs3

  real, dimension(:,:,:,:), allocatable :: bound3, scatter_buffs4

  character(len=70) :: infilename

  real :: dalpha

! =========================================================

  if (step == 1) then

  else if (step == 2) then

    amax = 256

    ! set up alpha array
    allocate(alpha(amax))
    dalpha = 2.0*pi/amax
    alpha(1) = dalpha/2.0
    do n = 2, amax
      alpha(n) = alpha(n-1) + dalpha
    enddo

    ! load data into ring array
    if (mype == 0) then

      ! allocate space for ring and buffer
      allocate(rings2(nvars,amax,6,jmax,kmax,2))
      allocate(scatter_buffs2(nvars,amax,6,js,ks,2*pey*pez))

      ! read ring data from file
      infilename = 'output/' // trim(prefix) // '.s1_' // 'boundary.dat'
      open(unit=4,file=infilename,status='old',form='unformatted')

      do y = 1, 2
        do k = 1, kmax
          do j = 1, jmax
            do i = 1, 6
              do n = 1, amax

                read(4) rings2(:,n,i,j,k,y)

              enddo
            enddo
          enddo
        enddo
      enddo


      ! increment through processors
      ii = 1
      do y = 1, 2
        do pek = 0, pez-1
          do pej = 0, pey-1

            ! and sort data into scatter_buff
            do k  = 1, ks
              kk  = pek*ks+k
              do j  = 1, js
                jj  = pej*js+j
                do i  = 1, 6
                  do n  = 1, amax
                    do nv = 1, nvars

                      scatter_buffs2(nv,n,i,j,k,ii) = rings2(nv,n,i,jj,kk,y)

                    enddo
                  enddo
                enddo
              enddo
            enddo

          ii = ii + 1
          enddo
        enddo
      enddo


      deallocate(rings2)
    endif

    ! scatter data
    scattersize = nvars*amax*6*js*ks
    ALLOCATE(myring(nvars,amax,6,js,ks))

    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    call MPI_SCATTER(scatter_buffs2,  &
                     scattersize,		  &
                     VH1_DATATYPE,	  &
                     myring,			    &
                     scattersize,		  &
                     VH1_DATATYPE,	  &
                     0,				        &
                     MPI_COMM_WORLD,	&
                     mpierr)
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    if (mype == 0) deallocate(scatter_buffs2)
  else ! ------------------------------------------------------------

    amax = 512

    ! set up alpha array
    allocate(alpha(amax))
    dalpha = (x2max/x1min)**(1.0/float(amax-1)) - 1.0 
    alpha(1) = x1min
    do n = 2, amax
      alpha(n) = alpha(n-1)*(1.0+dalpha)
    enddo

    ! load data into ring array
    if (mype == 0) then

      ! allocate space for ring and buffer
      ALLOCATE(rings3        (nvars,amax,jmax,kmax,    2))
      allocate(scatter_buffs3(nvars,amax,js,ks,pey*pez*2))

      ! read ring data from file
      infilename = 'output/' // trim(prefix) // '.s2_' // 'boundary.dat'
      open(unit=4,file=infilename,status='old',form='unformatted')

      do y = 1, 2
        do k = 1, kmax
          do j = 1, jmax
            do n = 1, amax
              read(4) rings3(:,n,j,k,y)
            enddo
          enddo
        enddo
      enddo
      close(4)


      ! increment through processors
      ii = 1
      do y = 1, 2
        do pek = 0, pez-1
          do pej = 0, pey-1

            ! and sort data into scatter_buff
            do k  = 1, ks
              kk  = pek*ks+k
              do j  = 1, js
                jj  = pej*js+j
                do n  = 1, amax
                  do nv = 1, 6

                    scatter_buffs3(nv,n,j,k,ii) = rings3(nv,n,jj,kk,y)

                  enddo
                enddo
              enddo
            enddo

          ii = ii + 1
          enddo
        enddo
      enddo


      deallocate(rings3)
    endif

    ! scatter data
    scattersize =  nvars*amax*js*ks
    ALLOCATE(myray(nvars,amax,js,ks))

    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    call MPI_SCATTER(scatter_buffs3,	&
                     scattersize,		  &
                     VH1_DATATYPE,	  &
                     myray,			      &
                     scattersize,		  &
                     VH1_DATATYPE,	  &
                     0,				        &
                     MPI_COMM_WORLD,	&
                     mpierr)
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    if (mype == 0) deallocate(scatter_buffs3)



    ! outer bound
    if (mype == 0) then

      allocate(bound3        (2,jmax,kmax,    2))
      allocate(scatter_buffs4(2,js,ks,pey*pez*2))

      infilename = 'output/' // trim(prefix) // '.s2_avgbound.dat'
      open(unit=4, file=infilename, status='old', form='unformatted')

      do y = 1, 2
        do k = 1, kmax
          do j = 1, jmax
            do nv = 1, 2
              read(4) bound3(nv,j,k,y)
            enddo
          enddo
        enddo
      enddo
      close(4)


      ii = 1
      do y = 1, 2
        do pek = 0, pez-1
          do pej = 0, pey-1

            do k = 1, ks
              kk = pek*ks+k
              do j = 1, js
                jj = pej*js+j
                do nv = 1, 2
                  scatter_buffs4(nv,j,k,ii) = bound3(nv,jj,kk,y)
                enddo
              enddo
            enddo

          ii = ii + 1
          enddo
        enddo
      enddo

      deallocate(bound3)
    endif

    ! scatter data
    scattersize = 2*js*ks
    ALLOCATE(mybound(2,js,ks))

    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    call MPI_SCATTER(scatter_buffs4,  &
                     scattersize,     &
                     VH1_DATATYPE,    &
                     mybound,         &
                     scattersize,     &
                     VH1_DATATYPE,    &
                     0,               &
                     MPI_COMM_WORLD,  &
                     mpierr)
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    if (mype == 0) deallocate(scatter_buffs4)

  endif


return
end

