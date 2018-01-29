subroutine read_icons
! ===============================================
! load initial conditions into memory, then scale
! ===============================================

   use global
   use zone

implicit none

   include 'mpif.h'

   integer, parameter                  :: nchrs = 54
   integer, parameter                  :: nints = 1
   integer, parameter                  :: nrels = 19

   character     (len=nchrs)           :: icons_chrs
   integer, dimension(nints)           :: icons_ints
   real,    dimension(nrels)           :: icons_rels

   character(len=50)                   :: buffer, label
   integer                             :: line, pos, ios
   integer                             :: mpierr

! ===============================================


   ! load icons file
   if (jcol == 0) then
      line = 0
      ios  = 0
      open(unit=15, file='icons', status='old', form='formatted')
         do while (ios == 0)
            read(15, '(A)', iostat=ios) buffer
            if (ios == 0) then
               line = line + 1
               pos  = scan(buffer, ' ')
               label = buffer(1:pos)
               buffer = buffer(pos+1:)

               pos  = scan(buffer, '=')
               buffer = buffer(pos+1:)

               select case (label)
                  case ('prefix')
                     read(buffer, *, iostat=ios) prefix
                  case ('step')
                     read(buffer, *, iostat=ios) step
                  case ('gam')
                     read(buffer, *, iostat=ios) gam
                  case ('T_flr')
                     read(buffer, *, iostat=ios) T_flr
                  case ('Rstar')
                     read(buffer, *, iostat=ios) radP    !## change me
                  case ('x1min')
                     read(buffer, *, iostat=ios) x1min
                  case ('x1max')
                     read(buffer, *, iostat=ios) x1max
                  case ('x2max')
                     read(buffer, *, iostat=ios) x2max
                  case ('GMP')
                     read(buffer, *, iostat=ios) GMP
                  case ('TmpP')
                     read(buffer, *, iostat=ios) TmpP
                  case ('dmP')
                     read(buffer, *, iostat=ios) dmP
                  case ('GMS')
                     read(buffer, *, iostat=ios) GMS
                  case ('sep')
                     read(buffer, *, iostat=ios) sep
                  case ('radlaw')
                     read(buffer, *, iostat=ios) radlaw
                  case ('rad1')
                     read(buffer, *, iostat=ios) rad1
                  case ('rad2')
                     read(buffer, *, iostat=ios) rad2
                  case ('rad3')
                     read(buffer, *, iostat=ios) rad3
                  case ('kap1')
                     read(buffer, *, iostat=ios) kap1
                  case ('kap2')
                     read(buffer, *, iostat=ios) kap2
                  case ('kap3')
                     read(buffer, *, iostat=ios) kap3
                  case ('dmscal')
                     read(buffer, *, iostat=ios) dmscal
                  case ('x0min')
                     read(buffer, *, iostat=ios) x0min
                  case default   ! do nothing
               end select

            endif
         enddo
      close(15)

      ! buffer icons
      icons_rels(01) = gam
      icons_rels(02) = T_flr
      icons_rels(03) = radP
      icons_rels(04) = x1min
      icons_rels(05) = x1max
      icons_rels(06) = x2max
      icons_rels(07) = GMP
      icons_rels(08) = TmpP
      icons_rels(09) = dmP
      icons_rels(10) = GMS
      icons_rels(11) = sep
      icons_rels(12) = rad1
      icons_rels(13) = rad2
      icons_rels(14) = rad3
      icons_rels(15) = kap1
      icons_rels(16) = kap2
      icons_rels(17) = kap3
      icons_rels(18) = dmscal
      icons_rels(19) = x0min

      icons_chrs(1:50)  = prefix
      icons_chrs(51:54) = radlaw

      icons_ints(01) = step
   endif


   call MPI_BCAST( icons_rels, nrels, VH1_DATATYPE, 0, MPI_COMM_ROW, mpierr )
      gam = icons_rels(01)
      T_flr = icons_rels(02)
      radP  = icons_rels(03) 
      x1min = icons_rels(04) 
      x1max = icons_rels(05)
      x2max = icons_rels(06)
      GMP   = icons_rels(07)
      TmpP  = icons_rels(08)
      dmP   = icons_rels(09)
      GMS   = icons_rels(10)
      sep   = icons_rels(11)
      rad1  = icons_rels(12)
      rad2  = icons_rels(13)
      rad3  = icons_rels(14)
      kap1  = icons_rels(15)
      kap2  = icons_rels(16)
      kap3  = icons_rels(17)
      dmscal= icons_rels(18)
      x0min = icons_rels(19)


   call MPI_BCAST( icons_chrs, nchrs, MPI_CHARACTER, 0, MPI_COMM_ROW, mpierr )
      prefix = icons_chrs(1:50)
      radlaw = icons_chrs(51:54)


   call MPI_BCAST( icons_ints, nints, MPI_INTEGER, 0, MPI_COMM_ROW, mpierr )
      step = icons_ints(01)


! ===============================================
return 
end