subroutine boundaryR
! =========================================================
! interpolate ring data based on alpha=omega*time
! =========================================================

  use global
  use zone

implicit none

  ! ...

  integer :: i,j,k
  integer :: nbelow, nabove
  real :: dalpha, c1, c2, arot, rad, dfrac
  real, dimension(6) :: c3, c4
  integer, dimension(6) :: mabove, mbelow

! =========================================================

  if (step == 1) then
  
  else if (step == 2) then !-----------------------------------------
    ! determine interpolation coefficients

    arot = omega*time
    do while (arot > 2.0*pi)
      arot =  arot - 2.0*pi
    enddo

    !arot = 2.0*pi-arot

    dalpha = alpha(2) - alpha(1)

    nbelow = (arot-alpha(1))/dalpha + 1
    nabove = nbelow                 + 1

    ! account for periodic data
    if (nabove == amax+1) then
      nabove = 1
      c1 = (alpha(nabove)-arot+2.0*pi)/dalpha
      c2 = (arot-alpha(nbelow))/dalpha
    else if (nbelow == 0) then
      nbelow = amax
      nabove = 1
      c1 = (alpha(nabove)-arot)/dalpha
      c2 = (arot-alpha(nbelow)+2.0*pi)/dalpha
    else
      c1 = (alpha(nabove)-arot)/dalpha
      c2 = (arot-alpha(nbelow))/dalpha
    endif

    ! load boundary conditions from ring
    do k = 1, ks
      do j = 1, js
        do i = 1, 6

          rin(i,j,k) = myring(1,nbelow,i,j,k)*c1 + myring(1,nabove,i,j,k)*c2
          pin(i,j,k) = myring(2,nbelow,i,j,k)*c1 + myring(2,nabove,i,j,k)*c2
          cin(i,j,k) = myring(3,nbelow,i,j,k)*c1 + myring(3,nabove,i,j,k)*c2
          uin(i,j,k) = myring(4,nbelow,i,j,k)*c1 + myring(4,nabove,i,j,k)*c2
          vin(i,j,k) = myring(5,nbelow,i,j,k)*c1 + myring(5,nabove,i,j,k)*c2
          win(i,j,k) = myring(6,nbelow,i,j,k)*c1 + myring(6,nabove,i,j,k)*c2

        enddo
      enddo
    enddo
  else ! ------------------------------------------------------------

    ! prepare outer boundary based on current rmax
    do i = 1, 6
      rad = zxc(imax) + zdx(imax)*i

      if (rad < alpha(amax)) then
        nabove = 2
        do while (alpha(nabove) < rad)
          nabove = nabove + 1
        enddo
        nbelow = nabove - 1

        c1 = (alpha(nabove)-rad)/(alpha(nabove)-alpha(nbelow))
        c2 = (rad-alpha(nbelow))/(alpha(nabove)-alpha(nbelow))

        do k = 1, ks
          do j = 1, js
            rin(i,j,k) = myray(1,nbelow,j,k)*c1 + myray(1,nabove,j,k)*c2
            pin(i,j,k) = myray(2,nbelow,j,k)*c1 + myray(2,nabove,j,k)*c2
            cin(i,j,k) = 0.0
            uin(i,j,k) = 0.0
            vin(i,j,k) = 0.0
            win(i,j,k) = 0.0
          enddo
        enddo
      else
        dfrac = (alpha(amax) / rad)**2

        do k = 1, ks
          do j = 1, js
            rin(i,j,k) = mybound(1,j,k)*dfrac
            pin(i,j,k) = mybound(2,j,k)*dfrac
            cin(i,j,k) = 0.0
            uin(i,j,k) = 0.0
            vin(i,j,k) = 0.0
            win(i,j,k) = 0.0
          enddo
        enddo
      endif
    enddo
  endif

return
end
