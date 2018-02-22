subroutine radiate

! use semi-implicit method to calculate the energy loss in each zone due to radiative cooling and heating
! dTemp = .5 * anrm*zrho*(Lambda(oldtemp) + Lambda(newtemp))
! iterate a few times, replacing netmin with tmpnew + dTemp

! after remap? call ppmlr->radiate in sweepx
! use truncated power law, need to produce a uniform cooling time
! agb star paper - what did they use

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS

REAL, DIMENSION(maxsweep) :: tmpold, tmpnew, qloss, xi, anrm, ff, fl, tmpf, tmpl, dtmp
REAL :: const1, const2, hdt, Tic, compton, xheat, cooling, angle, swap, hpi
REAL :: alpha, const1inv, df, temp, Tmax, Tmin, Tx, txm, xionmax, xionmin, xm2
INTEGER :: mithat


real :: tempw, avgmass, boltzmann

!-----------------------------------------------------------------------

Tic    = 2.9e99  !  this is hard-coded into some of the constants
tempw   = 1.0e-99  !  code this in

avgmass   = 1.67e-24
boltzmann = 1.67e-16

! what is avgmass?

const2 = gamm/avgmass/boltzmann*dt
const1 = avgmass/boltzmann
const1inv = 1. / const1
Tx     = 4. * Tic
!txm    = 2. * xlum * avgmass
hpi    = 0.33 * pi

Tmax     = Tic
Tmin     = tempw 
!xionmin  = 1.0
!xionmax  = 1.0e+05

! PV=nkT
! P = rho/mp*k*T
! T = C * P / rho -> C = mp/k


!if (.false.) then
!! Calculate xi 
!xm2 = xlum * avgmass * 2.0
!do n = nmin-3, nmax+3
!  xi(n) = min(max(xm2/(r(n)*radns(n)**2), xionmin), xionmax)
!
!  alpha = acos((xnstar-xa0(n)*stheta*cosphi)/radns(n))
!  if((alpha < subtend) .and. (abs(phi) > hpi)) xi(n) = xionmin   ! in the shadow of the primary
!
!  tmpold(n) = const1 * p(n) / r(n)				!!!! old temperature
!  tmpold(n) = min(max(tmpold(n),Tmin),Tmax)			!!!! old temperature
!
!  anrm  (n) = const2 * r(n) 
!enddo
!endif

do n = nmin-3, nmax+3
   xi    (n) = 1.0      ! ignore for now
   tmpold(n) = const1 * p(n) / r(n)
   tmpold(n) = min(max(tmpold(n),Tmin),Tmax)
   anrm  (n) = const2 * r(n)
enddo


! First guess at heating/cooling rate, then iterate using 
! Secant Method from Numerical Recipes

do n = nmin, nmax

     tmpl(n) = tmpold(n)
     temp    = tmpl(n)

     compton = 0.!3.56e-35 * xi(n) * (Tic - temp)
     xheat   = 0.!1.29e-29 * sqrt(sqrt(xi(n))/temp) * (Tx - temp)
     cooling = 2e-52!1.70e-18 * exp(-1.3e5/temp)/(xi(n)*sqrt(temp)) + 1.00e-24 + 3.30e-27 * sqrt(temp)
     qloss(n)= anrm(n) * (compton + xheat - cooling)   

     fl  (n) = - qloss(n)
     tmpf(n) = tmpl(n) + qloss(n)
     tmpf(n) = min(max(tmpf(n),Tmin),Tmax)
     temp    = tmpf(n)

     compton = 0.!3.56e-35*xi(n) * (Tic - temp)
     xheat   = 0.!1.29e-29*sqrt(sqrt(xi(n))/temp) * (Tx - temp)
     cooling = 2e-52!1.70e-18 * exp(-1.3e5/temp)/(xi(n)*sqrt(temp)) + 1.00e-24 + 3.30e-27 * sqrt(temp)
     qloss(n)= anrm(n) * (compton + xheat - cooling)

     ff  (n) = tmpf(n) - tmpold(n) - qloss(n)

enddo

do n = nmin, nmax
  if(abs(fl(n)) < abs(ff(n))) then
   swap    = tmpl(n)
   tmpl(n) = tmpf(n)
   tmpf(n) = swap
   swap    = fl(n)
   fl(n)   = ff(n)
   ff(n)   = swap
  endif
enddo

do n = nmin, nmax

    ! REPEAT THE FOLLOWING UNTIL CONVERGED  ######################
    do mithat = 1, 5

      df      = ff(n) - fl(n)
      if (df==0.0) df = small   
      dtmp(n) = (tmpl(n)-tmpf(n))*ff(n)/df
      tmpl(n) = tmpf(n)
      fl  (n) = ff(n)
      tmpf(n) = tmpf(n) + dtmp(n)						!!!! final temp
      tmpf(n) = min(max(tmpf(n),Tmin),Tmax)				!!!! final temp
      temp    = tmpf(n)							!!!! final temp

      compton = 0.!3.56e-35*xi(n) * (Tic - temp)
      xheat   = 0.!1.29e-29*sqrt(sqrt(xi(n))/temp) * (Tx - temp)
      cooling = 2e-52!1.70e-18 * exp(-1.3e5/temp)/(xi(n)*sqrt(temp)) + 1.00e-24 + 3.30e-27 * sqrt(temp)
      qloss(n)= anrm(n) * (compton + xheat - cooling)   		!!!! heat loss

      ff  (n) = tmpf(n) - tmpold(n) - qloss(n)				!!!!  ???

    enddo
    ! #############################################################

     tmpnew(n) = tmpold(n) + qloss(n)					!!!! new temperature
     tmpnew(n) = min(max(tmpnew(n),Tmin),Tmax)				!!!! new temperature
     p(n) = tmpnew(n) * r(n) * const1inv					!!!! new pressure
     !s(n) = qloss(n)	

enddo

! 1.70e-18*exp(-1.3e5/temp)*0.1/(xi(n) * sqrt(temp)) + 1.00e-24*0.1 + 3.30e-27 * sqrt(temp)
return
end