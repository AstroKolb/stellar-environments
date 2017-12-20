
! ===================================================================
! 3-D Perlin Noise Generator
! ===================================================================



! ===================================================================


subroutine initialize_noise

   use noise

   implicit none

   integer                                       :: i, j, rn
   real                                          :: amplitude
   logical                                       :: test
! -------------------------------------------------------------------

   ! randomization
   call srand(seed)

   pp   = 0
   perm = 0

   maxamplitude = 0.0
   amplitude    = 1.0

   i = 1
   do while (i < tiledim+1)
      rn = tiledim * rand() + 1

      test = .true.
      do j = 1, tiledim				! loop through each element of perm
         if (rn == perm(j)) test = .false.
      enddo

      if (test) then
         perm(i) = rn
         i = i + 1
      endif
   enddo
   
   do i = 1, tiledim
      pp(i)      = perm(i)
      pp(i+tiledim) = pp(i)
   enddo

   ! determine max amplitude
   do i = 1, octaves
      maxamplitude = maxamplitude + amplitude
      amplitude    = amplitude  * persistence
   enddo

   tilesize = float(tile) / repeats
   tilesize = tilesize    / tiledim

return
end

! ===================================================================


function fade(t) result(r)

   real, intent(in)                              :: t
   real                                          :: r
! -------------------------------------------------------------------

   r = t * t * t * (t * (t * 6 - 15) + 10)

end function fade

! ===================================================================


!function lerp(t, a, b) result(r)
function lerp(a, b, t) result(r)

   real, intent(in)                              :: t, a, b
   real                                          :: r
! -------------------------------------------------------------------

   r = a + t * (b - a)

end function lerp

! ===================================================================


function grad(hash, x, y, z) result(r)

   integer, intent(in)                           :: hash
   real,    intent(in)                           :: x, y, z
   real                                          :: r

   integer                                       :: h
   real                                          :: u, v, f, s
! -------------------------------------------------------------------

   h = iand(hash, 15)

   if (h < 8) then
      u = x
   else
      u = y
   endif

   if (h < 4) then
      v = y
   else if ((h == 12) .or. (h == 14)) then
      v = x
   else
      v = z
   endif

   if (iand(h, 1) == 0) then
      f =  u
   else
      f = -u
   endif

   if (iand(h, 2) == 0) then
      s =  v
   else
      s = -v
   endif

   r = f + s

end function grad

! ===================================================================


function get_noise(x, y, z) result(r)

   use noise

   real, intent(in)                              :: x, y, z
   real                                          :: r

   real                                          :: lerp, grad, fade
   
   real                                          :: i, j, k
   real                                          :: u, v, w

   integer                                       :: xx, yy, zz
   integer                                       :: A, AA, AB, B, BA, BB
   integer                                       :: AAA, ABA, AAB, ABB
   integer                                       :: BAA, BBA, BAB, BBB

   real                                          :: x1, x2, y1, y2
! -------------------------------------------------------------------

   xx = iand(int(x), tiledim-1)		! xi
   yy = iand(int(y), tiledim-1)		! yi
   zz = iand(int(z), tiledim-1)		! zi

   i  = x - int(x)		! xf
   j  = y - int(y)		! yf
   k  = z - int(z)		! zf

   u  = fade(i)
   v  = fade(j)
   w  = fade(k)

   AAA = pp(pp(pp(xx  )+yy  )+zz  )
   ABA = pp(pp(pp(xx  )+yy+1)+zz  )
   AAB = pp(pp(pp(xx  )+yy  )+zz+1)
   ABB = pp(pp(pp(xx  )+yy+1)+zz+1)
   BAA = pp(pp(pp(xx+1)+yy  )+zz  )
   BBA = pp(pp(pp(xx+1)+yy+1)+zz  )
   BAB = pp(pp(pp(xx+1)+yy  )+zz+1)
   BBB = pp(pp(pp(xx+1)+yy+1)+zz+1)


   x1  = lerp( grad(AAA,i  ,j  ,k  ), &
               grad(BAA,i-1,j  ,k  ), &
               u                      )
   x2  = lerp( grad(ABA,i  ,j-1,k  ), &
               grad(BBA,i-1,j-1,k  ), &
               u                      )

   y1  = lerp(x1,x2,v)


   x1  = lerp( grad(AAB,i  ,j  ,k-1), &
               grad(BAB,i-1,j  ,k-1), &
               u                      )
   x2  = lerp( grad(ABB,i  ,j-1,k-1), &
               grad(BBB,i-1,j-1,k-1), &
               u                      )

   y2  = lerp(x1,x2,v)


   r   = lerp(y1,y2,w)


   

end function get_noise

! ===================================================================


function gen_noise(x, y, z) result(r)

   use noise

   real, intent(in)                              :: x, y, z
   real                                          :: r

   real                                          :: fade, lerp, grad, get_noise

   real                                          :: amplitude
   real                                          :: frequency
   real                                          :: sc, color, grey

   integer                                       :: i, j, k
! -------------------------------------------------------------------

   frequency = 1.0
   amplitude = 1.0
   color     = 0.0
   sc        = float(tile) / tilesize

   grey      = 0.0

   do i = 1, octaves
      sc = sc * frequency
      grey = get_noise(sc*x/x_lscale+1000, sc*y/y_lscale+1000, sc*z/z_lscale+1000)
      grey = (1.0 + grey) / 2.0
      grey = grey * amplitude

      color = color + grey
      frequency = frequency * 2.0
      amplitude = amplitude * persistence
   enddo

   r = color / maxamplitude

end function gen_noise
