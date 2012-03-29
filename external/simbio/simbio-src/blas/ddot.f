c@**********************************************************************
c                                                                      *
c     function ddot (n,sx,incx,sy,incy)                                *
c                                                                      *
c     forms the dot product of two vectors.                            *
c     uses unrolled loops for increments equal to one.                 *
c     jack dongarra, linpack, 3/11/78                                  *
c                                                                      *
c***********************************************************************
      function ddot (n,sx,incx,sy,incy)
c
      integer i,incx,incy,ix,iy,m,mp1,n
      double precision sx(*),sy(*),ddot,stemp
c
      stemp = 0.d0
      if (n .ge. 1) then
         if (incx .ne. 1 .or. incy .ne. 1) then
c
c...........code for unequal increments or equal increments
c...........not equal to 1
c
            ix = 1
            iy = 1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            if (incy .lt. 0) iy = (-n+1)*incy + 1
            do i = 1,n
               stemp = stemp + sx(ix)*sy(iy)
               ix = ix + incx
               iy = iy + incy
            end do
         else
c
c...........code for both increments equal to 1
c
c...........clean-up loop
c
            m = mod(n,5)
            do i = 1,m
               stemp = stemp + sx(i)*sy(i)
            end do
            if (n .ge. 5) then
               mp1 = m + 1
               do i = mp1, n, 5
                  stemp =   stemp 
     &                    + sx(i  )*sy(i  )
     &                    + sx(i+1)*sy(i+1)
     &                    + sx(i+2)*sy(i+2)
     &                    + sx(i+3)*sy(i+3)
     &                    + sx(i+4)*sy(i+4)
               end do
            end if
         end if
      end if
      ddot = stemp
      return
      end
