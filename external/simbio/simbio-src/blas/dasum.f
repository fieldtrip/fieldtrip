c@**********************************************************************
c                                                                      *
c     function dasum (n,sx,incx)                                       *
c                                                                      *
c     takes the sum of the absolute values.                            *
c     uses unrolled loops for increment equal to one.                  *
c     jack dongarra, linpack, 3/11/78                                  *
c     modified to correct problem with negative increments, 9/29/88.   *
c                                                                      *
c***********************************************************************
      function dasum (n,sx,incx)

      integer n,incx,ix,i,m,mp1
      double precision sx(*),dasum,stemp

      stemp = 0.
      if (n .ge. 1) then
         if (incx .ne. 1) then

c...........code for increment not equal to 1

            ix = 1
            if( incx .lt. 0) ix = (-n+1)*incx + 1
            do i = 1,n
               stemp = stemp + abs(sx(ix))
            ix = ix + incx
            end do
         else

c...........code for increment equal to 1
c...........clean-up loop

            m = mod(n,6)
            do i = 1,m
               stemp = stemp + abs(sx(i))
            end do
            if (n .ge. 6) then
               mp1 = m + 1
               do i = mp1, n, 6
                  stemp = stemp
     &                    + abs(sx(i  )) + abs(sx(i+1))
     &                    + abs(sx(i+2)) + abs(sx(i+3))
     &                    + abs(sx(i+4)) + abs(sx(i+5))
               end do
            end if
         end if
      end if
      dasum = stemp
      return
      end
