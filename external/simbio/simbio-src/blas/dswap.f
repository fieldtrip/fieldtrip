c@**********************************************************************
c                                                                      *
c     subroutine dswap (n,sx,incx,sy,incy)                             *
c                                                                      *
c     interchanges two vektors.                                        *
c     uses unrolled loops for increments equal to 1.                   *
c     jack dongarra, linpack, 3/11/78                                  *
c                                                                      *
c***********************************************************************
      subroutine dswap (n,sx,incx,sy,incy)
c
      integer i,incx,incy,ix,iy,m,mp1,n
      double precision sx(*),sy(*),stemp
c
      if (n .ge. 1) then
         if (incx .ne. 1 .or. incy .ne. 1) then
c
c...........code for unequal increments or increments not equal to 1
c
            ix = 1
            iy = 1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            if (incy .lt. 0) iy = (-n+1)*incy + 1
            do i = 1, n
               stemp  = sx(ix)
               sx(ix) = sy(iy)
               sy(iy) = stemp
               ix = ix + incx
               iy = iy + incy
            end do
         else
c
c...........code for both increments equal to 1
c
c...........clean-up loop
c
            m = mod(n,3)
            do i = 1, m
               stemp = sx(i)
               sx(i) = sy(i)
               sy(i) = stemp
            end do
            if (n .ge. 3) then
               mp1 = m + 1
               do i = mp1, n, 3
                  stemp = sx(i)
                  sx(i) = sy(i)
                  sy(i) = stemp
                  stemp = sx(i + 1)
                  sx(i + 1) = sy(i + 1)
                  sy(i + 1) = stemp
                  stemp     = sx(i + 2)
                  sx(i + 2) = sy(i + 2)
                  sy(i + 2) = stemp
               end do
            end if
         end if
      end if
      return
      end
