c@**********************************************************************
c                                                                      *
c     subroutine dcopy (n,sx,incx,sy,incy)                             *
c                                                                      *
c     copies a vektor x, to a vector y.                                *
c     uses unrolled loops for increments equal to 1.                   *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      subroutine dcopy (n,sx,incx,sy,incy)
c
      integer i,incx,incy,ix,iy,m,mp1,n
      double precision  sx(*),sy(*)
c
      if (n .ge. 1) then
         if (incx .ne. 1 .or. incy .ne. 1) then
c
c...........code for unequal increments or equal increments
c...........not equal to 1
c
            ix=1
            iy=1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            if (incy .lt. 0) iy = (-n+1)*incy + 1
            do i = 1, n
               sy(iy) = sx(ix)
               ix = ix + incx
               iy = iy + incy
            end do
         else
c
c...........code for both increments equal to 1
c
c...........clean-up loop
c
            m = mod(n,7)
            do i = 1,m
               sy(i) = sx(i)
            end do
            if (n .ge. 7) then
               mp1 = m + 1
               do i = mp1, n, 7
                  sy(i  ) = sx(i  )
                  sy(i+1) = sx(i+1)
                  sy(i+2) = sx(i+2)
                  sy(i+3) = sx(i+3)
                  sy(i+4) = sx(i+4)
                  sy(i+5) = sx(i+5)
                  sy(i+6) = sx(i+6)
               end do
            end if
         end if
      end if
      return
      end
