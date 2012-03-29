c***********************************************************************
c                                                                      *
c     subroutine drot (n,sx,incx,sy,incy,c,s)                          *
c                                                                      *
c     applies a plane rotation.                                        *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      subroutine drot (n,sx,incx,sy,incy,c,s)

      double precision sx(*),sy(*),stemp,c,s
      integer i,incx,incy,ix,iy,n

      if (n .ge. 1) then
         if (incx .ne. 1 .or. incy .ne. 1) then

c........code for unequal increments or equal increments not equal
c........to 1

            ix = 1
            iy = 1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            if (incy .lt. 0) iy = (-n+1)*incy + 1
            do i = 1, n
              stemp = c*sx(ix) + s*sy(iy)
              sy(iy) = c*sy(iy) - s*sx(ix)
              sx(ix) = stemp
              ix = ix + incx
              iy = iy + incy
            end do
         else

c...........code for both increments equal to 1

            do i = 1, n
              stemp = c*sx(i) + s*sy(i)
              sy(i) = c*sy(i) - s*sx(i)
              sx(i) = stemp
            end do
         end if
      end if
      return
      end
