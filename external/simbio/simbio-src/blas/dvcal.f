c@**********************************************************************
c                                                                      *
c     subroutine dvcal (n, sa, sx, incx, sy, incy)                     *
c                                                                      *
c     multiply a vector sx by a scalar sa and store the result in      *
c     another vector sy,  sy = sa * sx.                                *
c                                                                      *
c***********************************************************************
      subroutine dvcal (n, sa, sx, incx, sy, incy)
c
      integer n,incx,incy,ix,iy,i,m,mp1
      double precision sx(*),sy(*),sa
c
      if (n .gt. 0) then
         if (incx .ne. 1 .or. incy .ne. 1) then
c
c...........code for unequal increments or equal
c...........increments not equal to 1
c
            ix = 1
            iy = 1
            do i=1, n
               sy(iy) = sa*sx(ix)
               ix = ix + incx
               iy = iy + incy
            end do
         else
c
c...........code for both increments equal to 1
c
            m = mod(n,4)
c
c...........clean-up loop
c
            do i=1, m
               sy(i) = sa*sx(i)
            end do
            if (n .ge. 4) then
               mp1 = m + 1
               do i=mp1, n, 4
                  sy(i) = sa*sx(i)
                  sy(i+1) = sa*sx(i+1)
                  sy(i+2) = sa*sx(i+2)
                  sy(i+3) = sa*sx(i+3)
               end do
            end if
         end if
      end if
      return
      end
