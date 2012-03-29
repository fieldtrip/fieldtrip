c@**********************************************************************
c                                                                      *
c     subroutine daxpy (n,da,dx,incx,dy,incy)                          *
c                                                                      *
c     constant times a vector plus a vector (dy(i)=dy(i)+sa*dx(i)).    *
c     uses unrolled loop for increments equal to one.                  *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      subroutine daxpy (n,sa,sx,incx,sy,incy)
c
      integer i,incx,incy,ix,iy,m,mp1,n
      double precision sx(*),sy(*),sa
c
      if (n .ge. 1 .and. sa .ne. 0.d0) then
         if (incx .ne. 1 .or. incy .ne. 1) then
c
c...........code for nonequal increments or equal increments
c...........not equal to 1
c
            ix = 1
            iy = 1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            if (incy .lt. 0) iy = (-n+1)*incy + 1
            if (sa .eq. 1.d0) then
               do i = 1, n
                  sy(iy) = sy(iy) + sx(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
            else
               do i = 1, n
                  sy(iy) = sy(iy) + sa*sx(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
            end if
         else
c
c...........code for both increments equal to 1
c
c...........clean-up loop
c
            m = mod(n,4)
            if (sa .eq. 1.d0) then
               do i = 1,m
                  sy(i) = sy(i) + sx(i)
               end do
               if (n .ge. 4) then
                  mp1 = m + 1
                  do i = mp1, n, 4
                     sy(i  ) = sy(i  ) + sx(i  )
                     sy(i+1) = sy(i+1) + sx(i+1)
                     sy(i+2) = sy(i+2) + sx(i+2)
                     sy(i+3) = sy(i+3) + sx(i+3)
                  end do
               end if
            else
               do i = 1,m
                  sy(i) = sy(i) + sa*sx(i)
               end do
               if (n .ge. 4) then
                  mp1 = m + 1
                  do i = mp1, n, 4
                     sy(i  ) = sy(i  ) + sa*sx(i  )
                     sy(i+1) = sy(i+1) + sa*sx(i+1)
                     sy(i+2) = sy(i+2) + sa*sx(i+2)
                     sy(i+3) = sy(i+3) + sa*sx(i+3)
                  end do
               end if
            end if
         end if
      end if
      return
      end
