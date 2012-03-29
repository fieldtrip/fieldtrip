c@**********************************************************************
c                                                                      *
c     subroutine iadd (n, sa, sx, incx)                                *
c                                                                      *
c     add a scalar da to each component of a vector                    *
c     uses unrolled loops for increment equal to one.                  *
c     jack dongarra, linpack, 3/11/78                                  *
c                                                                      *
c***********************************************************************
      subroutine iadd (n, sa, sx, incx)

      integer n,incx,nincx,i,m,mp1
      integer sx(*),sa

      if (n .gt. 0) then
         if (incx .ne. 1) then

c...........code for increment not equal to 1

            nincx = n*incx
            do i=1, nincx, incx
               sx(i) = sa + sx(i)
            end do
         else

c...........code for increment equal to 1

            m = mod(n,5)

c...........clean-up loop

            do i=1, m
               sx(i) = sa + sx(i)
            end do
            if (n .ge. 5) then
               mp1 = m + 1
               do i=mp1, n, 5
                  sx(i  ) = sa + sx(i  )
                  sx(i+1) = sa + sx(i+1)
                  sx(i+2) = sa + sx(i+2)
                  sx(i+3) = sa + sx(i+3)
                  sx(i+4) = sa + sx(i+4)
               end do
            end if
         end if
      end if
      return
      end
