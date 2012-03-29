c@**********************************************************************
c                                                                      *
c     subroutine dadd (n, sa, sx, incx)                                *
c                                                                      *
c     add a scalar da to each component of a vector                    *
c     uses unrolled loops for increment equal to one.                  *
c     jack dongarra, linpack, 3/11/78                                  *
c     modified to correct problem with negative increments, 9/29/88.   *
c                                                                      *
c***********************************************************************
      subroutine dadd (n, sa, sx, incx)
c
      integer n,incx,nincx,i,m,mp1
      double precision  sx(*),sa
c
      if (n .ge. 1) then
         if (incx .ne. 1) then
c
c...........code for increment not equal to 1
c
            nincx = n*incx
            do i=1, nincx, incx
               sx(i) = sa + sx(i)
            end do
         else
c
c...........code for increment equal to 1
c
            m = mod(n,5)
c
c...........clean-up loop
c
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
