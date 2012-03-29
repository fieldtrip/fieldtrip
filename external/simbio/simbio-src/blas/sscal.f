c@**********************************************************************
c                                                                      *
c     subroutine sscal (n,sa,sx,incx)                                  *
c                                                                      *
c     scales a vektor by a constant.                                   *
c     uses unrolled loops for increment equal to 1.                    *
c     jack dongarra, linpack, 3/11/78.                                 *
c     modified to correct problem with negative increments, 9/29/88.   *
c                                                                      *
c***********************************************************************
      subroutine sscal (n,sa,sx,incx)

      integer i,incx,m,mp1,n,ix
      real    sx(*),sa

      if (n .ge. 1) then
         if (incx .ne. 1) then

c...........code for increment not equal to 1

            ix = 1 
            if (incx .lt. 0) ix = (-n+1)*incx + 1 
            do i = 1, n 
              sx(ix) = sa*sx(ix)
              ix = ix + incx 
            end do
         else

c...........code for increment equal to 1
c...........clean-up loop

            m = mod(n,5)
            do i = 1,m
               sx(i) = sa*sx(i)
            end do
            if (n .ge. 5) then
               mp1 = m + 1
               do i = mp1, n, 5
                  sx(i  ) = sa*sx(i  )
                  sx(i+1) = sa*sx(i+1)
                  sx(i+2) = sa*sx(i+2)
                  sx(i+3) = sa*sx(i+3)
                  sx(i+4) = sa*sx(i+4)
               end do
            end if
         end if
      end if
      return
      end
