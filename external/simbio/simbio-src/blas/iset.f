c@**********************************************************************
c                                                                      *
c     subroutine iset (n,sa,sx,incx)                                   *
c                                                                      *
c     constant times a vector                                          *
c     uses unrolled loop for increments equal to one.                  *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      subroutine iset (n,sa,sx,incx)
c
      integer i,incx,n,m,mp1,ix
      integer sx(*),sa
c
      if (n .ge. 1) then
         if (incx .ne. 1) then
c
c...........code for increment not equal to 1
c
            ix = 1
            if (incx .lt. 0) ix = (-n+1)*incx + 1
            do i = 1, n
               sx(ix) = sa
               ix = ix + incx
            end do
         else
c
c...........code for both increments equal to 1
c
c...........clean-up loop
c
            m = mod(n,4)
            do i = 1, m
               sx(i) = sa
            end do
            if (n .ge. 4) then
               mp1 = m + 1
               do i = mp1, n, 4
                  sx(i)   = sa
                  sx(i+1) = sa
                  sx(i+2) = sa
                  sx(i+3) = sa
               end do
            end if
         end if
      end if
      return
      end
