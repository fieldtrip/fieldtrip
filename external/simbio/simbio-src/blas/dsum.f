c@**********************************************************************
c                                                                      *
c     function dsum (n,sx,incx)                                        *
c                                                                      *
c     takes the sum of the values.                                     *
c     uses unrolled loops for increment equal to one.                  *
c     jack dongarra, linpack, 3/11/78                                  *
c                                                                      *
c***********************************************************************
      function dsum (n,sx,incx)
c
      integer i,incx,nincx,n,m,mp1
      double precision sx(*),dsum,stemp
c
      stemp = 0.d0
      if (n .ge. 1) then
         if (incx .ne. 1) then
c
c...........code for increment not equal to 1
c
            nincx = n*incx
            do i = 1,nincx,incx
               stemp = stemp + sx(i)
            end do
         else
c
c...........code for increment equal to 1
c
c...........clean-up loop
c
            m = mod(n,6)
            do i = 1,m
               stemp = stemp + sx(i)
            end do
            if (n .ge. 6) then
               mp1 = m + 1
               do i = mp1, n, 6
                  stemp =   stemp 
     &                    + sx(i  ) + sx(i+1)
     &                    + sx(i+2) + sx(i+3)
     &                    + sx(i+4) + sx(i+5)
               end do
            end if
         end if
      end if
      dsum = stemp
      return
      end
