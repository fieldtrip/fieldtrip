c@**********************************************************************
c                                                                      *
c     function idmax (n,sx,incx)                                       *
c                                                                      *
c     finds the index of element having max. value.                    *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      function idmax (n,sx,incx)
c
      integer idmax,i,incx,ix,n
      double precision sx(*),smax
c
      idmax = 0
      if (n .ge. 1) then
         idmax = 1
         if (n .ge. 2) then
            if (incx .ne. 1) then
c
c..............code for increment not equal to 1
c
               ix = 1
               smax = sx(1)
               ix = ix + incx
               do i = 2, n
                  if (sx(ix) .gt. smax) then
                     idmax = i
                     smax = sx(ix)
                  end if
                  ix = ix + incx
               end do
            else
c
c..............code for increment equal to 1
c
               smax = sx(1)
               do i = 2, n
                  if (sx(i) .gt. smax) then
                     idmax = i
                     smax = sx(i)
                  end if
               end do
            end if
         end if
      end if
      return
      end
