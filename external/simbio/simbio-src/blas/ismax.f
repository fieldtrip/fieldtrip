c@**********************************************************************
c                                                                      *
c     function ismax (n,sx,incx)                                       *
c                                                                      *
c     finds the index of element having max. value.                    *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      function ismax (n,sx,incx)
c
      integer ismax,i,incx,ix,n
      real    sx(*),smax
c
      ismax = 0
      if (n .ge. 1) then
         ismax = 1
         if (n .ge. 2) then
            if (incx .ne. 1) then
c
c..............code for increment not equal to 1
c
               ix = 1
               smax = sx(1)
               ix = ix + incx
               do i = 2,n
                  if (sx(ix) .gt. smax) then
                     ismax = i
                     smax = sx(ix)
                  end if
                  ix = ix +incx
               end do
            else
c
c..............code for increment equal to 1
c
               smax = sx(1)
               do i = 2,n
                  if (sx(i) .gt. smax) then
                     ismax = i
                     smax = sx(i)
                  end if
               end do
            end if
         end if
      end if
      return
      end
