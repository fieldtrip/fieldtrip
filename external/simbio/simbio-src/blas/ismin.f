c@**********************************************************************
c                                                                      *
c     function ismin (n,sx,incx)                                       *
c                                                                      *
c     finds the index of element having min. value.                    *
c     jack dongarra, linpack, 3/11/78.                                 *
c                                                                      *
c***********************************************************************
      function ismin (n,sx,incx)
c
      integer ismin,i,incx,ix,n
      real    sx(*),smin
c
      ismin = 0
      if (n .ge. 1) then
         ismin = 1
         if (n .ge. 2) then
            if (incx .ne. 1) then
c
c..............code for increment not equal to 1
c
               ix = 1
               smin = sx(1)
               ix = ix + incx
               do i = 2,n
                  if (sx(ix) .lt. smin) then
                     ismin = i
                     smin = sx(ix)
                  end if
                  ix = ix +incx
               end do
            else
c
c..............code for increment equal to 1
c
               smin = sx(1)
               do i = 2,n
                  if (sx(i) .lt. smin) then
                     ismin = i
                     smin = sx(i)
                  end if
               end do
            end if
         end if
      end if
      return
      end
