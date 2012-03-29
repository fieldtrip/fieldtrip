c@**********************************************************************
c                                                                      *
c     function isamax (n,sx,incx)                                      *
c                                                                      *
c     finds the index of element having max. absolute value.           *
c     jack dongarra, linpack, 3/11/78.                                 *
c     modified to correct problem with negative increments, 9/29/88.   *
c                                                                      *
c***********************************************************************
      function isamax (n,sx,incx)

      integer isamax,i,incx,ix,n
      real    sx(*),smax

      isamax = 0
      if (n .ge. 1) then
         isamax = 1
         if (n .ge. 2) then
            if (incx .ne. 1) then

c..............code for increment not equal to 1

               ix = 1
               if (incx .lt. 0) ix = (-n+1)*incx + 1
               smax = abs(sx(ix))
               ix = ix + incx
               do i = 2,n
                  if (abs(sx(ix)) .gt. smax) then
                     isamax = i
                     smax = abs(sx(ix))
                  end if
                  ix = ix +incx
               end do
            else

c..............code for increment equal to 1

               smax = abs(sx(1))
               do i = 2,n
                  if (abs(sx(i)) .gt. smax) then
                     isamax = i
                     smax = abs(sx(i))
                  end if
               end do
            end if
         end if
      end if
      return
      end
