c@**********************************************************************
c                                                                      *
c     subroutine drotg (sa,sb,c,s)                                     *
c                                                                      *
c     construct givens plane rotation.                                 *
c     jack dongarra, linpack, 3/11/78.                                 *
c                    modified 9/27/86.                                 *
c                                                                      *
c***********************************************************************
      subroutine drotg (sa,sb,c,s)

      double precision sa,sb,c,s,roe,scale,r,z

      roe = sb
      if (abs(sa) .gt. abs(sb)) roe = sa
      scale = abs(sa) + abs(sb)
      if (scale .eq. 0.d0) then
         c = 1.d0
         s = 0.d0
         r = 0.d0
      else
         r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
         r = sign(1.d0,roe)*r
         c = sa/r
         s = sb/r
      end if
      z = s
      if (abs(c) .gt. 0.d0 .and. abs(c) .le. s) z = 1.d0/c
      sa = r
      sb = z
      return
      end
