c***********************************************************************
c                                                                      *
c     subroutine srotg (sa,sb,c,s)                                     *
c                                                                      *
c     construct givens plane rotation.                                 *
c     jack dongarra, linpack, 3/11/78.                                 *
c                    modified 9/27/86.                                 *
c                                                                      *
c***********************************************************************
      subroutine srotg (sa,sb,c,s)

      real sa,sb,c,s,roe,scale,r,z

      roe = sb
      if (abs(sa) .gt. abs(sb)) roe = sa
      scale = abs(sa) + abs(sb)
      if(scale .eq. 0.) then
         c = 1.
         s = 0.
         r = 0.
      else
         r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
         r = sign(1.,roe)*r
         c = sa/r
         s = sb/r
      end if
      z = s
      if (abs(c) .gt. 0. .and. abs(c) .le. s) z = 1./c
      sa = r
      sb = z
      return
      end
