!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     numutl.f :
!                              -------------------
!     begin                : Mai 2000
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!    NeuroFEM license:
!    =================
!    Copyright (c) 2007 by 
!    Dr.Carsten Wolters, Dr.Alfred Anwander,  
!    Prof.Dr.H.Buchner, Prof.Dr.G.Knoll, Dr.Adrian Rienaecker, 
!    Rainer Beckmann. 
!
!    Permission is hereby granted, free of charge, to any person 
!    obtaining a copy of the NeuroFEM software and associated 
!    documentation files (the "Software"), to deal in the Software 
!    without restrictions, including without limitations the rights to 
!    use, copy, modify, merge, publish, distribute, sublicense, 
!    and/or sell copies of the Software, and to permit persons to 
!    whom the Software is furnished to do so, subject to the following 
!    conditions:
!
!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
!    OR OTHER DEALINGS IN THE SOFTWARE. THE ENTIRE RISK AS TO THE 
!    QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH YOU.
!
!    The above copyright notice, permission and warranty notices 
!    shall be included in all copies of the Software. A copyright 
!    notice of the contributing authors and the above permission 
!    and warranty notices shall be included in all source code 
!    files of the Software. 
!
!---->---1---------2---------3---------4---------5---------6---------7--<

      function detang(y,x)

      implicit double precision (a-h,o-z)
      parameter (pi=3.141592653589793d0)
      parameter (z0=0.d0,z2=2.d0)
c
      if( (y.eq.z0).and.(x.eq.z0) ) then
         detang=z0
         return
      end if
c
      detang=atan2(y,x)
      if (detang.lt.z0) detang=detang+z2*pi
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine iswap(i1,i2)
         idum=i2
         i2=i1
         i1=idum
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function ranknu(idum)
c
c     produziert gleichverteilte zufallszahlen nach D.E. KNUTH
c     aus "numerical recipes". PORTIERBAR!!!
c     "scheinbar" double precision
c
      double precision ranknu
      parameter ( mbig = 1000000000, mseed = 161803398,
     &            mz   = 0         , fac   = 1./mbig    )
c
      dimension ma( 55 )
      save ma, inext, inextp, iff
      data iff / 0 /
c
      if ( idum .lt. 0  .or.  iff .eq. 0 ) then
         iff      = 1
         mj       = mseed - iabs( idum )
         mj       = mod( mj, mbig )
         ma( 55 ) = mj
         mk       = 1
         do 10 i  = 1,54
            ii       = mod( 21*i, 55 )
            ma( ii ) = mk
            mk       = mj - mk
            if ( mk .lt. mz ) mk = mk + mbig
            mj = ma( ii )
 10      continue
         do 30 k = 1,4
            do 20 i = 1,55
               ma( i ) = ma( i ) - ma(  1 + mod( i+30, 55 )  )
               if( ma( i ) .lt. mz ) ma( i ) = ma( i ) + mbig
 20         continue
 30      continue
         inext = 0
         inextp = 31
         idum = 1
      end if
      inext = inext + 1
      if( inext .eq. 56 ) inext = 1
      inextp = inextp + 1
      if( inextp .eq. 56 ) inextp = 1
      mj = ma( inext ) - ma( inextp )
      if( mj .lt. mz ) mj = mj + mbig
      ma( inext ) = mj
      ranknu = dble(mj * fac)
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine averag(number,avrges,valnow,avenow,avenum,nitges)
      implicit double precision (a-h,o-z)
      dimension avrges(number)
c
      iposit = mod(nitges,number)
      if(iposit.eq.0) iposit=number
      aveout = avrges(iposit)
      avenow = (dble(nitges-1)*avenow+valnow)/dble(nitges)
      avenum = ( dble(nitges)*avenow - dble(nitges-number)*aveout)
     &         / dble(number)
      avrges(iposit)=avenow
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function valint(f2,f1,x2,x1,x)
      implicit double precision (a-h,o-z)
c
      valint=f1+(f2-f1)*(x-x1)/(x2-x1)
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
c---->---1---------2---------3---------4---------5---------6---------7--<
c---->---1---------2---------3---------4---------5---------6---------7--<
c Shell-Sort, sortiert ia in aufsteigender Reihenfolge
c Numerical Recipes, 2.Ed, S.323
c
      subroutine iasort(ia,n)
      dimension ia(n)
      inc=1
    1 continue
      inc=3*inc+1
      if(inc.le.n) goto 1
    2 continue
      inc=inc/3
      do i=inc+1,n
         k=ia(i)
         j=i
    3    continue
         if(ia(j-inc).gt.k)then
            ia(j)=ia(j-inc)
            j=j-inc
            if(j.le.inc) goto 4
            goto 3
         end if
    4    continue
         ia(j)=k
      end do
      if(inc.gt.1) goto 2
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine difquo(vecto1,vecto2,differ,npoin)
      implicit double precision (a-h,o-z)
      dimension vecto1(npoin),vecto2(npoin)
c
      faktor=1.d0/differ
      do 10 i=1,npoin
         vecto1(i)=(vecto1(i)-vecto2(i))*faktor
   10 continue
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine mat4x4(array1,array2,array3,n)
      implicit double precision (a-h,o-z)
      dimension array1(n,n),array2(n,n),array3(n,n)
c
      do 10 i=1,n
         do 20 j=1,n
            dum=z0
            do 30 k=1,n
               dum=dum+array1(i,k)*array2(k,j)
   30       continue
            array3(i,j)=dum
   20    continue
   10 continue
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
