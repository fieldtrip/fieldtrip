c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine partch(gstif,ipodia,indexj,lenge,npoin,gstifc,ierr)
!
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
c
c
c     partch:
c     =======
c
c     partielle cholesky-zerlegung der skalierten matrix gstif
c     multiplikation der aussendiagonalelemente von gstif mit
c     1/(1+alpha) liefert die matrix gstifc, falls die zerlegung
c     moeglich ist; ( falls nicht, ierr=4, tip: versuche die
c     zerlegung mit rilu.f )
c
c variable  typ  art  dimension   erklaerung
c---------------------------------------------------------------------
c gstif      d    i    (lenge)    matrix
c ipodia     i    i    (npoin+1)  zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i    (lenge)    zeiger fuer spalte der elemente aus
c                                 gstif
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c gstifc     d    i    (lenge)    matrix
c ierr       i    o               fehlerparameter
c                                 0: alles klar
c                                 4: cholesky zerlegung nicht moeglich
c=====================================================================
      implicit double precision (a-h,o-z)
      dimension gstif(lenge),gstifc(lenge)
      dimension ipodia(npoin+1),indexj(lenge)
c
      ierr=4
      alpha=0.d1
c
c multiplizieren der nichtdiagonalelemente mit rfak
c
 5    rfak=1.d0/(1.d0+alpha)
      gstifc(1)=gstif(1)
      do 20 i=2,npoin
         jlow=ipodia(i-1)+1
         jup=ipodia(i)-1
         gstifc(ipodia(i))=gstif(ipodia(i))
         if (jlow .gt. jup) goto 20
         do 10 j=jlow,jup
            gstifc(j)=gstif(j)*rfak
 10      continue
 20   continue
c
c eigentliche cholesky-zerlegung
c
      do 70 i=2,npoin
         jlow=ipodia(i-1)+1
         jup=ipodia(i)-1
         if (jlow .gt. jup) goto 60
         do 50 j=jlow,jup
            gstifc(j)=gstifc(j)/gstifc(ipodia(indexj(j)))
            klow=j+1
            kup=ipodia(i)
            do 40 k=klow,kup
               l=indexj(k)
               ilow=ipodia(l-1)+1
               iup=ipodia(l)
               do 30 igstifc=ilow,iup
                  if (indexj(igstifc) .gt. indexj(j)) goto 40
                  if (indexj(igstifc) .lt. indexj(j)) goto 30
                  gstifc(k)=gstifc(k)-gstifc(j)*gstifc(igstifc)
                  goto 40
 30            continue
 40         continue
 50      continue
 60      id=ipodia(i)
         if (gstifc(id) .lt. 1.d-6)  then
            if (alpha .eq. 0.d1) alpha=0.5d-6
c$$$ adrian, 6.11.'95            alpha=alpha+alpha
c$$$ adrian, 6.11.'95            if (alpha .gt. 1.d-3) return
            alpha=alpha*4.d0
            if (alpha .gt. 1.d0) return
            goto 5
         end if
         gstifc(id)=dsqrt(gstifc(id))
 70   continue
      ierr=0
      return
      end
 
 
 
 
 
 
 
 
 
 
