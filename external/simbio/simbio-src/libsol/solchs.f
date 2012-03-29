c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine solchs(gstih,vecx,vecb,ipodia,indexj,lenge,npoin)
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
c     solchs:
c     ========
c     
c     Loest die Gleichung H*HT*x=b ( hier: gstih*vecx=vecb ) !
c     Wobei H eine linke untere Dreiecksmatrix ist.
c     H ist abgelegt in gstih !
c     
c
c
c variable  typ  art  dimension   erklaerung
c---------------------------------------------------------------------
c gstih      d    i    (lenge)    matrix
c vecx       d    o    (npoin)    loesungsvector
c vecb       d    i    (npoin)    konstantenvector (rechte Seite)
c ipodia     i    i    (npoin)    zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i    (lenge)    zeiger fuer spalte der elemente aus 
c                                 gstif
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c
c
c=====================================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension gstih(lenge),vecx(npoin),vecb(npoin)
      dimension ipodia(npoin),indexj(lenge)
c
      vecx(1)=vecb(1)
      do 30 k=2,npoin
         s=0.d0
         jlow=ipodia(k-1)+1
         jup=ipodia(k)-1
         if (jlow .gt. jup) goto 20
         do 10 j=jlow,jup
            s=s+gstih(j)*vecx(indexj(j))
 10      continue
 20      vecx(k)=(vecb(k) - s) / gstih(ipodia(k))
 30   continue
      do 50 kh=2,npoin
         k=npoin+2-kh
         vecx(k)=vecx(k)/gstih(ipodia(k))
         jlow=ipodia(k-1)+1
         jup=ipodia(k)-1
         if (jlow .gt. jup) goto 50
         do 40 j=jlow,jup
            vecx(indexj(j))=vecx(indexj(j))-gstih(j)*vecx(k)
 40      continue
 50   continue
      return
      end

         


    

