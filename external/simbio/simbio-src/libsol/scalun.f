      subroutine scalun(gstif,vecx,vecb,ipodia,indexj,dkond,lenge,npoin)
      implicit double precision (a-h,o-z)
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
c#
c     scalun:
c     =======
c
c     skalierung der steifigkeitsmatrix,
c     der rechten seite und der loesung
c     fuer UNsymmetrische matrizen!!!
c
c variable  typ  art  dimension   erklaerung
c-----------------------------------------------------------------------
c dkond      d    o   (npoin)     vektor mit skalierungsfaktoren
c gstif      d   i/o  (lenge)     matrix
c ipodia     i    i   (npoin+1)   zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i   (lenge)     zeiger fuer spalte der elemente aus gstif
c vecb       d   i/o  (npoin)     konstantenvektor (rechte seite)
c vecx       d   i/o  (npoin)     loesungsvektor
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c
c
c
c#
c
c=======================================================================
c
      parameter (z1=1.d0)
      dimension ipodia(npoin+1),indexj(lenge)
      dimension vecb(npoin),gstif(lenge),dkond(npoin),vecx(npoin)
c
      do 10 i=1,npoin
        dkond(i)=z1/sqrt(gstif(ipodia(i)))
        vecb(i)=vecb(i)*dkond(i)
        vecx(i)=vecx(i)/dkond(i)
   10 continue
c
      do 20 i=2,npoin+1
         im1=i-1
         do 30 j=ipodia(im1),ipodia(i)-1
            gstif(j)=gstif(j)*dkond(im1)*dkond(indexj(j))
 30      continue
 20   continue
c
      return
      end
c=======================================================================
