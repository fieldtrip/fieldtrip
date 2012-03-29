      subroutine unscal(gstif,vecx,vecb,ipodia,indexj,dkond,lenge,
     &                  npoin)
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
c#
c     unscal:
c     =======
c
c     skalierung der steifigkeitsmatrix rueckgaengig machen,
c     der rechten seite und der loesung
c     fuer symmetrische matrizen!!!
c
c
c variable  typ  art  dimension   erklaerung
c-----------------------------------------------------------------------
c dkond      d    i   (npoin)     vektor mit skalierungsfaktoren
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
      implicit double precision (a-h,o-z)
      dimension ipodia(npoin),indexj(lenge)
      dimension gstif(lenge),dkond(npoin),vecb(npoin),vecx(npoin)
c
      gstif(1) = gstif(1) / ( dkond(1)*dkond(indexj(1)) )
      vecb(1)=vecb(1)/dkond(1)
      vecx(1)=vecx(1)*dkond(1)
      do 20 ij=2,npoin
         vecb(ij)=vecb(ij)/dkond(ij)
         vecx(ij)=vecx(ij)*dkond(ij)
         do 30 j=ipodia(ij-1)+1,ipodia(ij)
            gstif(j)=gstif(j)/(dkond(ij)*dkond(indexj(j)))
 30      continue
 20   continue
c
      return
      end
c======================================================================
