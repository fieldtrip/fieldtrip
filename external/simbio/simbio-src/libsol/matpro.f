      subroutine matpro(gstif,vector,dprod,ipodia,indexj,lenge,npoin)
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
c
c     matpro:
c     =======
c
c     berechnet ein matrixprodukt gstif*vector=dprod
c     in kompakter speicherung (!!! fuer symmetrische matrizen !!!)
c
c
c
c
c variable  typ  art  dimension   erklaerung
c-----------------------------------------------------------------------
c dprod      d    o   (npoin)     matrixprodukt
c gstif      d    i   (lenge)     matrix
c ipodia     i    i   (npoin+1)   zeiger fuer position der diagonalelemente in gstif
c indexj     i    i   (lenge)     zeiger fuer spalte der elemente aus gstif
c vector     d    i   (npoin)     vektor
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c
c
c
c#
c=======================================================================
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      dimension gstif(lenge),vector(npoin),dprod(npoin)
      dimension ipodia(npoin),indexj(lenge)
c
      dprod(1)=gstif(1)*vector(1)
      ianf=2
      do 10 i=2,npoin
         vi=vector(i)
         id=ipodia(i)
         di=gstif(id)*vi
         do 20 k=ianf,id-1
            j=indexj(k)
            gj = gstif(k)
            di = di + gj * vector(j)
            dprod(j) = dprod(j) + gj * vi
 20      continue
         dprod(i)=di
         ianf=id+1
 10   continue
c
      return
      end
c====================================================================
 
