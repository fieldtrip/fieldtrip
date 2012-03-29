c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine solvlu(gstilu,vecx,vecb,vecy,ipodia,indexj,
     &                                                 lenge,npoin)
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
c     solvelu:
c     ========
c     
c     Loest die Gleichung L*U*x=b ( hier: gstilu*vecx=vecb ) !
c     Wobei L eine linke untere Dreiecksmatrix und
c           U eine rechte obere Dreiecksmatrix ist.
c     L und U sind abgelegt in gstilu, wobei die Diagonalelemente
c     von u den Wert 1 haben und in gstilu nicht abgelegt sind!
c     
c
c
c variable  typ  art  dimension   erklaerung
c---------------------------------------------------------------------
c gstilu     d    i    (lenge)    matrix
c vecx       d    o    (npoin)    loesungsvector
c vecb       d    i    (npoin)    konstantenvector (rechte Seite)
c ipodia     i    i    (npoin+1)  zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i    (lenge)    zeiger fuer spalte der elemente aus 
c                                 gstif
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c vecy       d    o    (npoin)    hilfsvector
c
c=====================================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension gstilu(lenge),vecx(npoin),vecb(npoin),vecy(npoin)
      dimension ipodia(npoin+1),indexj(lenge)
c
c L*y=b
c
      call solvel(gstilu,vecy,vecb,ipodia,indexj,lenge,npoin)
c
c U*x=y
c
      call solveu(gstilu,vecx,vecy,ipodia,indexj,lenge,npoin)
      return
      end

         


    

