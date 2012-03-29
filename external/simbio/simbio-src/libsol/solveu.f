c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine solveu(gstilu,vecx,vecb,ipodia,indexj,lenge,npoin)
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
c     solveu:
c     =======
c     
c     Loest die Gleichung U*x=b ( hier: gstilu*vecx=vecb ) !
c     Wobei U eine rechte obere Dreiecksmatrix ist
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
c 
c=====================================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension gstilu(lenge),vecx(npoin),vecb(npoin)
      dimension ipodia(npoin+1),indexj(lenge)
c
      do 10 i=1,npoin
         vecx(i)=vecb(i)
 10   continue
c
      do 20 i=npoin-1,1,-1
         do 30 l=ipodia(i+1)-1,ipodia(i)+1,-1
            if (indexj(l) .gt. i) then
               vecx(i)=vecx(i)-gstilu(l)*vecx(indexj(l))
            else
               goto 20
            end if
 30      continue
 20   continue
      return
      end


 
 

