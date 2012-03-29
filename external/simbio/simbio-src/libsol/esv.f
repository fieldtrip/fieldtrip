      subroutine esv(gstif,vecx,vecb,ipodia,indexj,lenge,
     &               npoin,dprod,resi,omega,tol,iter)
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
c     esv:
c     =======
c
c     subroutine zur Loesung der Gleichung A*x+b=0 mit dem Einzelschritt-
c     verfahren nach Gauss-Seidel mit Ueberrelaxation (SOR - Verfahren!)
c     omega kann nur 0<omega<2 gewaehlt werden!
c
c     nur fuer positiv definite matrizen geeignet!
c     ABER
c     kompakte speicherung wie fuer UNsymmetrische matrix notwendig!
c     => (siehe: schwarz, methode der finiten elemente,teubner)
c
c
c
c variable  typ  art  dimension   erklaerung
c-----------------------------------------------------------------------
c gstif      d    i   (lenge)     matrix
c ipodia     i    i   (npoin+1)   zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i   (lenge)     zeiger fuer spalte der elemente aus gstif
c vecb       d    i   (npoin)     konstantenvektor (rechte seite)
c vecx       d   i/o  (npoin)     loesungsvektor
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c tol        d    i               abbruchgrenze (euklidische norm)
c iter       i    o               zahl der iterationen
c omega      d    i               relaxationsfaktor
c
c#
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0,z1=1.d0,z2=2.d0)
      dimension gstif(lenge),ipodia(npoin+1),indexj(lenge),resi(npoin),
     &          vecx(npoin),dprod(npoin),vecb(npoin)
      iter = z0
 555  iter = iter + 1
      do 20 i=1,npoin
         dummy=vecx(i)
         wert=z0
         do 10 j=ipodia(i)+1,ipodia(i+1)-1
            wert=wert+gstif(j)*vecx(indexj(j))
 10      continue
         wert=wert+vecb(i)
         vecx(i)= - omega*wert/gstif(ipodia(i))+(z1-omega)*dummy
 20   continue
c---[A*x=ax!]
      call matpro(gstif,vecx,dprod,ipodia,indexj,lenge,npoin)
      rnorm=z0
      do 40 i=1,npoin
         resi(i)=dprod(i)+vecb(i)
         rnorm=rnorm+resi(i)*resi(i)
 40   continue
      rnorm=sqrt(rnorm)
      if (rnorm.gt.tol) goto 555
      return
      end
c=======================================================================
