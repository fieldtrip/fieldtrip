      subroutine cgstab (gstif,vecx,vecb,ipodia,indexj,lenge,npoin,ax,
     &                   vecp,vecr,rdach,vecs,vect,vecv,tol,iter)
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
c     cgstab:
c     =======
c
c     Subroutine zur Loesung der Gleichung A*x + b = 0, wobei A eine un-
c     symmetrische Matrix ist. Der verwendete Algorithmus ist der
c     Bi-CGSTAB (bi-cgstability) von H.A. Van der Vorst (siam 3.1992)
c     kompakte speicherung fuer UNsymmetrische matrizen verwenden
c               => (siehe: schwarz, methode der finiten elemente)
c
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
c
c
c#
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0,z1=1.d0,z2=2.d0)
      dimension vecs(npoin),vect(npoin),vecv(npoin),gstif(lenge),
     &          rdach(npoin),vecp(npoin),vecx(npoin),vecb(npoin),
     &          vecr(npoin),ax(npoin)
      dimension ipodia(npoin+1),indexj(lenge)
c=======================================================================
c                 INITIALISIERUNG
      iter=0
c=======================================================================
c---START, 0-TER SCHRITT
c
c---[A*x=ax!]
      call axunsy(gstif,vecx,ax,ipodia,indexj,lenge,npoin)
      do 10 i=1,npoin
         vecr(i)=vecb(i)+ax(i)
c---bestimmen von rdach (rdach=vecr ist deutlich besser als zufallsvektor
c---rdach!)
         rdach(i)=vecr(i)
 10   continue
c
 19   delta0=z1
      alfa=z1
      omega=z1
c
      do 20 i=1,npoin
         vecv(i)=z0
         vecp(i)=z0
 20   continue
c
c=======================================================================
c                 FOR n=1........nmax
c=======================================================================
c
 5    iter=iter+1
      delta1=z0
      do 90 i=1,npoin
         delta1=delta1+rdach(i)*vecr(i)
 90   continue
c
c---Nullabfrage eruebrigt sich fuer omega und delta!
c
      beta=(delta1/delta0)*(alfa/omega)
      do 40 i=1,npoin
         vecp(i)=vecr(i)+beta*(vecp(i)-omega*vecv(i))
 40   continue
c
      call axunsy (gstif,vecp,ax,ipodia,indexj,lenge,npoin)
      alfa=z0
      do 50 i=1,npoin
         vecv(i)=ax(i)
      alfa=alfa+rdach(i)*vecv(i)
 50   continue
c
      alfa=delta1/alfa
      delta0=delta1
      do 60 i=1,npoin
         vecs(i)=vecr(i)-alfa*vecv(i)
 60   continue
c
      call axunsy (gstif,vecs,ax,ipodia,indexj,lenge,npoin)
      do 70 i=1,npoin
         vect(i)=ax(i)
 70   continue
c
      omega1=z0
      omega2=z0
      do 72 i=1,npoin
         omega1=omega1+vect(i)*vecs(i)
         omega2=omega2+vect(i)*vect(i)
 72   continue
c
      omega=omega1/omega2
      do 80 i=1,npoin
         vecx(i)=vecx(i)-alfa*vecp(i)-omega*vecs(i)
         vecr(i)=vecs(i)-omega*vect(i)
 80   continue
c
c---residuum bestimmen
c---[A*x=ax!]
      rnorm=z0
      do 120 i=1,npoin
         rnorm=rnorm+vecr(i)*vecr(i)
 120  continue
      rnorm=sqrt(rnorm)
c
c---abbruchkriterium
c
      if (rnorm.lt.tol) then
         return
      end if
c-----------------------------------------------------------------------
      goto 5
      end
