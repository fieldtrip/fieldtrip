      subroutine prcgst (gstif,vecx,vecb,ipodia,indexj,lenge,npoin,
     &                      vecp,vecr,rdach,vecs,vect,vecv,tol,maxit,
     &                      gstilu,vecy,vecz,vecth,vecsh)
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
c     precgstab:
c     ==========
c
c     Subroutine zur Loesung der Gleichung A*x+b=0,
c     die den vorkonditionierten Algorithmus
c     Bi-CGSTAB (bi-cgstability) von H.A. Van der Vorst (siam 13.1992)
c     mit kompakter Speicherung fuer UNsymmetrische matrizen verwendet
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
c indexj     i    i   (lenge)     zeiger fuer spalte der
c                                 elemente aus gstif
c vecb       d    i   (npoin)     konstantenvektor (rechte seite)
c vecx       d   i/o  (npoin)     loesungsvektor
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c tol        d    i               abbruchgrenze (euklidische norm)
c maxit      i    i               maximale anzahl der iterationen
c gstilu     d    i   (lenge)     vorkonditionierungsmatrix
c
c
c vecth      d        (npoin)     hilfsvector fuer viele anwendungen
c=======================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (z0=0.d0,z1=1.d0,z2=2.d0)
      dimension vecs(npoin),vect(npoin),vecv(npoin),gstif(lenge),
     &          rdach(npoin),vecp(npoin),vecx(npoin),vecb(npoin),
     &          vecr(npoin),gstilu(lenge),vecy(npoin),
     &          vecz(npoin),vecth(npoin),vecsh(npoin)
      dimension ipodia(npoin+1),indexj(lenge)
c=======================================================================
c                 INITIALISIERUNG
      iter=0
c=======================================================================
c---START, 0-TER SCHRITT
c
c---[A*x=vecth!]
      call axunsy(gstif,vecx,vecth,ipodia,indexj,lenge,npoin)
      do 10 i=1,npoin
         vecr(i)=vecb(i)+vecth(i)
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
      call solvlu (gstilu,vecy,vecp,vecth,ipodia,indexj,lenge,npoin)
      call axunsy (gstif,vecy,vecv,ipodia,indexj,lenge,npoin)
      alfa=z0
      do 50 i=1,npoin
         alfa=alfa+rdach(i)*vecv(i)
 50   continue
c
      alfa=delta1/alfa
      delta0=delta1
      do 60 i=1,npoin
         vecs(i)=vecr(i)-alfa*vecv(i)
 60   continue
c
      call solvlu (gstilu,vecz,vecs,vecth,ipodia,indexj,lenge,npoin)
      call axunsy (gstif,vecz,vect,ipodia,indexj,lenge,npoin)
c
      call solvel(gstilu,vecth,vect,ipodia,indexj,lenge,npoin)
      call solvel(gstilu,vecsh,vecs,ipodia,indexj,lenge,npoin)
c
      omega1=z0
      omega2=z0
      do 72 i=1,npoin
         omega1=omega1+vecth(i)*vecsh(i)
         omega2=omega2+vecth(i)*vecth(i)
 72   continue
c
      omega=omega1/omega2
      do 80 i=1,npoin
         vecx(i)=vecx(i)-alfa*vecy(i)-omega*vecz(i)
         vecr(i)=vecs(i)-omega*vect(i)
 80   continue
c
c---residuum bestimmen
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
      else if (iter .gt. maxit) then
         ierr=4
         return
      end if
c-----------------------------------------------------------------------
      goto 5
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
