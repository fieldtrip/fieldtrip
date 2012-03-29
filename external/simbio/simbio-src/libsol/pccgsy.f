c !  $2 04.05.2007  Seok Lew: change to energy norm as PCG stop condition

      subroutine pccgsy(gstif,vecx,vecb,ipodia,indexj,lenge,npoin,
     &                  dprod,vecr,vecp,
     &                  resi,tol,maxit,ierr,gstifc,vecrho)
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
c     prekograd:
c     ==========
c
c     vorkonditionierter gleichungsloeser fuer A*x+b=0;
c     nur positiv definite matrizen
c     nach loesungsvorschrift von hestenes und stiefel (matrizen-
c     numerik von schwarz-rutishauser-stiefel)
c     kompakte speicherung fuer symmetrische matrizen verwenden
c             => (siehe: schwarz, methode der finiten elemente)
c
c
c
c variable  typ  art  dimension   erklaerung
c-----------------------------------------------------------------------
c gstif      d    i   (lenge)     matrix
c ipodia     i    i   (npoin+1)   zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i   (lenge)     zeiger fuer spalte der elemente
c                                 aus gstif
c vecb       d    i   (npoin)     konstantenvektor (rechte seite)
c vecx       d   i/o  (npoin)     loesungsvektor
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c tol        d    i               abbruchgrenze (euklidische norm)
c maxit      i    o               maximale anzahl der iterationen
c resi       d    o               euklidische norm des restvektors
c ierr       i    o               fehlerparameter
c                                 0: alles klar
c                                 1: anzahl iterationen groesser als npoin
c gstifc     d    i   (lenge)     vorkonditionierungsmatrix
c
c
c
c=======================================================================
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      dimension ipodia(npoin),indexj(lenge)
      dimension gstif(lenge),vecb(npoin),vecx(npoin),vecr(npoin),
     &          dprod(npoin),vecp(npoin),vecrho(npoin),gstifc(lenge)
c
c---initialisierung
      drn=z0
      ierr=0
      iter=0
      dnormb=z0
      do 5 i=1,npoin
         dnormb=dnormb+vecb(i)*vecb(i)
 5    continue
      dnormb=sqrt(dnormb)

c
c---0-ter schritt
      resi = z0
      call matpro(gstif,vecx,dprod,ipodia,indexj,lenge,npoin)
      do 10 i=1,npoin
         vecr(i) =  dprod(i)+vecb(i)
         resi  = resi  + vecr(i) *    vecr(i)
 10   continue
c---- Seok: initial \gamma for calculating energy norm
      resi0=resi
c-- for L2 norm    resi = sqrt(resi)
c-- for relative L2 norm    resi=resi/dnormb
c---- Seok: tol must be squared for PCG stop condition when energy norm is used 
      if (resi.lt.tol*tol) then
         return
      end if
      call solchs(gstifc,vecrho,vecr,ipodia,indexj,lenge,npoin)
      do 15 i=1,npoin
         vecp(i) =  -vecrho(i)
 15   continue
 
c
c---allgemeines schema ab dem 1-ten schritt
 20   iter=iter+1
c
      dra = drn
      drn = z0
      do 30 i=1,npoin
         drn = drn + vecr(i) * vecrho(i)
 30   continue
      if(iter.ge.2)then
         de = drn / dra
         do 40 i=1,npoin
            vecp(i) = vecp(i) * de - vecrho(i)
 40      continue
      else
         do 50 i=1,npoin
            vecp(i) = - vecrho(i)
 50      continue
      end if
c
      call matpro(gstif,vecp,dprod,ipodia,indexj,lenge,npoin)
c
      dnen = z0
      do 80 i=1,npoin
         dnen = dnen + vecp(i) * dprod(i)
 80   continue
      dq = drn / dnen
c
      resi = z0
      do 90 i=1,npoin
         vecx(i) = vecx(i) + dq    *    vecp(i)
         vecr(i) = vecr(i) + dq    * dprod(i)
         resi  = resi  + vecr(i) *    vecr(i)
   90 continue
c---- Seok: current \gamma divided by the initial \gamma for calculating energy norm
      resi=resi/resi0
c-- for L2 nrom     resi = sqrt(resi)
c-- for relative L2     resi = resi/dnormb
c---- Seok: tol must be squared for PCG stop condition when energy norm is used 
      if (resi.lt.tol*tol) then
         ierr=-iter
         return
      else if (iter.ge.maxit) then
         ierr=1
         return
      end if
         call solchs(gstifc,vecrho,vecr,ipodia,indexj,lenge,npoin)
         goto 20
c
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
