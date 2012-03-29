      subroutine ptfqmr(gstif,vecx,vecb,ipodia,indexj,lenge,npoin,ax,
     &                 astriy,vecd,rdach,vecv,vecw,vecy,tol,maxit,ierr,
     &                 gstilu,astrix,vecxst,vecbst,veco,vecp)
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
c     prefreund:
c     ==========
c
c
c     Subroutine zur Loesung einer Gleichung A*x + b = 0 .
c
c     (A ist eine unsymmetrische Matrix). Der verwendete Algorithmus 
c     ist der vorkonditionierte
c     Transpose Free Quasi Minimal Residual Algorithm von 
c     Roland W.Freund
c     kompakte speicherung fuer UNsymmetrische matrizen verwenden
c               => (siehe: schwarz, methode der finiten elemente)
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
c maxit      i    i               maximale anzahl der iterationen
c ierr       i    o               fehlerparameter 
c                                 0: alles klar
c                                 2: vorgegebene toleranz kann
c                                    nicht erreicht
c                                    werden, automatisch bis zur groesst-
c                                    moeglichen genauigkeit berechnet
c gstilu     d    i   (lenge)     vorkonditionierungsmatrix
c 
c
c---pruefe ob erwartete genauigkeit erreicht!!!!!!!!!!!!!!!!!
c=======================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      real ransol
      parameter (z0=0.d0,z1=1.d0,z2=2.d0,dd5=1.d-2)
      dimension vecy(npoin),vecw(npoin),vecv(npoin),gstif(lenge),
     &          rdach(npoin),vecd(npoin),vecx(npoin),vecb(npoin),
     &          astriy(npoin),astrix(npoin),vecxst(npoin),
     &          gstilu(lenge),veco(npoin),vecp(npoin),ax(npoin)
      dimension vecbst(npoin),ipodia(npoin+1),indexj(lenge)
c
c=======================================================================
c                 INITIALISIERUNG
c=======================================================================
c---[x'=U*x]
      call axunsu(gstilu,vecx,vecxst,ipodia,indexj,lenge,npoin)
c
      tolin=tol
      rnorm=1.d10
      iter=0
c=======================================================================
c---START, 0-TER SCHRITT
c=======================================================================
c---[A*x=ax]
      call axunsy(gstif,vecx,ax,ipodia,indexj,lenge,npoin)
c---[A'x'=astrix=L^(-1)*A*x]
      call solvel(gstilu,astrix,ax,ipodia,indexj,lenge,npoin)
c---[b'=L^(-1)*b]
      call solvel(gstilu,vecbst,vecb,ipodia,indexj,lenge,npoin)
      do 10 i=1,npoin
         vecy(i)=vecbst(i)+astrix(i)
         vecw(i)=vecy(i)
 10   continue
c=======================================================================
c---[ A'* y = astriy ]
c=======================================================================
c---[o=U^(-1)*y]
      call solveu(gstilu,veco,vecy,ipodia,indexj,lenge,npoin)
c---[p=A*o=A*U^(-1)*y]
      call axunsy(gstif,veco,vecp,ipodia,indexj,lenge,npoin)
c---[astriy=L^(-1)*p=L^(-1)*A*U^(-1)*y]
      call solvel(gstilu,astriy,vecp,ipodia,indexj,lenge,npoin)
c-----------------------------------------------------------------------
      tau=z0
      do 20 i=1,npoin
         vecv(i)=astriy(i)            
         tau=tau+vecy(i)*vecy(i)
         vecd(i)=z0
 20   continue            
      tau=sqrt(tau)
      teta=z0
      eta=z0
      r0=tau
c-----------------------------------------------------------------------
c---bestimmen von rdach mit ransol.f
c---Komponenten von rdach mit dem gleichen Vorzeichen wie vecv(i)
c---versehen!    
c-----------------------------------------------------------------------
      iseed=-3
      vbetra=z0
      rdbetr=z0
 70   do 30 i=1,npoin
         vbetra=vbetra+(vecv(i)*vecv(i))
         rdach(i)=dble(ransol(iseed))
         rdach(i)=sign(rdach(i),vecv(i))
         rdbetr=rdbetr+(rdach(i)*rdach(i))
 30   continue
      vbetra=sqrt(vbetra)
      rdbetr=sqrt(rdbetr)
c-----------------------------------------------------------------------
c--- ueberpruefen, ob rdach senkrecht zu vecv
c-----------------------------------------------------------------------
      contro=z0
      do 50 i=1,npoin
         contro=contro+(rdach(i)*vecv(i))
 50   continue
      contro=contro/(vbetra*rdbetr)
      conbet=dabs(contro)
      if (conbet.gt.dd5) goto 60 
      goto 70
c-----------------------------------------------------------------------
 60   delta=z0
      do 80 i=1,npoin
         delta=delta+rdach(i)*vecy(i)
 80   continue
c-----------------------------------------------------------------------
c                 FOR n=1...nmax
c-----------------------------------------------------------------------
 5    iter=iter+1
      sigma=z0
      do 90 i=1,npoin
         sigma=sigma+rdach(i)*vecv(i)
 90   continue
c      
      alfa=delta/sigma
c=======================================================================
      wbetra=z0
      dd1=z1/alfa
      dd2=teta*teta*eta*dd1
      do 100 i=1,npoin
         vecd(i)=vecy(i)+dd2*vecd(i)
         vecw(i)=vecw(i)-alfa*astriy(i)            
         wbetra=wbetra+vecw(i)*vecw(i)
100   continue
c      
      wbetra=sqrt(wbetra)
      teta=wbetra/tau
      c=z1/(sqrt(z1+teta*teta))
      tau=tau*teta*c
      eta=c*c*alfa
c      
      do 110 i=1,npoin
         vecxst(i)=vecxst(i)-eta*vecd(i)                        
110   continue
c-----------------------------------------------------------------------
c---abbruchkriterium
c-----------------------------------------------------------------------
      if (tau*(sqrt(dble(iter)+z1)).lt.tolin*r0) goto 6
c=======================================================================
c---zweiter teil des tfqmr von freund    
c=======================================================================
      do 130 i=1,npoin
         vecy(i)=vecy(i)-alfa*vecv(i)
 130     continue
c=======================================================================
c---[ A'* y = astriy ]
c=======================================================================
c---[o=U^(-1)*y]
      call solveu(gstilu,veco,vecy,ipodia,indexj,lenge,npoin)
c---[p=A*o=A*U^(-1)*y]
      call axunsy(gstif,veco,vecp,ipodia,indexj,lenge,npoin)
c---[astriy=L^(-1)*p=L^(-1)*A*U^(-1)*y]
      call solvel(gstilu,astriy,vecp,ipodia,indexj,lenge,npoin)
c-----------------------------------------------------------------------
      wbetra=z0
      dd3=z1/alfa
      dd4=teta*teta*eta*dd3
      do 140 i=1,npoin
         vecd(i)=vecy(i)+dd4*vecd(i)
         vecw(i)=vecw(i)-alfa*astriy(i)            
         wbetra=wbetra+vecw(i)*vecw(i)
140   continue
c      
      wbetra=sqrt(wbetra)
      teta=wbetra/tau
      c=z1/(sqrt(z1+teta*teta))
      tau=tau*teta*c
      eta=c*c*alfa
c      
      do 150 i=1,npoin
         vecxst(i)=vecxst(i)-eta*vecd(i)                        
150   continue
c-----------------------------------------------------------------------
c---abbruchkriterium
c-----------------------------------------------------------------------
      if (tau*(sqrt(dble(iter)+z1)).lt.tolin*r0) goto 6
c=======================================================================
c---letzter teil des tfqmr von freund
c=======================================================================
      beta=delta
c      
      delta=z0
      do 160 i=1,npoin
         delta=delta+rdach(i)*vecw(i)
160   continue
      delbet=dabs(delta)
      if (beta.eq.z0) ierr=3
      beta=delta/beta
c      
      do 170 i=1,npoin
         vecv(i)=beta*(astriy(i)+beta*vecv(i))                  
         vecy(i)=vecw(i)+beta*vecy(i)
 170  continue
c=======================================================================
c---[ A'* y = astriy ]
c=======================================================================
c---[o=U^(-1)*y]
      call solveu(gstilu,veco,vecy,ipodia,indexj,lenge,npoin)
c---[p=A*o=A*U^(-1)*y]
      call axunsy(gstif,veco,vecp,ipodia,indexj,lenge,npoin)
c---[astriy=L^(-1)*p=L^(-1)*A*U^(-1)*y]
      call solvel(gstilu,astriy,vecp,ipodia,indexj,lenge,npoin)
c-----------------------------------------------------------------------
      do 180 i=1,npoin 
         vecv(i)=vecv(i)+astriy(i)
 180  continue
c
c-----------------------------------------------------------------------
c---abbruchkriterium ( notbremse bei divergenter iteration )
c-----------------------------------------------------------------------
      if (iter.gt.maxit) goto 6
c
      goto 5
c---Ende des freund Algorithmus!
c=======================================================================
c---[x=U^(-1)*x']
 6    call solveu(gstilu,vecx,vecxst,ipodia,indexj,lenge,npoin)
      return
      end







