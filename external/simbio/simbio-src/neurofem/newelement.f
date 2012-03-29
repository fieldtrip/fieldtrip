!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     newelement.f : 
!                              -------------------
!     begin                : Mai 2000
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

      subroutine newelement(knoflf,necflf,itygeo,necgeo,knogeo,xyzgeo,
     &                  pergeo,ipkt,isoll,itet, nelgeo, npogeo, nelflf,
     &                  maxgeo, eegpos, neeg, naniso, perele)
c 
      implicit double precision ( a-h,o-z )
c
      parameter( maxelk = 10000 )
c maxelk: max. Anzahl der Elektroden
      parameter( pi = 3.1415926535897931, z0 = 0.d0, z1 = 1.d0 , 
     &           z2 = 2.d0 )
      parameter( bignum=.888d88 )
      parameter( ndmgeo=3 )
      parameter( wuzeps=3.d-8 )
c
      logical isoll, itet
c
      dimension knoflf(maxgeo,nelflf),necflf(nelflf)
      dimension knogeo(maxgeo,nelgeo+neeg),necgeo(nelgeo+neeg)
      dimension xyzgeo(npogeo+neeg,ndmgeo),itygeo(nelgeo+neeg)
      dimension pergeo(nelgeo+neeg,naniso)
c
      dimension eegpos(ndmgeo,neeg)
c
      dimension ifix(maxelk)
      dimension p1(3),p2(3),p3(3),ps(3)
c
      if( isoll ) then
         write(*     ,900)' moving the nodes next to the electrodes'
      else
         write(*     ,900)' searching the nodes next to the electrodes'
      end if
c
c initialization of variables
      kn1 = 0
      ipkt = 0
      dissum = 0.0
      elxmax = -bignum
      elymax = -bignum
      elzmax = -bignum
      elxmin = bignum
      elymin = bignum
      elzmin = bignum
c
c---reading electrode coordinates
c
c   10 read(*,*,err=9999,end=99) elekx, eleky, elekz
   10 if( kn1 .ge. neeg ) goto 99
      kn1 = kn1 + 1
      elekx= eegpos(1,kn1)
      eleky= eegpos(2,kn1)
      elekz= eegpos(3,kn1)
      if( elekx .gt. elxmax ) elxmax = elekx
      if( eleky .gt. elymax ) elymax = eleky
      if( elekz .gt. elzmax ) elzmax = elekz
      if( elekx .lt. elxmin ) elxmin = elekx
      if( eleky .lt. elymin ) elymin = eleky
      if( elekz .lt. elzmin ) elzmin = elekz
c
c---search of next surface node
      dismin = bignum
      knomin = 0
      do ielflf = 1,nelflf
         do necken = 1, necflf( ielflf )
            kn2 = knoflf( necken, ielflf )
            if( kn2 .ne. 0 ) then
c..............Vergleich Fixknoten und Oberflaechenknoten
               xdiff = xyzgeo(kn2,1)-elekx
               ydiff = xyzgeo(kn2,2)-eleky
               zdiff = xyzgeo(kn2,3)-elekz
               dis   = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
               if( dis .lt. dismin ) then
                  knomin = kn2
                  dismin = dis
                  if( dismin .eq. z0 ) goto 30
               end if
            end if
         end do
      end do
   30 dismin = sqrt(dismin)
      dissum = dissum + dismin
      ipkt = ipkt + 1
      ifix(ipkt) = knomin
      write(*,800)' electrode ',kn1,
     &     ': Position  X Y Z: ', elekx, eleky, elekz
      write(*,810)'clos. node / pos / distance ',knomin,
     &     (xyzgeo(knomin,k),k=1,ndmgeo),dismin
c
      if (isoll) then
c........bewege den gefundenen Knoten auf den Fixpunkt
         xyzgeo(knomin,1) = elekx
         xyzgeo(knomin,2) = eleky
         xyzgeo(knomin,3) = elekz
      else if (.not.itet) then
c........Staebchenbildung
         if (dismin.lt.wuzeps) then
            elekx = elekx + wuzeps
            eleky = eleky + wuzeps
            elekz = elekz + wuzeps
            write(*,807)' new electrode pos ',elekx, eleky, elekz
         end if
         xyzgeo(npogeo+ipkt,1) = elekx
         xyzgeo(npogeo+ipkt,2) = eleky
         xyzgeo(npogeo+ipkt,3) = elekz
         knogeo(1,nelgeo+ipkt) = knomin
         knogeo(2,nelgeo+ipkt) = npogeo+ipkt
         necgeo(  nelgeo+ipkt) = 2
         itygeo(  nelgeo+ipkt) = 301
         do idiman = 1,naniso
            pergeo(  nelgeo+ipkt,idiman) = perele
         end do
         write(*,810)' new element / permeability ',
     &               nelgeo+ipkt,perele
      else
         xyzgeo(npogeo+ipkt,1) = elekx
         xyzgeo(npogeo+ipkt,2) = eleky
         xyzgeo(npogeo+ipkt,3) = elekz
c
c Suche zu einem Knoten ein Oberflaechenelement:
c 1. Projektion dieses Knotens auf die Ebene und 
c    Bestimmung des (Dreiecks)Flaechenwertes
c 2. Auswahl eines Obfl.Elementes und anschliessender Tertraederbildung 
c  
         wminmn = -bignum
         wminps = bignum
         minmin = 0
         minps =  0
         do ielflf = 1,nelflf
            knofnd = 0
            do necken = 1,necflf(ielflf)
               if( knoflf(necken,ielflf) .eq. knomin ) knofnd = 1
            end do
            if( knofnd .eq. 1 ) then
               p1(1)=xyzgeo(knoflf(1,ielflf),1)
               p2(1)=xyzgeo(knoflf(2,ielflf),1)
               p3(1)=xyzgeo(knoflf(3,ielflf),1)
               p1(2)=xyzgeo(knoflf(1,ielflf),2)
               p2(2)=xyzgeo(knoflf(2,ielflf),2)
               p3(2)=xyzgeo(knoflf(3,ielflf),2)
               p1(3)=xyzgeo(knoflf(1,ielflf),3)
               p2(3)=xyzgeo(knoflf(2,ielflf),3)
               p3(3)=xyzgeo(knoflf(3,ielflf),3)
               ps(1)=elekx
               ps(2)=eleky
               ps(3)=elekz
c
c berechne die Dreieckskoordinaten der Projektion und den Abstand
c von Ebene und betrachtetem Knotenpunkt
c
               call flawrt( p1, p2, p3, ps, wl1, wl2, wl3, abst )
               wmin = min(wl1,wl2,wl3)
               if( wmin .gt. 0. ) then
                  if( wmin .lt. wminps ) then  
                     wminps = wmin
                     minps = ielflf
                     abstps = abst
                  end if
               else if( wmin .gt. wminmn ) then
                  wminmn = wmin
                  minmin = ielflf
                  abstmn = abst
               end if
            end if
         end do
c
c     Elementbildung
         minabs = 0
         if( minps .gt. 0 .and. minps .ne. bignum ) then
            write(*     ,910)
     &           ' next (POS.) surface element / dist.',minps,abstps
            write(*     ,920)' surface nodes ',
     &           (knoflf(i,minps),i=1,ndmgeo)
            minabs = minps
         else
            write(*     ,910)
     &           ' next (NEG.) surface element / dist.',minmin,abstmn
            write(*     ,920)' surface nodes ',
     &           (knoflf(i,minmin),i=1,ndmgeo)
            minabs = minmin
         end if
c........Tetraederbildung
         if( minabs .eq. 0 ) then
            write( nterm, *) 'Error at making tetrahedrons.'
            ierr=1
c           call qfclal
            goto 9999
         end if
         knogeo(1,nelgeo+ipkt) = knoflf(1,minabs)
         knogeo(2,nelgeo+ipkt) = knoflf(3,minabs)
         knogeo(3,nelgeo+ipkt) = knoflf(2,minabs)
         knogeo(4,nelgeo+ipkt) = npogeo+ipkt
         necgeo(  nelgeo+ipkt) = 4
         itygeo(  nelgeo+ipkt) = 303
         do idiman = 1,naniso
             pergeo(  nelgeo+ipkt,idiman) = perele
         end do
      end if
c     
      goto 10
c
c...end of file reached, write results to disk and screen
c      
   99 write(*     ,930)' sum of all distances: ', dissum
c
      if (.not.isoll) knomin=npogeo+ipkt
  
c         write(*,900)' writing Dirchlet node file '
c         write(iunit9,'(7x,12i6)' ) knomin
c
c         write(*,900)' writing boundary condition file ',
c         write(iunit9,'(2x,i6,a,e12.5)') knomin,': ',z0
c
c         write(*,900)' writing electrode nodes file ',
c         write(iunit9,'(7x,12i6)')(i,i=npogeo+1,npogeo+ipkt-1)
c
c         write(*,900)'writing EIPP nodes file ', 
c         write(iunit3,'(7x,12i6)')(ifix(i),i=1,ipkt)

c      write(*     ,950)' reference electrode No./node/pos.',
c     &     kn1, knomin, elekx, eleky, elekz

      write(*     ,900)' range of electrodes positions'
      write(*     ,940) elxmin, ' < x < ',  elxmax
      write(*     ,940) elymin, ' < y < ',  elymax
      write(*     ,940) elzmin, ' < z < ',  elzmax
c     
      RETURN
c
  800 format(a10,i4,1x,a,3f10.3)
  805 format(16x,3f10.3)
  807 format(a,3f15.9)
  810 format(a,i6,3f10.3,f10.5)
  820 format(79('='))
  900 format(2a)
  910 format(a,i5,f7.3)
  920 format(a,3i6)
  930 format(a,f8.3)
  940 format(1x,f6.2,a,f6.2)
  950 format(a,2i7,2x,3(1x,f6.2))
c     
 9999 STOP
      END
c======================================================================
c----------------------------------------------------------------------
c
c  Subroutine FLAEWE
c             ======
c
c Diese Routine ermittelt einen Kennwert fuer ein Flaechen-Elemente,
c der bestimmt wie genau die Flaeche durch vorgegebenen Support-Points
c verlaeuft
c
c EINGABE:
c nflae                 max. Anzahl Flaechen
c mknpfl                max. Anzahl Knoten pro Flaeche
c knofl(mknpel,nelem)   Knoten-Flaechen-Zuordnung
c p1,p2,p3 = Koordinaten des Flaechenelementes
c ps       = Knotenkoordinaten des Messpunktes
c
c AUSGABE:
c
c HILFSFELDER:
c
c LOKAL:
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine flawrt( p1, p2, p3, ps, wl1, wl2, wl3, abst )
c
      implicit double precision ( a-h,o-z )
c
      parameter( pi = 3.1415926535897931, z0 = 0.d0, z1 = 1.d0 , 
     &           z2 = 2.d0 )

      dimension
     & p1(3),  p2(3),  p3(3),  pp(3), ps(3),
     & v12(3), v23(3), v31(3), vn(3), vp(3), vnp(3)

c Flaechen-Vektoren & Normale ermitteln

      do iri=1,3
         v12(iri)=p2(iri)-p1(iri)
         v23(iri)=p3(iri)-p2(iri)
         v31(iri)=p1(iri)-p3(iri)
      end do
      call crosp(v12,v23,vn)
      d=-scalp(vn,p1)

      bn2=scalp(vn,vn)
      bn=sqrt(bn2)
      edbn2=1./bn2
      a=bn*0.5
      eda=1./a

c      amax=spabst
c      edamax=1./amax

c Flaechenwert ermitteln

      flw=0.
      nnp=0
      t=edbn2*( scalp(ps,vn) + d )
      abst=abs(t*bn)
         
      pp(1)=ps(1)-t*vn(1)
      pp(2)=ps(2)-t*vn(2)
      pp(3)=ps(3)-t*vn(3)

      vp(1)=pp(1)-p3(1)
      vp(2)=pp(2)-p3(2)
      vp(3)=pp(3)-p3(3)
      call crosp(v23,vp,vnp)
      a1=0.5*sqrt(scalp(vnp,vnp))*sign(z1,scalp(vn,vnp))

      vp(1)=pp(1)-p1(1)
      vp(2)=pp(2)-p1(2)
      vp(3)=pp(3)-p1(3)
      call crosp(v31,vp,vnp)
      a2=0.5*sqrt(scalp(vnp,vnp))*sign(z1,scalp(vn,vnp))
         
      vp(1)=pp(1)-p2(1)
      vp(2)=pp(2)-p2(2)
      vp(3)=pp(3)-p2(3)
      call crosp(v12,vp,vnp)
      a3=0.5*sqrt(scalp(vnp,vnp))*sign(z1,scalp(vn,vnp))
         
      wl1=a1*eda
      wl2=a2*eda
      wl3=a3*eda
      wlmin=min(wl1,wl2,wl3)
      wlmax=max(wl1,wl2,wl3)
      wsum=wl1+wl2+wl3
      if(wsum.lt.0.999.or.wsum.gt.1.001)then
         write(*,7000) wsum
 7000    format(1x,'E-flaewe: Summe 3ecks-Koord. =',g12.5)
      end if
c      if(wlmin.gt.wlkrit)then
c            wlamax=max(wlamax,wlmax)
c            wlamin=min(wlamin,wlmin)
c            nnp=nnp+1
cc            df=wlmin*(amax-abst)
c            flw=flw+df
c         end if
c$$$      end if
c$$$      end do
c      if(nnp.ne.0) flw=flw/real(nnp)
c      flw=3.*flw*edamax

      RETURN
      END
c----------------------------------------------------------------------
c Kreuzprodunkt zweier Vektoren
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine crosp(a,b,v)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),v(3)
      v(1)=a(2)*b(3)-a(3)*b(2)
      v(2)=a(3)*b(1)-a(1)*b(3)
      v(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

c----------------------------------------------------------------------
c Skalarprodukt zweier Vektoren
c---->---1---------2---------3---------4---------5---------6---------7--<
      function scalp(a,b)
      implicit double precision(a-h,o-z)
      dimension a(3),b(3)
      scalp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
c======================================================================
