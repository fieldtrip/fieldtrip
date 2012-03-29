!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem7.f :          Subroutines related to the treatment of the 
!                            surface 
!                            Regularized inverse solutions
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

      subroutine detnrm(knoflf,ityflf,necflf,ieldum,                    &
     &                  xyzflf,vcnflf)
      include 'neurofem.inc'
!
      dimension knoflf(mxkn2d,nelflf),ityflf(nelflf),necflf(nelflf),    &
     &          ieldum(nelflf)
      dimension xyzflf(npoflf,3),vcnflf(npoflf,3)
!
      ndmgeo=3
!---zero relevant vectors
      call dset(ndmgeo*npoflf,z0,vcnflf,1)
!
!---test surface structure for integrity [(c) Richard Schoenen, Oct. 94]
!---(all normal vectors shall be defined in the same way)
      if (lnorch) then
         if (ndmgeo.eq.3) then
            cpuchk=qscput(4,0,ierr)
            call norchk(nelflf,mxkn2d,knoflf,ityflf,ieldum)
            write(*,930) qscput(4,1,ierr)
         else
            cpuchk=qscput(4,0,ierr)
            call norch2(nelflf,mxkn2d,knoflf,ityflf,ieldum)
            write(*,930) qscput(4,1,ierr)
         end if
      end if
!
!---determine normal vectors for elements and nodes
      spaang=z0
      do 10 ielflf=1,nelflf
         call eleflf(knoflf,ityflf,necflf,                              &
     &        xyzflf,vcnflf,ielflf,spalok)
         spaang=spaang+spalok
   10 continue
!
      write(*,910) spaang*z180/pi,720/(4-ndmgeo),360/(4-ndmgeo),0
!
!---scale nodal normal vectors (in this loop nodes are encountered more than once!!!
!---this is no problem, however since scaling of unit vectors yields unit vectors)
      do 20 ielflf=1,nelflf
         do 30 iecke=1,necflf(ielflf)
            kno=knoflf(iecke,ielflf)
            valnrm=z0
            do 40 idim=1,ndmgeo
               valnrm=valnrm+vcnflf(kno,idim)*vcnflf(kno,idim)
   40       continue
!---prevent sqrt(zero) for nonregular meshes
            if (valnrm.gt.0) then
               valnrm=z1/sqrt(valnrm)
            end if
            do 50 idim=1,ndmgeo
               vcnflf(kno,idim)=vcnflf(kno,idim)*valnrm
   50       continue
   30    continue
   20 continue
!
      return
  900 format(1x,'surface number ',i6,' normal vector ',3(g12.5,1x))
  910 format(1x,'Space angle for the origin <--> surface ',g12.5,/,     &
     &          'NOTE: If the surface is closed, the space angle',/,    &
     &          '   should have the following values: ',/,              &
     &          '   --> ',i3,' if the origin is  inside the surface',/, &
     &          '   --> ',i3,' if the origin is      on the surface',/, &
     &          '   --> ',i3,' if the origin is outside the surface')
  930 format(1x,'CPU-Time for surface structure check ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine eleflf(knoflf,ityflf,necflf,                           &
     &                  xyzflf,vcnflf,                                  &
     &                  ielflf,spalok)
      include 'neurofem.inc'
      dimension knoflf(mxkn2d,nelflf),ityflf(nelflf),necflf(nelflf)
      dimension lokkno(32)
      dimension xyzflf(npoflf,3),vcnflf(npoflf,3)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn),      &
     &          adjunc(3),xyzglo(3),rjacob(3,3),vcnlok(mxkn2d,3)
!
      ndmgeo=3
!---initialize everything
      call iset(32,0,lokkno,1)
      spalok=z0
      gamma=z0
!
!---turn to local node numbers
      do 10 i=1,necflf(ielflf)
        kno=knoflf(i,ielflf)
        lokkno(i)=kno
        do 15 idim=1,ndmgeo
           xlk(idim,i) = xyzflf(kno,idim)
           vcnlok(i,idim)=z0
           fag(idim,i)=z0
   15   continue
   10 continue
      call itpanz(ityflf(ielflf),necflf(ielflf),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'eleflf 1'
!
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
         call itpdat(ityflf(ielflf),intgrd,iint,dinpkt,ierr)
         if (ierr.ne.0) stop 'eleflf 2'
         call fofagn(ityflf(ielflf),dinpkt(1),dinpkt(2),dinpkt(3),      &
     &        lokkno(1),xlk,fwe,fag,detj,ierr)
         if (ierr.ne.0) stop 'eleflf 3'
!
!---determine global coordinates at gaussian point
         call dset(ndmgeo,z0,xyzglo,1)
         do 90 iecke=1,necflf(ielflf)
            do 100 idim=1,ndmgeo
               xyzglo(idim)=xyzglo(idim)+fwe(iecke)*xlk(idim,iecke)
  100       continue
   90    continue
!
!---determine length of the vector from the origin to the gaussian point
         veclen=z0
         do 110 idim=1,ndmgeo
            veclen=veclen+xyzglo(idim)*xyzglo(idim)
  110    continue
         veclen=sqrt(veclen)
!
!---determine components of the normal vector r_{n}
         call getadj(ityflf(ielflf),adjunc(1),adjunc(2),adjunc(3),ierr)
         if(ierr.ne.0) stop 'eleflf 4'
!
!---scale r_{n} to unit length
         absrno=z0
         do 70 idim=1,ndmgeo
            absrno=absrno+adjunc(  idim)*adjunc(  idim)
   70    continue
         absrno=z1/sqrt(absrno)
         do 120 idim=1,ndmgeo
            adjunc(  idim)=adjunc(  idim)*absrno
  120    continue
!
!---determine scalar product between normal and radius vectors
         cosang=z0
         do 140 idim=1,ndmgeo
            cosang=cosang+xyzglo(idim)*adjunc(idim)
  140    continue
!
!---normalize cosine (normal vector is already of unit length)
         cosang=cosang/veclen
!
!---evaluate space angle integrand at gaussian point
         spawin=cosang/(veclen**(ndmgeo-1))
!
         xfak=dinpkt(4)*detj
         do 30 iecke=1,necflf(ielflf)
            kno=knoflf(iecke,ielflf)
            gamma=gamma+xfak*fwe(iecke)
            spalok=spalok+xfak*spawin*fwe(iecke)
            do 40 idmgeo=1,ndmgeo
               vcnlok(iecke ,idmgeo)=vcnlok(iecke ,idmgeo)+             &
     &              adjunc(  idmgeo)*xfak*fwe(iecke)
   40       continue
   30    continue
!
   20 continue
!
!---distribute normal vector components to nodes
      do 150 iecke=1,necflf(ielflf)
         kno=knoflf(iecke,ielflf)
         do 160 idmgeo=1,ndmgeo
            vcnflf(kno,idmgeo)=vcnflf(kno,idmgeo)+vcnlok(iecke,idmgeo)
  160    continue
  150 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!  Subroutine NORCHK
!             ======
!
!  (c) Richard Schoenen, G. Knoll, Okt 1994
!
! EINGABE:
! nelem                 Anzahl der Elemente
! mknpel                max. Anzahl Knoten pro Element
! iknel(mknpel,nelem)   Knoten-Element-Zuordnung
! ielty       (nelem)   Elementtyp
!
! HILFSFELDER:
!
! LOKAL:
! list1         ( 5  )  Hilfs-Liste fuer Kantenknoten
! list2         ( 5  )  Hilfs-Liste fuer Kantenknoten
! ngef1         ( 4  )  Zaehler fuer gefundene Kanten
! unterstuetzt werden die Elemente vom Typ 302 und 312
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine norchk(nelem ,mknpel, iknel , ielty, ielrf)
      dimension iknel(mknpel,nelem),ielty(nelem),ielrf(nelem),          &
     &          list1(5),list2(5),ngef1(4)
!
      nknty(ityp)=(ityp-302)/10+3
!
      write(*,7000)
 7000 format(1x,'Testing the surface structure ...')
!
      do i=1,nelem
        ielrf(i)=i
      end do
      nflae=1
      ndreh=0
      neben=0
!
      ianf=1
      iend=1
!
! Schleife ueber Ebenen
!
   10 continue
!
! Schleife ueber Elemente der Ebene
!
      iesav=iend
      do ir1=ianf,iend
         ie1=ielrf(ir1)
         nk1=nknty(ielty(ie1))
         call icopy(nk1,iknel(1,ie1),1,list1,1)
         list1(nk1+1)=list1(1)
         call iset(nk1,0,ngef1,1)
!
         do ir2=iend+1,nelem
            ie2=ielrf(ir2)
            nk2=nknty(ielty(ie2))
            call icopy(nk2,iknel(1,ie2),1,list2,1)
            list2(nk2+1)=list2(1)
!
! gemeinsame Knoten zaehlen
!
            ngem=0
            do i1=1,nk1
               k1=list1(i1)
               do i2=1,nk2
                  if(k1.eq.list2(i2)) ngem=ngem+1
               end do
            end do
!
! nur bei 2 gemeinsamen Knoten kann eine gemeinsame Kante existieren
!
            if(ngem.gt.1)then
               if(ngem.gt.2)then
                  write(*,7004) ie1,ie2,ngem
 7004             format(/,                                             &
     &            '  Flaeche nicht regulaer, ',i1,' (>2) gemeinsame',/, &
     &            '  Knoten zw. den Elementen ',2i7)
               end if
!
               igef=1
               do ik1=1,nk1
                  ka1=list1(ik1  )
                  ke1=list1(ik1+1)
                  do ik2=1,nk2
                     ka2=list2(ik2  )
                     ke2=list2(ik2+1)
                     if( (ka1.eq.ke2.and.ke1.eq.ka2).or.                &
     &                   (ka1.eq.ka2.and.ke1.eq.ke2)   ) goto 20
                  end do
               end do
               igef=0
!
   20          continue
               if(igef.eq.1)then
                  if(ka1.eq.ka2.and.ke1.eq.ke2)then
!
! Kante falsch rum --> Element rumdrehen
! !! Achtung !!
!    Hier werden nur die Hauptknoten vertauscht
!
                     do ik=1,nk2/2
                        ii=nk2+1-ik
                        idum=iknel(ik,ie2)
                        iknel(ik,ie2)=iknel(ii,ie2)
                        iknel(ii,ie2)=idum
                     end do
                     ndreh=ndreh+1
                  end if
!
                  ngef1(ik1)=ngef1(ik1)+1
!
                  iend=iend+1
                  idum        = ielrf(iend)
                  ielrf(iend) = ielrf(ir2 )
                  ielrf(ir2 ) = idum
               end if
!
            end if
!
         end do
!
! Test, ob Kante in mehr als 2 Elementen vorkommt
!
         do ik1=1,nk1
            if(ngef1(ik1).gt.1)then
               write(*,7001) ie1,list1(ik1),list1(ik1+1)
 7001          format(1x,/,1x,'  Kante in mehr als 2 Elementen',/,      &
     &                1x,'  Element:',i7,', Knoten:',I7,',',I7)
            end if
         end do
      end do
!
! Test, ob alle Elemente abgearbeitet
!
      neben=neben+1
      if(iend.lt.nelem)then
         if(iend.eq.iesav)then
            nflae=nflae+1
            iend=iend+1
         end if
         ianf=iesav+1
         goto 10
      end if
!
!
      write(*,7002) nflae,ndreh,neben
 7002 format(/,1x,'  Number of independent Surface structures :',i7,/,  &
     &       1x,  '  Number of re-oriented surface elements   :',i7,/,  &
     &       1x,  '  Number of internal ordering stages       :',i7,/)
!
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!  Subroutine NORCH2   (2D-Version)
!             ======
!
!  (c) Richard Schoenen, G. Knoll, Okt 1994
!
! EINGABE:
! nelem                 Anzahl der Elemente
! mknpel                max. Anzahl Knoten pro Element
! iknel(mknpel,nelem)   Knoten-Element-Zuordnung
! ielty       (nelem)   Elementtyp
!
! HILFSFELDER:
!
! LOKAL:
! ngef1         ( 2  )  Zaehler fuer gefundene Knoten
!
! unterstuetzt werden nur Elemente vom Typ 201
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine norch2(nelem ,mknpel, iknel , ielty, ielrf)
      dimension iknel(mknpel,nelem),ielty(nelem),ielrf(nelem),          &
     &          ngef1(2)
!
      write(*,7000)
 7000 format(1x,'Testing the surface structure ...')
!
      do i=1,nelem
        ielrf(i)=i
        if(ielty(i).ne.201) stop'S-cauchy.norch2: falscher Element-Typ'
      end do
      nflae=1
      ndreh=0
      neben=0
!
      ianf=1
      iend=1
!
! Schleife ueber Ebenen
!
   10 continue
!
! Schleife ueber Elemente der Ebene
!
      iesav=iend
      do ir1=ianf,iend
         ie1=ielrf(ir1)
         ka1=iknel(1,ie1)
         ke1=iknel(2,ie1)
         ngef1(1)=0
         ngef1(2)=0
!
         do ir2=iend+1,nelem
            ie2=ielrf(ir2)
            ka2=iknel(1,ie2)
            ke2=iknel(2,ie2)
!
! gemeinsame Knoten pruefen
!
            nga1=0
            nge1=0
            if(ka1.eq.ka2.or.ka1.eq.ke2) nga1=nga1+1
            if(ke1.eq.ka2.or.ke1.eq.ke2) nge1=nge1+1
            ngem=nga1+nge1
!
            if(ngem.gt.0)then
               if(ngem.gt.1)then
                  write(*,7004) ie1,ie2,ngem
 7004             format(/,                                             &
     &            '  Rand nicht regulaer, ',i1,' (>1) gemeinsame',/,    &
     &            '  Knoten zw. den Rand-Elementen ',2i7)
               end if
!
               if(ka2.eq.ka1.or.ke2.eq.ke1)then
!
! Stab falsch rum --> Element rumdrehen
! !! Achtung !!
!    Hier werden nur die Hauptknoten vertauscht
!
                  iknel(1,ie2)=ke2
                  iknel(2,ie2)=ka2
                  ndreh=ndreh+1
               end if
               ngef1(1)=ngef1(1)+nga1
               ngef1(2)=ngef1(2)+nge1
!
               iend=iend+1
               idum        = ielrf(iend)
               ielrf(iend) = ielrf(ir2 )
               ielrf(ir2 ) = idum
            end if
!
         end do
!
! Test, ob Knoten in mehr als 2 Elementen vorkommt
!
         do ik=1,2
            if(ngef1(ik).gt.1)then
               write(*,7001) iknel(ik,ie1)
 7001          format(/,'  Knoten in mehr als 2 Elementen',/,           &
     &                  '  Knoten:',i7,',',i7)
            end if
         end do
      end do
!
! Test, ob alle Rand-Elemente abgearbeitet
!
      neben=neben+1
      if(iend.lt.nelem)then
         if(iend.eq.iesav)then
            nflae=nflae+1
            iend=iend+1
         end if
         ianf=iesav+1
         goto 10
      end if
!
!
      write(*,7002) nflae,ndreh,neben
 7002 format(/,1x,'  Number of independent Surface structures :',i7,/,  &
     &       1x,  '  Number of re-oriented surface elements   :',i7,/,  &
     &       1x,  '  Number of internal ordering stages       :',i7,/)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine nodadj(neigeo,xyzflf,xyzgeo)
      include 'neurofem.inc'
      dimension neigeo(npoflf)
      dimension xyzflf(npoflf,3),xyzgeo(npogeo,3)
      dimension rflf(3)
!
      save icount
      data icount /0/
      ndmgeo=3
!
!---locates most adjacent volume node for each surface node
      if (icount.eq.0) then
         write(*,910) 
      end if
      if (invway.eq.2) then
         if (icount.eq.0) then
            call percen(0,npoflf)
         end if
         dstmax=z0
         do 10 ipoflf=1,npoflf
            do 20 k=1,ndmgeo
               rflf(k)=xyzflf(ipoflf,k)
   20       continue
            rmin=bignum
            do 30 ipogeo=1,npogeo
               dstflf=z0
               do 40 k=1,ndmgeo
                  dum=xyzgeo(ipogeo,k)-rflf(k)
                  dstflf=dstflf+dum*dum
   40          continue
!%%%               dstflf=sqrt(dstflf)
               if (dstflf.lt.rmin) then
                  neigeo(ipoflf)=ipogeo
                  rmin=dstflf
               end if
   30       continue
            dstmax=max(dstmax,rmin)
            if (icount.eq.0) then
               call percen(ipoflf,npoflf)
            end if
   10    continue
!%%% Wurzelziehen reicht einmal!
         dstmax=sqrt(dstmax)
         if (dstmax.gt.10) then
            write(*,900) dstmax
         end if
         else if (invway.eq.1) then
            write(*,920) 
            do 50 i=1,npoflf
               neigeo(i)=i
   50       continue
         end if
      icount=1
!
      return
  900 format(' Max. distance to a volume node ',g12.5)
  910 format(' Finding neighbours for all influence geometry nodes: ')
  920 format(' Shortcut search if influence equals volume geometry')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!
!
!
