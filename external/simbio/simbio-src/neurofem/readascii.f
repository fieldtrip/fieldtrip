!  $1	21.03.2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     Library for ascii file I/O
!         No usage of FORTRAN 'include' Statements
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

!------------SUBROUTINES ZUM LESEN
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine inigeo(filnam,npogeo,nelgeo,ndmgeo)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character zeil*80,filnam*80
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'error in inigeo'
      end if
 10   read (iunit,'(a)',end=100) zeil
      call qctolo(zeil)
!
!---suche steuerkarte
      isteu=index(zeil,'steuerkarte')
      if(isteu.gt.0) then
 20      read(iunit,900,end=100) zeil(1:80)
         call qctolo(zeil)
         iend=index(zeil,'eoi - ')
         if(iend.gt.0) then
            goto 100
         end if
!
!---suche nach stichworten und einlesen
         ikno=index(zeil,'anzahl der knoten')
         if(ikno.gt.0) then
            read(zeil,910) npogeo
         end if
         iel=index(zeil,'anzahl der elemente')
         if(iel.gt.0) then
            read(zeil,910) nelgeo
         end if
         idim=index(zeil,'geometr. struktur - dimension')
         if(idim.gt.0) then
            read(zeil,910) ndmgeo
         end if
         goto 20
      end if
      goto 10
!
 100  call qfclos(iunit,0)
!
      return
  900 format(a)
  910 format(bn,34x,i20)
  920 format(bn,34x,e12.0)
 1000 stop 'error while reading in inigeo'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reageo(itys3d,knos3d,necs3d,xyzs3d,
     &     nels3d,mxkn3d,npos3d,ndms3d,filnam)
      implicit double precision (a-h,o-z)
      character zeil*80,filnam*80
      integer qfreeu
      dimension itys3d(nels3d),knos3d(mxkn3d,nels3d),necs3d(nels3d)
      dimension  xyzs3d(npos3d,ndms3d)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'error in reageo'
      end if
      ikart=0
 10   read (iunit,'(a)',end=30) zeil
      call qctolo(zeil)
      ikoor=index(zeil,'boi - koordinatenkarte')
      if (ikoor.gt.0) then
         ikart=ikart+1
         if (ndms3d.eq.2) then
            read(iunit,900) ((xyzs3d(i,j),j=1,ndms3d),i=1,npos3d)
         else if (ndms3d.eq.3) then
            read(iunit,910) ((xyzs3d(i,j),j=1,ndms3d),i=1,npos3d)
         else
            stop 'error in reageo: ndms3d'
         end if
      end if
      ielkno=index(zeil,'boi - elementknotenkarte')
      if (ielkno.gt.0) then
         ikart=ikart+1
         do 20 iel=1,nels3d
            if (iel.eq.nels3d-1) then
               idum=0
            end if
            call reaelm(itys3d(iel),knos3d(1,iel),iunit,iel,zeil,
     &                  mxkn3d)
            necs3d(iel)=ncorn(knos3d(1,iel),mxkn3d)
 20      continue
      end if
      if (ikart.ge.2) goto 30
      goto 10
 30   call qfclos(iunit,0)
!
      return
 900  format (bn,3(2x,2g12.0))
 910  format (bn,2x,3g12.0,4x,3g12.0)
 9999 stop 'error reageo'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaelm(ieltyp,loknzu,iunit,iel,zeil,maxkno)
      implicit double precision (a-h,o-z)
      character zeil*80
      dimension loknzu(maxkno)
!
      if (iel.eq.1) then
         read(iunit,'(a)') zeil
         call qctolo(zeil)
      end if
      idop=index(zeil,':')
      if (idop.gt.0) then
!        read(zeil,900) ieltyp,(loknzu(i),i=1,min(maxkno,12))
!         inumb=qclen(zeil)-7
         inumb=LEN(zeil)-7
         inumb=inumb/6
         inumb=min(maxkno,inumb)
         read(zeil,900) ieltyp,(loknzu(i),i=1,min(inumb,12))
         read (iunit,'(a)') zeil
         call qctolo(zeil)
         idop=index(zeil,':')
         iend=index(zeil,'eoi - ')
         if (idop.gt.0.or.iend.gt.0) then
            return
         else
!
!---1. fortsetzungszeile
            read(zeil,910) (loknzu(j),j=13,min(maxkno,24))
            read (iunit,'(a)') zeil
            call qctolo(zeil)
            idop=index(zeil,':')
            iend=index(zeil,'eoi - ')
            if (idop.gt.0.or.iend.gt.0) then
               return
            else
!
!---2. fortsetzungszeile
               read(zeil,910) (loknzu(j),j=25,min(maxkno,32))
               return
            end if
         end if
      else
         return
      end if
!
      return
 900  format(bn,2x,i3,2x,12(i6,:))
 910  format(bn,      7x,12(i6,:))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reapar(itys3d,knos3d,necs3d,xyzs3d,ipart,
     &     nels3d,mxkn3d,npos3d,ndms3d,filnam)
      implicit double precision (a-h,o-z)
      character zeil*80,filnam*80
      integer qfreeu
      dimension itys3d(nels3d),knos3d(mxkn3d,nels3d),necs3d(nels3d)
      dimension  xyzs3d(npos3d,ndms3d),ipart(nels3d) 
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'error in reageo'
      end if
      ikart=0
 10   read (iunit,'(a)',end=30) zeil
      call qctolo(zeil)
      ikoor=index(zeil,'boi - koordinatenkarte')
      if (ikoor.gt.0) then
         ikart=ikart+1
         if (ndms3d.eq.2) then
            read(iunit,900) ((xyzs3d(i,j),j=1,ndms3d),i=1,npos3d)
         else if (ndms3d.eq.3) then
            read(iunit,910) ((xyzs3d(i,j),j=1,ndms3d),i=1,npos3d)
         else
            stop 'error in reageo: ndms3d'
         end if
      end if
      ielkno=index(zeil,'boi - elementknotenkarte')
      if (ielkno.gt.0) then
         ikart=ikart+1
         do 20 iel=1,nels3d
            if (iel.eq.nels3d-1) then
               idum=0
            end if
            call elmpar(itys3d(iel),knos3d(1,iel),ipart(iel),iunit,iel,
     &                  zeil,mxkn3d)
            necs3d(iel)=ncorn(knos3d(1,iel),mxkn3d)
 20      continue
      end if
      if (ikart.ge.2) goto 30
      goto 10
 30   call qfclos(iunit,0)
!
      return
 900  format (bn,3(2x,2g12.0))
 910  format (bn,2x,3g12.0,4x,3g12.0)
 9999 stop 'error reageo'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elmpar(ieltyp,loknzu,ipart,iunit,iel,zeil,maxkno)
      implicit double precision (a-h,o-z)
      character zeil*80
      dimension loknzu(maxkno)
	  integer qclen
!
      if (iel.eq.1) then
         read(iunit,'(a)') zeil
         call qctolo(zeil)
      end if
      idop=index(zeil,':')
      if (idop.gt.0) then
!         read(zeil,900) ieltyp,(loknzu(i),i=1,min(maxkno,12)),ipart
         inumb=qclen(zeil)-7
		 inumb=inumb/6 -1
      	 if (inumb.gt.1) then
            read(zeil,900) ieltyp,(loknzu(i),i=1,inumb),ipart
         else 
            read(zeil,900) ieltyp,(loknzu(i),i=1,min(maxkno,12))
         end if
         read (iunit,'(a)') zeil
         call qctolo(zeil)
         idop=index(zeil,':')
         iend=index(zeil,'eoi - ')
         if (idop.gt.0.or.iend.gt.0) then
            return
         end if
      else
         return
      end if
!
      return
 900  format(bn,2x,i3,2x,12(i6,:))
 910  format(bn,      7x,12(i6,:))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine iniflf(filnam,nelflf)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character*80 filnam,zeil
!
      nelflf=0
      iflstr=0
      iflend=0
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'iniflf: error'
      end if
 10   read (iunit,'(a)',end=30) zeil
      call qctolo(zeil)
!
!---finish
      idum=index(zeil,'eoi - flaechenknotenkarte')
      if (idum.gt.0) then
         iflend=1
         goto 30
      end if
!
!---count surfaces
      if(iflstr.gt.0) nelflf=nelflf+1
!
!--start surface counting
      idum=index(zeil,'boi - flaechenknotenkarte')
      if (idum.gt.0) then
         iflstr=1
      end if
      goto 10
 30   call qfclos(iunit,0)
      if (iflstr.eq.0) write(*,900)
      if (iflend.eq.0) write(*,910)
!
      return
 900  format('iniflf: boi - flaechenknotenkarte not found')
 910  format('iniflf: eoi - flaechenknotenkarte not found')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaflf(filnam,knoflf,ityflf,necflf,maxkno,nelflf)
      implicit double precision (a-h,o-z)
      character*80 zeil,filnam
      integer qfreeu
      dimension knoflf(maxkno,nelflf),ityflf(nelflf),necflf(nelflf)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'reaflf: error'
      end if
 10   read (iunit,'(a)',end=30) zeil
      call qctolo(zeil)
      idum=index(zeil,'boi - flaechenknotenkarte')
      if (idum.gt.0) then
         do 20 iel=1,nelflf
            call reaelm(ityflf(iel),knoflf(1,iel),iunit,iel,zeil,
     &                  maxkno)
            necflf(iel)=ncorn(knoflf(1,iel),maxkno)
            if (index(zeil,'eoi - ').gt.0) then
               nflae=iel
               goto 30
            end if
 20      continue
         if (index(zeil,'eoi - ').le.0) stop 'to many surfaces'
      end if
      goto 10
 30   call qfclos(iunit,0)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaknw(filnam,markkn,werkno,npoin,iknowe,numfre)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      integer qfreeu,qclen
      character*80 zeil, filnam
      dimension markkn(npoin,numfre),iread(6)
      dimension werkno(npoin,numfre),wread(4)
!
!--initialize arrays
      call iset(npoin*numfre, 0,markkn,1)
      call dset(npoin*numfre,z0,werkno,1)
!
      markfl=0
      iknowe=0
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'error in reaknw'
      end if
 10   read (iunit,'(a)',end=100) zeil
      call qctolo(zeil)
      if(index(zeil,'boi - knotenwertefile').gt.0) markfl=1
      if(index(zeil,'eoi - knotenwertefile').gt.0) goto 100
!
      if ( index(zeil,'boi - skalar').gt.0 ) then
!%%%
         if (numfre.ne.1) then
            write(*,950) ' reaknw: the  file ',filnam(1:qclen(filnam)),
     &                   ' is scalar, but',numfre,
     &                   ' degrees of freedom at each node are needed.'
            stop 'reaknw'
         end if
!%%%
 20      read (iunit,'(a)',end=100) zeil
         call qctolo(zeil)
         if ( index(zeil,'eoi - skalar').gt.0 ) goto 10
         call iset(3,0,iread,1)
         read(zeil,910) (iread(i),wread(i),i=1,3)
         do 30 i=1,3
            if(iread(i).gt.0) then
               iknowe=iknowe+1
               if (iknowe.gt.npoin*numfre) then
                  write(*,940) 'reaknw: too many entries in file: ',
     &                          filnam(1:qclen(filnam))
                  stop 'reaknw'
               end if
               markkn(iread(i),1)=1
               werkno(iread(i),1)=wread(i)
            end if
 30      continue
         goto 20
      end if
      if ( index(zeil,'boi - 2d-vektor').gt.0 ) then
!%%%
         if (numfre.ne.2) then
            write(*,950) ' reaknw: the  file ',filnam(1:qclen(filnam)),
     &                   ' is 2-dimensional, but',numfre,
     &                   ' degrees of freedom at each node are needed.'
            stop 'reaknw'
         end if
!%%%
  400    read (iunit,'(a)',end=100) zeil
         call qctolo(zeil)
         if ( index(zeil,'eoi - 2d-vektor').gt.0 ) goto 10
         kn2=0
         read(zeil,920) (iread(i),i=1,3), (wread(j),j=1,2),
     &                  (iread(i),i=4,6), (wread(j),j=3,4)
         if (iread(1).ne.0) then
            if (iread(2).ne.0) iknowe=iknowe+1
            if (iread(3).ne.0) iknowe=iknowe+1
            if (iknowe.gt.npoin*numfre) then
               write(*,940) 'reaknw: too many entries in file: ',
     &              filnam(1:qclen(filnam))
               stop 'reaknw'
            end if
            kno=iread(1)
            if (iread(2).ne.0) markkn(kno,1)=1
            if (iread(3).ne.0) markkn(kno,2)=1
            werkno(kno,1)=wread(1)
            werkno(kno,2)=wread(2)
         end if
         if (iread(4).ne.0) then
            if (iread(5).ne.0) iknowe=iknowe+1
            if (iread(6).ne.0) iknowe=iknowe+1
            if (iknowe.gt.npoin*numfre) then
               write(*,940) 'reaknw: too many entries in file: ',
     &              filnam(1:qclen(filnam))
               stop 'reaknw'
            end if
            kno=iread(4)
            if (iread(5).ne.0) markkn(kno,1)=1
            if (iread(6).ne.0) markkn(kno,2)=1
            werkno(kno,1)=wread(3)
            werkno(kno,2)=wread(4)
         end if
         goto 400
      end if
      if ( index(zeil,'boi - 3d-vektor').gt.0 ) then
!%%%
         if (numfre.ne.3) then
            write(*,950) ' reaknw: the  file ',filnam(1:qclen(filnam)),
     &                   ' is 3-dimensional, but',numfre,
     &                   ' degrees of freedom at each node are needed.'
            stop 'reaknw'
         end if
!%%%
  200    read (iunit,'(a)',end=100) zeil
         call qctolo(zeil)
         if ( index(zeil,'eoi - 3d-vektor').gt.0 ) goto 10
         read(zeil,930) kno,(markkn(kno,i),i=1,3),(werkno(kno,i),i=1,3)
         do 210 k=1,3
            if (markkn(kno,k).ne.0) then
               iknowe=iknowe+1
               if (iknowe.gt.npoin*numfre) then
                  write(*,940) 'reaknw: too many entries in file: ',
     &                          filnam(1:qclen(filnam))
                  stop 'reaknw'
               end if
               markkn(kno,k)=1
            end if
  210    continue
         goto 200
      end if
      goto 10
 100  call qfclos(iunit,0)
!
      if (markfl.ne.1) then
         write(*,900)
         stop
      end if
      return
  900 format('keyword for node-value file is missing!')
  910 format(bn,3(2x,i6,2x,f12.0))
  920 format(bn,2(2x,i6,2x,2i1,2f12.0))
  930 format(bn,2x,i6,2x,3i1,2x,3f12.0)
  940 format(a,a)
  950 format(a,a,/,a,i2,a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!------------SUBROUTINES ZUM SCHREIBEN
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine georad (  ityp, necke, knoel, xyz, wert,
     &     filnam,namlen,nelem,npoin,maxkno,radius,ndim)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character  filnam*(*)
!
      dimension  ityp(nelem), necke(nelem), knoel(maxkno,nelem)
      dimension  xyz(npoin,ndim-1),wert(npoin)
!
      wmax=abs(wert(idamax(npoin,wert,1)))
      fdum=radius/(wmax*3.d0)
      iunit = qfreeu()
      call qfopse(iunit,filnam(1:namlen),'un','fo',ierr)
      write(iunit,900) npoin,nelem,ndim
      if (ndim.eq.3) then
         write (iunit,910) ( ( radius+fdum*wert(i) )*cos(xyz(i,1)),
     &        ( radius+fdum*wert(i) )*sin(xyz(i,1)),
     &          radius     *  xyz(i,2) ,i=1,npoin )
      else
         write (iunit,810) ( ( radius+fdum*wert(i) )*cos(xyz(i,1)),
     &        ( radius+fdum*wert(i) )*sin(xyz(i,1)),i=1,npoin )
      end if
      write (iunit,920) 
!
      do 10 ielem=1,nelem
         write(iunit,930) ityp(ielem)+100,
     &        (knoel(i,ielem),i=1,necke(ielem))
 10   continue
      write(iunit,940)
      call qfclos(iunit,0)
!
      return
 900  format ('BOI - GEOMETRIEFILE',/,78('='),/,78('='),/,
     &        'BOI - STEUERKARTE',/,
     &        '  ANZAHL DER KNOTEN             : ',I6,/,
     &        '  ANZAHL DER ELEMENTE           : ',I6,/,
     &        '  GEOMETR. STRUKTUR - DIMENSION :      ',i1,/,
     &        'EOI - STEUERKARTE',/,78('='),/,78('='),/,
     &        'BOI - KOORDINATENKARTE' )
 910  format (bn,2x,3d12.5,4x,3d12.5)
 810  format (bn,3(2x,2d12.5))
 920  format ('EOI - KOORDINATENKARTE',/,78('='),/,78('='),/,
     &        'BOI - ELEMENTKNOTENKARTE' )
 930  format (bn,2x,i3,': ',12i6)
 940  format ('EOI - ELEMENTKNOTENKARTE',/,78('='),/,78('='),/,
     &        'EOI - GEOMETRIEFILE' )
 9999 stop 'georad'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wrigeo(filnam,namlen,ityp,necke,knoel,xyz,
     &     ndim,nelem,npoin,maxkno)
      implicit double precision (a-h,o-z)
      integer qfreeu 
      character filnam*80
      dimension ityp(nelem),necke(nelem),knoel(maxkno,nelem)
      dimension xyz(npoin,ndim)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         stop 'wrigeo'
      end if
!
!---schreibe steuerkarte
      write(iunit,900) npoin,nelem,ndim
!
!---schreibe koordinatenkarte
      if (ndim.eq.3) then
         write(iunit,910) ((xyz(i,j),j=1,3 ),i=1,npoin)
      else if (ndim.eq.2) then
         write(iunit,810) ((xyz(i,j),j=1,2 ),i=1,npoin)
      else
         stop 'error of -ndim-, subroutine wrigeo'
      end if
      write(iunit,920)
!
!--schreibe elementknotenzuordnung  
      do 10 iel=1,nelem
         write(iunit,930) ityp(iel),(knoel(i,iel),i=1,necke(iel))
 10   continue
      write(iunit,940)
      call qfclos(iunit,0)
!
 900  format ('BOI - GEOMETRIEFILE',/,78('='),/,78('='),/,
     &        'BOI - STEUERKARTE',/,
     &        '  ANZAHL DER KNOTEN             : ',I6,/,
     &        '  ANZAHL DER ELEMENTE           : ',I6,/,
     &        '  GEOMETR. STRUKTUR - DIMENSION :      ',i1,/,
     &        'EOI - STEUERKARTE',/,78('='),/,78('='),/,
     &        'BOI - KOORDINATENKARTE' )
 910  format (bn,2x,3d12.5,4x,3d12.5)
 810  format (bn,3(2x,2d12.5))
 920  format ('EOI - KOORDINATENKARTE',/,78('='),/,78('='),/,
     &        'BOI - ELEMENTKNOTENKARTE' )
 930  format (bn,2x,i3,': ',12i6,:,/,7x,12i6,:,/,7x,12i6)
 940  format ('EOI - ELEMENTKNOTENKARTE',/,78('='),/,78('='),/,
     &        'EOI - GEOMETRIEFILE' )
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!@modul-creerg
!@======================================================================
! Erstellen eines Ergebnisfiles
!=======================================================================
      subroutine creerg(fname,npoin,nelem,ndim,aktion,nerge)
      implicit double precision (a-h,o-z)
      parameter(nwepre=500)
      character*(*) fname,aktion
      integer qfreeu,qfrecl
      dimension ibuff(1000),  kactre(3,1)
 
!.....Recordlaenge des Ergebnisplotfiles bestimmen
 
      lrerpl=qfrecl(0,nwepre,0)

!.....Steuerdaten generieren
!.....MODERG, IDRERG immer gleich 0 (nur fuer Festigkeitsberechnungen)

      moderg=0
      idrerg=0

!.....Anzahl der Aktionen, immer 1 wenn HPPL oder TEHD

      nactio=1

!.....Anzahl der Rekords im Ergebnisfile

      nrerpl=2+nerge*((npoin-1)/nwepre+1)
 
!.....Binaeren Save-File oeffnen
 
      nterpl=qfreeu()
      call qfsize(nrerpl)
      call qfcrdi(nterpl,fname,'un',lrerpl,ierror)
      if (ierror.ne.0)  then
         write(*,8001) fname
         stop 'creerg: error'
      endif
 
!.....Steuerdaten sichern
 
      ibuff(1)=npoin
      ibuff(2)=nelem
      ibuff(3)=ndim
      ibuff(4)=moderg
      ibuff(5)=idrerg
      ibuff(6)=nactio
      ibuff(7)=nrerpl
      ibuff(8)=ndim
      write(nterpl,rec=1) ibuff
 
!.....Kopfdaten generieren und sichern
!.....Nummer der Aktion 12=TEHD

      if (aktion.eq.'hppl') kactre(1,1)=10
      if (aktion.eq.'tehd') kactre(1,1)=12
      if (aktion.eq.'test') kactre(1,1)=13
!
!.....Record, wo die Ergebnisse beginnen

      kactre(2,1)=3

!.....Anzahl der Ergebnisse 

      kactre(3,1)=nerge
      ipoib=0
      do 10 i=1,nactio
         ibuff(ipoib+1)=kactre(1,i)
         ibuff(ipoib+2)=kactre(2,i)
         ibuff(ipoib+3)=kactre(3,i)
         ipoib=ipoib+3
   10 continue
      if (kactre(1,1).eq.13) then
         write(nterpl,rec=2) (ibuff(i),i=1,3),(' ',i=1,13*20)
      else
         write(nterpl,rec=2) ibuff
      endif
      call qfclos(nterpl,0)
      return
 8001 format(1x,'error while writing the result file ',a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!@modul-wrispa
!@======================================================================
! Wegschreiben einer Ergebnisspalte
!=======================================================================
      subroutine wrispa(fname,erge,npoin,ispalt,text)            
      implicit double precision (a-h,o-z)
      parameter(nwepre=500)
      character*(*) fname,text,fotext*20
      dimension erge(npoin),fotext(13)
      integer qfreeu,qfrecl
      logical lnull

!.....ISPALT kleiner als Null -> Spalte mit Nullen fuellen

      if (ispalt.lt.0) then
         ispalt=-ispalt
         lnull=.true.
      else
         lnull=.false.
      endif

!.....Recordlaenge ermitteln

      lrerpl=qfrecl(0,nwepre,0)

!.....Anzahl der Records pro Spalte des Ergebnisplotfiles ermitteln

      nrepwe=(npoin-1)/nwepre+1

!.....Startrecord der gewuenschten Spalte ermitteln

      ireact=3+(ispalt-1)*nrepwe

!.....Ergebnisdatei oeffnen

      ntlese=qfreeu()
      call qfopdi(ntlese,fname,'ol','un',lrerpl,ierror)
      if (ierror.ne.0) goto 9990

!.....NPOIN lesen (zur Kontrolle)

      read(ntlese,rec=1,err=9990) npo
      if (npo.ne.npoin) goto 9990

!.....Zahl der Spalten lesen (zur Kontrolle)

      read(ntlese,rec=2,err=9990) iactio,istart,nspalt
      if (ispalt.gt.nspalt) goto 9990

!.....Text wegschreiben

      if (iactio.eq.13) then
         read(ntlese,rec=2,err=9990) 
     *        iactio,istart,idum,(fotext(i),i=1,nspalt)
         fotext(ispalt)=text
         write(ntlese,rec=2,err=9990)
     *     iactio,istart,nspalt,(fotext(i),i=1,nspalt)
      endif

!.....Werte wegschreiben

      nloop=npoin/nwepre+min( mod(npoin,nwepre),1 )
      do 100 iloop=1,nloop
         i1=(iloop-1)*nwepre+1
         i2=min(i1+nwepre-1,npoin)
         if (lnull) then
            write(ntlese,rec=ireact,err=9990) (0.d0,i=i1,i2)
         else
            write(ntlese,rec=ireact,err=9990) (erge(i),i=i1,i2)
         endif
         ireact=ireact+1
 100  continue
      call qfclos(ntlese,0)
      return
 9990 stop 'wrispa: error while writing of a result column '
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!@modul-wrierg
!@======================================================================
! Eine Spalte in einen EIPP-Ergebnisplotfile schreiben
!
! Eingabeparameter      Beschreibung
!  nterpl                File-Pointer
!  erge(npoin)           Ergebniswerte
!  ispalt                Spalte die gelesen werden soll
!  ierge                 =1 --> Spalte wegschreiben, 
!                        =0 --> Spalte mit Nullen fuellen
!  npoin
!=======================================================================
      subroutine wrierg(nterpl,erge,ispalt,ierge,npoin)            
      implicit double precision (a-h,o-z)
      parameter(nwepre=500)
      dimension buffer(nwepre), erge(npoin) 
      integer qfrecl
 
!.....Recordlaenge des Ergebnisplotfiles bestimmen
 
      lrerpl=qfrecl(0,nwepre,0)

!.....Die ersten zwei Records beinhalten Kopfdaten (NPOIN,NELEM,..)

      irebeg=3

!.....Anzahl der Records fuer das Ergebnisfeld ermitteln

      nrepwe=(npoin-1)/nwepre+1

!.....Start-Record festlegen

      irakt=irebeg+(ispalt-1)*nrepwe

!.....Zaehler initialisieren

      ianf=1
      nloop=npoin/nwepre
      nrest=mod(npoin,nwepre)

!.....Ergebnisfeld wegschreiben

      if (ierge.eq.0) call dset(nwepre,0.0d0,buffer,1)
      do 50 iloop=1,nloop
         if (ierge.eq.1) call dcopy(nwepre,erge(ianf),1,buffer,1)
!$$$         call writra(nterpl,irakt,buffer,lrerpl,*9990)
         write(nterpl,rec=irakt) buffer
         irakt=irakt+1
         ianf=ianf+nwepre
 50   continue

!.....Rest des Ergebnisfeldes wegschreiben, der kein ganzes Record fuellt

      if (nrest.ne.0) then
         if (ierge.eq.1) call dcopy(nrest,erge(ianf),1,buffer,1)
!$$$         call writra(nterpl,irakt,buffer,lrerpl,*9990)
         write(nterpl,rec=irakt) buffer
         irakt=irakt+1
      end if
      return
 9990 write(*,9991)
 9991 format(1x,'error while writing of a result file ')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine eiperg( vecerg,vecto1,
     &                   vecto2,vecto3,vecto4,vecto5,
     &                   vecto6,vecto7,vecto8,vecto9,
     &                   vect10,vect11,vect12,vect13,
     &                   npof2d,nelf2d,nglodm,nspalt,
     &                   filnam)
      implicit double precision (a-h,o-z)
      parameter (nwepre=500,z1=1.d0)
      integer qfreeu,qfrecl
      character filnam*80
      dimension  ibuff(2*nwepre)
      dimension buffer(nwepre)
      dimension vecerg(npof2d),vecto1(npof2d),vecto2(npof2d),
     &          vecto3(npof2d),vecto4(npof2d),vecto5(npof2d),
     &          vecto6(npof2d),vecto7(npof2d),vecto8(npof2d),
     &          vecto9(npof2d),vect10(npof2d),vect11(npof2d),
     &          vect12(npof2d),vect13(npof2d)
!
      if (nspalt.eq.0) goto 1000
      iunit=qfreeu()
      irerpl=qfrecl(0,nwepre,0)
      call qfsize(npof2d)
      call qfcrdi(iunit,filnam,'un',irerpl,ierr)
      if (ierr.ne.0) stop 'error output 1'

      irerpl=qfrecl(0,nwepre,0)
      moderg = 0
      idrerg = 0
      nactio = 1
      nrerpl = 2+nspalt*((npof2d-1)/nwepre+1)
!
!---kopfzeile
      ibuff(1) = npof2d
      ibuff(2) = nelf2d
      ibuff(3) = nglodm
      ibuff(4) = moderg
      ibuff(5) = idrerg
      ibuff(6) = nactio
      ibuff(7) = nrerpl
      ibuff(8) = nglodm
!      call qfwrdi(iunit,1,ibuff,irerpl,ierr)
      write(iunit,rec=1) (ibuff(i),i=1,8)
      if (ierr.gt.0) stop 'eiperg 1'

!
!---zweite zeile
      ibuff(1) = 12
      ibuff(2) = 3
      ibuff(3) = nspalt
      write(iunit,rec=2) (ibuff(i),i=1,3)
!      call qfwrdi(iunit,2,ibuff,irerpl,ierr)
      if (ierr.gt.0) stop 'eiperg 2'
!
      call dcopy(npof2d,vecto1,1,vecerg,1)
      call wriold(iunit,1,1,3,vecerg,buffer,npof2d)
      if (nspalt.eq.1) goto 900
!
      call dcopy(npof2d,vecto2,1,vecerg,1)
      call wriold(iunit,1,2,3,vecerg,buffer,npof2d)
      if (nspalt.eq.2) goto 900
!
      call dcopy(npof2d,vecto3,1,vecerg,1)
      call wriold(iunit,1,3,3,vecerg,buffer,npof2d)
      if (nspalt.eq.3) goto 900
!
      call dcopy(npof2d,vecto4,1,vecerg,1)
      call wriold(iunit,1,4,3,vecerg,buffer,npof2d)
      if (nspalt.eq.4) goto 900
!
      call dcopy(npof2d,vecto5,1,vecerg,1)
      call wriold(iunit,1,5,3,vecerg,buffer,npof2d)
      if (nspalt.eq.5) goto 900
!
      call dcopy(npof2d,vecto6,1,vecerg,1)
      call wriold(iunit,1,6,3,vecerg,buffer,npof2d)
      if (nspalt.eq.6) goto 900
!
      call dcopy(npof2d,vecto7,1,vecerg,1)
      call wriold(iunit,1,7,3,vecerg,buffer,npof2d)
      if (nspalt.eq.7) goto 900
!
      call dcopy(npof2d,vecto8,1,vecerg,1)
      call wriold(iunit,1,8,3,vecerg,buffer,npof2d)
      if (nspalt.eq.8) goto 900
!
      call dcopy(npof2d,vecto9,1,vecerg,1)
      call wriold(iunit,1,9,3,vecerg,buffer,npof2d)
      if (nspalt.eq.9) goto 900
!
      call dcopy(npof2d,vect10,1,vecerg,1)
      call wriold(iunit,1,10,3,vecerg,buffer,npof2d)
      if (nspalt.eq.10) goto 900
!
      call dcopy(npof2d,vect11,1,vecerg,1)
      call wriold(iunit,1,11,3,vecerg,buffer,npof2d)
      if (nspalt.eq.11) goto 900
!
      call dcopy(npof2d,vect12,1,vecerg,1)
      call wriold(iunit,1,12,3,vecerg,buffer,npof2d)
      if (nspalt.eq.12) goto 900
!
      call dcopy(npof2d,vect13,1,vecerg,1)
      call wriold(iunit,1,13,3,vecerg,buffer,npof2d)
!
 900  call qfclos(iunit,0)
 1000 return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wriold(ntwrit,iactio,iweact,irebeg,erge2d,buffer,
     &                  npof2d)
      implicit double precision (a-h,o-z)
      parameter(nwepre=500)
      integer qfrecl
      dimension buffer(nwepre), erge2d(npof2d)
!
      irerpl=qfrecl(0,nwepre,0)
      nrepwe=(npof2d-1)/nwepre+1
      ireact=irebeg+(iweact-1)*nrepwe
      ianf=0
      nloop=npof2d/nwepre
      do 30 iloop=1,nloop
         do 20 iwepre=1,nwepre
            buffer(iwepre)=erge2d(ianf+iwepre)
 20      continue
         call qfwrdi(ntwrit,ireact,buffer,irerpl,ierr)
         if (ierr.ne.0) stop 'wriold 1'
         ireact=ireact+1
         ianf=ianf+nwepre
 30   continue
      nr=npof2d-ianf
      if (nr.ne.0) then
         do 40 iwepre=1,nr
            buffer(iwepre)=erge2d(ianf+iwepre)
 40      continue
         call qfwrdi(ntwrit,ireact,buffer,irerpl,ierr)
         if (ierr.ne.0) stop 'wriold 2'
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wriknw(ifrgeo,vecwri,npogeo,numfre,filnam)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      character*80 filnam
      integer qfreeu
      dimension ifrgeo(npogeo,numfre),ibuff(4),kno(2)
      dimension vecwri(npogeo,numfre),wbuff(4)
!
!---write nodal-values-file
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,900) 'error filnam in wriknw'
         stop
      end if
      write (iunit,900) 
      if (numfre.eq.1) then
         write (iunit,920)
         icount=0
         do 10 i=1,npogeo
            if (ifrgeo(i,1).gt.0) icount=icount+1
   10    continue
!     
         idum=0
         inode=0
   50    inode=inode+1
         if(inode.gt.npogeo) goto 222
         if (ifrgeo(inode,1).gt.0) then
            icount=icount-1
            idum=idum+1
            ibuff(idum)=inode
            wbuff(idum)=vecwri(inode,1)
            if (idum.eq.3)  then
               write(iunit,980) (ibuff(j),wbuff(j),j=1,idum)
               idum=0
            end if
         end if
         goto 50
  222    write(iunit,980) (ibuff(j),wbuff(j),j=1,idum)
         write (iunit,930) 
      else if (numfre.eq.2) then
         write (iunit,940)
         idum=0
         iwri=0
         do 400 i=1,npogeo
            ifre=0
            do 390 j=1,numfre
               ifre=ifre+ifrgeo(i,j)
  390       continue
            if (ifre.gt.0) then
               idum=idum+1
               kno(idum)=i
               do 410 k=1,numfre
                  iwri=iwri+1
                  if (ifrgeo(i,k).gt.0) then
                     ibuff(iwri)=k
                     wbuff(iwri)=vecwri(i,k)
                  else 
                     ibuff(iwri)=0
                     wbuff(iwri)=z0
                  end if
  410          continue
               if (idum.eq.2) then
                  write (iunit,990) kno(1),ibuff(1),ibuff(2),
     &                                     wbuff(1),wbuff(2),
     &                              kno(2),ibuff(3),ibuff(4),
     &                                     wbuff(3),wbuff(4)
                  idum=0
                  iwri=0
               end if
            end if
  400    continue
         if (idum.gt.0) then
            write (iunit,990) kno(1),(ibuff(j),j=1,2),
     &           (wbuff(k),k=1,2)
         end if
         write (iunit,950)
      else if (numfre.eq.3) then
         write (iunit,960)
         do 500 i=1,npogeo
            idum=0
            do 510 k=1,numfre
               if (ifrgeo(i,k).eq.1) then
                  idum=idum+1
                  ibuff(k)=k
                  wbuff(k)=vecwri(i,k)
               else 
                  ibuff(k)=0
                  wbuff(k)=0.d0
               end if
  510       continue
            if (idum.gt.0) then
               write (iunit,995) i,(ibuff(j),j=1,3),(wbuff(k),k=1,3)
            end if
  500    continue
         write (iunit,970)
      else
         stop 'wriknw: degree of freedom not implemented'
      end if
      write (iunit,910)
      call qfclos(iunit,0)
!
      return
  900 format('BOI - KNOTENWERTEFILE')
  910 format('EOI - KNOTENWERTEFILE')
  920 format('BOI - SKALAR')
  930 format('EOI - SKALAR')
  940 format('BOI - 2D-VEKTOR')
  950 format('EOI - 2D-VEKTOR')
  960 format('BOI - 3D-VEKTOR')
  970 format('EOI - 3D-VEKTOR')
  980 format(bn,3(2x,i6,2x,g12.5,:))
  990 format(bn,2(2x,i6,2x,2i1,2g12.5,:))
  995 format(bn,2x,i6,2x,3i1,2x,3g12.5)
      end 
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reanod(filnam,marked,npoin,nfound)
!
!---liest Knoten von {Knotenfiles,Flaechenfiles,Elementefiles} ein
!   und markiert die entsprechenden Knoten im Vektor "marked"
      implicit double precision (a-h,o-z)
      integer qfreeu
      character*80 filnam,zeil
      dimension marked(npoin),kno(12)
!
      nfound=0
      call iset(npoin,0,marked,1)
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.gt.0) then
         write(*,'(a)') 'error while opening ',filnam,' in reanod'
         stop
      end if
!
!---fileart feststellen
 10   read(iunit,'(a)',end=9998) zeil
      call qctolo(zeil)
      iflae=index(zeil,'boi - flaechenfile')
      ielem=index(zeil,'boi - elementefile')
      iknof=index(zeil,'boi - knotenfile')
      if( (iknof+ielem+iflae).le.0 ) goto 10
!
!---weiter wenn fileart klar, entsprechende karte suchen
 20   read(iunit,'(a)',end=9998) zeil
      call qctolo(zeil)
      if (iflae.gt.0) ikart=index(zeil,'boi - flaechenknotenkarte')
      if (iknof.gt.0) ikart=index(zeil,'boi - knotennummernkarte')
      if (ielem.gt.0) ikart=index(zeil,'boi - elementnummernkarte')
      if (ikart.le.0) goto 20
!
!---wenn richtige karte gefunden, karte bis 'eoi -' abarbeiten
 30   read(iunit,'(a)') zeil
      call qctolo(zeil)
      if (index(zeil,'eoi -').gt.0) goto 300
      read (zeil,'(bn,7x,12i6)',end=9998) (kno(i),i=1,12)
      do 40 i=1,12
         node=kno(i)
         if (node.ne.0) then
            if (marked(node).ne.1) then
               marked(node)=1
               nfound=nfound+1
            end if
         end if
 40   continue
      goto 30
 300  call qfclos(iunit,0)
!
      return
 9998 stop 'error reanod'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!@======================================================================
! 3D-GNUPLOT-File wegschreiben. Es werden drei Spalten weggeschrieben,
! wobei nach den ersten zwei Spalten sortiert wird. Der Feldzugriff
! auf die drei Spalten <x1>, <x2> und <x3> erfolgt ueber ein Pointer-
! Feld <ipoin>.
!
! Eingabeparameter      Beschreibung
!  fname                 Dateiname fuer den zu erstellenden GNUPLOT-File
!  x1(npoin)             Erste Spalte fuer den GNUPLOT-File
!  x2(npoin)             Zweite Spalte fuer den GNUPLOT-File
!  x3(npoin)             Dritte Spalte fuer den GNUPLOT-File
!  ipoin(nanz)           Pointerfeld
!  npoin                 Dimension der Felder <x1>, <x2> und <x3>
!  nanz                  Dimension des Pointerfeldes
!  tolsp2                Toleranz fuer die Blockeinteilung 
!
! Ausgabeparameter      Beschreibung
!  ipoin(nanz)           siehe oben
!=======================================================================
      subroutine gnupl3(fname,x1,x2,x3,ipoin,npoin,nanz,tolsp2)
      implicit double precision (a-h,o-z)
      parameter(tol=1.0d-2)
      dimension x1(npoin),x2(npoin),x3(npoin),ipoin(nanz)
      character fname*(*)
      integer qfreeu

!.....Wenn ipoin(1) Null ist, dann muss das Pointerfeld initialisiert 
!.....werden

      if (ipoin(1).le.0) then
         do 100 i=1,nanz
            ipoin(i)=i
 100     continue
      endif
      call gnsort(x1,x2,ipoin,nanz,tolsp2)
      iunit=qfreeu()
      call qfopse( iunit,fname,'un','fo',ierr)
      if (ierr.gt.0) goto 9990
      do 200 i=1,nanz
         write(iunit,9100) x1(ipoin(i)),x2(ipoin(i)),x3(ipoin(i))
         if (i.ne.nanz) then
            diff=abs(x1(ipoin(i))-x1(ipoin(i+1)))
            if (diff.gt.tol) write(iunit,9200)
         endif
 200  continue
      call qfclos(iunit,0)
      return
 9100 format(f12.5,' ',f12.5,' ',g12.5)
 9200 format()
 9990 stop 'error while creating a gnuplot file !'
      end
!@======================================================================
! Sortieralgorithmus fuer 3D-GNUPLOT-Files. Ein 3D-GNUPLOT-File besteht
! aus zwei Koordinatenspalten und einer Wertespalte. Der Zufriff auf
! die Koordinatenwerte erfolgt ueber ein Pointerfeld <ipoin>, in dem
! nach dem Sortiervorgang die neue Reihenfolge der Indizes fuer die
! Koordinatenfelder <feld1> und <feld2> steht. Sortiert wird erst
! nach der ersten Spalte (<feld1>), und bei gleichen Werten im <feld1>
! wird nach <feld2> sortiert. 
!
! Eingabeparameter      Beschreibung
!  feld1(*)              Erste Spalte des GNUPLOT-Files
!  feld2(*)              Zweite Spalte des GNUPLOT-Files
!  ipoin(n)              Pointerfeld fuer den Zugriff auf die Felder
!  n                     Dimension des Pointerfeldes
!  tol                   Toleranz fuer das Bilden von Bloecken
!
! Ausgabeparameter      Beschreibung
!  ipoin(*)              siehe oben
!===================================================================
      subroutine gnsort(feld1,feld2,ipoin,n,tol)
      implicit double precision (a-h,o-z)
      dimension feld1(*),feld2(*),ipoin(n)

!.....Pointerfeld nach aufsteigender erster Spalte sortieren

      call busort(feld1,ipoin,1,n)

!.....Bloecke mit gleicher erster Spalte nach der zweiten Spalte 
!.....sortieren

      istart=1
      do 100 i=2,n
         diff=abs(feld1(ipoin(i))-feld1(ipoin(istart)))
         if (diff.gt.tol) then
            call busort(feld2,ipoin,istart,i-1)
            istart=i
         endif
 100  continue

!.....Letzten Block sortieren

      call busort(feld2,ipoin,istart,n)
      return
      end
!@======================================================================
! Bubble-Sort Algorithmus. Das Feld <ipoin> ist ein Pointer-Feld, ueber
! das auf die Elemente des zu sortierenden Feldes <feld> zugegriffen
! wird. Somit koennen auch Teilbereiche eines Feldes sortierwerden.
!
! Eingabeparameter      Beschreibung
!  feld(*)               Zahlenwerte, nach denen sortiert wird
!  ipoin(*)              Pointerfeld, sortiert wird feld(ipoin(n1..n2))
!  n1                    Erster Index im Pointerfeld
!  n2                    Letzter Index im Pointerfeld
!
! Ausgabeparameter      Beschreibung
!  ipoin(*)              siehe oben
!=======================================================================
      subroutine busort(feld,ipoin,n1,n2)
      implicit double precision (a-h,o-z)
      dimension feld(*),ipoin(*)
      logical lwechs

!.....Vorwaerts oder rueckwaerts sortieren ?

      if (n1.le.n2) then
         istep=1
      else
         istep=-1
      endif

!.....Sortier-Loop

 100  continue
      lwechs=.false.
      do 200 i=n1+istep,n2,istep
         if (feld(ipoin(i-1)).gt.feld(ipoin(i))) then
            ihilf=ipoin(i)
            ipoin(i)=ipoin(i-1)
            ipoin(i-1)=ihilf
            lwechs=.true.
         endif
 200  continue
      if (lwechs) goto 100
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine iniint(filnam,npof2d,nelf2d,ndmf2d,diaf2d,widf2d)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character zeil*80,filnam*80
!---liest die steuerkarte des interfacefiles
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.ne. 0) stop 'error iniint'
 10   read(iunit,900,end=800) zeil(1:80)
      call qctolo(zeil)
!
!---suche steuerkarte
      isteu=index(zeil,'steuerkarte')
        if(isteu.gt.0) then
   20     read(iunit,900,end=800) zeil(1:80)
          call qctolo(zeil)
          iend=index(zeil,'eoi - ')
          if(iend.gt.0) then
            goto 800
          end if
!
!---suche nach stichworten und einlesen
          ikno=index(zeil,'anzahl der knoten')
          if(ikno.gt.0) then
            read(zeil,910) npof2d
          end if
          iel=index(zeil,'anzahl der elemente')
          if(iel.gt.0) then
            read(zeil,910) nelf2d
          end if
          idim=index(zeil,'geometr. struktur - dimension')
          if(idim.gt.0) then
            read(zeil,910) ndmf2d
          end if
          idim=index(zeil,'lagerdurchmesser')
          if(idim.gt.0) then
            read(zeil,920) diaf2d
          end if
          idim=index(zeil,'lagerbreite')
          if(idim.gt.0) then
            read(zeil,920) widf2d
          end if
          goto 20
          end if
      goto 10
!
  800 call qfclos(iunit,0)
!
      return
  900 format(a)
  910 format(bn,34x,i20)
  920 format(bn,34x,e12.0)
 1000 stop 'error while reading -iniint-'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reaint(knof2d,ktrf23,ityf2d,ielc2d,necf2d,
     &                  xyzf3d,xyzf2d,npof2d,nelf2d,ndmf2d,
     &                  filnam,diaf2d)
      implicit double precision (a-h,o-z)
      integer qfreeu
      character zeil*80,filnam*80
      parameter(bignum=1.d100,z2=2.d0)
      parameter(mxkn2d=12)
      parameter (pi=3.141592653589793d0)
      dimension knof2d(mxkn2d,nelf2d),ktrf23(npof2d),ityf2d(nelf2d),
     &          ielc2d(       nelf2d),necf2d(nelf2d)
      dimension xyzf3d(npof2d,ndmf2d),xyzf2d(npof2d,ndmf2d-1)
      dimension xyzglo(mxkn2d,3)
!
      iunit=qfreeu()
      call qfopse(iunit,filnam,'un','fo',ierr)
      if (ierr.ne. 0) stop 'error reaint'
!
!---liest den restlichen interfacefile
   10 read(iunit,900,end=200) zeil
      call qctolo(zeil)
      if(index(zeil,'boi - ').le.0) goto 10
!
!---koordinatenkarte
      if(index(zeil,'koordinatenkarte').gt.0) then
         if (ndmf2d-1.eq.1) then
            read(iunit,810,end=9998) 
     &      ((xyzf3d(i,j),j=1,ndmf2d), i=1,npof2d)
         else if (ndmf2d-1.eq.2) then
            read(iunit,910,end=9998) 
     &      ((xyzf3d(i,j),j=1,ndmf2d), i=1,npof2d)
         end if
         read(iunit,900,end=9998)  zeil
         call qctolo(zeil)
         if(index(zeil,'eoi - ').le.0) goto 9998
!
!---flaechenkarte
      else if(index(zeil,'flaechenkarte').gt.0) then
         do 20 iel=1,nelf2d
            read(iunit,920,end=9998) ityf2d(iel),
     &           (knof2d(j,iel), j=1,mxkn2d)
!
!---die elementdimension wird um einen Grad erniedrigt
            ityf2d(iel)=ityf2d(iel)-100
   20    continue
         read(iunit,900,end=9998) zeil
         call qctolo(zeil)
         if(index(zeil,'eoi - ').le.0) goto 9998
!
!---knotentransformationskarte
      else if(index(zeil,'knotentransf.-karte').gt.0) then
         read(iunit,930,end=9998) (ktrf23(i),i=1,npof2d)
         read(iunit,900,end=9998) zeil
         call qctolo(zeil)
         if(index(zeil,'eoi - ').le.0) goto 9998
      else
      goto 10
      end if
      goto 10
!
!---das feld necf2d(nelf2d) besetzen
 200  do 30 iel=1,nelf2d
         necf2d(iel)=ncorn(knof2d(1,iel),mxkn2d)
   30 continue
!
!---netzkoordinaten dimensionslos machen
      do 40 i=1,npof2d
         phi=detang(xyzf3d(i,2),xyzf3d(i,1) )
         xyzf2d(i,1)=phi
         if (ndmf2d.eq.3) xyzf2d(i,2)=xyzf3d(i,3)*z2/diaf2d
 40   continue
!
      do 50 iel = 1, nelf2d
         kno=knof2d(1,iel)
         xmax=-bignum
         xmin= bignum
         do 60 iecke=1,necf2d(iel)
            kno=knof2d(iecke,iel)
            xmax=max( xmax,xyzf2d(kno,1) )
            xmin=min( xmin,xyzf2d(kno,1) )
            do 70 idim=1,ndmf2d-1
               xyzglo(iecke,idim) =  xyzf2d(kno,idim)
 70         continue
 60      continue
!
!---kennung fuer elemente die ueber den ursprung gehen
        if(abs(xmax-xmin).gt.pi) then
           ielc2d(iel)=1
           do 80 iecke=1,necf2d(iel)
              if(xyzglo(iecke,1).lt.pi) xyzglo(iecke,1)
     &                                = xyzglo(iecke,1) + z2*pi
 80        continue
        end if
        call elchck(knof2d(1,iel),xyzglo,ityf2d(iel),mxkn2d,ierr)
 50   continue
!
      call qfclos(iunit,0)
      return
!
 900  format (a)
 810  format (bn,3(2x,2d12.0))
 910  format (bn,2x,3f12.0,4x,3f12.0)
 920  format (bn,2x,i3,2x,12i6)
 930  format (bn,9(2x,i6))
 1000 stop '****error:reaint****'
 9998 stop 'data mismatch reaint'
      end
!---->---1---------2---------3---------4---------5---------6---------7--<

!@modul-flaanz
!@======================================================================
! Test wie viele Flaechen in der Alfa-Teta-U-Karte stehen, um dann
! die Felder entsprechend dimensionieren zu koennen.
!
! Eingabeparameter      Beschreibung
!  fname                 Dateiname des EITED-Inputfiles
!
! Ausgabeparameter      Beschreibung
!  nflae                 Anzahl gefundener Flaechen
!=======================================================================
      subroutine flaanz(nflae,fname)
      implicit double precision(a-h,o-z)
      character*80 zeil,fname*(*)
      logical lalfa
      integer qfreeu

      nflae=0
      lalfa=.false.
      iunit=qfreeu()
      call qfopse(iunit,fname,'ol','fo',ierror)
      if (ierror.ne.0) then
         write(*,8001) fname
         stop 'flaanz: error '
      endif

!.....Anfang des Lese-Loops

 100  continue
      read(iunit,8002,end=1000) zeil(1:80)
      call qctolo(zeil)

!.....Alfa-Teta-U-Karte gefunden ?

      if (index(zeil,'boi - alfa-teta-u-karte').gt.0) then
         lalfa=.true.
         goto 100
      endif

!.....Ende der Alfa-Teta-U-Karte gefunden ?

      if (index(zeil,'eoi - alfa-teta-u-karte').gt.0) goto 1000
      if (index(zeil,'eoi - inputfile').gt.0) goto 1000

!.....Wenn File-Pointer innerhalb der Alfa-Teta-U-Karte, dann lesen

      if (lalfa) then
         if (index(zeil,'flaechen:').gt.0) then
 200        read(iunit,8002) zeil(1:80)
            call qctolo(zeil)
            if (zeil.eq.' ') goto 200
            if (index(zeil,'e').gt.0) goto 100
            if (index(zeil,':').gt.0) nflae=nflae+1
            goto 200
         endif
      endif
      goto 100

 1000 call qfclos(iunit,0)
      return
 8001 format(1x,'error while opening the thd-file ',a)
 8002 format(a)
      end
!@modul-raltuk
!@======================================================================
! Alfa-Teta-U-Karte einlesen
!
! Eingabeparameter      Beschreibung
!  iunit                 File-Pointer
!  ityp
!  knoel
!  arand(nflae)          Alfa-Werte, arand(1)=Default-Alfa
!  trand(nflae)          Teta-U-Werte, trand(1)=Default-Teta-U
!  nflae                 Anzahl der Flaechen in der Alfa-Teta-U-Karte
!  maxkno
!
! Ausgabeparameter      Beschreibung
!  arand                 Alfa-Werte
!  trand                 Teta-U-Werte
!=======================================================================
      subroutine raltuk(filnam,ityp,knoel,necke,
     &     arand,trand,nflae,maxkno)
      implicit double precision(a-h,o-z)
      character*80 zeil,filnam
      logical lende
      integer qfreeu
      dimension knoel(maxkno,nflae),trand(nflae),
     *          arand(nflae),
     *          ityp(nflae),necke(nflae)

      iunit=qfreeu()
      call qfopse(iunit,filnam,'ol','fo',ierror)
      if (ierror.ne.0) then
         write(*,8001) fname
         stop 'raltuk: error'
      endif
!
!---vorspulen ...
 10   read(iunit,'(a)',end=9991) zeil
      call qctolo(zeil)
      if(index(zeil,'boi - alfa-teta-u-karte').gt.0) goto 20
      goto 10
 20   continue
!
      call iset(maxkno*nflae,0,knoel,1)

!.....Ende der Alfa-Teta-U-Karte gefunden --> lende=.true.

      lende=.false.

!.....Die Default-Werte stehen jeweils im ersten Feldelement 

      alfad=arand(1)
      tetaud=trand(1)

!.....Index fuer die Randflaechenfelder <arand> und <trand> nullen

      iflae=0
 100  continue
      if (lende) goto 8000

!.....Anfang des aktuellen Bereichs aus der Alfa-Teta-U-Karte

      istart=iflae+1

!.....Default-Werte fuer Alfa und Teta-U fuer den aktuellen Bereich

      al=alfad
      tu=tetaud
 200  continue
      read(iunit,9000,end=9991) zeil(1:80)
      call qctolo(zeil)
 201  continue
      if (zeil.eq.' ') goto 200
      if (index(zeil,'eoi - alfa-teta-u-karte').gt.0) lende=.true.

!.....Neuer Bereich gefunden

      if (index(zeil,'neuer bereich').gt.0 .or. lende) then
         do 300 i=istart,iflae
            arand(i)=al
            trand(i)=tu
 300     continue         
         goto 100
      endif
      if (index(zeil,'default - teta').gt.0) then
         read(zeil,9010) tu
      else if (index(zeil,'default - alfa').gt.0) then
         read(zeil,9020) al
      else if (index(zeil,'flaechen').gt.0) then
 800     continue
         read(iunit,9000,end=9991) zeil(1:80)
         call qctolo(zeil)
         if (zeil.eq.' ') goto 800
         if (index(zeil,'e').gt.0) goto 201
         iflae=iflae+1
         read(zeil,9030) ityp(iflae),(knoel(j,iflae),j=1,maxkno)
         necke(iflae)=ncorn(knoel(1,iflae),maxkno)
         goto 800
      endif
      goto 200
 8000 call qfclos(iunit,0)
      return
 8001 format(1x,'error while opening the thd-file ',a)
 9000 format(a)
 9010 format(bn,20x,e12.0)
 9020 format(bn,18x,e12.0)
 9030 format(bn,2x,i3,2x,12i6)
 9991 stop 'raltuk: corrupt alfa-teta-u-card'
      end

