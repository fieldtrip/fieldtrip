!  $2 04/06/2003  Anwander A.  removed check for equal area for usage in Jena
!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem12.f :      	 Subroutines related to the MEG Treatment
!                              -------------------
!     begin                : Mai 2000
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!    NeuroFEM license:
!    =================
!    Copyright (c) 2007 by 
!    Dr.Carsten Wolters, Dr.Alfred Anwander,  
!    Prof.Dr.H.Buchner, Prof.Dr.G.Knoll, Dr.Adrian Rienaecker, 
!    Rainer Beckmann, Robert Pohlmeier. 
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
      subroutine inimeg(filakt)
!     #################################
!     # initializes meg configuration #
!     #################################
      include 'neurofem.inc'
!
      iunit=qfreeu()
      call qfopse(iunit,filakt,'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,900) 'fehler filakt in inimeg'
         stop
      end if
!
      read (iunit,910) nummeg
      read (iunit,920) numtyp
      read (iunit,925) numsys
!
      melmeg=0
      do 10 i=1,numtyp
         read(iunit,930) num,ianf,iend,dum
         if ((num.lt.1).or.(num.gt.numtyp)) then
            write(*,940)
            stop 'in inimeg'
         end if
         if ((ianf.lt.1).or.(ianf.ge.iend)) then
            write(*,950)
            stop 'in inimeg'
         end if
         if (iend.gt.nelmeg) then
            write(*,950)
            stop 'in inimeg'
         end if
!---determine max elements=f(typ)
         num=iend-ianf+1
         if (num.gt.melmeg) melmeg=num
   10 continue
!
      if (melmeg.eq.0) stop 'error in inimeg'
      call qfclos(iunit,0)
!
      return
  900 format(a)
  910 format(24x,i5)
  920 format(24x,i5)
  925 format(24x,i5)
  930 format(24x,3i5,1x,f12.0)
  940 format('not allowed typ number ',i5)
  950 format('not allowed element number ',i5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reameg(filakt,megtyp,itypmg,factyp)
!     ###########################
!     # reads meg configuration #
!     ###########################
      include 'neurofem.inc'
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2),factyp(numtyp)
!
      iunit=qfreeu()
      call qfopse(iunit,filakt,'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,900) 'fehler filakt in reameg'
         stop
      end if
!
      read (iunit,910) nummeg
      read (iunit,920) numtyp
      read (iunit,925) numsys
!     
      call iset(numtyp*(mxloop+1),0,megtyp,1)
      do 10 i=1,numtyp
         read(iunit,930) num,ianf,iend,factyp(i)
         if ((num.lt.1).or.(num.gt.numtyp)) then
            write(*,940)
            stop 'in reameg'
         end if
         if ((ianf.lt.1).or.(ianf.ge.iend)) then
            write(*,950)
            stop 'in reameg'
         end if
         if (iend.gt.nelmeg) then
            write(*,950)
            stop 'in reameg'
         end if
         megtyp(num,1)=ianf
         megtyp(num,2)=iend
   10 continue
      do 40 i=1,numtyp
         if (megtyp(i,1).eq.0) then
            write(*,980) i
            stop 'in reameg'
         end if
   40 continue
!
      call iset(nummeg*2,0,itypmg,1)
      do 20 i=1,nummeg
         read(iunit,935) inum,ityp,isys
         if ((inum.lt.1).or.(inum.gt.nummeg)) then
            write(*,940) i
            stop 'in reameg'
         end if
         if ((ityp.lt.1).or.(ityp.gt.numtyp)) then
            write(*,960) i
            stop 'in reameg'
         end if
         if ((isys.lt.1).or.(isys.gt.numsys)) then
            write(*,980) i
            stop 'in reameg'
         end if
         itypmg(inum,1)=ityp
         itypmg(inum,2)=isys
   20 continue
      do 50 isys=1,numsys
         do 30 i=1,nummeg
            if (itypmg(i,isys).eq.0) then
               write(*,970) i
               stop 'in reameg'
            end if
   30    continue
   50 continue
!
      call qfclos(iunit,0)
!
      return
  900 format(a)
  910 format(24x,i5)
  920 format(24x,i5)
  925 format(24x,i5)
  930 format(24x,3i5,1x,f12.0)
  935 format(20x,3i5)
  940 format('not allowed typ number ',i5)
  950 format('not allowed element number ',i5)
  960 format('not allowed typ number ',i5)
  970 format('typ of meg-channel ',i5,' not defined')
  980 format('meg typ ',i5,' not defined')
  990 format('meg system ',i5,' not defined')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine gettra(itypmg,trahed,trameg)
!     #######################################
!     # determining transformation matrices #
!     #######################################
      include 'neurofem.inc'
      dimension itypmg(nummeg,2)
      dimension trahed(4,4,numsys),trameg(4,4,nummeg)
      dimension array3(4,4)
!      
!      call reatra(numsys,trahed,filnam(41))
!      call reatra(nummeg,trameg,filnam(37))
      do 10 inum=1,nummeg
!$$$         ityp=itypmg(inum,1)
         isys=itypmg(inum,2)
         call dgemm('N','N',4,4,4,z1,trahed(1,1,isys),4,                &
     &               trameg(1,1,inum),4,z0,array3,4)
         call dcopy(4*4,array3,1,trameg(1,1,inum),1)
   10 continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine reatra(numtra,tramat,filakt)
!     #################################
!     # reads transformation matrices #
!     #################################
      include 'neurofem.inc'
      dimension tramat(4,4,numtra)
!
      call dset(numtra*16,z0,tramat,1)
!
      iunit=qfreeu()
      call qfopse(iunit,filakt,'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,'(2a)') ' error in reatra because of file: ',         &
     &        filakt
         stop
      end if
!
      if (numtra.lt.1) then         
         write(*,940)
         stop
      end if
      do 10 i=1,numtra
         read (iunit,910) num
         if ((num.lt.1).or.(num.gt.numtra)) then
            write(*,930) num
            stop
         end if
         read(iunit,920) ((tramat(j,k,num),k=1,4),j=1,4)
   10 continue
!      
      call qfclos(iunit,0)
!
      return
  900 format(a)
  910 format(30x,i5)
  920 format(1x,4f12.0)
  930 format('illegal number of transformation matrix: ',i5)
  940 format('not allowed number of transformation matrices')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine gampre(knomeg,itymeg,necmeg,xyzmeg,intgam,fstgam,      &
     &           xyzmec,trameg,megtyp,itypmg,factyp)
!    ##################################################
!    # preprocessing of magneto-/gradiometer geometry #
!    ##################################################
      include 'neurofem.inc'
      dimension knomeg(mxkn1d,nelmeg),itymeg(nelmeg),necmeg(nelmeg),    &
     &          loknod(32)
      dimension xyzmeg(npomeg,ndmmeg),gamvec(3),xyzmec(npomeg,ndmmeg)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn),      &
     &          rjac(3,3)
      dimension intgam(melmeg,nummeg)
      dimension fstgam(mvagam,mgagam,melmeg,nummeg)
      dimension trameg(4,4,nummeg)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2),factyp(numtyp)
!
!---zero arrays
      call iset(melmeg*nummeg,0,intgam,1)
      call dset(mvagam*mgagam*melmeg*nummeg,z0,fstgam,1)
      call iset(32,0,loknod(1),1)
!
      do 600 iakmeg=1,nummeg
!
!---begin/end of magnetometer elements
         iaktyp=itypmg(iakmeg,1)
         ianf=megtyp(iaktyp,1)
         if (mxloop.gt.1) then
            iend=megtyp(iaktyp,3)
            if (iend.eq.0) iend=megtyp(iaktyp,2)
         else
            iend=megtyp(iaktyp,2)
         end if
!
!---factor to compute fT (SImm to SI!,1*fT=1d-15*Vs/m^2) 
         call flaemm(knomeg,itymeg,necmeg,xyzmeg,ianf,iend,flaech)
         facfla=ftconv*pi4mue*flconv/flaech*factyp(iaktyp)
         iend=megtyp(iaktyp,2)
!
!---mapping of geometrie
         call mapmeg(trameg,xyzmeg,xyzmec,iakmeg)
!
!---loop over all elements of actuell typ
         do 500 iel=ianf,iend
!
!---turn to local node numbers
            do 10 i=1,necmeg(iel)
               kno=knomeg(i,iel)
               loknod(i)=kno
               do 15 idim=1,ndmmeg
                  xlk(idim,i) = xyzmec(kno,idim)
   15          continue
   10       continue
!
!---determine number of gaussian points
            call itpanz(itymeg(iel),necmeg(iel),intgrd,ninpkt,ierr)
            if (ierr.ne.0) stop 'gampre 1'
            if (ninpkt.gt.mgagam) stop 'enlarge mgagam'
            intgam(iel-ianf+1,iakmeg)=ninpkt
!
!---loop over all gaussian points
            do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
               call itpdat(itymeg(iel),intgrd,iint,dinpkt,ierr)
               if (ierr.ne.0) stop 'gampre 2'
               call fofagn(itymeg(iel),dinpkt(1),dinpkt(2),dinpkt(3),   &
     &              loknod,xlk,fwe,fag,detj,ierr)
               if (ierr.ne.0) stop 'gampre 3'
!     
!---determine coordinates of gaussian point
               do 40 k=1,ndmmeg
                  gamvec(k)=z0
                  do 50 i=1,necmeg(iel)
                     gamvec(k)=gamvec(k)+xlk(k,i)*fwe(i)
   50             continue
   40          continue
!
!---rjac(1,1..3)= {dx/dxi, dy/dxi, dz/dxi}^T
               call getjac(rjac)
!
!---assign fstgam               
               d4=dinpkt(4)*facfla
               ifst=iel-ianf+1
               fstgam(1,iint,ifst,iakmeg)=rjac(1,1)*d4
               fstgam(2,iint,ifst,iakmeg)=rjac(1,2)*d4
               fstgam(3,iint,ifst,iakmeg)=rjac(1,3)*d4
               fstgam(4,iint,ifst,iakmeg)=gamvec(1)
               fstgam(5,iint,ifst,iakmeg)=gamvec(2)
               fstgam(6,iint,ifst,iakmeg)=gamvec(3)
!
   20       continue
            call iset(necmeg(iel),0,loknod(1),1)
  500       continue
  600    continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine flaemm(knomeg,itymeg,necmeg,xyzmeg,ianf,iend,flaech)
!     #####################################################
!     # determines the area of a loop in a z=const plane! #
!     #####################################################
      include 'neurofem.inc'
      dimension knomeg(mxkn1d,nelmeg),itymeg(nelmeg),necmeg(nelmeg),    &
     &          loknod(32)
      dimension xyzmeg(npomeg,ndmmeg)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn),      &
     &          rjac(3,3)
!      
!---zero arrays
      call iset(32,0,loknod(1),1)
      flaech=z0
!
!---loop over elements of the loop
      do 20 iel=ianf,iend
!
!---turn to local node numbers
         do 30 i=1,necmeg(iel)
            kno=knomeg(i,iel)
            loknod(i)=kno
            do 40 idim=1,ndmmeg
               xlk(idim,i) = xyzmeg(kno,idim)
   40       continue
   30    continue
!
!---determine number of gaussian points
         call itpanz(itymeg(iel),necmeg(iel),intgrd,ninpkt,ierr)
         if (ierr.ne.0) stop 'flaemm 1'
!
!---loop over all gaussian points
         do 50 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
            call itpdat(itymeg(iel),intgrd,iint,dinpkt,ierr)
            if (ierr.ne.0) stop 'flaemm 2'
            call fofagn(itymeg(iel),dinpkt(1),dinpkt(2),dinpkt(3),      &
     &                  loknod,xlk,fwe,fag,detj,ierr)
            if (ierr.ne.0) stop 'flaemm 3'
!
!---rjac(1,1..3)= {dx/dxi, dy/dxi, dz/dxi}^T
            call getjac(rjac)
!
            xkoord=z0
            do 60 i=1,necmeg(iel)
               xkoord=xkoord+xlk(1,i)*fwe(i)
   60       continue
!
!---rot{0,x,0}^T={0,0,1}^T 
            flaech=flaech+xkoord*rjac(1,2)*dinpkt(4)
!     
   50    continue
         call iset(necmeg(iel),0,loknod(1),1)
   20 continue
!     
      if (flaech.lt.z0) flaech=-flaech
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine mapmeg(trameg,xyzmeg,xyzmec,iakmeg)
!     ############################################
!     # mapping of flux-sensor to given position #
!     ############################################
      include 'neurofem.inc'
      dimension trameg(4,4,nummeg),                                     &
     &          xyzmeg(npomeg,ndmmeg),xyzmec(npomeg,ndmmeg)
      dimension xvek(4)
!
!---transformation of coordinates
      do 10 i=1,npomeg
         do 20 j=1,3
            xvek(j)=xyzmeg(i,j)
   20    continue
         xvek(4)=z1
         do 30 j=1,3
            dum=z0
            do 40 k=1,4
               dum=dum+trameg(j,k,iakmeg)*xvek(k)
   40       continue
            xyzmec(i,j)=dum
   30    continue
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine intmeg(knogeo,itygeo,necgeo,sekmeg,xyzgeo,pergeo,      &
     &                  intgam,fstgam,itypmg,megtyp)
!     #################################################################
!     # determines integration matrix for secondary flux contribution #
!     #################################################################
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo)
      dimension sekmeg(npogeo,nummeg),xyzgeo(npogeo,ndmgeo)
      dimension pergeo(numper,naniso)
      dimension intgam(melmeg,nummeg)
      dimension fstgam(mvagam,mgagam,melmeg,nummeg)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2)
!
      write(*,910)
!
!---read secondary flux matrix, if present and not (lead field matrix)
      if (qfilda(filnam(40)).and.(.not.linmat)) then
         write(*,940) filnam(40)(1:qclen(filnam(40)))
         iunit=qfreeu()
         call qfopse(iunit,filnam(40),'un','un',ierr)
         if (ierr.gt.0) then
            write (*,950) 'error filnam(40) in intmeg'
            stop
         end if
         read(iunit) num1
         read(iunit) num2
         if (num1.eq.nummeg.and.num2.eq.npogeo) then
            read(iunit)((sekmeg(i,j),i=1,npogeo),j=1,nummeg)
            call qfclos(iunit,0)
            return
         else
            call qfclos(iunit,0) 
         end if
      end if
!
!---compile secondary meg matrix
      call dset(npogeo*nummeg,z0,sekmeg,1)
      cpuvol=qscput(9,0,ierr)
      call percen(0,nelgeo)
      do 10 iel=1,nelgeo
         call sekele(knogeo,itygeo,necgeo,xyzgeo,pergeo,                &
     &               intgam,fstgam,itypmg,megtyp,                       &
     &               sekmeg,iel)
         call percen(iel,nelgeo)
   10 continue
      write(*,900) qscput(9,1,ierr)
      write(*,920)
!---rauskommentieren unterbindet das Schreiben der Sekundaermatrix
      write(*,930) filnam(40)(1:qclen(filnam(40)))
      iunit=qfreeu()
      call qfopse(iunit,filnam(40),'un','un',ierr)
      if (ierr.gt.0) then
         write (*,950) 'error filnam(40) in intmeg'
         stop
      end if
      write(iunit) nummeg
      write(iunit) npogeo
      write(iunit)((sekmeg(i,j),i=1,npogeo),j=1,nummeg)

      call qfclos(iunit,0)
!
      return
  900 format (' CPU--time for secondary meg integration: ',g12.5)
  910 format (' Preparing magnetic secondary flux integration matrix')
  920 format (' Finished  magnetic secondary flux integration matrix')
  930 format (' Storing secondary flux matrix   on file ',a)
  940 format (' Reading secondary flux matrix from file ',a)
  950 format (a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine sekele(knogeo,itygeo,necgeo,xyzgeo,pergeo,             &
     &                  intgam,fstgam,itypmg,megtyp,                    &
     &                  sekmeg,iel)
!     ######################################
!     # contribution of one volume element #
!     ######################################
!
!---calculation of element matrices
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo),    &
     &          loknod(32)
      dimension pergeo(numper,naniso),xyzgeo(npogeo,ndmgeo)
      dimension perlok(mxkn3d,6),pergau(6)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn),      &
     &          volvec(3),valint(3),fagcpy(3)
      dimension intgam(melmeg,nummeg)
      dimension fstgam(mvagam,mgagam,melmeg,nummeg)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2)
      dimension sekmeg(npogeo,nummeg)
!
!---turn to local node numbers
      call iset(32,0,loknod(1),1)
      do 10 i=1,necgeo(iel)
         kno       = knogeo(i,iel)
         xlk(1,i)  = xyzgeo(kno,1)
         xlk(2,i)  = xyzgeo(kno,2)
         xlk(3,i)  = xyzgeo(kno,3)
         loknod(i) = kno
!
!---element related properties
         idum=iel
!
         do 15 j=1,naniso
            perlok(i,j) = pergeo(idum,j)
   15    continue
   10 continue
!
!---determine number of gaussian points
      call itpanz(itygeo(iel),necgeo(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'sekele 1'
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
         call itpdat(itygeo(iel),intgrd,iint,dinpkt,ierr)
         if (ierr.ne.0) stop 'sekele 2'
         call fofagn(itygeo(iel),dinpkt(1),dinpkt(2),dinpkt(3),         &
     &               loknod,xlk,fwe,fag,detj,ierr)
         if (ierr.ne.0) stop 'sekele 3'
!
!---determine values at gaussian points
         volvec(1) = z0
         volvec(2) = z0
         volvec(3) = z0
         call dset(naniso,z0,pergau,1)
         do 40 i = 1, necgeo(iel)
            volvec(1) = volvec(1) + xlk(1,i) * fwe(i)
            volvec(2) = volvec(2) + xlk(2,i) * fwe(i)
            volvec(3) = volvec(3) + xlk(3,i) * fwe(i)
            do 50 j=1,naniso
               pergau(j) = pergau(j) + perlok(i,j) * fwe(i)
   50       continue
   40    continue
         xfak = detj * dinpkt(4)
!
!---loop over meg-sensors
         do 100 iakmeg=1,nummeg            
!
!---integration over meg loop
            call gamint(intgam,fstgam,volvec,itypmg,                    &
     &                  megtyp,valint,iakmeg,xfak)
!
!---loop over all nodal points and summation of contributions
            do 30 i=1,necgeo(iel)
!
!---multiply with electric conductivity
               if (laniso) then
                  fagcpy(1) = pergau(1)*fag(1,i)                        &
     &                      + pergau(4)*fag(2,i)+pergau(6)*fag(3,i)
                  fagcpy(2) = pergau(2)*fag(2,i)                        &
     &                      + pergau(5)*fag(3,i)+pergau(4)*fag(1,i) 
                  fagcpy(3) = pergau(3)*fag(3,i)                        &
     &                      + pergau(6)*fag(1,i)+pergau(5)*fag(2,i) 
               else
                  fagcpy(1) = pergau(1)*fag(1,i)
                  fagcpy(2) = pergau(1)*fag(2,i)
                  fagcpy(3) = pergau(1)*fag(3,i)
               end if
               dum = fagcpy(1)*valint(1)                                &
     &             + fagcpy(2)*valint(2)                                &
     &             + fagcpy(3)*valint(3)
               iglo=knogeo(i,iel)
               if (iglo.gt.0) then
                  sekmeg(iglo,iakmeg)=sekmeg(iglo,iakmeg)+dum
               end if
   30       continue
  100    continue
   20 continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine gamint(intgam,fstgam,volvec,itypmg,                    &
     &                  megtyp,valint,iakmeg,perfak)
!     ##################################
!     # integration over a flux-sensor #
!     ##################################
      include 'neurofem.inc'
      dimension intgam(melmeg,nummeg),itypmg(nummeg,2),                 &
     &          megtyp(numtyp,mxloop+1)
      dimension fstgam(mvagam,mgagam,melmeg,nummeg)
      dimension volvec(3),valint(3)
!
      valint(1)=z0
      valint(2)=z0
      valint(3)=z0
      iaktyp=itypmg(iakmeg,1)
      iend=megtyp(iaktyp,2)-megtyp(iaktyp,1)+1      
      do 20 iel=1,iend
         do 30 iint=1,intgam(iel,iakmeg)
            dif=fstgam(4,iint,iel,iakmeg)-volvec(1)
            dt=dif*dif
            dif=fstgam(5,iint,iel,iakmeg)-volvec(2)
            dt=dt+dif*dif
            dif=fstgam(6,iint,iel,iakmeg)-volvec(3)
            dt=dt+dif*dif
!---dt is now squared distance between gaussian and integration point
            dt=perfak/sqrt(dt)
            valint(1)=valint(1)+dt*fstgam(1,iint,iel,iakmeg)
            valint(2)=valint(2)+dt*fstgam(2,iint,iel,iakmeg)
            valint(3)=valint(3)+dt*fstgam(3,iint,iel,iakmeg)
   30    continue
   20 continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine primeg(intgam,fstgam,itypmg,megtyp,                    &
     &                  xyzflf,dipole,potmeg,iknflf)
!     #######################################
!     # determines primary flux contibution #
!     #######################################
      include 'neurofem.inc'
      dimension intgam(melmeg,nummeg)
      dimension fstgam(mvagam,mgagam,melmeg,nummeg)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2)
      dimension xyzflf(npoflf,ndmflf),dipole(npoflf,ndmflf),            &
     &          potmeg(nummeg)
      dimension volvec(3),valint(3)
!
      do 10 k=1,3
         volvec(k)=xyzflf(iknflf,k)
   10 continue
      do 20 iakmeg=1,nummeg
         call gamint(intgam,fstgam,volvec,itypmg,                       &
     &               megtyp,valint,iakmeg,z1)
         scalar=z0
         do 30 k=1,3
            scalar=scalar+valint(k)*dipole(iknflf,k)
   30    continue
!
!---primary magnetic flux
         potmeg(iakmeg)=potmeg(iakmeg)+scalar
   20 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine megana(xyzmeg,xyzmec,trameg,knomeg,itypmg,megtyp,      &
     &                  dipole,mrkdip,xyzflf,anameg,potmeg,factyp)
!###############################################################
!# analytical reference solution, use sensors with small area! #
!###############################################################
      include 'neurofem.inc'
      dimension xyzmeg(npomeg,ndmmeg),xyzmec(npomeg,ndmmeg)
      dimension mrkdip(npoflf,ndmflf),dipole(npoflf,ndmflf)
      dimension xyzflf(npoflf,ndmflf)
      dimension anameg(nummeg),potmeg(nummeg)
      dimension trameg(4,4,nummeg)
      dimension knomeg(mxkn1d,nelmeg)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2),factyp(numtyp)
      dimension dipoll(3),ortdip(3),ortflu(3),fluxde(3)
      dimension veknor(3)
!
!---anlytical solution only possible if ndmflf=3
      if (ndmflf.ne.3) then
         write(*,'(a)') 'meg analytical solution not possible!'
         write(*,'(a)') '      ndmflf=3 is necessary!         '
         return
      end if
!
!---loop over all sensors
      do 10 iakmeg=1,nummeg
         call mapmeg(trameg,xyzmeg,xyzmec,iakmeg)
         anameg(iakmeg)=z0
         iaktyp=itypmg(iakmeg,1)
!     
!---find  number of loops
         nloop=1
         do 20 iloop=2,mxloop
            if (megtyp(iaktyp,iloop+1).gt.0) nloop=nloop+1
   20    continue
!
!---loop over loops
         ianf=megtyp(iaktyp,1)
         do 30 iloop=1,nloop
            if (iloop.eq.nloop) then
               iend=megtyp(iaktyp,2)
            else
               iend=megtyp(iaktyp,iloop+1)            
            end if
            call geoana(xyzmec,knomeg,trameg,factyp,                    &
     &                  ortflu,veknor,factor,                           &
     &                  ianf,iend,iakmeg,iaktyp)
!
!---loop over dipoles
            do 40 idip=1,npoflf
               mrkdum=0
               do 50 j=1,3
                  mrkdum=mrkdum+mrkdip(idip,j)
   50          continue
               if (mrkdum.gt.0) then
                  do 60 j=1,3
                     dipoll(j)=dipole(idip,j)
                     ortdip(j)=xyzflf(idip,j)
   60             continue
!
!---flux through one loop
                  call dset(3,z0,fluxde,1)
                  call fluana(dipoll,ortdip,ortflu,fluxde)
                  du2=z0
                  do 80 j=1,3
                     du2=du2+fluxde(j)*veknor(j)*factor
   80             continue
                  anameg(iakmeg)=anameg(iakmeg)+du2
               end if
   40       continue
!     
            ianf=iend+1
   30    continue
   10 continue
!
!---comparison of solutions
      write(*,900)
      write(*,910)
      euknrm=z0
      anaana=z0
      anapot=z0
      potpot=z0
      do 300 i=1,nummeg
         dif=anameg(i)-potmeg(i)
         anaana=anaana+anameg(i)*anameg(i)
         anapot=anapot+anameg(i)*potmeg(i)
         potpot=potpot+potmeg(i)*potmeg(i)
         write(*,920) i,anameg(i),potmeg(i),dif
         euknrm=euknrm+dif*dif
  300 continue
      euknrm=sqrt(euknrm)
      write(*,930) euknrm
      corcof=sqrt(anaana)*sqrt(potpot)
      if (corcof.ne.z0) then
         corcof=anapot/corcof
         write(*,940) corcof
      else
         write(*,950)
      end if
      return
  900 format(' Comparing MEG sensor fluxes [fT]!')
  910 format('      #   analytical   numerical    difference')
  920 format(1x,i6,1x,e12.5,1x,e12.5,1x,e12.5)
  930 format(' MEG (ana-num) Euclidean norm is: ',e12.5)
  940 format(' MEG (ana-num)  Correlation   is: ',e12.5)
  950 format(' MEG (ana-num)  Correlation   is:  not defined')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine fluana(dipoll,ortdip,ortflu,fluxde)
!     #################################################
!     # determines magnetic flux density analytically #
!     #    outside a spherical symmetric conductor    #
!     #################################################
      include 'neurofem.inc'
      dimension dipoll(3),ortdip(3),ortflu(3),fluxde(3)
      dimension prokrz(3),disvek(3)
!
!---formula derivated by SARVAS (Phys. Med. Biol., Vol 32, No 1, 11-22)
!
!---dipole cross location of dipole
      prokrz(1)=dipoll(2)*ortdip(3)-dipoll(3)*ortdip(2)
      prokrz(2)=dipoll(3)*ortdip(1)-dipoll(1)*ortdip(3)
      prokrz(3)=dipoll(1)*ortdip(2)-dipoll(2)*ortdip(1)
!
!---vector from dipole to point of interest
      do 10 i=1,3
         disvek(i)=ortflu(i)-ortdip(i)
   10 continue
!
!---some scalar products
      ortfl2=z0
      disve2=z0
      proskm=z0
      proskr=z0
      prospa=z0
      do 20 i=1,3
         ortfl2=ortfl2+ortflu(i)*ortflu(i)
         disve2=disve2+disvek(i)*disvek(i)
         proskm=proskm+ortflu(i)*ortdip(i)
         proskr=proskr+ortflu(i)*disvek(i)
         prospa=prospa+prokrz(i)*ortflu(i)
   20 continue
!
!---distances
      ortfl1=sqrt(ortfl2)
      disve1=sqrt(disve2)
      compdi=proskr/disve1
!
!---terms for volume current
      divide=(ortfl2+ortfl1*disve1-proskm)*disve1
      divide=z1/divide
      fact1=z2*ortfl1+z2*disve1+disve2/ortfl1+compdi
      fact2=z2*ortfl1+disve1+compdi
      fact1=fact1*divide*prospa
      fact2=fact2*divide*prospa
!
!---factor to compute fT (SImm to SI!, fT=(10**-9)*uVs/m^2)
      facfla=pi4mue*divide*ftconv*flconv
!
!---flux density
      do 30 i=1,3
         fluxde(i)=facfla*(prokrz(i)-fact1*ortflu(i)+fact2*ortdip(i))
   30 continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine premeg(xyzmeg,knomeg,itymeg,necmeg,megtyp)
!     ###############################################
!     # preprocessing and checking of meg geometrie #
!     ###############################################
      include 'neurofem.inc'
      dimension xyzmeg(npomeg,ndmmeg)
      dimension knomeg(mxkn1d,nelmeg),itymeg(nelmeg),necmeg(nelmeg)
      dimension megtyp(numtyp,mxloop+1)
!
      write(*,900) ' Checking and sorting meg-geometrie '
!
!---checking type of elements, 301 is the code of barelements
      do 100 iaktyp=1,numtyp
         ielanf=megtyp(iaktyp,1)
         ielend=megtyp(iaktyp,2)
         do 200 iel=ielanf,ielend
            if (itymeg(iel).ne.301) then
               write(*,900)                                             &
     &      ' Magneto-/gradiometer elements have to be bar elements!'
               stop 'in srtmeg'
            end if
  200    continue
  100 continue
!
!---sorting the elements and checking the loops
      do 10 iaktyp=1,numtyp
         ielanf=megtyp(iaktyp,1)
         ielend=megtyp(iaktyp,2)
         do 20 iloop=1,mxloop
            if (ielanf.lt.ielend) then
               call findlp(knomeg,necmeg,ielanf,ielend,iendlp)
               if (iendlp.lt.ielend) then
                  if (iloop.lt.mxloop) then
                     megtyp(iaktyp,iloop+2)=iendlp
                  else
                     write(*,900)                                       &
     &               ' Increase allowed number of loops !'
                     write(*,900) ' parameter mxloop in > cauchy.inc < '
                     stop 'in srtmeg'
                  end if                  
               end if
               ielanf=iendlp+1
            else if (ielanf.eq.ielend) then
               write(*,900)                                             &
     &         ' A loop has to consist of at least two elements!'
               stop 'in srtmeg'
            else if (iloop.lt.mxloop) then
               megtyp(iaktyp,iloop+2)=0
            end if
   20    continue
   10 continue
!
!---all loops are restricted to have constant z-value 
      do 110 iaktyp=1,numtyp
         nloops=0
         do 120 iloop=1,mxloop
            if (megtyp(iaktyp,iloop+1).gt.0) nloops=nloops+1
  120    continue
         ianf=megtyp(iaktyp,1)
         do 130 iloop=1,nloops
            if (iloop.eq.nloops) then
               iend=megtyp(iaktyp,2)
            else
               iend=megtyp(iaktyp,iloop+2)            
            end if
            n=0
            sum1=z0
            sum2=z0
            do 140 iel=ianf,iend
               do 150 ikno=1,necmeg(iel)
                  kno=knomeg(ikno,iel)
                  zval=xyzmeg(kno,3)
                  n=n+1
                  sum1=sum1+zval
                  sum2=sum2+zval*zval
  150          continue
  140       continue
            test=abs(dble(n)*sum2-sum1*sum1)
            if (test.gt.tol) then
               write(*,900)                                             &
     &              ' The z-value of each loop has to be constant!'
               stop 'in premeg'
            end if
            ianf=iend+1
  130    continue
  110 continue
!
!---the loops of each typ have to have equal area
      do 210 iaktyp=1,numtyp
         n=1
         m=0
         fl1=z0
         fl2=z0
         do 220 iloop=2,mxloop
            if (megtyp(iaktyp,iloop+1).gt.0) n=n+1
  220    continue
         ianf=megtyp(iaktyp,1)
         do 230 iloop=1,n
            if (iloop.eq.n) then
               iend=megtyp(iaktyp,2)
            else
               iend=megtyp(iaktyp,iloop+1)            
            end if
            call flaemm(knomeg,itymeg,necmeg,xyzmeg,ianf,iend,flaech)
            m=m+1
            fl1=fl1+flaech
            fl2=fl2+flaech*flaech
            ianf=iend+1
  230    continue
         test=abs(dble(m)*fl2-fl1*fl1)
         if (test.gt.tol) then
            write(*,900)                                                &
     &    ' The loops of each sensor typ have to have equal area!'
!            stop 'in premeg'
         end if
  210 continue
!
      return
  900 format (a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine findlp(knomeg,necmeg,ielanf,ielend,iendlp)
!     ###############################
!     # finds loops of meg elements #
!     ###############################
      include 'neurofem.inc'
      dimension knomeg(mxkn1d,nelmeg),necmeg(nelmeg)
!
      knoanf=knomeg(1,ielanf)
      knoend=knomeg(2,ielanf)
      iendlp=ielanf+1
   30 iel=iendlp
!---find next element of loop
   40 if (knoend.ne.knomeg(1,iel)) then
         iel=iel+1
         if (iel.gt.ielend) then
            write(*,900) ' A loop has to be closed!'
            stop 'in findlp'
         end if
         goto 40
      end if
!---next element found, new start node, swapp elements
      knoend=knomeg(2,iel)
      call swapel(knomeg,necmeg,iel,iendlp)
!---if loop not yet closed
      if (knoend.ne.knoanf) then
         iendlp=iendlp+1
         if (iendlp.gt.ielend) then
            write(*,900) ' A loop has to be closed!'
            stop 'in findlp'
         end if
         goto 30
      end if
!
      return
  900 format(a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine swapel(knomeg,necmeg,i1,i2)
!     ###########################
!     # swapps two meg elements #
!     ###########################
      include 'neurofem.inc'
      dimension knomeg(mxkn1d,nelmeg),necmeg(nelmeg)
!
      if (i1.ne.i2) then
         idum=necmeg(i1)
         necmeg(i1)=necmeg(2)
         necmeg(i2)=idum
         do 50 j=1,mxkn1d
            idum=knomeg(j,i1)
            knomeg(j,i1)=knomeg(j,i2)
            knomeg(j,i2)=idum
   50    continue
      end if
!
      return
  900 format(a)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine geoana(xyzmec,knomeg,trameg,factyp,                    &
     &                  ortflu,veknor,factor,ianf,iend,iakmeg,iaktyp)
!     ###################################################
!     # determines midpoint and orientation of the loop #
!     ###################################################
      include 'neurofem.inc'
      dimension xyzmec(npomeg,ndmmeg),trameg(4,4,nummeg)
      dimension knomeg(mxkn1d,nelmeg),factyp(numtyp)
      dimension ortflu(3),prokrz(3),dum(3),veknor(3)  
!     
!---midpoint and normal vector of the loop
      call dset(3,z0,ortflu,1)
      do 10 i=ianf,iend
         kno=knomeg(1,i)
         do 20 j=1,3
            ortflu(j)=ortflu(j)+xyzmec(kno,j)
   20    continue
   10 continue
      dummy=z1/dble(iend-ianf+1)
      do 30 j=1,3
         ortflu(j)=ortflu(j)*dummy
         veknor(j)=trameg(j,3,iakmeg)
   30 continue
!
!---orientation of the loop
      ntest1=knomeg(1,ianf)
      ntest2=knomeg(2,ianf)
      do 40 j=1,3
         dum(j)=xyzmec(ntest2,j)-xyzmec(ntest1,j)
   40 continue
      prokrz(1)=dum(2)*veknor(3)-dum(3)*veknor(2)
      prokrz(2)=dum(3)*veknor(1)-dum(1)*veknor(3)
      prokrz(3)=dum(1)*veknor(2)-dum(2)*veknor(1)
      ska=z0
      do 50 j=1,3
         dummy=xyzmec(ianf,j)-ortflu(j)
         ska=ska+dummy*prokrz(j)
   50 continue
      if (ska.gt.z0) then
         factor=factyp(iaktyp)
      else
         factor=-factyp(iaktyp)
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
