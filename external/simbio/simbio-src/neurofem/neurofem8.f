!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem8.f :      	 calculation of the power, used in each 
!                                compartment
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

      subroutine powdiv(knogeo,itygeo,necgeo,xyzgeo,pergeo,vecsol)
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo)
      dimension xyzgeo(npogeo,ndmgeo),pergeo(numper,naniso)
      dimension vecsol(npogeo)
!
      call dset(maxlei,z0,volule,1)
      call dset(maxlei,z0,powele,1)
!
      cpuvol=qscput(9,0,ierr)
      call percen(0,nelgeo)
      powtot=z0
      voltot=z0
      idone=0
      do 20 j=1,numlei
         do 10 iel=ianfle(j),iendle(j)
            idone=idone+1
            call elepow(knogeo,itygeo,necgeo,xyzgeo,                    &
     &           pergeo,vecsol,vollok,powlok,iel )
            volule(j)=volule(j)+vollok
            powele(j)=powele(j)+powlok
            call percen(idone,nelgeo)
   10    continue
         voltot=voltot+volule(j)
         powtot=powtot+powele(j)
   20 continue
      write(*,900) qscput(9,1,ierr)
!
      write(*,920) voltot,powtot
      do 30 i=1,numlei
         write(*,910) i,volule(i),powele(i),                            &
     &   volule(i)/voltot,powele(i)/powtot
   30 continue
!
 1000 return
  900 format (' CPU--time for current-compartment search: ',g12.5)
  910 format (' comp. ',i2,' vol ',e12.5,' pow ',e12.5,                 &
     &        ' v-r ',f12.5,' p-r ',f12.5)
  920 format (' Volume-total ',g12.5,' Power total ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elepow(knogeo,itygeo,necgeo,xyzgeo,                    &
     &                  pergeo,vecsol,vollok,powlok,iel )
!
!---calculation of element matrices
      include 'neurofem.inc'
      dimension knogeo(mxkn3d,nelgeo),itygeo(nelgeo),necgeo(nelgeo),    &
     &          loknod(32)
      dimension pergeo(numper,naniso),xyzgeo(npogeo,ndmgeo)
      dimension vecsol(npogeo)
      dimension perlok(mxkn3d,6),potlok(mxkn3d),pergau(6)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn)
      dimension dum(3),dum2(3)
!
      vollok=z0
      powlok=z0
!
!---turn to local node numbers
      call iset(32,0,loknod(1),1)
      do 10 i=1,necgeo(iel)
         kno=knogeo(i,iel)
         loknod(i)=kno
         do 12 idim=1,ndmgeo
            xlk(idim,i) = xyzgeo(kno,idim)
   12    continue
         potlok(i) = vecsol(kno)
!
!---element related properties
         idum=iel
!
         do 15 j=1,naniso
            perlok(i,j) = pergeo(idum,j)
   15    continue
!
   10 continue
!
!---determine number of gaussian points
      call itpanz(itygeo(iel),necgeo(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'elepow 1'
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
         call itpdat(itygeo(iel),intgrd,iint,dinpkt,ierr)
         if (ierr.ne.0) stop 'elepow 2'
         call fofagn(itygeo(iel),dinpkt(1),dinpkt(2),dinpkt(3),         &
     &               loknod,xlk,fwe,fag,detj,ierr)
         if (ierr.ne.0) stop 'elepow 3'
!
!---determine values at gaussian points
         call dset(naniso,z0,pergau,1)
         call dset(3     ,z0,dum   ,1)
         do 40 i = 1, necgeo(iel)
            do 45 j=1,naniso
               pergau(j) = pergau(j) + perlok(i,j) * fwe(i)
   45       continue
            do 50 idim=1,ndmgeo
               dum(idim)=dum(idim)+potlok(i)*fag(idim,i)
   50       continue
   40    continue
!
!---dum  contains the local potentential gradient
!---dum2 gets the value of the local current flow
         if (laniso) then
         else
            dum2(1)=pergau(1)*dum(1)
            dum2(2)=pergau(1)*dum(2)
            dum2(3)=pergau(1)*dum(3)
         end if
         potgra = z0
         do 60 idim=1,ndmgeo
            potgra=potgra+dum2(idim)*dum(idim)
   60    continue
         xfak   =   detj * dinpkt(4)
         powlok = powlok + potgra*xfak
         vollok = vollok + xfak
!
   20 continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!
!
!
!
!
!
!
