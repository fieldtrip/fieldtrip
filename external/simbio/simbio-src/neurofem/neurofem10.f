!  $2 21/03/2003  Anwander A.  released for SIMBIO
!  $1 30/06/00aa subroutins which uses work.inc moved to cauchy.f

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem10.f :      	 Subroutines related to the eigenvalue problem
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
      subroutine reduce(ipogeo,indgeo,syscpy,valeig,veceig,             &
     &     orthog,dprodg,redinv,volmat,dkondg,ortsys)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo)
      dimension syscpy(lengeo),valeig(nvala,nustei),volmat(lengeo),     &
     &          veceig(npogeo,nvala),orthog(nmodes,2),dprodg(npogeo),   &
     &          redinv(nmodes),dkondg(npogeo),ortsys(nmodes,2)
!
      call dset(nmodes*2,z0,orthog,1)
!
!---sound check for the eigenvectors \phi^i_r i=1,nvala r=1,npogeo:
!   1) orthogonality of eigenvectors \phi^i_r with respect to 
!       a) the volume integration matrix (mass matrix, generalized eigenproblem)
!          \phi^i_r      m_{rs} \phi^j_s = \delta^{ij}
!       b) the kronecker    delta matrix (unit matrix, special     eigenproblem)
!          \phi^i_r \delta_{rs} \phi^j_s = \delta^{ij}
!   2) eigenvectors must form eigenvalues with stiffness matrix
!          \phi^i_r      k_{rs} \phi^j_s = \lambda^{ij}
!
      ortflt=z0
      write(*,920)
      call percen(0,nmodes)
      do 10 i=1,nmodes
         itru=i
         if (lconsi) then
            call matpro(volmat,veceig(1,itru),dprodg,                   &
     &           ipogeo,indgeo,lengeo,npogeo)
            scapro=ddot(npogeo,dprodg,1,veceig(1,itru),1)
            orthog(i,1)=orthog(i,1)+scapro
         else
            scapro=ddot(npogeo,veceig(1,itru),1,veceig(1,itru),1)
            orthog(i,1)=orthog(i,1)+scapro
         end if
         call matpro(syscpy,veceig(1,itru),dkondg,                      &
     &        ipogeo,indgeo,lengeo,npogeo)
         scapro=ddot(npogeo,dkondg,1,veceig(1,itru),1)
         ortsys(i,1)=ortsys(i,1)+scapro
         do 20 j=i+1,nmodes
            jtru=j
            if (lconsi) then
               scapro=ddot(npogeo,dprodg,1,veceig(1,jtru),1)
               orthog(i,2)=orthog(i,2)+z2*abs(scapro)
            else
               scapro=ddot(npogeo,veceig(1,itru),1,veceig(1,jtru),1)
               orthog(i,2)=orthog(i,2)+z2*abs(scapro)
            end if
            scapro=ddot(npogeo,dkondg,1,veceig(1,jtru),1)
            ortsys(i,2)=ortsys(i,2)+z2*abs(scapro)
   20    continue
         call percen(i,nmodes)
   10 continue
      scaort=z1/sqrt(dble(nmodes))
      ortflt=dnrm2(nmodes,orthog(1,2),1)*scaort
      ortdia=dnrm2(nmodes,orthog(1,1),1)*scaort
      write(*,930) ortdia
      write(*,940) ortflt
!
!---determine inverse of the reduced system matrix
      write(*,950)
      eigflt=z0
      ortflt=z0
      do 40 i=1,nmodes
         itru=i
         diff=ortsys(i,1)-valeig(i,1)
         eigflt=eigflt+diff*diff
         ortflt=ortflt+ortsys(i,2)*ortsys(i,2)
         redinv(i)=z1/valeig(i,1)
   40 continue
      eigflt=sqrt(eigflt)*scaort
      ortflt=sqrt(ortflt)*scaort
      write(*,960) eigflt
      write(*,970) ortflt
!
      return
  920 format(' Basic orthogonality check for the Eigenvectors')
  930 format(' Average  diagonal entry in ortho-matrix ',g12.5,' (?=1)')
  940 format(' Av.  non-diagonal entry in ortho-matrix ',g12.5,' (?=0)')
  950 format(' Basic validity check for the eigenvectors ',/,           &
     &       ' Eigenvectors must diagonalize the system matrix ')
  960 format(' sysmat:    eigenvalue-fault is ',e12.5)
  970 format(' sysmat: orthogonality-fault is ',e12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$23/07/03aa removed: undefined reference to `iovect_' `op_' `dnlaso_'
!
!  Implementation of the LASO2 Package  (c) D.S. Scott, 1984
!
!!---->---1---------2---------3---------4---------5---------6---------7--<
!	   subroutine eigval(indeig,valeig,veceig,wrkeig)
!	   include 'neurofem.inc'
!	   external op,iovect
!	   dimension indeig(nvala)
!	   dimension valeig(nvala,nustei),veceig(npogeo,nvala),wrkeig(nuwrk)
!!
!	   if (.not.liomem) then
!		  iunio=qfreeu()
!		  call qfopse(iunio,'lanczos.tmp','un','un',ierr)
!		  if (ierr.gt.0) stop 'lanczos.tmp'
!	   end if
!!
!	   nperm=0
!	   write(*,920)
!	   write(*,930) nvala,nblock,maxj,nperm
!	   write(*,935) nfig,nuwrk,maxop,liomem
!	   call dset(npogeo*nblock,z0,wrkeig,1)
!	   cpueig=qscput(4,0,ierr)
!	   call dnlaso(op, iovect,npogeo, nval, nfig, nperm,				 &
!	  &   nvala,valeig,npogeo,veceig, nblock, maxop, maxj,wrkeig,		 &
!	  &   indeig, ierr)
!	   cpueig=qscput(4,1,ierr)
!	   write(*,900) indeig(1)
!	   write(*,910) cpueig
!	   call decode(ierr)
!	   if (ierr.gt.0) then
!		  stop 'laso error'
!	   end if
!	   write(*,940) nperm
!	   nwriln=min(nperm,nvala)
!!
!	   if (.not.liomem) then
!		  call qfclos(iunio,-1)
!	   end if
!!
!	   return
!  900 format(' number of calls to OP ',i8)
!  910 format(' CPU-time for the eigenvalues: ',g12.5,' [s]')
!  920 format(' Parameters used for the eigenvalue analysis: ')
!  930 format(' nvala: ',i6,' nblock: ',i6,'  maxj: ',i6,' nperm: ',i6)
!  935 format('  nfig: ',i6,'  nuwrk: ',i6,' maxop: ',i6,' iomem: ',l1)
!  940 format(' LASO determined ',i6,' eigenvalues')
!	   end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine decode(ierr)
!
      write(*,800) ierr
      if (ierr.eq.0) then
         write(*,1000)
      else if (ierr.gt.0) then
         if(ierr/512.gt.0) write(*,900)
         ierr=mod(ierr,512)
         if(ierr/256.gt.0) write(*,910)
         ierr=mod(ierr,256)
         if(ierr/128.gt.0) write(*,920)
         ierr=mod(ierr,128)
         if(ierr/64.gt.0) write(*,930)
         ierr=mod(ierr,64)
         if(ierr/32.gt.0) write(*,940)
         ierr=mod(ierr,32)
         if(ierr/16.gt.0) write(*,950)
         ierr=mod(ierr,16)
         if(ierr/8.gt.0) write(*,960)
         ierr=mod(ierr,8)
         if(ierr/4.gt.0) write(*,970)
         ierr=mod(ierr,4)
         if(ierr/2.gt.0) write(*,980)
         ierr=mod(ierr,2)
         if(ierr/1.gt.0) write(*,990)
      else
         if(ierr.eq.-1) write(*,1010)
         if(ierr.eq.-2) write(*,1020)
         if(ierr.eq.-8) write(*,1030)
      end if
!
      return
  800 format(' dnlaso returned with error code ',i6)
  900 format(' nblock < 1')
  910 format(' abs(nval) > maxj/2')
  920 format(' abs(nval) > maxop')
  930 format(' abs(nval) > nvala')
  940 format(' abs(nval) < max(1, nperm)')
  950 format(' maxj < 6*nblock')
  960 format(' nperm < 0')
  970 format(' npogeo < n')
  980 format(' nfig < 0')
  990 format('n < 6*nblock')
 1000 format(' No error occurred, eigenvectors are OK')
 1010 format(' nperm > 0')
 1020 format(' more than maxop calls to op')
 1030 format(' disastrous loss of orthogonality! ',/,                   &
     &       ' check routine op or iovect',/,                           &
     &       ' or fiddle with nblock,maxj,maxop ')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$22/07/00aa copy of thease functions to cauchy.f to be complied with f90
!      subroutine op(n,m,p,q)
!      include 'neurofem.inc'
!      include 'work.inc'
!      character*8 chatot
!      dimension p(n,m),q(n,m)
!      save icouop,initia,ntotal
!      data icouop,initia,ntotal /0,0,0/
!
!      ntotal=ntotal+1
!      if (initia.eq.0) then
!         write(*,890)
!         write(*,900)
!         initia=1
!      end if
!      if (icouop.eq.0) then
!         call qdtext(' >',1)
!      end if
!      icouop=icouop+1
!      call dumyop(n,m,p,q,ip_ipogeo,ip_indgeo,dp_syscpy)
!      if (mod(icouop,2).eq.0.and.icouop.gt.0) call qdtext('.',1)
!      if (icouop.eq.100) then
!         write(chatot,'(i8)') ntotal
!         call qdtext('<',1)
!         call qdtext(chatot,1)
!         write(*,910)
!         icouop=0
!      end if
!
!      return
!  890 format(' Number of calls to routine OP indicated by dots')
!  900 format(' 0---10---20---30---40---50---60---70---80---90--100')
!  910 format(a)
!  920 format(' Total Number of calls to OP ',i8)
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dumyop(n,m,p,q,ipogeo,indgeo,syscpy)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo)
      dimension syscpy(lengeo)
      dimension p(n,m),q(n,m)
!
      do 10 i=1,m
         call matpro(syscpy,p(1,i),q(1,i),ipogeo,indgeo,lengeo,npogeo)
   10 continue
!
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$22/07/00aa copy of thease functions to cauchy.f to be complied with f90
!      subroutine iovect(n,m,q,j,k)
!      include 'neurofem.inc'
!      include 'work.inc'
!      dimension q(n,m)
!
!      if (liomem) then
!         call iovmem(n,m,q,j,k,dp_qeivec)
!      else
!         call iovdsk(n,m,q,j,k,iunio)
!      end if
!
!      return
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine iovmem(n,m,q,j,k,qvec)
      include 'neurofem.inc'
      dimension q(n,m),qvec(npogeo,maxj)
!
      if(k .eq. 1) go to 30
      do 20 l = 1,m
         l1 = j - m + l
         do 10 i = 1,n
            qvec(i,l1) = q(i,l)
   10    continue
   20 continue
      return
!
   30 do 50 l = 1,m
         l1 = j - m + l
         do 40 i = 1,n
            q(i,l) = qvec(i,l1)
   40    continue
   50 continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine iovdsk(n,m,q,j,k,iunio)
      integer n,m,j,k,i,l,iunio
      double precision q(n,m)
      if(j .eq. m) rewind iunio
      if(k .eq. 0) write(iunio) ((q(i,l), i=1,n), l=1,m)
      if(k .eq. 1)  read(iunio) ((q(i,l), i=1,n), l=1,m)
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$23/07/03aa removed this functions : undefined reference to `lanz_'
!
!     Implementation of the LANZ Package
!
!     LANZ is (c) M.T. Jones, M.L. Patrick, 1990
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!	   subroutine lnzcal(ipogeo,indgeo,lanrel,lanrow,lancol,			 &
!	  & 				 mascol,masrow,sysmat,sysmas,veclnz,			 &
!	  & 				 vallnz,eigess,syslan,volmat)
!	   include 'neurofem.inc'
!	   dimension ipogeo(npogeo),indgeo(lengeo),lanrel(lenval),  		 &
!	  & 		 lanrow(lenrow),lancol(lencol), 						 &
!	  & 		 mascol(lencol),masrow(lenrow)
!	   dimension veclnz(npogeo,lenval),vallnz(lenval,2),				 &
!	  & 		 eigess(npogeo),sysmat(lengeo),sysmas(lenlan),  		 &
!	  & 		 syslan(lenlan),volmat(lengeo)
!!
!	   call dset(npogeo*lenval,z0,veclnz,1)
!	   call dset(lenval*2,z0,vallnz,1)
!!
!	   call parlan(lanrel,veclnz,eigess)
!!
!	   call wandel(npogeo,lengeo,sysmat,indgeo,ipogeo,lancol(lengeo+1),  &
!	  & 	lanrow(npogeo+2),syslan,lancol(1),lanrow(1))
!	   if (lconsi) then
!		  call wandel(npogeo,lengeo,volmat,indgeo,ipogeo,				 &
!	  & 	   mascol(lengeo+1),masrow(npogeo+2),sysmas,				 &
!	  & 	   mascol(1),masrow(1))
!	   else
!		  call masmat(masrow,mascol,sysmas)
!	   end if
!	   call lanz(ldylan,ipar,rpar,vallnz(1,1),vallnz(1,2),  			 &
!	  & 	veclnz,lanrel,syslan,lanrow,lancol,sysmas,  				 &
!	  & 	masrow,mascol,eigess)
!!
!	   if (ipar(5).eq.0) then
!		  write(*,900)
!	   else
!		  write(*,910) ipar(5)
!	   end if
!	   write(*,920) ipar(4)
!	   write(*,930) ipar(6)
!	   nwriln=min(ipar(6),nvala)
!!
!	   call chklan(ipogeo,indgeo,lanrel,sysmat,volmat,vallnz,veclnz,	 &
!	  & 		  syslan,sysmas)
!!
!	   return
!  900 format(' LANZ calculation OK')
!  910 format(' LANZ reports error #',i8)
!  920 format(' LANZ took ',i8,' steps in this calculation')
!  930 format(' LANZ determined ',i8,' eigenvalues')
!	   end
!---->---1---------2---------3---------4---------5---------6---------7--<
	   subroutine parlan(lanrel,veclnz,eigess)
	   include 'neurofem.inc'
      dimension lanrel(lenval)
      dimension veclnz(npogeo,lenval),eigess(npogeo)
!
!---set parameters according to the description of LANZ
!---ldylan=leading dimension of the eigenvector array y
      ldylan=npogeo
!---order of the matrices K and M
      ipar( 1)=npogeo
!---maximum number of eigenpairs that can be stored (>nvala)
      ipar( 2)=lenval
!---number of desired eigenvalues
      ipar( 3)=nvala
!---maximum number of Lanczos steps (default=3*nvala)
      ipar( 4)=maxop
!---number of known eigenpairs
      ipar( 6)=0
      call dset(npogeo,z0,eigess,1)
!---debugging level (must be 0)
      ipar( 7)=0
!---type of problem: 0 K positive definite and M positive semi-definite
      ipar( 8)=0
!---inertia calculations (0=no inertia check)
      ipar( 9)=inerti
!---printed report (0=none)
      ipar(10)=0
!---maximum number of steps per shift (default=50)
      ipar(11)=max(50,nstpsh)
!---number/range of eigenvalues (0=number)
      ipar(12)=0
!---storage format of K (must be 0)
      ipar(13)=0
!---storage format of M (must be 0)
      ipar(14)=0
!---reserved for future use (level of blas--unrolling)
      ipar(15)=0
!---factorization type
!---0 banded LDL or Bunch--Kaufman
!---1 banded LDL
!---2 sparse LDL and sparse Bunch--Kaufman
!---3 sparse LDL
      ipar(16)=3
!---dynamic shifting (0=on)
      ipar(17)=0
!---no initial guess
      ipar(18)=0
!---leading index of array y
      ipar(19)=npogeo
!---for ipar(16)=2 number of eigenvalues less than the largest searched for
      ipar(20)=nvala-1
!
!---array rpar
!---shift to be searched around
      rpar(1)=dshift
!---relative accuracy
      rpar(2)=relacc
!---left end of the range (for ipar(12)=1)
      rpar(3)=z0
!---right end of the range (for ipar(12)=1)
      rpar(4)=z1
!---storage factor if ipar(16)=0
      rpar(5)=1.1d0
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine masmat(masrow,mascol,sysmas)
      include 'neurofem.inc'
      dimension masrow(lenrow),mascol(lencol)
      dimension sysmas(lenlan)
!
      do 10 i=1,npogeo
         masrow(i)=1
         mascol(i)=1
         sysmas(i)=z1
         sysmas(i+npogeo)=z0
   10 continue
      masrow(npogeo+1)=1
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     subroutine   w a n d e l
!                  ===========
!     written by Richard Schoenen (c)
!
! Matrix in Lanz-Format wandeln
!
!  iinp            : Input-Unit fuer IME-Matrix
!  iout            : Output-Unit fuer LANZ-Matrix
!  ipodia(nfrhg)   : Position der Diagonalelemente in der Matrix IME
!  indexj(laenge)  : Spalte fuer jedes Matrixelement IME
!  gesmat(laenge)  : Matrix in kompakter Speicherung IME
!  irow(nfrhg+1)   : Zeilen-Pointer LANZ
!  icol(laenge)    : Spalten-Indizes LANZ
!  a(laenge)       : Matrix in kompakter Speicherung LANZ
!
! HILFSGROESSEN:
!  iend(nfrhg)     : Zeigt auf das zuletzt gefundene Element einer Spalte
!  ivor(laenge)    : Zeigt auf das vorher gefundene Element gleicher Spalte
!
! Die Felder iend und ivor ermoeglichen das schnelle Durchlaufen einer
! Spalte der kompakt gespeicherten Matrix
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wandel(nfrhg,laenge,                                   &
     &           gesmat,indexj,ipodia,ivor,iend,a,icol,irow)
      implicit double precision (a-h,o-z)
      dimension ipodia(nfrhg),indexj(laenge),ivor(laenge),iend(nfrhg),  &
     &          icol(laenge),irow(nfrhg+1)
      dimension gesmat(laenge),a(laenge)
!
! Matrix im IME binaer-Format einlesen
!
! Felder nullen
!
      do i=1,nfrhg
         iend(i)=0
      end do
!
! Index-Feld zum Spaltendurchlauf bestimmen
! Durch die rueckwaerts laufende Schleife zeigt iend(j) schliesslich
! auf die erste Zeile, die Spalte j enthaelt, und der Pointer ivor(ip)
! zeigt fuer jedes Matrixelement auf das naechste der gleichen Spalte
!
      do i=nfrhg,1,-1
         ipa=1
         if(i.gt.1) ipa=ipodia(i-1)+1
         ipe=ipodia(i)
         do ip=ipa,ipe
            j=indexj(ip)
            ivor(ip)=iend(j)
            iend(j)=ip
         end do
      end do
!
! in Lanz-Format wandeln (oberer rechter Teil)
! zuerst die Diagonalelemente
!
      do i=1,nfrhg
        a(i)=gesmat(ipodia(i))
      end do
!
! jetzt die Nichtdiagonalelemente der rechten oberen Matrix
! wegen der Symmetrie entsprechen die Spalten von gesmat den
! Zeilen von a
!
      iplanz=0
      irow(1)=iplanz+1
      do jime=1,nfrhg
         ioff=jime
         ip=iend(jime)
         ip=ivor(ip)   ! Diagonalelement ueberspringen
   10    continue
         if(ip.ne.0)then
            if(indexj(ip).ne.jime) stop 'Fehler 1 in WANDEL (intern)'
!
! weiteres Element der Spalte gefunden
! --> Zeile ermitteln & Element merken
!
            do i=ioff,nfrhg
               if(ipodia(i).ge.ip) goto 15
            end do
            stop'Fehler 2 in WANDEL (intern)'
   15       continue
            ioff=i
!
            iplanz=iplanz+1
            a(nfrhg+iplanz)=gesmat(ip)
            icol(iplanz)=i
!
            ip=ivor(ip)
            goto 10
         end if
         irow(jime+1)=iplanz+1
      end do
      if(iplanz+nfrhg.ne.laenge) stop'Fehler 3 in WANDEL (intern)'
!
      call iset (laenge,0,ivor,1)
      call iset ( nfrhg,0,iend,1)
!
      return
   90 continue
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chklan(ipogeo,indgeo,lanrel,sysmat,volmat,             &
     &           vallnz,veclnz,dprodg,dkondg)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),lanrel(lenval)
      dimension sysmat(lengeo),volmat(lengeo),                          &
     &          veclnz(npogeo,lenval),vallnz(lenval,2),                 &
     &          dprodg(npogeo),dkondg(npogeo)
!
      if (lconsi) then
         write(*,800) 'general'
      else
         write(*,800) 'special'
      end if
!
      ntrumo=nwriln
      do 10 k=1,nwriln
         val=vallnz(k,1)
         call matpro(sysmat, veclnz(1,lanrel(k)),dprodg,                &
     &        ipogeo,indgeo,lengeo,npogeo)
         if (lconsi) then
            call matpro(volmat, veclnz(1,lanrel(k)),dkondg,             &
     &           ipogeo,indgeo,lengeo,npogeo)
         else
            call dcopy(npogeo,veclnz(1,lanrel(k)),1,dkondg,1)
         end if
         call daxpy (npogeo,-val,dkondg,1,dprodg,1)
         residu=dnrm2(npogeo,dprodg,1)
         write(*,900) k,val,residu
         if (val.lt.tolsol.and.residu.gt.tol) then
            ntrumo=ntrumo-1
         end if
   10 continue
      if (nwriln.ne.ntrumo) then
         write(*,910) nwriln,ntrumo
         nwriln=ntrumo
      end if
!
      return
  800 format(' Checking residuals of the ',a,' eigenvalue problem')
  900 format(' Eigenvalue number ',i8,' is ',e12.5,' accuracy ',e12.5)
  910 format(' Number of eigenvalues reduced from ',i8,' to ',i8)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function lnzsrt(lanrel,nvala,i)
      integer lnzsrt,i,lanrel(nvala)
!
      lnzsrt=lanrel(i)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine svdinf(nodinf,flxmat,vtrval,sngval,usvval,             &
     &                  wrksvd,weiinv,covinv,itask)
      include 'neurofem.inc'
      dimension nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          vtrval(nevlpo,ninfpo*ndmflx),                           &
     &          sngval(nevlpo),usvval(nevlpo,nevlpo),                   &
     &          wrksvd(nwrksv),weiinv(npoinv,ndmflx),                   &
     &          covinv(nevlpo)
!
!---store flxmat into twodimensional array vtrval
      write(*,900)
      do 10 k=1,ndmflx
         do 20 j=1,ninfpo
            jmod=(j-1)*ndmflx+k
            jinv=nodinf(j)
            we=weiinv(jinv,k)
            do 30 i=1,nevlpo
               co=sqrt(covinv(i))
               vtrval(i,jmod)=flxmat(i,j,k)*we*co
   30       continue
   20    continue
   10 continue
!
      dum=qscput(1,0,ierr)
      chasvu='S'
      chasvv='O'
      call dgesvd(chasvu,chasvv,nevlpo,ninfpo*ndmflx,vtrval,nevlpo,     &
     &            sngval,usvval,nevlpo,vtrval,nevlpo,wrksvd,nwrksv,     &
     &     infsvd)
!
      if (infsvd.eq.0.and.itask.gt.0) then
         write(*,910)
         write(*,915) nwrksv,int(wrksvd(1)+f2)
         write(*,990) qscput(1,1,ierr)
      else if (infsvd.gt.0) then
         write(*,920) 
         write(*,930) infsvd 
         if (infsvd.lt.0) write(*,940) -infsvd 
         if (infsvd.gt.0) write(*,950)  infsvd 
         stop 'dgesvd'
      end if
!
      if (itask.gt.0) then
         write(*,960)
         difmax=z0
         do 40 ievl=1,nevlpo
            do 50 infl=1,ninfpo
               do 60 k=1,ndmflx
                  tesinf=flxmat(ievl,infl,k)
                  ipoj=(infl-1)*ndmflx+k
                  do 80 i=1,nevlpo
                     wrksvd(i)=sngval(i)*vtrval(i,ipoj)
   80             continue
                  sum=z0
                  do 70 i=1,nevlpo
                     sum=sum+usvval(ievl,i)*wrksvd(i)
   70             continue
                  difmax=max(abs(sum-tesinf),difmax)
   60          continue
   50       continue
   40    continue
!
         write(*,965) difmax
         if (difmax.lt.wuzeps) then
            write(*,970) '     OK'
         else
            write(*,970) ' not OK'
         end if
         write(*,980)
      end if
!
      return
  900 format(' Preparing singular value decomposition')
  910 format(' Singular value decomposition successfull ')
  915 format(' work array size for SVD was',i10,/,                      &
     &       ' optimum work array size  is',i10)
  920 format(' Singular value decomposition failed (!) ')
  930 format(' Error code of SVD is ',i10)
  940 format(' Argument # ',i10,' had an illegal value')
  950 format(' ',i10,' off-diagonal elements did not converge to zero')
  960 format(' Checking SVD by comparison to influence matrix')
  965 format(' Maximum difference between SVD and FLXMAT is ',g12.5)
  970 format(' SVD was checked ',a)
  980 format(' Singular value decomposition finished')
  990 format(' CPU--time for SVD was ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine svdvec(nodinf,vtrval,dprflf,ivcakt,idmakt)
      include 'neurofem.inc'
      dimension nodinf(ninfpo)
      dimension vtrval(nevlpo,ninfpo*ndmflx),dprflf(npoflf)
!
      call dset(npoflf,z0,dprflf,1)
      dprmax=z0
      do 10 i=1,ninfpo
         iflpoi=nodinf(i)
         ipos=(i-1)*ndmflx+idmakt
         dprflf(iflpoi)=vtrval(ivcakt,ipos)
         dprmax=max(abs(dprflf(iflpoi)),dprmax)
   10 continue
      if (dprmax.gt.z0) then
         dprsca=z1/dprmax
      else
         dprsca=z1
      end if
      call dscal(npoflf,dprsca,dprflf,1)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine wrisvd(usvval,sngval,vtrval,icall)
      include 'neurofem.inc'
      dimension usvval(nevlpo,nevlpo),sngval(nevlpo),                   &
     &          vtrval(nevlpo,ninfpo*ndmflx)
!
      iunfil=qfreeu()
      call qfopse(iunfil,filnam(33),'un','un',ierr)
      if (ierr.gt.0) stop 'qfopse filnam(33)'
      if (icall.eq.1) then
         write(iunfil) ((usvval(j,i),j=1,nevlpo),i=1,nevlpo)
         write(iunfil)  (sngval(j)  ,j=1,nevlpo)
         write(iunfil) ((vtrval(j,i),j=1,nevlpo),i=1,ninfpo*ndmflx)
         iungnu=qfreeu()
         ipunkt=index(filnam(33),'.')
         zeil(1:ipunkt)=filnam(33)(1:ipunkt)
         zeil(ipunkt+1:ipunkt+3)='gnu'
         zeil(ipunkt+4:80)=' '
         call qfopse(iungnu,zeil,'un','fo',ierr)
         if (ierr.gt.0) stop 'qfopse filnam(33).gnu'
         write(iungnu,900)
         write(iungnu,910) (j,sngval(j)  ,j=1,nevlpo)
         call qfclos(iungnu,0)
      else
         read(iunfil) ((usvval(j,i),j=1,nevlpo),i=1,nevlpo)
         read(iunfil)  (sngval(j)  ,j=1,nevlpo)
         read(iunfil) ((vtrval(j,i),j=1,nevlpo),i=1,ninfpo*ndmflx)
      end if
      call qfclos(iunfil,0)
!
      return
  900 format('# no singulaerwert')
  910 format(1x,i10,1x,g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!
