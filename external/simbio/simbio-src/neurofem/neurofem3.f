!  $2 21/03/2003  Anwander A.  released for SIMBIO
!  $1 30/06/00aa multiple return in do loops removed

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem3.f :      	 Regularized inverse solutions
!     	                        in volumes or on surfaces in 3D
!     		                FOCUSS algorithm
!     		                truncated CG regularization
!     				depth weighting
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

      subroutine invers( nodevl,nodinf,invbnd,                          &
     &     mrkinv,knoinv,ityinv,necinv,                                 &
     &     ipoinv,indinv,                                               &
     &     stfinv,vcbinv,flxmat,                                        &
     &     potmes,bndinv,solinv,                                        &
     &     vcrinv,vcpinv,dprdin,dkndin,                                 &
     &     vechey,vechew,vechev,scainv,                                 &
     &     weiinv,dipinv,dirmat,dirrhs,                                 &
     &     vcnflf,srcinv,xyzinv,covinv,                                 &
     &     solold,xivinv,pcoinv,xicinv,xtcinv)
!
      include 'neurofem.inc'
!
      dimension nodevl(nevlpo),nodinf(ninfpo),invbnd(npoinv),           &
     &          mrkinv(npoinv,ndminv),knoinv(mxknin,nelinv),            &
     &          ityinv(nelinv),necinv(nelinv),                          &
     &          ipoinv(npoinv),indinv(leninv)
      dimension stfinv(leninv,ndmflx),vcbinv(npoinv,ndmflx),            &
     &          flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          potmes(nevlpo),bndinv(npoinv),                          &
     &          solinv(npoinv,ndmflx)
      dimension vcrinv(npoinv,ndmflx),vcpinv(npoinv,ndmflx),            &
     &          dprdin(npoinv,ndmflx),dkndin(npoinv,ndmflx),            &
     &          vechey(nevlpo),vechew(npoinv,ndmflx),                   &
     &          vechev(npoinv,ndmflx),scainv(npoinv,ndmflx),            &
     &          weiinv(npoinv,ndmflx)
      dimension dipinv(npoinv,ndminv),dirmat(ndimat),dirrhs(ndirhs),    &
     &          vcnflf(npoflf,ndmflf),srcinv(npogeo),                   &
     &          xyzinv(npoinv,ndminv),covinv(nevlpo),                   &
     &          solold(npoinv,ndmflx)
      dimension xivinv(npoinv,ndmflx),pcoinv(npoinv,ndmflx),            &
     &          xicinv(npoinv,ndmflx),xtcinv(npoinv,ndmflx)
!
      xlmdn=wuzeps
      xlmup=z100*z10
!
!---build smoothing operator for tikhonov--phillips 
      call invsur(knoinv,ityinv,ipoinv,indinv,necinv,                   &
     &            stfinv,xyzinv)
!
      call dset(npoinv*ndmflx,z1,weiinv,1)
      if (lfocus) call dset(npoinv*ndmflx,z0,solold,1)
      nfocus=0
  444 call matrhs(invbnd,nodinf,flxmat,solinv,vcbinv,potmes,            &
     &            weiinv,covinv)
!
      if (lkonve) then
         iunco=qfreeu()
         call qfopse(iunco,filnam(19),'un','fo',ierr)
         if (ierr.gt.0) then
            write (*,900) 'fehler filnam(19) in invers'
            stop
         else 
            write(iunco,910)
         end if
      end if
!
      if (lchisl) write(*,*)' desired chi-square value ', chival
  700 write(*,*)
      write(*,*)' actual value of lambda: ',xlmbda
!
      if (ldirsl) then
!
!---compute a direct (CHOLESKY) solution of the inverse problem
!
         call dset(ndimat,z0,dirmat,1)
         call dset(ndirhs,z0,dirrhs,1)
!---compile system of equations (packed storage for lower half of the matrix)
         call dirsys(ipoinv,indinv,nodevl,nodinf,invbnd,                &
     &        stfinv,vcbinv,flxmat,dirmat,dirrhs,covinv)
!---factor the matrix and estimate its condition number
         call dppco(dirmat,ndirhs,rcond,dprdin,info)
         if (info.ne.0) write(*,'(a,i8)') 'dirmat singular? ',info
         write(*,1040) z1/rcond
!---solve
         call dppsl(dirmat,ndirhs,dirrhs)
         ikenn=0
      else
!
!---CG solution of the inverse problem
!
!---set start solution for cg algorithm
         do 4000 i=1,ninfpo
            iglo=nodinf(i)
            do 4010 k=1,ndmflx
               if (invbnd(iglo).gt.0) then
                  solinv(iglo,k)=z0
               else 
                  solinv(iglo,k)=z0
               end if
 4010       continue
 4000    continue
!---depth--weight factors (scainv)
         call depsca(ipoinv,invbnd,nodinf,indinv,                       &
     &        flxmat,scainv,stfinv,solinv,vcbinv,                       &
     &        weiinv,covinv)
!---solve
         if (lcglin) then
            call cgrinv(ipoinv,indinv,invbnd,nodinf,                    &
     &           stfinv,solinv,vcbinv,flxmat,                           &
     &           dprdin,vcrinv,vcpinv,vechey,                           &
     &           vechew,vechev,scainv,weiinv,                           &
     &           covinv,potmes,resi,  iter,icount,  ierr)
            write(*,1010) resi
         else
            eps=1.d-10
            nout=100
            if (ninfpo.lt.1000) nout=20
            if (ninfpo.lt.100) nout=10
            if (ninfpo.lt.10) nout=1
            itmax=max(ninfpo*ndmflx,10)
            do 100 iover=1,5
               call frprmn(solinv,npoinv*ndmflx,tolinv,  iter,          &
     &              fret,vcrinv,vcpinv,xivinv, itmax,                   &
     &              eps,pcoinv,xicinv,xtcinv,wuzeps,nout)
               if (iter.lt.itmax) goto 200
  100       continue
  200       continue
         end if
!
!---check if data term is within the expected range
         if (lkonve) then
            call funchk(solinv,dat,xl1,xl2,xlc,ent,fun)
            write(iunco,920) xlmbda,dat,xl2,xl1,ent,xlc,fun
!            write(    *,930) xlmbda,dat,xl2,xl1,ent,xlc,fun
         end if
!
         call chichk(nodinf,invbnd,nodevl,flxmat,covinv,solinv,         &
     &               potmes,vechey,weiinv,ispru)
         if (ispru.lt.0) goto 700
!$$$11.9.97         if ((llamit.or.lkonve).and.ispru.lt.0) goto 700
         if (lkonve) then
            call qfclos(iunco,0)
         end if
!
!---check convergence of the focuss procedure
         if (lfocus) then
            call focuss(solinv,weiinv,solold,nfocus,ispru)
            if (ispru.lt.0) goto 444
         end if
         ikenn=1
      end if
!
!---plot data and different model terms for different norms
      call funchk(solinv,dat,xl1,xl2,xlc,ent,fun)
      write(*,930) xlmbda,dat,xl2,xl1,ent,xlc,fun
!
!---transfer the results to global numbering
      call resglo(nodinf,mrkinv,                                        &
     &            dirrhs,dipinv,solinv,                                 &
     &            srcinv,vcnflf,ikenn)
!
      return
  900 format(a)
  910 format('#  lambda       chi-square   l2-norm      l1-norm',       &
     &       '      entropy      l2-c         fun') 
  920 format(7(1x,e12.5,:))
  930 format(' l ',e11.5,' d ',e11.5,' l2 ',e11.5,                      &
     &      ' l1 ',e11.5,' e ',e11.5,' lc ',e11.5,' f ',e11.5)
 1010 format(' Conjugate Gradient Inverse Residual: ',g12.5)
 1040 format(' Condition number ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine invsur(knoinv,ityinv,ipoinv,indinv,necinv,             &
     &                  stfinv,xyzinv)
      include 'neurofem.inc'
!
      dimension knoinv(mxknin,nelinv),ityinv(nelinv),indinv(leninv),    &
     &          necinv(nelinv),ipoinv(npoinv)
      dimension stfinv(leninv,ndmflx),xyzinv(npoinv,ndminv)
      dimension syslok(mxkn3d,mxkn3d)
!
      call dset(leninv*ndmflx,z0,stfinv,1)
!
!---if kronecker delta matrix (iunit matrix) is requested
      if (iglatt.eq.0) then
         do 200 i=1,npoinv
            do 210 k=1,ndmflx
               stfinv(ipoinv(i),k)=z1
  210       continue
  200    continue
         trasmo=dble(npoinv)
      else
!
!---if regularization of values or gradients is requested
         do 10 iel=1,nelinv
            call elesur(knoinv,ityinv,necinv,                           &
     &           xyzinv,syslok, iel )
            do 20 i=1,necinv(iel)
               iglo=knoinv(i,iel)
               do 30 j=1,necinv(iel)
                  jglo=knoinv(j,iel)
                  if (jglo.le.iglo) then
                     if (iglo.gt.1) then
                        ia=ipoinv(iglo-1)+1
                        in=ipoinv(iglo)-ipoinv(iglo-1)
                        iposi=ia+ipoj(jglo,indinv(ia),in)
                     else
                        iposi=1
                     end if
                     do 40 k=1,ndmflx
                        stfinv(iposi,k) = stfinv(iposi,k) + syslok(i,j)
   40                continue
                  end if
   30          continue
   20       continue
   10    continue
         trasmo=z0
         do 50 i=1,npoinv
            id=ipoinv(i)
            do 55 k=1,ndmflx
               trasmo=trasmo+stfinv(id,k)
   55       continue
   50    continue
      end if
      write(*,900) trasmo
!
      return
  900 format( ' Trace of the smoothing operator is: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elesur(knoinv,ityinv,necinv,                           &
     &     xyzinv,syslok,iel)
      include 'neurofem.inc'
      dimension knoinv(mxknin,nelinv),ityinv(nelinv),necinv(nelinv),    &
     &          lokkno(32)
      dimension xyzinv(npoinv,ndminv)
      dimension syslok(mxkn3d,mxkn3d)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn)
!
      call iset(32,0,lokkno,1)
!
!---turn to local node numbers
      do 10 i=1,necinv(iel)
        kno=knoinv(i,iel)
        lokkno(i)=kno
        do 15 idim=1,ndminv
           xlk(idim,i) = xyzinv(kno,idim)
 15     continue
!     
!---zero element stiffness matrix and right hand side
         fag(2,i) = z0
         fag(3,i) = z0
         do 10 j=1,necinv(iel)
            syslok(j,i) = z0
   10 continue
!
!---determine number of gaussian points
      call itpanz(ityinv(iel),necinv(iel),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'eleinv 1'
!
!---loop over all gaussian points
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
      call itpdat(ityinv(iel),intgrd,iint,dinpkt,ierr)
      if (ierr.ne.0) stop 'eleinv 2'
      call fofagn(ityinv(iel),dinpkt(1),dinpkt(2),dinpkt(3),            &
     &            lokkno(1),xlk,fwe,fag,detj,ierr)
      if (ierr.ne.0) stop 'eleinv 3'
!
      xfak =   detj * dinpkt(4)
!
!---loop over all nodal points and summation of contributions
      do 30 i=1,necinv(iel)
!$$$30/06/00aa multible return not possible in f90
!         do 30 j=1,necinv(iel)
         do 35 j=1,necinv(iel)
            if (iglatt.eq.2) then
               sysdum=z0
!$$$               do 50 k=1,ndmgeo-1
               do 50 k=1,ndminv
                  sysdum = sysdum + fag(k,i)*fag(k,j)
 50            continue
               syslok(j,i) = syslok(j,i) + sysdum * xfak
            else
               syslok(j,i) = syslok(j,i) + fwe(i)*fwe(j) * xfak
            end if
 35      continue
 30      continue
 20   continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dirsys(ipoinv,indinv,nodevl,nodinf,invbnd,             &
     &                  stfinv,vcbinv,flxmat,dirmat,dirrhs,             &
     &                  covinv)
      include 'neurofem.inc'
      dimension ipoinv(npoinv),indinv(leninv),nodevl(nevlpo),           &
     &          nodinf(ninfpo),invbnd(npoinv)
      dimension stfinv(leninv,ndmflx),vcbinv(npoinv,ndmflx),            &
     &          dirmat(ndimat),dirrhs(ndirhs),                          &
     &          flxmat(nevlpo,ninfpo,ndmflx),covinv(nevlpo)
!
!---zero vectors
      call dset(ndimat,z0,dirmat,1)
      call dset(ndirhs,z0,dirrhs,1)
!
!---apply homogeneous boundary conditions
!---by setting each diagonal element z1
!---(elements without homogeneous boundary conditions will be
!---overwritten later)
      do 400 i=1,npoinv
         if (invbnd(i).lt.1) goto 400
         do 410 k=1,ndmflx
            ipodof=(i-1)*ndmflx+k
            ipomat=ipodof*(ipodof+1)/2
            dirmat(ipomat)=z1
  410    continue
  400 continue
!
!---right hand side
      do 10 i=1,npoinv
         do 20 k=1,ndmflx
            ipodof=(i-1)*ndmflx+k
            dirrhs(ipodof)=-vcbinv(i,k)
   20    continue
   10 continue
!
!---system matrix for normal equations ( E_{iu}^v E_{ir}^s )
      do 30 j=1,ninfpo
         iu=nodinf(j)
         do 40 kv=1,ndmflx
            do 50 l=1,ninfpo
               ir=nodinf(l)
               do 60 ks=1,ndmflx
                  ipodof=(iu-1)*ndmflx+kv
                  jpodof=(ir-1)*ndmflx+ks
                  if (jpodof.gt.ipodof) goto 60
                  ipomat=ipodof*(ipodof-1)/2+jpodof
                  sum=z0
                  do 70 ievl=1,nevlpo
                     sum=sum+flxmat(ievl,j,kv)*flxmat(ievl,l,ks)        &
     &                      *covinv(ievl)
   70             continue
                  dirmat(ipomat) = sum
   60          continue
   50       continue
   40    continue
   30 continue
!     
!---smoothing operator
      do 330 k=1,ndmflx
         ipodof=k
         ipomat=(ipodof)*(ipodof+1)/2
         dirmat(ipomat)=dirmat(ipomat)+xlmbda*stfinv(1,k)
  330 continue
      do 300 i=2,npoinv
         do 310 j=ipoinv(i-1)+1,ipoinv(i)
            jglo=indinv(j)
            do 320 k=1,ndmflx
               ipodof=(i-1)*ndmflx+k
               ipopak=ipodof*(ipodof-1)/2
               jpodof=(jglo-1)*ndmflx+k
               ipomat=ipopak+jpodof
               dirmat(ipomat)=dirmat(ipomat)+xlmbda*stfinv(j,k)
  320       continue
  310    continue
  300 continue
!
!---sysmat will be checked
      do 500 i=1,npoinv
         do 510 k=1,ndmflx
            ipodof=(i-1)*ndmflx+k
            ipomat=ipodof*(ipodof+1)/2
            if (dirmat(ipomat).lt.wuzeps) then
               write (*,900) ipodof
            end if
            if (dirmat(ipomat).eq.z1.and.dirrhs(ipodof).ne.z0)then
               write(*,910) ipodof
            end if
  510    continue
  500 continue
!
      return
  900 format('Dipmat problems with small pivot at ',i8)
  910 format('Homogeneous boundary condition problem at ',i8)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine focuss(solinv,weiinv,solold,nfocus,ispru)
      include 'neurofem.inc'
      dimension solinv(npoinv,ndmflx),weiinv(npoinv,ndmflx),            &
     &          solold(npoinv,ndmflx)
      save tolfoc,tolmax
      data tolfoc,tolmax /1.d-6,0.9/
!
      if (nfocus.eq.0) then
!         open(55,file='focus.3d')
         write(*,920) 
         xlmbda=z0
      end if
      nfocus=nfocus+1
      iactiv=0
      vecmax=z0
      do 510 i=1,npoinv
         dum=z0
         do 500 k=1,ndmflx
            solinv(i,k)=solinv(i,k)*weiinv(i,k)
            dum=dum+solinv(i,k)*solinv(i,k)
  500    continue
         dum=sqrt(dum)
         vecmax=max(vecmax,dum)
         if (dum.gt.z0) iactiv=iactiv+1
!         write(55,910) dble(nfocus),dble(i),dum
  510 continue
!      write(55,910)
!
      tolakt=max(tolfoc*vecmax,tolfoc)
      change=z0
      do 610 i=1,npoinv
         dum=z0
         do 600 k=1,ndmflx
            dum=dum+solinv(i,k)*solinv(i,k)
            dif=solold(i,k)-solinv(i,k)
            change=change+abs(dif)
            solold(i,k)=solinv(i,k)
  600    continue
         dum=sqrt(dum)
         if (dum.lt.tolakt) then
            do 800 k=1,ndmflx
               weiinv(i,k)=z0
  800       continue
         else
            do 810 k=1,ndmflx
               weiinv(i,k)=weiinv(i,k)*solinv(i,k)
  810       continue
         end if
  610 continue
      change=change/dble(npoinv*ndmflx)
      ispru=-1
      write(*,900) iactiv,change
      if (iactiv.le.numdip.or.change.lt.tolakt.or.chisqr.gt.chimax) then
         if (chisqr.gt.chimax) then
            call dcopy(npoinv*ndmflx,solold,1,solinv,1)
         end if
         ispru=1
      end if
!     
      return
  900 format(' Active sources: ',i6,' average change: ',e12.5)
  910 format(3(1x,e12.5))
  920 format(' FOCUSS procedure redefines lambda inverse to zero')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chichk(nodinf,invbnd,nodevl,flxmat,covinv,solinv,      &
     &           potmes,potcal,weiinv,ispru)
      include 'neurofem.inc'
      dimension nodinf(ninfpo),invbnd(npoinv),nodevl(nevlpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),covinv(nevlpo),            &
     &          solinv(npoinv,ndmflx),potmes(nevlpo),potcal(nevlpo),    &
     &          weiinv(npoinv,ndmflx)
!
      ispru=1
      do 10 i=1,nevlpo
         potcal(i)=z0
         do 20 j=1,ninfpo
            jglo=nodinf(j)
            if (invbnd(jglo).lt.1) then
               do 30 k=1,ndmflx
                  fx=flxmat(i,j,k)
                  sv=solinv(jglo,k)*weiinv(jglo,k)
                  potcal(i)=potcal(i)+fx*sv
   30          continue
            end if
   20    continue
   10 continue
!
!      iunpo=qfreeu()
!      call qfopse(iunpo,'pot-vgl.gnu','un','fo',ierr)
!      if (ierr.gt.0) then
!         write (*,900) 'error pot-vgl.gnu'
!         stop
!      end if
!      write(iunpo,950)
      chisqr=z0
      do 40 i=1,nevlpo
         node=nodevl(i)
         dif=potcal(i)-potmes(i)
         co=z1/sqrt(covinv(i))
         po=potmes(i)
!         write(iunpo,960) i,node,po,po+co,po-co,potcal(i),dif
         chisqr=chisqr+dif*covinv(i)*dif
   40 continue
!      call qfclos(iunpo,0)
!
! -- put lamdba-iteration variance on 10% of -chival-             
      chivar = chival * .1d0
      chimin = chival - chivar
      chimax = chival + chivar
      chimid = chival
!
! -- statistical range of chi-square defined by -nevlpo-
!     chivar=sqrt(z2*nevlpo)
!     chimin=dble(nevlpo)-chivar
!     chimax=dble(nevlpo)+chivar
!     chimid=dble(nevlpo)
!
      if (llamit.or.lchisl) then
         if (chisqr.gt.chimin.and.chisqr.lt.chimax) then
            write(*,880)
            ispru=1
            write(*,920) xlmbda
         else
            write(*,890) 
            ispru=-1
            if (chisqr.lt.chimin) then
               xlmdn=xlmbda
               xlmdif=f2*(xlmup-xlmdn)
               xlmbda=xlmbda+xlmdif
            else if (chisqr.gt.chimax) then
               xlmup=xlmbda
               xlmdif=f2*(xlmup-xlmdn)
               xlmbda=xlmbda-xlmdif
            end if
            write(*,910) xlmdn,xlmbda,xlmup
         end if
      else if (lkonve) then
         if (xlmbin.eq.z1) then
            ispru=1
         else if (xlmbin.gt.z1) then
            if (xlmbda.lt.xlmbfi) then
               xlmbda=xlmbda*xlmbin
               ispru=-1
            else
               ispru=1
            end if
         else if (xlmbin.gt.z0) then
            if (xlmbda.gt.xlmbfi) then
               xlmbda=xlmbda*xlmbin
               ispru=-1
            else
               ispru=1
            end if
         else
            ispru=1
         end if
      end if
!
! -- statistical and desired range of chi-square defined by -nevlpo-/-chival-
      if(.not.lkonve) then
         if (lchisl) then
            write(*,900) chisqr,dble(nevlpo)-sqrt(z2*nevlpo),chimin,    &
     &                          dble(nevlpo)+sqrt(z2*nevlpo),chimax,    &
     &                          dble(nevlpo)                ,chimid
         else
            write(*,905) chisqr,dble(nevlpo)-sqrt(z2*nevlpo),           &
     &                          dble(nevlpo)+sqrt(z2*nevlpo),           &
     &                          dble(nevlpo)
         end if
      end if
!
      RETURN
  880 format(' Chi-square confidence range     reached')
  890 format(' Chi-square confidence range not reached')
  900 format('        Chi-square-data term is: ',g12.5,/,               &
     &       33x,'standard     desired-value',/                         &
     &       ' it should be in the range from: ',g12.5,2x,g16.5,/       &
     &       '                             to: ',g12.5,2x,g16.5,/       &
     &       '         with expectation value: ',g12.5,2x,g16.5   )
  905 format('        Chi-square-data term is: ',g12.5,/                &
     &       ' it should be in the range from: ',g12.5,/                &
     &       '                             to: ',g12.5,/                &
     &       '         with expectation value: ',g12.5  )
!$$$  900 format('        Chi-square-data term is: ',g12.5,/,
!$$$     &       ' it should be in the range from: ',g12.5,/
!$$$     &       '                             to: ',g12.5,/,
!$$$     &       '         with expectation value: ',g12.5 )
  910 format(' lam-min ',e12.6,' lam-akt ',e12.6,' lam-max ',e12.6)
  920 format(' final   lambda value is: ',g12.5)
  950 format('# i node pot-mes pot-mes+co pot-mes-co pot-cal dif')
  960 format(2(1x,i6),5(1x,e12.5))
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine george(nodinf,invbnd,flxmat,vtrval,sngval,usvval,      &
     &                  wrksvd,weiinv,covinv,solinv,potmes,             &
     &                  cutsol,solold,vechey)
      include 'neurofem.inc'
      dimension nodinf(ninfpo),invbnd(npoinv)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          vtrval(nevlpo,ninfpo*ndmflx),                           &
     &          sngval(nevlpo),usvval(nevlpo,nevlpo),                   &
     &          wrksvd(nwrksv),weiinv(npoinv,ndmflx),                   &
     &          covinv(nevlpo),solinv(npoinv,ndmflx),                   &
     &          potmes(nevlpo,numtim),cutsol(ninfpo*ndmflx),            &
     &          solold(npoinv,ndmflx),vechey(nevlpo)
!
      call dset(npoinv*ndmflx,z0,solold,1)
      deforg=z0
      do 10 i=1,nevlpo
         deforg=deforg+covinv(i)*potmes(i,1)*potmes(i,1)
   10 continue
      nfocus=0
   20 defect=deforg
      call svdinf(nodinf,flxmat,vtrval,sngval,usvval,wrksvd,weiinv,     &
     &                  covinv,0)
      call solsvd(nodinf,vtrval,sngval,usvval,weiinv,solinv,            &
     &                  covinv,potmes,cutsol)
!
      do 30 i=1,ninfpo
         iglo=nodinf(i)
         do 40 k=1,ndmflx
            j=(i-1)*ndmflx+k
            solinv(iglo,k)=cutsol(j)
   40    continue
   30 continue
      call focuss(solinv,weiinv,solold,nfocus,ispru)
      if (ispru.lt.0) goto 20
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine solsvd(nodinf,vtrval,sngval,usvval,weiinv,solinv,      &
     &                  covinv,potmes,cutsol)
      include 'neurofem.inc'
      dimension nodinf(ninfpo)
      dimension vtrval(nevlpo,ninfpo*ndmflx),sngval(nevlpo),            &
     &          usvval(nevlpo,nevlpo),weiinv(npoinv,ndmflx),            &
     &          solinv(npoinv,ndmflx),covinv(nevlpo),                   &
     &          potmes(nevlpo,numtim),cutsol(ninfpo*ndmflx)
!
      call dset(ninfpo*ndmflx,z0,cutsol,1)
      do 10 i=1,nevlpo
         if (llucut) then
            if (defect.lt.rpostl*rpostl) goto 20
            if (sngval(i).lt.z1) goto 20
         end if
         alphan=z0
         do 30 j=1,nevlpo
            alphan=alphan+usvval(j,i)*potmes(j,1)*sqrt(covinv(j))
   30    continue
         faktor=alphan/sngval(i)
         do 40 j=1,ninfpo*ndmflx
            cutsol(j)=cutsol(j)+faktor*vtrval(i,j)
   40    continue
         defect=defect-alphan*alphan
   10 continue
   20 continue
      write(*,900) i,defect
!
      return
  900 format(' Truncating SVD from i=',i8,' defect: ',e12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!
!
!
!
