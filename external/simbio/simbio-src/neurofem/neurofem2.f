!  $5 04.05.2007  Seok Lew: change to energy norm as PCG stop condition
!  $4 21/03/2003  Anwander A.  released for SIMBIO
!  $3 25.04.2001  aa Since ipvtdi(maxnei)-Parameter in subroutine "diplod" 
!                     infmat, forwrd, rhsblu is never used, it was deleted.
!  $2 12/21/00    cw parameter "srcgeo(npogeo)" was deleted in the calls 
!                    of "subroutine diplod"
!  $1 11/07/00    aa full matrix constructed and write out the equation system

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem2.f : Solution of the forward problem
!     ---> Lead field (Influence) matrices
!     ---> Reference solutions
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

      subroutine infmat(ipogeo,indgeo,mrkvec,                           &
     &                  ipoasy,indasy,ibndgo,inflor,                    &
     &                  ievlor,neigeo,intgam,itypmg,                    &
     &                  megtyp,knoobf,knosur,itysur,necsur,             &
     &                  xyzgeo,difvec,dipmat,diprhs,                    &
     &                  dipsol,vecbgo,mrkdip,dipole,                    &
     &                  sysmat,rbdgeo,vecsol,vecrgo,                    &
     &                  vecpgo,dprodg,dkondg,gstifc,                    &
     &                  vecrho,flxmat,volmat,srcgeo,                    &
     &                  dipreg,surmat,alfmat,srcsur,                    &
     &                  vcbamb,vcnflf,vecref,syscpy,                    &
     &                  xyzflf,veceig,redinv,fstgam,                    &
     &                  potmeg,sekmeg,potsng,sysrng)
!
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkvec(maxnei),           &
     &          mrkdip(npoflf,ndmflf),                                  &
     &          ipoasy(npogeo+1),indasy(lenasy),                        &
     &          ibndgo(npogeo),inflor(ninfpo),ievlor(nevlpo),           &
     &          neigeo(npoflf),knoobf(npogeo)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2)
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
      dimension xyzgeo(npogeo,ndmgeo),difvec(ndmgeo,maxnei),            &
     &          dipmat(ndimmt),diprhs(maxnei),                          &
     &          dipsol(ndmgeo*norder),vecbgo(npogeo),                   &
     &          dipole(npoflf,ndmflf),sysmat(lengeo),                   &
     &          sysrng(lengeo),                                         &
     &          rbdgeo(npogeo),vecsol(npogeo),vecrgo(npogeo),           &
     &          vecpgo(npogeo),dprodg(npogeo),dkondg(npogeo),           &
     &          gstifc(lengeo),vecrho(npogeo),vecref(npogeo),           &
     &          flxmat(nevlpo,ninfpo,ndmflx),syscpy(lengeo),            &
     &          volmat(lengeo),srcgeo(npogeo),dipreg(maxnei)
      dimension surmat(lengeo),alfmat(lengeo),srcsur(npogeo),           &
     &          vcbamb(npogeo),vcnflf(npoflf,ndmflf),                   &
     &          xyzflf(npoflf,ndmflf),potsng(npogeo)
      dimension veceig(npogeo,nvala),redinv(nmodes)
      dimension intgam(nelmeg,nummeg)
      dimension fstgam(mvagam,mgagam,nelmeg,nummeg),potmeg(nummeg),     &
     &          sekmeg(npogeo,nummeg)
!
      cpumat=qscput(1,0,ierr)
      if (ldipol) then
         call iset(npoflf*ndmflf, 0,mrkdip,1)
         call dset(npoflf*ndmflf,z0,dipole,1)
      else
         call dset(npogeo,z0,srcgeo,1)
         call dset(npogeo,z0,srcsur,1)
      end if
!
!--determine "vecref" (remember: lead field matrices are gradient operators!)
      call dcopy(lengeo,sysmat,1,syscpy,1)
      iknflf = -1
      call forwrd(ipogeo,indgeo,mrkvec,                                 &
     &     ipoasy,indasy,sysasy,ibndgo,mrkdip,neigeo,intgam,            &
     &     itypmg,megtyp,knoobf,knosur,itysur,necsur,                   &
     &     xyzgeo,difvec,dipmat,diprhs,                                 &
     &     dipsol,vecbgo,dipole,potsng,                                 &
     &     syscpy,rbdgeo,vecref,vecrgo,                                 &
     &     vecpgo,dprodg,dkondg,gstifc,                                 &
     &     vecrho,volmat,srcgeo,dipreg,                                 &
     &     surmat,alfmat,srcsur,vcbamb,                                 &
     &     xyzflf,veceig,redinv,sysrng,fstgam,                          &
     &     potmeg,sekmeg,iknflf)
!
      call percen(0,ninfpo)
      do 10 inf=1,ninfpo
         iknflf=inflor(inf)
         cpupoi=qscput(2,0,ierr)
         do 20 idim=1,ndmflx
            call inflod(mrkdip,dipole,srcgeo,srcsur,vcnflf,             &
     &           iknflf,idim,1)
            call dcopy(lengeo,sysmat,1,syscpy,1)
            call forwrd(ipogeo,indgeo,mrkvec,                           &
     &           ipoasy,indasy,sysasy,ibndgo,mrkdip,neigeo,intgam,      &
     &           itypmg,megtyp,knoobf,knosur,itysur,necsur,             &
     &           xyzgeo,difvec,dipmat,diprhs,                           &
     &           dipsol,vecbgo,dipole,potsng,                           &
     &           syscpy,rbdgeo,vecsol,vecrgo,                           &
     &           vecpgo,dprodg,dkondg,gstifc,                           &
     &           vecrho,volmat,srcgeo,dipreg,                           &
     &           surmat,alfmat,srcsur,vcbamb,                           &
     &           xyzflf,veceig,redinv,sysrng,fstgam,                    &
     &           potmeg,sekmeg,iknflf)
            do 30 jkno=1,numeeg
               knoglo=ievlor(jkno)
               flxmat(jkno,inf,idim)= vecsol(knoglo)-vecref(knoglo)
   30       continue
            if (logmeg) then
               do 40 j=1,nummeg
                  jpos=numeeg+j
                  flxmat(jpos,inf,idim)=potmeg(j)
   40          continue
            end if
            call inflod(mrkdip,dipole,srcgeo,srcsur,vcnflf,             &
     &           iknflf,idim,-1)
   20    continue
         call percen(inf,ninfpo)
   10 continue
      write(*,970) qscput(1,1,ierr)
!
      return
 970  format(1x,' Influence matrix completed.',                         &
     &          ' CPU-Time needed: ',e12.5,' [s]')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dpinck(ipogeo,indgeo,mrkvec,                           &
     &                  ipoasy,indasy,ibndgo,inflno,                    &
     &                  ievlno,nodevl,nodinf,neigeo,                    &
     &                  intgam,itypmg,megtyp,knoobf,                    &
     &                  knosur,itysur,necsur,                           &
     &                  xyzgeo,difvec,dipmat,diprhs,                    &
     &                  dipsol,vecbgo,mrkdip,dipole,                    &
     &                  sysmat,rbdgeo,vecsol,vecrgo,                    &
     &                  vecpgo,dprodg,dkondg,gstifc,                    &
     &                  vecrho,flxmat,volmat,srcgeo,                    &
     &                  dipreg,potinf,vecref,syscpy,                    &
     &                  surmat,alfmat,srcsur,vcbamb,                    &
     &                  vcnflf,xyzflf,veceig,redinv,                    &
     &                  fstgam,potmeg,sekmeg,potsng)
!
!---Check influence matrix solution versus forward solution
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkvec(maxnei),           &
     &          mrkdip(npoflf,ndmflf),                                  &
     &          ipoasy(npogeo+1),indasy(lenasy),knoobf(npogeo),         &
     &          ibndgo(npogeo),inflno(npogeo),ievlno(npogeo),           &
     &          nodevl(nevlpo),nodinf(ninfpo),neigeo(npoflf)
      dimension megtyp(numtyp,mxloop+1),itypmg(nummeg,2)
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
      dimension xyzgeo(npogeo,ndmgeo),difvec(ndmgeo,maxnei),            &
     &          dipmat(ndimmt),diprhs(maxnei),                          &
     &          dipsol(ndmgeo*norder),vecbgo(npogeo),                   &
     &          dipole(npoflf,ndmflf),sysmat(lengeo),                   &
     &          rbdgeo(npogeo),vecsol(npogeo),vecrgo(npogeo),           &
     &          vecpgo(npogeo),dprodg(npogeo),dkondg(npogeo),           &
     &          gstifc(lengeo),vecrho(npogeo),                          &
     &          flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          volmat(lengeo),srcgeo(npogeo),dipreg(maxnei),           &
     &          potinf(nevlpo),vecref(npogeo),syscpy(lengeo)
      dimension surmat(lengeo),alfmat(lengeo),srcsur(npogeo),           &
     &          vcbamb(npogeo),vcnflf(npoflf,ndmflf),                   &
     &          xyzflf(npoflf,ndmflf),potsng(npogeo)
      dimension veceig(npogeo,nvala),redinv(nmodes)
      dimension intgam(nelmeg,nummeg)
      dimension fstgam(mvagam,mgagam,nelmeg,nummeg),potmeg(nummeg),     &
     &          sekmeg(npogeo,nummeg)
!
      data idum /-1/
!
      if (ldipol) then
         call iset(npoflf*ndmflf, 0,mrkdip,1)
         call dset(npoflf*ndmflf,z0,dipole,1)
      else
         call dset(npogeo,z0,srcgeo,1)
         call dset(npogeo,z0,srcsur,1)
      end if
!
!---select random dipole load for 10 nodes
      write(*,'(a)') 'Selection of random load for nodes: '
      numran=min(ninfpo-1,10)
      do 10 i=1,numran
  222    iload=1+int( dble( ninfpo-1 ) * ran2(idum) )
         iflf=nodinf(iload)
!
!---if node already selected ---> goto 222
         if (ldipol) then
            if(mrkdip(iflf,1).gt.0) goto 222
         else
            if (invway.eq.1) then
               if(abs(srcgeo(iflf)).gt.tolsol) goto 222
            else if (invway.eq.2) then
               if(abs(srcsur(iflf)).gt.tolsol) goto 222
            end if
         end if
         write(zeil,'(i8)') iflf 
         call qdtext(zeil(1:8),1)
         if(ldipol) then
            if (lnorco) then
               srcval=ranknu(idum)
               do 13 k=1,ndmflf
                  mrkdip(iflf,k)=1
                  dipole(iflf,k)=srcval*vcnflf(iflf,k)
   13          continue
            else
               do 15 k=1,ndmflx
                  mrkdip(iflf,k)=1
                  dipole(iflf,k)=ranknu(idum)
   15          continue
            end if
         else
            if (invway.eq.1) srcgeo(iflf)=ranknu(idum)
            if (invway.eq.2) srcsur(iflf)=ranknu(idum)
         end if
   10 continue
      write(*,'(a)')
      write(*,'(a)') '  Random load selected'
!     
!---determine regular forward solution
      call dcopy(lengeo,sysmat,1,syscpy,1)
      iknflf = 0
      call forwrd(ipogeo,indgeo,mrkvec,                                 &
     &     ipoasy,indasy,sysasy,ibndgo,mrkdip,                          &
     &     neigeo,intgam,itypmg,megtyp,knoobf,                          &
     &     knosur,itysur,necsur,                                        &
     &     xyzgeo,difvec,dipmat,diprhs,                                 &
     &     dipsol,vecbgo,dipole,potsng,                                 &
     &     syscpy,rbdgeo,vecsol,vecrgo,                                 &
     &     vecpgo,dprodg,dkondg,gstifc,                                 &
     &     vecrho,volmat,srcgeo,dipreg,                                 &
     &     surmat,alfmat,srcsur,vcbamb,                                 &
     &     xyzflf,veceig,redinv,sysrng,fstgam,                          &
     &     potmeg,sekmeg,iknflf)
!
!---determine solution with influence matrix
      do 20 ievl=1,nevlpo
         sum=z0
         do 30 infl=1,ninfpo
            iflinf=nodinf(infl)
            if (ldipol) then
               if (lnorco) then
                  vallod=z0
                  do 35 idim=1,ndmflf
                     vallod=vallod+dipole(iflinf,idim)*                 &
     &                             dipole(iflinf,idim)
   35             continue
                  vallod=sqrt(vallod)
                  sum=sum+flxmat(ievl,infl,ndmflx)*vallod
               else
                  do 40 idim=1,ndmflf
                     sum=sum+flxmat(ievl,infl,idim)*dipole(iflinf,idim)
   40             continue
               end if
            else
               if(invway.eq.1)sum=sum+flxmat(ievl,infl,1)*srcgeo(iflinf)
               if(invway.eq.2)sum=sum+flxmat(ievl,infl,1)*srcsur(iflinf)
            end if
   30    continue
         if (ievl.le.numeeg) then
            potinf(ievl)=sum+vecref(ievl)
         else
            potinf(ievl)=sum
         end if
   20 continue
!
!---compare the solutions
      funeeg=z0
      do 50 ievl=1,numeeg
         iglo=nodevl(ievl)
         differ=potinf(ievl)-vecsol(iglo)
         funeeg=funeeg+differ*differ
!         if (abs(differ).gt.tol) then
            write(*,900) iglo,vecsol(iglo),potinf(ievl)
!         end if
   50 continue
      funeeg=sqrt(funeeg)
      funmeg=z0
      if (logmeg) then
         do 60 ievl=1,nummeg
            differ=potinf(ievl+numeeg)-potmeg(ievl)
            funmeg=funmeg+differ*differ
!            if (abs(differ).gt.tolmeg) then
               write(*,905) ievl,potmeg(ievl),potinf(ievl+numeeg)
!            end if
   60    continue
         funmeg=sqrt(funmeg)
      end if
      write(*,910) funeeg,funmeg
!
      return
  900 format(' Node ',i6,' Fwd ',g12.5,' Infl ',g12.5)
  905 format(' MEG  ',i6,' Fwd ',g12.5,' Infl ',g12.5)
  910 format(' Euclidean Difference  Fwd-Infl: ',/,                     &
     &       ' EEG-norm is: ',g12.5,/,                                  &
     &       ' MEG-norm is: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine forwrd(ipogeo,indgeo,mrkvec,                           &
     &                  ipoasy,indasy,sysasy,ibndgo,mrkdip,             &
     &                  neigeo,intgam,itypmg,megtyp,                    &
     &                  knoobf,knosur,itysur,necsur,                    &
     &                  xyzgeo,difvec,dipmat,diprhs,                    &
     &                  dipsol,vecbgo,dipole,potsng,                    &
     &                  sysmat,rbdgeo,vecsol,vecrgo,                    &
     &                  vecpgo,dprodg,dkondg,gstifc,                    &
     &                  vecrho,volmat,srcgeo,dipreg,                    &
     &                  surmat,alfmat,srcsur,vcbamb,                    &
     &                  xyzflf,veceig,redinv,sysrng,                    &
     &                  fstgam,potmeg,sekmeg,                           &
     &                  iknflf)
!
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkvec(maxnei),           &
     &          mrkdip(npoflf,ndmflf),                                  &
     &          ipoasy(npogeo+1),indasy(lenasy),                        &
     &          ibndgo(npogeo),neigeo(npoflf),knoobf(npogeo)
      dimension xyzgeo(npogeo,ndmgeo),difvec(ndmgeo,maxnei),            &
     &          dipmat(ndimmt),diprhs(maxnei),                          &
     &          dipsol(ndmgeo*norder),vecbgo(npogeo),                   &
     &          dipole(npoflf,ndmflf),sysmat(lengeo),                   &
     &          sysrng(lengeo),sysasy(lenasy),                          &
     &          rbdgeo(npogeo),vecsol(npogeo),vecrgo(npogeo),           &
     &          vecpgo(npogeo),dprodg(npogeo),dkondg(npogeo),           &
     &          gstifc(lengeo),vecrho(npogeo),volmat(lengeo),           &
     &          srcgeo(npogeo),dipreg(maxnei),potsng(npogeo)
      dimension surmat(lengeo),alfmat(lengeo),srcsur(npogeo),           &
     &          vcbamb(npogeo),xyzflf(npoflf,ndmflf)
      dimension veceig(npogeo,nvala),redinv(nmodes)
      dimension intgam(nelmeg,nummeg),megtyp(numtyp,mxloop+1),          &
     &          itypmg(nummeg,2)
      dimension fstgam(mvagam,mgagam,nelmeg,nummeg),potmeg(nummeg),     &
     &          sekmeg(npogeo,nummeg)
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
!
      cpulod=qscput(4,0,ierr)
      call dset(npogeo,z0,vecbgo,1)
      if (logmeg) call dset(nummeg,z0,potmeg,1)
      if (lrango) call dset(npogeo,z0,potsng,1)
!
!---if dipoles shall be modelled
      if (ldipol) then
         potint=z0
!
!---determine desired dipole moment load combination for normal operation
         if (iknflf.eq.0) then
            do 10 iknflf=1,npoflf
               icount=0
               do 30 k=1,ndmflf
                  icount=icount+mrkdip(iknflf,k)
                  if (mrkdip(iknflf,k).eq.0) dipole(iknflf,k)=z0
   30          continue
               if (icount.gt.0) then
                  if (lrango) then
                     call avefey(knosur,itysur,necsur,                  &
     &                    xyzgeo,xyzflf,dipole,potsng,                  &
     &                    potint,iknflf)
                  else
                     call soldip(dipole,dipsol,iknflf)
!     
!---determine load vector for system of equations
                     call diplod(ipoasy,indasy,mrkvec,                  &
     &                    ipogeo,indgeo,neigeo,                         &
     &                    xyzgeo,difvec,dipmat,diprhs,                  &
     &                    dipsol,vecbgo,dipreg,                         &
     &                    xyzflf,iknflf,numnei)
                  end if
                  if (logmeg) call primeg(intgam,fstgam,itypmg,megtyp,  &
     &                 xyzflf,dipole,potmeg,iknflf)
               end if
   10       continue
         else if (iknflf.gt.0) then
!
!---determine load vector for system of equations for influence matrices (faster)
!---first order
            if (lrango) then
               call avefey(knosur,itysur,necsur,                        &
     &              xyzgeo,xyzflf,dipole,potsng,                        &
     &              potint,iknflf)
            else
               call soldip(dipole,dipsol,iknflf)
               call diplod(ipoasy,indasy,mrkvec,                        &
     &              ipogeo,indgeo,neigeo,                               &
     &              xyzgeo,difvec,dipmat,diprhs,                        &
     &              dipsol,vecbgo,dipreg,                               &
     &              xyzflf,iknflf,numnei)
            end if
            if (logmeg) call primeg(intgam,fstgam,itypmg,megtyp,        &
     &           xyzflf,dipole,potmeg,iknflf)
         end if
         if (logsur) pomisg=potint/volume
      else
!
!---determine monopole loads 
         if (lnodlo) then
            call dcopy(npogeo,srcgeo,1,vecbgo,1)
            call daxpy(npogeo,z1,srcsur,1,vecbgo,1)
            call dscal(npogeo,-z1,vecbgo,1)
         else
            call matpro(volmat,srcgeo,vecbgo,ipogeo,indgeo,             &
     &           lengeo,npogeo)
            call dscal(npogeo,-z1,vecbgo,1)
!     
!---set Cauchy boundary conditions 
            call matpro(surmat,srcsur,dprodg,ipogeo,indgeo,             &
     &           lengeo,npogeo)
            call daxpy(npogeo,-z1,dprodg,1,vecbgo,1)
!---set Robbins (spring-type) boundary conditions (heat transfer conditions)
         end if
         if (lalpha) then
            call daxpy(lengeo, z1,alfmat,1,sysmat,1)
            call daxpy(npogeo,-z1,vcbamb,1,vecbgo,1)
         end if
      end if
!
      if (lrango) then
         call matpro(sysrng,potsng,dprodg,ipogeo,indgeo,                &
     &        lengeo,npogeo)
         do 301 i=1,npogeo
            if (knoobf(i).eq.0) then
               vecbgo(i)=vecbgo(i)+dprodg(i)
            end if
  301    continue
         call matpro(sysmat,potsng,dprodg,ipogeo,indgeo,                &
     &        lengeo,npogeo)
         do 300 i=1,npogeo
            if (knoobf(i).gt.0) then
               vecbgo(i)=vecbgo(i)+dprodg(i)
            end if
  300    continue
      end if
!
      cpulod=qscput(4,1,ierr)
      cpusol=qscput(4,0,ierr)
!
!---for modal operation
      if (lmocgs.or.lmofwd) then
         call dset(npogeo,z0,vecsol,1)
         do 200 i=1,nmodes
            potred=-ddot(npogeo,veceig(1,i),1,vecbgo,1)*redinv(i)
            call daxpy(npogeo,potred,veceig(1,i),1,vecsol,1)
  200    continue
      end if
!
!---for normal operation
      if (.not.lmofwd) then
!     
!---consider dirichlet boundary conditions
         call bndcon(ipogeo,indgeo,ibndgo,sysmat,                       &
     &        vecbgo,rbdgeo,    npogeo,   lengeo)
         avelod=dsum(npogeo,vecbgo,1)
!         write(*,5000) avelod
! 5000    format(' av-b ',g12.5)
!
!---solve system of equations
         call solver(ipogeo,indgeo,                                     &
     &        sysmat,vecbgo,vecsol,vecrgo,                              &
     &        vecpgo,dprodg,dkondg,                                     &
     &        gstifc,vecrho,                                            &
     &        npogeo,lengeo,tolsol,isolvr,resiun)
         cpusol=qscput(4,1,ierr)
         if (.not.linmat.and.ldipol) then
            write(*,900) cpulod/max(cpulod+cpusol,wuzeps)*z100
            write(*,920) resiun
         end if
      end if
!
!-----For AMG-Krylov-test-purpose: Write out the equation system. 
!-----vecrgo is only used as an auxiliary vector 
!-----C.Wolters, 24.5.2000 
      if ((isolvr.eq.3).or.(isolvr.eq.4)) then
         write(*,*) 'Writing out full equation system for fast solver: '
         call wriequ(filnam(43),filnam(44),ipogeo,indgeo,sysmat,vecsol, &
     &                vecbgo,vecrgo,ipoasy,indasy,sysasy)
      end if
!
!---if average potential shall be zero
      if (lmncor.or.logneu) then
         call matpro(volmat,potsng,dprodg,                              &
     &        ipogeo,indgeo,lengeo,npogeo)
         pomis2=dsum(npogeo,dprodg,1)/volume
         call matpro(volmat,vecsol,dprodg,                              &
     &        ipogeo,indgeo,lengeo,npogeo)
         pomico=dsum(npogeo,dprodg,1)/volume
         if(.not.logneu) write(*,910) pomico
         do 800 i=1,npogeo
            vecsol(i)=vecsol(i)-pomico
            if(lrango) then
               if (logsur) then
                  potsng(i)=potsng(i)-pomisg
               else
                  potsng(i)=potsng(i)-pomis2
               end if
            end if
  800    continue
         if(logsur) then
            potcor=pomico+pomisg
         else
            potcor=pomico+pomis2
         end if
      end if
!
      if (lrango) then
         call daxpy(npogeo,z1,potsng,1,vecsol,1)
      end if
      call matpro(volmat,vecsol,dprodg,                                 &
     &     ipogeo,indgeo,lengeo,npogeo)
      pomite=dsum(npogeo,dprodg,1)/volume
!
!--for meg-treatment (dipoles only!)
      if (logmeg.and.ldipol) then
         do 600 iakmeg=1,nummeg
!---secondary magnetic flux
            valsec=-ddot(npogeo,sekmeg(1,iakmeg),1,vecsol,1)
            potmeg(iakmeg)=potmeg(iakmeg)+valsec
  600    continue
      end if
!
      return
  900 format('Ratio of load determination/total solution time '         &
     &     ,g12.5,' [%]')
  910 format(' FEM-potentials corrected by the average value ',g12.5)
  920 format(' Residual in unscaled system: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!-----The following subroutine calculates the right-hand-side of the 
!-----FE equation system.
!-----The blurred Buchner97-dipole-model is used and it is assumed that 
!-----all the dipole locations are grid-nodes. 
!-----The "source-grid" can be 
!-----   a) FE-grid itself (xyzflf=xyzgeo, npoflf=npogeo, ndmflf=ndmgeo)
!-----   b) surface grid 
!-----The dipoles are marked on the grid through "mrkdip" (positions) 
!-----and "dipole" (strengths).
!-----The variable "mrkvec" will be used to implement a dipole as 
!-----monopole loads on neighboring nodes in "diplod" and "neibor".
      subroutine rhsblu(ipogeo,indgeo,mrkvec,                           &
     &                  ipoasy,indasy,mrkdip,neigeo,                    &
     &                  xyzgeo,difvec,dipmat,diprhs,                    &
     &                  dipsol,vecbgo,dipole,dipreg,                    &
     &                  xyzflf,iknflf,numnei)
!
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkvec(maxnei),           &
     &          mrkdip(npoflf,ndmflf),                                  &
     &          ipoasy(npogeo+1),indasy(lenasy),                        &
     &          neigeo(npoflf)
      dimension xyzgeo(npogeo,ndmgeo),difvec(ndmgeo,maxnei),            &
     &          dipmat(ndimmt),diprhs(maxnei),dipsol(ndmgeo*norder),    &
     &          vecbgo(npogeo),dipole(npoflf,ndmflf),dipreg(maxnei)
      dimension xyzflf(npoflf,ndmflf)
!
      call dset(npogeo,z0,vecbgo,1)
      potint=z0
!
      if (iknflf.eq.0) then
!--------In this case, the dipoles marked in "mrkdip" will be implemented 
!--------alltogether in the right-hand-side-vector for a source simulation
         do 10 iknflf=1,npoflf
            icount=0
            do 30 k=1,ndmflf
               icount=icount+mrkdip(iknflf,k)
               if (mrkdip(iknflf,k).eq.0) dipole(iknflf,k)=z0
   30       continue
            if (icount.gt.0) then
               call soldip(dipole,dipsol,iknflf)
!     
               call diplod(ipoasy,indasy,mrkvec,                        &
     &              ipogeo,indgeo,neigeo,                               &
     &              xyzgeo,difvec,dipmat,diprhs,                        &
     &              dipsol,vecbgo,dipreg,xyzflf,                        &
     &              iknflf,numnei)
            end if
   10    continue
      else if (iknflf.gt.0) then
!--------In this case, only one dipole "iknflf" will be implemented 
!--------in the right-hand-side-vector for the generation of an
!--------influence matrix
            call soldip(dipole,dipsol,iknflf)
            call diplod(ipoasy,indasy,mrkvec,                           &
     &           ipogeo,indgeo,neigeo,                                  &
     &           xyzgeo,difvec,dipmat,diprhs,                           &
     &           dipsol,vecbgo,dipreg,xyzflf,                           &
     &           iknflf,numnei)
      end if

      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!-----The following subroutine calculates the right-hand-side of the 
!-----FE equation system.
!-----The analytical von-Rango-dipole-model is used and it is assumed that 
!-----all the dipole locations are grid-nodes. 
!-----The "source-grid" can be 
!-----   a) FE-grid itself (xyzflf=xyzgeo, npoflf=npogeo, ndmflf=ndmgeo)
!-----   b) surface grid 
!-----The dipoles are marked on the grid through "mrkdip" (positions) 
!-----and "dipole" (strengths)
      subroutine rhsrng(ipogeo,indgeo,mrkdip,                           &
     &                  knoobf,knosur,itysur,necsur,                    &
     &                  xyzgeo,vecbgo,dipole,potsng,                    &
     &                  sysmat,dprodg,xyzflf,sysrng,                    &
     &                  iknflf)
!
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkdip(npoflf,ndmflf),    &
     &          knoobf(npogeo)
      dimension xyzgeo(npogeo,ndmgeo),vecbgo(npogeo),                   &
     &          dipole(npoflf,ndmflf),dprodg(npogeo),potsng(npogeo)
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
      dimension xyzflf(npoflf,ndmflf)
      dimension sysrng(lengeo),sysmat(lengeo)     
!
      call dset(npogeo,z0,vecbgo,1)
      call dset(npogeo,z0,potsng,1)
      potint=z0
!
         if (iknflf.eq.0) then
!--------In this case, the dipoles marked in "mrkdip" will be implemented 
!--------alltogether in the right-hand-side-vector for a source simulation
            do 10 iknflf=1,npoflf
               icount=0
               do 30 k=1,ndmflf
                  icount=icount+mrkdip(iknflf,k)
                  if (mrkdip(iknflf,k).eq.0) dipole(iknflf,k)=z0
   30          continue
               if (icount.gt.0) then
                  call avefey(knosur,itysur,necsur,                     &
     &                 xyzgeo,xyzflf,dipole,potsng,                     &
     &                 potint,iknflf)
               end if
   10       continue
         else if (iknflf.gt.0) then
!--------In this case, only one dipole "iknflf" will be implemented 
!--------in the right-hand-side-vector for the generation of an
!--------influence matrix
            call avefey(knosur,itysur,necsur,                           &
     &           xyzgeo,xyzflf,dipole,potsng,                           &
     &           potint,iknflf)
         end if
!
!-------- Roberts first implementation just works for spheres, not for 
!-------- general headmodels. It has to be replaced by a surface-integral.
!-------- Therefore, the following
         call matpro(sysrng,potsng,dprodg,ipogeo,indgeo,                &
     &        lengeo,npogeo)
         do 301 i=1,npogeo
               vecbgo(i)=vecbgo(i)+dprodg(i)
  301    continue
!-------- replaces the former implementation
!         call matpro(sysrng,potsng,dprodg,ipogeo,indgeo,                &
!     &        lengeo,npogeo)
!         do 301 i=1,npogeo
!            if (knoobf(i).eq.0) then
!               vecbgo(i)=vecbgo(i)+dprodg(i)
!            end if
!  301    continue
!         call matpro(sysmat,potsng,dprodg,ipogeo,indgeo,                &
!     &        lengeo,npogeo)
!         do 300 i=1,npogeo
!            if (knoobf(i).gt.0) then
!               vecbgo(i)=vecbgo(i)+dprodg(i)
!            end if
!  300    continue

      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine solver(ipodia,indexj,gstif,vecb,vecx,vecr,vecp,dprod,  &
     &                  dkond,gstifc,vecrho,npoin,lenge,tolsol,isolvr,  &
     &                  resiun)
      implicit double precision (a-h,o-z)
      parameter (z100=100.d0, z1=1.d0, z0=0.d0, mrepet=3)
!
      dimension ipodia(npoin),indexj(lenge)
!
      dimension  gstif(lenge),vecb(npoin),vecx(npoin),vecr(npoin),      &
     &           vecp(npoin),dprod(npoin),dkond(npoin),                 &
     &           gstifc(lenge),vecrho(npoin)
!
      save icount
      data icount /0/
!
      icount=icount+1
      nrepet=0
      iter=0;
      call scalen (gstif,vecx,vecb,ipodia,indexj,dkond,lenge,npoin)
!
      if((tolsol*tolsol).ge.1) goto 30       ! in case of singularity        
!    
 10   if (isolvr.eq.2) then
         cpucho=qscput(6,0,ier)
         call partch(gstif,ipodia,indexj,lenge,npoin,gstifc,ierr)
         cpucho=qscput(6,1,ier)
         if (ierr.ne.0) then
            write(*,920)
            isolvr=1
            goto 20
         end if
!
         call dscal (npoin,-z1,vecb,1)
         call solchs(gstifc,vecx,vecb,ipodia,indexj,lenge,npoin)         
         call dscal (npoin,-z1,vecb,1)
         call matpro(gstif,vecx,dprod,ipodia,indexj,lenge,npoin)
         call daxpy(npoin,z1,vecb,1,dprod,1)
         residu=dnrm2(npoin,dprod,1)
         if (icount.lt.2) write(*,910) 'Solver-partch: ',residu
!
         cpusol=qscput(5,0,ier)
         call pccgsy(gstif,vecx,vecb,ipodia,indexj,lenge,npoin,         &
     &        dprod,vecr,vecp,                                          &
     &        resi,tolsol,npoin,ierr,gstifc,vecrho)
         cpusol=qscput(5,1,ier)
         if (ierr.lt.0) then
            iter=abs(ierr)
         end if
      end if
 20   if (isolvr.eq.1) then
         cpusol=qscput(5,0,ier)
         call kograd (gstif,vecx,vecb,ipodia,indexj,lenge,npoin,        &
     &        dprod,vecr,vecp,resi,tolsol,iter,ierr)
         cpusol=qscput(5,1,ier)
      end if
 30   if (ierr.gt.0) then
         write(*,900) 'Warning: solver problems ',nrepet,'. Iteration'
         write(*,910) 'Solver-Fault: ',resi
         nrepet = nrepet + 1
         if (nrepet.lt.mrepet) then
            goto 10
         else
            stop 'solver not convergent'
         end if
      end if
      call unscal (gstif,vecx,vecb,ipodia,indexj,dkond,lenge,npoin)
      call matpro (gstif,vecx,dprod,ipodia,indexj,lenge,npoin)
      call daxpy  (npoin,z1,vecb,1,dprod,1)
      resiun=dnrm2(npoin,dprod,1)
!
!      if (icount.lt.2) then
         if (isolvr.eq.2) write(*,940) cpucho
          write(*,930) cpusol,iter,dble(iter)/dble(npoin)*z100
!      end if
!
      return
 900  format (a,i2,a)
 910  format (a,e23.16)
 920  format (' Partial cholesky decomposition not possible',           &
     &        ' changing solver method !!!!!')
 930  format ('. CPU-time for the solver: ',g12.5,/,                     &
     &        ' Number of steps for CG-solution: ',i12,/,               &
     &        ' Relative efficiency CG-solution: ',g12.5,' [%]')
 940  format ('. CPU-time for IC(0) partial decomposition: ',g12.5)
!
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine inrmom(xyzgeo,vecint)
      include 'neurofem.inc'
      dimension xyzgeo(npogeo,ndmgeo),vecint(npogeo)
!
      do 10 idim=1,ndmgeo
         xyzmom(idim)=z0
         do 20 i=1,npogeo
            radi=z0
            do 30 j=1,ndmgeo-1
               k=mod(idim+j-1,ndmgeo)+1
               dummy=xyzgeo(i,k)-xyzcgr(k)
               radi = radi+dummy*dummy
   30       continue
            xyzmom(idim) = xyzmom(idim) + radi * vecint(i)
   20    continue
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
