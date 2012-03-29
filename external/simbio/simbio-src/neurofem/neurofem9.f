!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem9.f :      	 Subroutines related to singular inverse 
!                                solutions
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

      subroutine anneal(nodinf,                                         &
     &                  mrkinv,ichold,ichnew,                           &
     &                  ichopt,infcho,nodevl,                           &
     &                  lapivo,                                         &
     &                  flxmat,potmes,potcal,                           &
     &                  anlmat,anlrhs,                                  &
     &                  dipinv,vcnflf,                                  &
     &                  workla,sarray,diptim,                           &
     &                  vtrval,sngval,covinv,                           &
     &                  devval,vecref)
      include 'neurofem.inc'
      dimension nodinf(ninfpo),nodevl(nevlpo),                          &
     &          mrkinv(npoflf,ndmflf)
      dimension dipinv(npoflf,ndmflf)
      dimension ichold(numdip),ichnew(numdip),ichopt(numdip),           &
     &          infcho(ninfpo),lapivo(numpvt)
      dimension flxmat(nevlpo,ninfpo,ndmflx),potmes(nevlpo,numtim),     &
     &          potcal(nevlpo,numtim),anlmat(nmanmt),                   &
     &          vcnflf(npoflf,ndmflf),sarray(numsar),                   &
     &          workla(lawork),anlrhs(nmanrs,numtim),                   &
     &          diptim(numsol,numtim),vecref(npogeo),                   &
     &          vtrval(numsol,numsol),sngval(numsol),                   &
     &          covinv(nevlpo),devval(npoflf)
!
!$$$rb SRCOUT
      dimension dipsed(100,3)
      logical srcout
      srcout = .true.
      call dset(300,z0,dipsed,1)
!$$$rb SRCOUT
!
      call iset(npoflf*ndmflf, 0,mrkinv,1)
      call dset(npoflf*ndmflf,z0,dipinv,1)
      call iset(ninfpo       , 0,infcho,1)
      call dset(numsol*numtim,z0,diptim,1)
      if (ldevsc) call dset(npoflf,z0,devval,1)
!
!---testing overdetermined system (two times is required)
      if ((2*numdip*ndmflx).gt.nevlpo) then
         write(*,*)
         write(*,*) 'Too many unknown components'
         write(*,*) 'Allowed unknowns ', nevlpo/(2*ndmflx)
         write(*,*) 'Chosen  unknowns ', numdip
         STOP
      end if
!
!---initialization of protocoll file
      iunit=qfreeu()
      call qfopse(iunit,filnam(26),'un','fo',ierr)
      if (ierr.gt.0) then
         write (*,900) 'error filnam(26) in anneal'
         stop
      end if
!
!---initialization of variables
      info=0
      iseed=-1
      funold=dble(ran2(iseed))
      funold=bignum
      funopt=bignum
      icount=0
      icoun2=0
      temper=tmpsta
      faksta=tmpfak**32
      isuccs=0
      thresh=wuzeps
      istart=0
      lblopt = 0
!
!---determine starting value of functional
      call select(ichold,ichnew,infcho,icount)
      call dset  (nmanrs*numtim,z0,anlrhs,1)
      call funcnl(ichold,anlrhs,flxmat,potmes,potcal,                   &
     &            covinv,anlrhs,funbez)
      write(*,1150) funbez
      if (temper.le.z0) temper=z10*funbez
!
!-----------------------------------------------------------------------
!                         START MAIN ITERATION
!-----------------------------------------------------------------------
      cpuanl=qscput(6,0,ierr)
  100 icount=icount+1
      icoun2=icoun2+1
!
!---randomly perturb the configuration
      call select(ichold,ichnew,infcho,icount)
!
!---build system of equations ... and solve it
      call anlsys(ichnew,flxmat,potmes,anlmat,anlrhs,covinv)
      if (lcholy) then
         call dppsv('u',nmanrs,numtim,anlmat,anlrhs,nmanrs,info)
      else if (lqrfak) then
         call dgels('n',nevlpo,numsol,numtim,anlmat,nevlpo,             &
     &        anlrhs,nevlpo,workla,lawork,info)
      else if (lqrfax) then
         call iset(numpvt,0,lapivo,1)
         call dgelsx(nevlpo,numsol,numtim,anlmat,nevlpo,                &
     &        anlrhs,nevlpo,lapivo,thresh,larank,workla,info)
      else if (lsngvl) then
         call cutsvd(anlmat,vtrval,sngval,potmes,anlrhs,workla,covinv)
      else 
         stop 'anneal: solver not specified'
      end if
      if (info.ne.0) then
         write(*,1080) info
         stop 'anneal: lapack'
      end if
!     
!---evaluate outcome
      if (lsngvl) then
         funnew=funcut
      else
         call funcnl(ichnew,anlrhs,flxmat,potmes,potcal,                &
     &        covinv,anlrhs,funnew)
      end if
!
!---in case of deviation scan for a single time point,
!---each solution is stored for further evalution
      if (ldevsc.and..not.ltrans) then
         j=nodinf(icount)
         do 700 k=1,ndmflx
            dipinv(j,k)=anlrhs(k,1)
  700    continue
         devval(j)=funnew
      end if
!
!---record optimum
      if (funnew.lt.funopt) then
         funopt=funnew
         call icopy(numdip,ichnew,1,ichopt,1)
         if (ltrans) call rectim(anlrhs,anlrhs,diptim)
         lblopt = 1
      end if
!
!---accepting configuration (Metropolis criterium)
      if (lmetro(funnew,funold).and..not.ldevsc) then
         isuccs=isuccs+1
         call icopy(numdip,ichnew,1,ichold,1)
         funold=funnew
      end if
!
! -- output
      if (.not.ldevsc) then
         if ((icoun2.eq.10*ninfpo).or.(isuccs.ge.ninfpo)) then
            relati=dble(isuccs)/dble(icoun2)
            icoun2=0
            isuccs=0
            cpustp=qscput(6,1,ierr)/dble(icount)
            dumold=(funbez-funold)/funbez*z100
            dumopt=(funbez-funopt)/funbez*z100
!
            if( srcout ) then
               call anlsys(ichopt,flxmat,potmes,anlmat,anlrhs,covinv)
               if (lcholy) then
                 call dppsv('u',nmanrs,numtim,anlmat,anlrhs,nmanrs,info)
               else if (lqrfak) then
                  call dgels('n',nevlpo,numsol,numtim,anlmat,nevlpo,    &
     &                 anlrhs,nevlpo,workla,lawork,info)
               else if (lqrfax) then
                  call iset(numpvt,0,lapivo,1)
                  call dgelsx(nevlpo,numsol,numtim,anlmat,nevlpo,       &
     &                 anlrhs,nevlpo,lapivo,thresh,larank,workla,info)
               else if (lsngvl) then
                  call cutsvd(anlmat,vtrval,sngval,potmes,              &
     &                        anlrhs,workla,covinv)
               end if
               if (info.ne.0) then
                  write(*,1080) info
                  stop 'anneal: lapack'
               end if
!
               if ( lblopt .eq. 1 ) then                     
                  write(*,900) icount,temper,dumold,                    &
     &                 relati,dumopt
                  write(iunit,1010) icount,temper,dumold,               &
     &                 relati,dumopt
                  do i=1,numdip
                     nodflf=nodinf(ichopt(i))
                     do l=1,ndmflf
                        dipsed(i,l)=z0
                     end do
                     dum = z0
                     if (lnorco) then
                        do l=1,ndmflf
                           dipsed(i,l)=anlrhs(i,1)*vcnflf(nodflf,l)
                           dum = dum + dipsed(i,l)*dipsed(i,l)
                        end do
                     else
                        do k=1,ndmflx
                           ipos=(i-1)*ndmflx+k
                           dipsed(i,k)=anlrhs(ipos,1)
                           dum = dum + dipsed(i,k)*dipsed(i,k)
                        end do
                     end if
                     iout = ndmflf
                     if( .not. ldipol ) iout = 1
                     write(*,961) nodflf,(dipsed(i,j),j=1,iout)
                     write(iunit,961) nodflf,(dipsed(i,j),j=1,iout)
                  end do
               end if
            end if
            if( .not. srcout .or. (srcout.and.lblopt.eq.0)) then
               write(*,900) icount,temper,dumold,                       &
     &              relati,dumopt,                                      &
     &              (nodinf(ichopt(k)),k=1,numdip)
               write(iunit,1010) icount,temper,dumold,                  &
     &              relati,dumopt,                                      &
     &              (nodinf(ichopt(k)),k=1,numdip)
            end if
            lblopt = 0
!
            if (istart.gt.0) then
               temper=temper*tmpfak
            else 
               if (relati.lt.0.90) then
                  temper=temper/faksta
                  istart=-1
               else if ((relati.gt.0.99).and.(istart.eq.0)) then
                  temper=temper*faksta
               else 
                  istart=1
               end if
            end if
            if (qfilda('stop_anneal_now')) lostop=.true.
         end if
      end if
!     
!---check convergence
      if (ldevsc.and.icount.eq.ninfpo) goto 2000
      if ((icount.le.maxtry).and.(funold.ge.tolinv)                     &
     &     .and.(.not.lostop)) goto 100
!      if (icount.gt.maxtry) then
!      else if (funold.lt.tolinv) then
!      else if (lostop) then
!      else
!         goto 100
!      end if
!
!-----------------------------------------------------------------------
! output selected dipoles
!      write(*,910) 
!      do 300 i=1,numdip
!         nodflf=nodinf(ichnew(i))
!         write(zeil,'(i8)') nodflf
!         call qdtext(zeil(1:8),1)
!  300 continue
!      write(*,920) funnew
!      write(*,'(/)')
!
!-----------------------------------------------------------------------
!                     END OF MAIN ITERATION
!-----------------------------------------------------------------------
!---record optimal solution
 2000 call anlsys(ichopt,flxmat,potmes,anlmat,anlrhs,covinv)
      if (lcholy) then
         call dppsv('u',nmanrs,numtim,anlmat,anlrhs,nmanrs,info) 
      else if (lqrfak) then
         call dgels('n',nevlpo,numsol,numtim,anlmat,nevlpo,             &
     &        anlrhs,nevlpo,workla,lawork,info)
      else if (lqrfax) then
         call iset(numpvt,0,lapivo,1)
         call dgelsx(nevlpo,numsol,numtim,anlmat,nevlpo,                &
     &        anlrhs,nevlpo,lapivo,thresh,larank,workla,info)
      else if (lsngvl) then
         call cutsvd(anlmat,vtrval,sngval,potmes,anlrhs,workla,covinv)
      end if
      if (info.ne.0) then
         write(*,1080) info
         stop 'anneal: lapack'
      end if
!
!---evaluate outcome
      call funcnl(ichopt,anlrhs,flxmat,potmes,                          &
     &     potcal,covinv,anlrhs,funop2)
      if (lsngvl) then
         if (funop2-funcut.gt.tol) then
            write(*,1140)
         end if
      end if
      if (abs(funop2-funopt).gt.tol) then
         write(*,930) funop2,funopt
      end if
!
      write(*,915) 
      do 315 i=1,numdip
         nodflf=nodinf(ichopt(i))
         write(zeil,'(i8)') nodflf
         call qdtext(zeil(1:8),1)
  315 continue
      write(*,920) funopt
      write(*,'(/)')
!
!$rb 28.8.97      if (.not.ltrans) then
         do 220 i=1,numdip
            nodflf=nodinf(ichopt(i))
            do 225 l=1,ndmflf
               dipinv(nodflf,l)=z0
  225       continue
            do 230 k=1,ndmflx
               if (lnorco) then
                  do 620 l=1,ndmflf
                     mrkinv(nodflf,l)=1
                     dipinv(nodflf,l)=anlrhs(i,1)*vcnflf(nodflf,l)
  620             continue
               else
                  ipos=(i-1)*ndmflx+k
                  mrkinv(nodflf,k)=1
                  dipinv(nodflf,k)=anlrhs(ipos,1)
               end if
  230       continue
            iout = ndmflf
            if( .not. ldipol ) iout = 1
            write(*,960) nodflf,(dipinv(nodflf,j),j=1,iout)
  220    continue
!$rb 28.8.97      else
!$rb 28.8.97         write(*,965) (ichopt(k),k=1,numdip)
!$rb 28.8.97         write(*,968) funopt
!$rb 28.8.97      end if
!
      if (ldevsc) then
         devmin=devval(nodinf(1))
         do 400 i=2,ninfpo
            iflf=nodinf(i)
            devmin=min((devval(iflf)+wuzeps),devmin)
  400    continue
         do 410 i=1,ninfpo
            iflf=nodinf(i)
            devval(iflf)=devmin/(devval(iflf)+wuzeps)*z100
  410    continue
      end if
!
      if (qfilda('stop_anneal_now')) then
         call qfdel('stop_anneal_now',ierr)
         if (ierr.gt.0) then 
            write(*,940) 
         end if
      end if
!
!$$$1.9.97      if (.not.ltrans) then
          do 800 i=1,nevlpo
            us = vecref(nodevl(i))
            if (i.le.numeeg) then
               write(*,970) nodevl(i),potmes(i,1)+us,potcal(i,1)+us
               write(iunit,1070) i,potmes(i,1)+us,potcal(i,1)+us,       &
     &                           nodevl(i)
            else
               iak=i-numeeg
               write(*,975) iak,potmes(i,1)+us,potcal(i,1)+us
               write(iunit,1075) i,potmes(i,1)+us,potcal(i,1)+us,iak
            end if
  800    continue
!$$$1.9.97      end if
!
      call qfclos(iunit,0)
!
      RETURN
  900 format('s: ',i9,' t:',e9.2,' ev:',f8.2,' x:',e9.2,                &
     &       ' ev_o:',f11.5,1x,100(1x,i10))
!$$$rb  910 format(' selected dipole nodes ')
  915 format(' optimum nodes ')
  920 format(' euclidean difference in potentials ',f12.5)
  930 format(' warning: funop2: ',e12.5,' funopt: ',e12.5)
  940 format(' warning: could not delete file: stop_anneal_now')
  950 format(' average (euclidean) potential: ',f12.5,/,                &
     &       ' error norm for the truncated svd: ',f12.5)
  960 format(' opt. node ',i8,' source: ',3(g12.5,1x,:))
  961 format(2x, 'opt. node ', i8,' source: ',3(g12.5,1x,:))
  965 format(' opt. nodes ',100(i8,:))
  968 format(' Funktional: ',e12.5)
  970 format(' EEG-node: ',i8,' measured: ',e12.5,' calculated: ',e12.5)
  975 format(' MEG-node: ',i8,' measured: ',e12.5,' calculated: ',e12.5)
 1000 format('# icount temper fun-old kappa funopt ')
!$$$rb 1010 format(1x,i8,4(1x,e12.5),100i5)
 1010 format(1x,i8,4(1x,e12.5),100i11)
 1050 format('# average (euclidean) potential: ',e12.5)
 1070 format(1x,i5,1x,e12.5,e12.5,i8)
 1075 format(1x,i5,1x,e12.5,e12.5,i8)
 1080 format(' lapack returned error code ',i8)
 1110 format('#',100(i8,:))
 1120 format(a)
 1130 format(i8,100(e12.5,:))
 1140 format(' Problem: funop2-funcut > tol! ')
 1150 format(' starting value of functional is: ',e12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine select(ichold,ichnew,infcho,icall)
      include 'neurofem.inc'
      dimension ichold(numdip),ichnew(numdip),infcho(ninfpo)
!
      if (ldevsc) then
         if (icall.gt.0) then
            ichold(1)=ichnew(1)
            ichnew(1)=icall
         else
            ichold(1)=1
         end if
         return
      end if
      if (icall.eq.0) then
         do 600 i=1,numdip
            ichold(i)=i
  600    continue
!
      else
         do 5 i=1,numdip
            ichnew(i)=ichold(i)
    5    continue
!
         ipos = 1 + int( dble(numdip)*dble(ran2(iseed)) )
  100    ival = 1 + int( dble(ninfpo)*dble(ran2(iseed)) )
!     
         do 10 i=1,numdip
            if (ichold(i).eq.ival) goto 100
   10    continue
!     
         ichnew(ipos)=ival
!
         do 50 i=1,numdip
            inf=ichnew(i)
            infcho(inf)=infcho(inf)+1
   50    continue
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine anlsys(ichose,flxmat,potmes,anlmat,anlrhs,covinv)
      include 'neurofem.inc'
      dimension ichose(numdip)
      dimension flxmat(nevlpo,ninfpo,ndmflx),potmes(nevlpo,numtim),     &
     &          anlmat(nmanmt),anlrhs(nmanrs,numtim),covinv(nevlpo)
!
      call dset(nmanmt,z0,anlmat,1)
      call dset(nmanrs*numtim,z0,anlrhs,1)
!
      if (lcholy) then
         do 10 ir=1,numdip
            i=ichose(ir)
            ip=(ir-1)*ndmflx
            do 20 is=1,ndmflx
               ipodof=ip+is
               ipodia=ipodof*(ipodof-1)/2
               do 30 iu=1,numdip
                  j=ichose(iu)
                  iq=(iu-1)*ndmflx
                  do 40 iv=1,ndmflx
                     jpodof=iq+iv
                     if (jpodof.gt.ipodof) goto 40
                     sum=z0
                     do 50 ievl=1,nevlpo
                        co=covinv(ievl)
                        sum=sum+flxmat(ievl,i,is)*flxmat(ievl,j,iv)*co
   50                continue
                     ipos=ipodia+jpodof
                     anlmat(ipos)=sum
   40             continue
   30          continue
               do 65 itim=1,numtim
                  sum=z0
                  do 60 ievl=1,nevlpo
                     co=covinv(ievl)
                     sum=sum+flxmat(ievl,i,is)*potmes(ievl,itim)*co
   60             continue
                  ipos=(ir-1)*ndmflx+is
                  anlrhs(ipos,itim)=sum
   65          continue
   20       continue
   10    continue
      else
         do 110 ir=1,numdip
            i=ichose(ir)   
            ip=(ir-1)*ndmflx                  
            do 120 is=1,ndmflx
               iq=ip+is    
               do 130 j=1,nevlpo
                  k=(iq-1)*nevlpo+j  
                  co=sqrt(covinv(j))
                  anlmat(k)=flxmat(j,i,is)*co
  130          continue
  120       continue
  110    continue
         do 140 i=1,numtim
            do 150 j=1,nevlpo
               co=sqrt(covinv(j))
               anlrhs(j,i)=potmes(j,i)*co
  150       continue
  140    continue
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine funcnl(ichose,anlrhs,flxmat,potmes,potcal,             &
     &           covinv,cutsol,funval)
      include 'neurofem.inc'
      dimension ichose(numdip)
      dimension flxmat(nevlpo,ninfpo,ndmflx),potmes(nevlpo,numtim),     &
     &          anlrhs(nmanrs,numtim),potcal(nevlpo,numtim),            &
     &          covinv(nevlpo),cutsol(numsol,numtim)
!
      funval=z0
      do 5 itim=1,numtim
         do 10 ievl=1,nevlpo
            potcal(ievl,itim)=z0
            do 20 i=1,numdip
               inf=ichose(i)
               ibase=(i-1)*ndmflx
               do 30 k=1,ndmflx
                  l=ibase+k
                  if (lsngvl) then
                     potcal(ievl,itim)=potcal(ievl,itim)+               &
     &                    flxmat(ievl,inf,k)*cutsol(l,itim)
                  else
                     potcal(ievl,itim)=potcal(ievl,itim)+               &
     &                    flxmat(ievl,inf,k)*anlrhs(l,itim)
                  end if
   30          continue
   20       continue
            dif = potcal(ievl,itim)-potmes(ievl,itim)
            funval=funval+dif*dif*covinv(ievl)
   10    continue
    5 continue
      funval=funval
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function lmetro(funnew,funold)
      include 'neurofem.inc'
!
      lmetro=.false.
      dif=funnew-funold
!
      if (dif.le.z0) then
         lmetro=.true.
      else
         ran=dble(ran2(iseed))
         val=exp(-dif/temper)
         if (ran.lt.val) lmetro=.true.
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!     random number generator from numerical recipes
!---->---1---------2---------3---------4---------5---------6---------7--<
      function ran2(idum)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,    &
     &ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,       &
     &ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
!
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
!
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cutsvd(usvval,vtrval,sngval,potmes,                    &
     &           cutsol,workla,covinv)
      include 'neurofem.inc'
      dimension usvval(nevlpo,numsol),vtrval(numsol,numsol),            &
     &          sngval(numsol),workla(lawork),potmes(nevlpo,numtim),    &
     &          cutsol(numsol,numtim),covinv(nevlpo)
!
      save icount
      save conmax
      data icount /0/
      data conmax /z0/
!
      chasvu='O'
      chasvv='S'
      call dgesvd(chasvu,chasvv,nevlpo,numsol,usvval,nevlpo,            &
     &            sngval,usvval,nevlpo,vtrval,numsol,workla,lawork,     &
     &     infsvd)
!
      if (infsvd.eq.0) then
         if (icount.eq.0) then
            write(*,910)
            write(*,915) lawork,int(workla(1)+f2)
         end if
      else
         write(*,920) 
         write(*,930) infsvd 
         if (infsvd.lt.0) write(*,940) -infsvd 
         if (infsvd.gt.0) write(*,950)  infsvd 
         stop 'dgesvd in cutsvd'
      end if
      icount=1
!
!---louis' algorithm for the truncated svd
      funcut=z0
      call dset(numsol*numtim,z0,cutsol,1)
      do 20 itim=1,numtim
         defect=z0
         do 5 i=1,nevlpo
            defect=defect+potmes(i,itim)*potmes(i,itim)*covinv(i)
    5    continue
         do 10 i=1,numsol
            if (defect.lt.rpostl*rpostl) goto 333
!%%%am 10.4.97 auskommentiert Robert Pohlmeier
!%%%            if (sngval(i).le.z1) goto 333
            alphan=z0
            do 30 j=1,nevlpo
               co=sqrt(covinv(j))
               alphan=alphan+usvval(j,i)*potmes(j,itim)*co
   30       continue
            faktor=alphan/sngval(i)
            do 40 j=1,numsol
               cutsol(j,itim)=cutsol(j,itim)+faktor*vtrval(i,j)
   40       continue
            defect=defect-alphan*alphan
   10    continue
  333    funcut=funcut+defect
   20 continue
!
      return
  910 format(' Singular value decomposition successfull ')
  915 format(' work array size for SVD was',i10,/,                      &
     &       ' optimum work array size  is',i10)
  920 format(' Singular value decomposition failed (!) ')
  930 format(' Error code of SVD is ',i10)
  940 format(' Argument # ',i10,' had an illegal value')
  950 format(' ',i10,' off-diagonal elements did not converge to zero')
  960 format(g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine rectim(anlrhs,cutsol,diptim)
      include 'neurofem.inc'
      dimension anlrhs(nmanrs,numtim),cutsol(numsol,numtim),            &
     &          diptim(numsol,numtim)
!
      if (lsngvl) then
         do 10 i=1,numtim
            do 20 j=1,numdip*ndmflx
               diptim(j,i)=cutsol(j,i)
   20       continue
   10    continue
      else
         do 30 i=1,numtim
            do 40 j=1,numdip*ndmflx
               diptim(j,i)=anlrhs(j,i)
   40       continue
   30    continue
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine restar(icount,funold,relati,funopt,funbez,             &
     &                  ichopt,ichold,nodinf,istart)
      include 'neurofem.inc'
!
      dimension nodinf(ninfpo)
      dimension ichold(numdip),ichopt(numdip)
!
      if (qfilda('continue_anneal_now')) then
         iunid=qfreeu()
         call qfopse(iunid,'continue_anneal_now','un','fo',ierr)
         if (ierr.gt.0) then
            write (*,'(a)') 'file-error by continue_anneal_now'
            stop 'continue_anneal_now'
         end if
         if (numdip.gt.100) then
            write (*,'(a)')                                             &
     &      'continueing annealing not possible, because numdip > 100'
            stop
         end if
         read(iunid,900) icount,temper,funold,relati,funopt,            &
     &                  (ichopt(k),k=1,numdip)
         call qfclos(iunid,0)
         do 1 k=1,numdip
            iglo=ichopt(k)
            ichk=0
            do 2 i=1,ninfpo
               if (iglo.eq.nodinf(i)) then
                  ichopt(k)=i
                  ichk=ichk+1
               end if
    2       continue
            if (ichk.ne.1) then
               write (*,'(a)') 'error at continue_anneal_now'
               stop
            end if               
    1    continue
         call icopy(numdip,ichopt,1,ichold,1)
         funopt=funopt*funbez
         funold=funopt
         istart=1
      end if
      return
  900 format('s: ',i9,' t:',e9.2,' rv:',f8.2,' x:',e9.2,                &
     &       ' rv_o:',f11.5,1x,100i5)
!
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
