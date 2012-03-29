!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem5.f :      	 inverse solver (conjugate gradients)
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

      subroutine cgrinv(ipoinv,indinv,invbnd,nodinf,                    &
     &                  stfinv,solinv,vcbinv,flxmat,                    &
     &                  dprdin,vcrinv,vcpinv,vechey,                    &
     &                  vechew,vechev,scainv,weiinv,                    &
     &                  covinv,potmes,  resi,  iter,                    &
     &                  irepet, ierr)
!
      include 'neurofem.inc'
!---choose desired accuracy: (see your machine's manuals for availability
!                             and definition of extended precision)
!-extended precision      real *16 resext,de,dra,drn,dnen,dq,q0
!-extended precision      parameter (q0=0.q0)
!-double   precision      real *8  resext,de,dra,drn,dnen,dq,q0
!-double   precision      parameter (q0=0.d0)
      real *8 resext,de,dra,drn,dnen,dq,q0
      parameter (q0=0.d0)
      dimension ipoinv(npoinv),indinv(leninv),invbnd(npoinv),           &
     &          nodinf(ninfpo)
      dimension stfinv(leninv,ndmflx),vcbinv(npoinv,ndmflx),            &
     &          solinv(npoinv,ndmflx),vcrinv(npoinv,ndmflx),            &
     &          dprdin(npoinv,ndmflx),vcpinv(npoinv,ndmflx),            &
     &          vechey(nevlpo),vechew(npoinv,ndmflx),                   &
     &          vechev(npoinv,ndmflx),scainv(npoinv,ndmflx),            &
     &          flxmat(nevlpo,ninfpo,ndmflx),weiinv(npoinv,ndmflx),     &
     &          covinv(nevlpo),potmes(nevlpo)
!
!---start
      cpustp = qscput(8,0,ierr)
      call matinv(invbnd,ipoinv,indinv,nodinf,                          &
     &            flxmat,solinv,vechey,vechew,                          &
     &            vechev,stfinv,dprdin,scainv,                          &
     &            weiinv,covinv)
!
      nout=10
      if (ninfpo.lt.10) nout=1
      drn=q0
      do 15 k=1,ndmflx
         do 10 i=1,npoinv
            drn=drn+vcbinv(i,k)*vcbinv(i,k)
            vcrinv(i,k)= dprdin(i,k)+vcbinv(i,k)
            vcpinv(i,k)=-vcrinv(i,k)
   10    continue
   15 continue
      res0=sqrt(dble(drn))
      res0=max(res0,z1)
      write(*,890) res0
!
!---iteration
      drn  = q0
      iter =  0
      res10= z0
  333 iter=iter+1
!
      dra = drn
      drn = q0
      do 25 k=1,ndmflx
         do 20 i=1,npoinv
            drn=drn+vcrinv(i,k)*vcrinv(i,k)
   20    continue
   25 continue
      resext=sqrt(drn)
      resi=dble(resext)
      res10=res10+resi
!
      cpustp = qscput(8,1,ierr)
      if (mod(iter,nout).eq.0) then
         res10=res10/z10
         write(*,900) iter,cpustp,res10
         res10=z0
      end if
      cpustp = qscput(8,0,ierr)
!
      if (resi.lt.tolinv*res0) then
         ierr=0
         goto 1000
      else if (iter.gt.ninfpo*ndmflx) then
         ierr=1
         goto 1000
      end if
!      
      if (iter.gt.1) then
         de = drn / dra
      else
         de = q0
      end if
!
      do 35 k=1,ndmflx
         do 30 i=1,npoinv
            vcpinv(i,k)=vcpinv(i,k)*de-vcrinv(i,k)
   30    continue
   35 continue
!
      call matinv(invbnd,ipoinv,indinv,nodinf,                          &
     &            flxmat,vcpinv,vechey,vechew,                          &
     &            vechev,stfinv,dprdin,scainv,                          &
     &            weiinv,covinv)
!
      dnen=q0
      do 45 k=1,ndmflx
         do 40 i=1,npoinv
            dnen=dnen+vcpinv(i,k)*dprdin(i,k)
   40    continue
   45 continue
      dq = drn / dnen
!
!---parameter dq yields minimum in present search direction
      do 55 k=1,ndmflx
         do 50 i=1,npoinv
            solinv(i,k) = solinv(i,k) + dq * vcpinv(i,k)
            vcrinv(i,k) = vcrinv(i,k) + dq * dprdin(i,k)
   50    continue
   55 continue
!
!%%%darf nur bei Rechnung mit Kovarianzmatrix wirken
      if (lcovar) then
!
!---louis' truncated conjugate gradient algorithm
         fkstop=funnrm(nodinf,invbnd,flxmat,solinv,                     &
     &        vcbinv,covinv,weiinv,potmes)
         if (fkstop.lt.rpostl) then
            write(*,910) iter
            goto 1000
         end if
      end if
!
      goto 333
!
 1000 return
  890 format(' Principal residual for the inverse CG solution: ',g12.5)
  900 format(' CG-iter: ',i6,' CG-Step: ',g12.5,' [s], Res10: ',g12.5)
  910 format(' Further CG Iterations make no sense, Niter: ',i6)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine matinv(invbnd,ipoinv,indinv,nodinf,                    &
     &                  flxmat,solinv,vechey,vechew,                    &
     &                  vechev,stfinv,dprdin,scainv,                    &
     &                  weiinv,covinv)
      include 'neurofem.inc'
!
      dimension invbnd(npoinv),ipoinv(npoinv),indinv(leninv),           &
     &          nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          solinv(npoinv,ndmflx),vechey(nevlpo),                   &
     &          vechew(npoinv,ndmflx),vechev(npoinv,ndmflx),            &
     &          stfinv(leninv,ndmflx),                                  &
     &          dprdin(npoinv,ndmflx),scainv(npoinv,ndmflx),            &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo)
!
      call dset(nevlpo,z0,vechey,1)
      do 20 j=1,ninfpo
         jglo=nodinf(j)
         do 130 k=1,ndmflx
            we=weiinv(jglo,k)
            vw=solinv(jglo,k)*we
            do 10 i=1,nevlpo
               fx=flxmat(i,j,k)*covinv(i)
               vechey(i)=vechey(i)+fx*vw
   10       continue
  130    continue
   20 continue
!
      do 35 k=1,ndmflx
         do 30 j=1,ninfpo
            jglo=nodinf(j)
            we=weiinv(jglo,k)
            sum=z0
            do 40 i=1,nevlpo
               sum = sum + flxmat(i,j,k)*vechey(i)*we
   40       continue
            vechew(jglo,k) = sum
   30    continue
   35 continue
!
      do 60 k=1,ndmflx
         call matpro(stfinv(1,k),solinv(1,k),dprdin(1,k),ipoinv,indinv, &
     &        leninv,npoinv)
!---please check! (4.9.1996)
         do 65 i=1,npoinv
            dprdin(i,k)=dprdin(i,k)*scainv(i,k)
   65    continue
!$$$         call daxpy(npoinv,z1,scainv(1,k),1,dprdin(1,k),1)
   60 continue
!
      denom=trasmo*tradep
      if (denom.lt.wuzeps) stop ' denom close to zero in matinv'
!$$      xlmakt=xlmbda*tradat/denom
      xlmakt=xlmbda
      do 70 k=1,ndmflx
         do 80 iglo=1,npoinv
            if (invbnd(iglo).lt.1) then
               dprdin(iglo,k)=vechew(iglo,k)+xlmakt*dprdin(iglo,k)
            else
               dprdin(iglo,k)=z0
            end if
   80    continue
   70 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine matrhs(invbnd,nodinf,flxmat,                           &
     &           solinv,vcbinv,potmes,weiinv,covinv)
      include 'neurofem.inc'
!
      dimension invbnd(npoinv),nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),solinv(npoinv,ndmflx),     &
     &          vcbinv(npoinv,ndmflx),potmes(nevlpo),                   &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo)
!
      call dset(npoinv*ndmflx,z0,vcbinv,1)
!
      tradat=z0
      do 5 k=1,ndmflx
         do 10 i=1,ninfpo
            iglo=nodinf(i)
            we=weiinv(iglo,k)
            sum=z0
            if (invbnd(iglo).lt.1) then
               do 20 j=1,nevlpo
                  sum = sum - flxmat(j,i,k)*potmes(j)                   &
     &                        *covinv(j)*we
                  tradat=tradat+flxmat(j,i,k)*flxmat(j,i,k)*covinv(j)
   20          continue
            end if
            vcbinv(iglo,k) = sum
   10    continue
    5 continue
      write(*,900) tradat
!
      return
  900 format(' Trace of the data term: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine depsca(ipoinv,invbnd,nodinf,indinv,                    &
     &     flxmat,scainv,stfinv,solinv,vcbinv,weiinv,covinv)
      include 'neurofem.inc'
      dimension invbnd(npoinv),nodinf(ninfpo),                          &
     &          ipoinv(npoinv),indinv(leninv)
      dimension flxmat(nevlpo,ninfpo,ndmflx),scainv(npoinv,ndmflx),     &
     &          stfinv(leninv,ndmflx),solinv(npoinv,ndmflx),            &
     &          vcbinv(npoinv,ndmflx),weiinv(npoinv,ndmflx),            &
     &          covinv(nevlpo)
!
!---preset depth scaling factors
      call dset(npoinv*ndmflx,z1,scainv,1)
      tradep=dble(npoinv*ndmflx)
      if (.not.ldepsc) return
!
!---determine depth scaling factors
!---see fuchs-wischmann-wagner ISBET Newsletter, 1994, (5) 8--11
      scamax=z0
      do 10 i=1,ninfpo
         iglo=nodinf(i)
         valmax=z0
         do 20 j=1,nevlpo
            sum=z0
            do 30 k=1,ndmflx
               sum = sum + flxmat(j,i,k)*flxmat(j,i,k)
   30       continue
            sum=sqrt(sum)
            valmax=max(valmax,sum)
   20    continue
         scainv(iglo,1)=valmax
         scamax=max(valmax,scamax)
   10 continue
      betsqr=scamax*scamax/(signoi*signoi)
!
      do 40 i=1,ninfpo
         iglo=nodinf(i)
         faktor=scainv(iglo,1)
         do 50 k=1,ndmflx
            scainv(iglo,k)=(faktor+betsqr/faktor)*(faktor+betsqr/faktor)
   50    continue
   40 continue
!
      scamax = -bignum
      scamin =  bignum
      tradep = z0
      do 60 i=1,npoinv
         val=scainv(i,1)
         do 65 k=1,ndmflx
            tradep=tradep+scainv(i,k)
   65    continue
         scamax=max(scamax,val)
         scamin=min(scamin,val)
   60 continue
      write(*,900) scamax
      write(*,910) scamin
      write(*,920) tradep
!
      return
  900 format(' Maximum scaling factor: ',g12.5)
  910 format(' Minimum scaling factor: ',g12.5)
  920 format(' Trace of the scaling operator ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function funnrm(nodinf,invbnd,flxmat,solinv,                      &
     &           vcbinv,covinv,weiinv,potmes)
      include 'neurofem.inc'
      dimension invbnd(npoinv),nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),solinv(npoinv,ndmflx),     &
     &          vcbinv(npoinv,ndmflx),potmes(nevlpo),                   &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo)
!
      funnrm=z0
      do 10 i=1,nevlpo
         sum=z0
         if (invbnd(i).gt.0) goto 40
         do 20 j=1,ninfpo
            jglo=nodinf(j)
            do 30 k=1,ndmflx
               sum=sum+flxmat(i,j,k)*solinv(jglo,k)
   30       continue
   20    continue
!         sum=(sum-potmes(i))*sqrt(covinv(i))
!         funnrm=funnrm+sum*sum
   40    sum=sum-potmes(i)
         funnrm=funnrm+sum*sum*covinv(i)
   10 continue
      funnrm=sqrt(funnrm)
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
