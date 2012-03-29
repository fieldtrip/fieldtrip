!  $2 21/03/2003  Anwander A.  released for SIMBIO
!  $1 30/06/00aa subroutins which uses work.inc moved to cauchy.f
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem11.f :      	 Subroutines related to nonlinear 
!                            regularization (entropy, L1-norm)
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

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     user supplied subroutines for evaluating functional and gradient
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$22/07/00aa copy of thease functions to cauchy.f to be complied with f90
!      function func(p)
!      include 'neurofem.inc'
!      include 'work.inc'
!      dimension p(npoinv,ndmflx)
!
!      func=funcm(ip_invbnd,ip_ipoinv,ip_indinv,ip_nodinf,               &
!     &           dp_flxmat,         p,dpvechey,dp_vechew,               &
!     &           dp_vechev,dp_stfinv,dp_scainv,                         &
!     &           dp_weiinv,dp_covinv,dp_dprdin,dp_potmes)
!
!      return
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!      subroutine dfunc(x,df)
!      include 'neurofem.inc'
!      include 'work.inc'
!      dimension x(npoinv,ndmflx),df(npoinv,ndmflx)
!
!$$$      call dfuncm(x,df,dp_sysmat,dp_vecbgo)
!      call dfuncm(ip_invbnd,ip_ipoinv,ip_indinv,ip_nodinf,              &
!     &            dp_flxmat,         x,dp_vechey,dp_vechew,             &
!     &            dp_vechev,dp_stfinv,        df,dp_scainv,             &
!     &            dp_weiinv,dp_covinv,dp_potmes) 
!
!      return
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!      subroutine funchk(p,dat,xl1,xl2,xlc,ent,fun)
!      include 'neurofem.inc'
!      include 'work.inc'
!      dimension p(npoinv,ndmflx)
!
!      call funchm(ip_invbnd,ip_ipoinv,ip_indinv,ip_nodinf,              &
!     &            dp_flxmat,         p,dp_vechey,dp_vechew,             &
!     &            dp_vechev,dp_stfinv,dp_scainv,                        &
!     &            dp_weiinv,dp_covinv,dp_dprdin,dp_potmes,              &
!     &            dat,xl1,xl2,xlc,ent,fun)
!
!      return
!      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function funcm(invbnd,ipoinv,indinv,nodinf,                       &
     &               flxmat,solinv,vechey,vechew,                       &
     &               vechev,stfinv,scainv,weiinv,                       &
     &               covinv,dprdin,potmes)
      include 'neurofem.inc'
      dimension invbnd(npoinv),ipoinv(npoinv),indinv(leninv),           &
     &          nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          solinv(npoinv,ndmflx),vechey(nevlpo),                   &
     &          vechew(npoinv,ndmflx),vechev(npoinv,ndmflx),            &
     &          stfinv(leninv,ndmflx),scainv(npoinv,ndmflx),            &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo),                   &
     &          dprdin(npoinv,ndmflx),potmes(nevlpo)
!
!---determine calculated potentials
      call dset(nevlpo,z0,vechey,1)
      do 20 j=1,ninfpo
         jglo=nodinf(j)
         do 130 k=1,ndmflx
            vw=solinv(jglo,k)
            do 10 i=1,nevlpo
               fx=flxmat(i,j,k)
               vechey(i)=vechey(i)+fx*vw
   10       continue
  130    continue
   20 continue
!
!---data part of the functional
      datafu=z0
      do 30 i=1,nevlpo
         dif=vechey(i)-potmes(i)
         co=covinv(i)
         datafu=datafu+dif*co*dif
   30 continue
!
!---regularization part of the functional
      if (lentro) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         gesint=z0
         do 200 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 210 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  210       continue
            streng=sqrt(streng)
            gesint=gesint+streng
  200    continue
         gesint=max(gesint,wuzeps)
         tiny=wuzeps*gesint
         regula=z0
         do 220 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 230 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  230       continue
            streng=sqrt(streng)/gesint
            if(streng.gt.z0) regula=regula-streng*log(streng)
  220    continue
      else if (ll1nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=z0
         do 400 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 410 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  410       continue
            streng=sqrt(streng)
            regula=regula+streng
  400    continue
      else if (ll2nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=z0
         do 300 i=1,ninfpo
            iglo=nodinf(i)
            do 310 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               regula=regula+so*so
  310       continue
  300    continue
      else
         regula=z0
         do 60 k=1,ndmflx
            call matpro(stfinv(1,k),solinv(1,k),dprdin(1,k),            &
     &           ipoinv,indinv,leninv,npoinv)
            do 65 i=1,npoinv
               dprdin(i,k)=dprdin(i,k)*scainv(i,k)
   65       continue
!$$$            call daxpy(npoinv,z1,scainv(1,k),1,dprdin(1,k),1)
            regula=regula+ddot(npoinv,dprdin(1,k),1,solinv(1,k),1)
   60    continue
      end if
!
      denom=trasmo*tradep
      if (denom.lt.wuzeps) stop ' denom close to zero in funcm'
!$$$      xlmakt=xlmbda*tradat/denom
      xlmakt=xlmbda
!
      funcm=datafu+xlmakt*regula
!$$$      write(*,900) datafu,regula,funcm
!
      return
  900 format(' Data: ',e12.5,' Regular: ',e12.5,' Sum: ',e12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dfuncm(invbnd,ipoinv,indinv,nodinf,                    &
     &                  flxmat,solinv,vechey,vechew,                    &
     &                  vechev,stfinv,dprdin,scainv,                    &
     &                  weiinv,covinv,potmes)
      include 'neurofem.inc'
!
      dimension invbnd(npoinv),ipoinv(npoinv),indinv(leninv),           &
     &          nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          solinv(npoinv,ndmflx),vechey(nevlpo),                   &
     &          vechew(npoinv,ndmflx),vechev(npoinv,ndmflx),            &
     &          stfinv(leninv,ndmflx),                                  &
     &          dprdin(npoinv,ndmflx),scainv(npoinv,ndmflx),            &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo),potmes(nevlpo)
!
!---data part of the gradient
      call dset(nevlpo,z0,vechey,1)
      do 20 j=1,ninfpo
         jglo=nodinf(j)
         do 130 k=1,ndmflx
            vw=solinv(jglo,k)
            do 10 i=1,nevlpo
               fx=flxmat(i,j,k)
               co=covinv(i)
               vechey(i)=vechey(i)+co*fx*vw
   10       continue
  130    continue
   20 continue
!
      do 33 i=1,nevlpo
         co=covinv(i)
         po=potmes(i)
         vechey(i)=vechey(i)-co*po
   33 continue
!
      do 35 k=1,ndmflx
         do 30 j=1,ninfpo
            jglo=nodinf(j)
!$$$            we=weiinv(jglo,k)
            sum=z0
            do 40 i=1,nevlpo
!$$$               sum = sum + flxmat(i,j,k)*vechey(i)*we
               sum = sum + flxmat(i,j,k)*vechey(i)
   40       continue
            vechew(jglo,k) = sum
   30    continue
   35 continue
!
! regularization part of the gradient
      if (lentro) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         gesint=z0
         do 200 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 210 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  210       continue
            streng=sqrt(streng)
            gesint=gesint+streng
  200    continue
         gesint=max(gesint,wuzeps)
         tiny=wuzeps*gesint
!
         func=z0
         do 250 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 260 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  260       continue
            streng=max(sqrt(streng),tiny)
            entlog=log(streng/gesint)
            func=func+entlog*streng/gesint
            do 240 k=1,ndmflx
               dprdin(iglo,k)=entlog*solinv(iglo,k)                     &
     &              *sqrt(scainv(iglo,k))/streng/gesint
  240       continue
  250    continue
!
         do 270 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 280 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  280       continue
            streng=max(sqrt(streng),tiny)
            fak=func/(gesint*streng)
            do 290 k=1,ndmflx
               dprdin(iglo,k)=-(dprdin(iglo,k)-solinv(iglo,k)           &
     &              *sqrt(scainv(iglo,k))*fak)
  290       continue
  270    continue
      else if (ll1nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         do 400 i=1,ninfpo
            iglo=nodinf(i)
            streng=z0
            do 420 k=1,ndmflx
               so=solinv(iglo,k)*sqrt(scainv(iglo,k))
               streng=streng+so*so
  420       continue
            if(streng.gt.z0) strinv=z1/sqrt(streng)
            do 410 k=1,ndmflx
               if (streng.gt.z0) then
                  dprdin(iglo,k)=solinv(iglo,k)*strinv                  &
     &                 *sqrt(scainv(iglo,k))
               else
                  dprdin(iglo,k)=z0
!$$$                  dprdin(iglo,k)=z1
               end if
  410       continue
  400    continue
      else if (ll2nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         do 300 i=1,ninfpo
            iglo=nodinf(i)
            do 310 k=1,ndmflx
               dprdin(iglo,k)=solinv(iglo,k)*scainv(iglo,k)
  310       continue
  300    continue
      else
         do 60 k=1,ndmflx
            call matpro(stfinv(1,k),solinv(1,k),dprdin(1,k),            &
     &           ipoinv,indinv,leninv,npoinv)
            do 65 i=1,npoinv
               dprdin(i,k)=dprdin(i,k)*scainv(i,k)
   65       continue
!$$$            call daxpy(npoinv,z1,scainv(1,k),1,dprdin(1,k),1)
   60    continue
      end if
!
      denom=trasmo*tradep
      if (denom.lt.wuzeps) stop ' denom close to zero in dfuncm'
!$$$      xlmakt=xlmbda*tradat/denom
      xlmakt=xlmbda
      do 70 k=1,ndmflx
         do 80 i=1,ninfpo
            iglo=nodinf(i)
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
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine funchm(invbnd,ipoinv,indinv,nodinf,                    &
     &               flxmat,solinv,vechey,vechew,                       &
     &               vechev,stfinv,scainv,weiinv,                       &
     &               covinv,dprdin,potmes,datafu,                       &
     &               xl1reg,xl2reg,xlcreg,entreg,funreg)
      include 'neurofem.inc'
      dimension invbnd(npoinv),ipoinv(npoinv),indinv(leninv),           &
     &          nodinf(ninfpo)
      dimension flxmat(nevlpo,ninfpo,ndmflx),                           &
     &          solinv(npoinv,ndmflx),vechey(nevlpo),                   &
     &          vechew(npoinv,ndmflx),vechev(npoinv,ndmflx),            &
     &          stfinv(leninv,ndmflx),scainv(npoinv,ndmflx),            &
     &          weiinv(npoinv,ndmflx),covinv(nevlpo),                   &
     &          dprdin(npoinv,ndmflx),potmes(nevlpo)
!
!---determine calculated potentials
      call dset(nevlpo,z0,vechey,1)
      do 20 j=1,ninfpo
         jglo=nodinf(j)
         do 130 k=1,ndmflx
            vw=solinv(jglo,k)
            do 10 i=1,nevlpo
               fx=flxmat(i,j,k)
               vechey(i)=vechey(i)+fx*vw
   10       continue
  130    continue
   20 continue
!
!---data part of the functional
      datafu=z0
      do 30 i=1,nevlpo
         dif=vechey(i)-potmes(i)
         co=covinv(i)
         datafu=datafu+dif*co*dif
   30 continue
!
!---regularization part of the functional
!
!---entropy
      gesint=z0
      do 200 i=1,ninfpo
         iglo=nodinf(i)
         streng=z0
         do 210 k=1,ndmflx
            so=solinv(iglo,k)*sqrt(scainv(iglo,k))
            streng=streng+so*so
  210    continue
         streng=sqrt(streng)
         gesint=gesint+streng
  200 continue
      gesint=max(gesint,wuzeps)
      tiny=wuzeps*gesint
      entreg=z0
      do 220 i=1,ninfpo
         iglo=nodinf(i)
         streng=z0
         do 230 k=1,ndmflx
            so=solinv(iglo,k)*sqrt(scainv(iglo,k))
            streng=streng+so*so
  230    continue
         streng=sqrt(streng)/gesint
         if(streng.gt.z0) entreg=entreg-streng*log(streng)
  220 continue
!---l1-norm
      xl1reg=z0
      do 400 i=1,ninfpo
         iglo=nodinf(i)
         streng=z0
         do 410 k=1,ndmflx
            so=solinv(iglo,k)*sqrt(scainv(iglo,k))
            streng=streng+so*so
  410    continue
         streng=sqrt(streng)
         xl1reg=xl1reg+streng
  400 continue
!---l2-norm
      xl2reg=z0
      do 300 i=1,ninfpo
         iglo=nodinf(i)
         do 310 k=1,ndmflx
            so=solinv(iglo,k)*sqrt(scainv(iglo,k))
            xl2reg=xl2reg+so*so
  310    continue
  300 continue
!---continuous l2-norm
      xlcreg=z0
      do 60 k=1,ndmflx
         call matpro(stfinv(1,k),solinv(1,k),dprdin(1,k),               &
     &        ipoinv,indinv,leninv,npoinv)
         do 65 i=1,npoinv
            dprdin(i,k)=dprdin(i,k)*scainv(i,k)
   65    continue
!$$$         call daxpy(npoinv,z1,scainv(1,k),1,dprdin(1,k),1)
         xlcreg=xlcreg+ddot(npoinv,dprdin(1,k),1,solinv(1,k),1)
   60 continue
!
      if (ll1nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=xl1reg
      else if (ll2nrm) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=xl2reg
      else if (lentro) then
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=entreg
      else 
         trasmo=dble(ninfpo*ndmflx)
         tradep=z1
         regula=xlcreg
      end if
      denom=trasmo*tradep
      if (denom.lt.wuzeps) stop ' denom close to zero in funcm'
!$$$      xlmakt=xlmbda*tradat/denom
      xlmakt=xlmbda
!
      funreg=datafu+xlmakt*regula
!$$$      write(*,900) datafu,regula,funcm
!
      return
  900 format(' Data: ',e12.5,' Regular: ',e12.5,' Sum: ',e12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!
!
!
!
!
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     These routines are modified versions of corresponding
!     "Numerical Recipes" routines
!
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine frprmn(p,n,ftol,iter,fret,g,h,xi,itmax,eps,            &
     &           pcom,xicom,xt,wuzeps,nout)
      integer iter,n,itmax
      double precision fret,ftol,p(n),eps,func,wuzeps,                  &
     &                 dat,xl1,xl2,xlc,ent,fun
!$$$      external func
!u    uses dfunc,func,linmin
      integer its,j
      double precision dgg,fp,gam,gg,g(n),h(n),xi(n)
      double precision pcom(n),xicom(n),xt(n)
!
      fp=func(p)
      call dfunc(p,xi)
      write(*,900) fp
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,itmax
        iter=its
        call linmin(p,xi,n,fret,pcom,xicom,xt,wuzeps)
        if(mod(iter,nout).eq.0) then
           call funchk(p,dat,xl1,xl2,xlc,ent,fun)
           write(*,920) iter,dat,xl1,xl2,xlc,ent,fun
        end if
        if(2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) goto 2000
        fp=func(p)
        call dfunc(p,xi)
        gg=0.d0
        dgg=0.d0
        do 12 j=1,n
          gg=gg+g(j)*g(j)
!---fletcher reeves         dgg=dgg+xi(j)*xi(j)
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.d0) goto 2000
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
      write(*,'(a)') 'frprmn maximum iterations exceeded'
 2000 continue
      write(*,930)
      call funchk(p,dat,xl1,xl2,xlc,ent,fun)
      write(*,920) iter,dat,xl1,xl2,xlc,ent,fun
!
      return
!
  900 format(' Initial functional value ',e12.5)
  910 format('         functional value ',e12.5,' iter ',i6)
  920 format(' i ',i5,' d ',e10.3,' l1 ',e10.3,' l2 ',e10.3,            &
     &                ' lc ',e10.3,' e ',e10.3,' f ',e10.3)
  930 format(' Final evaluation of the functional')
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine linmin(p,xi,n,fret,pcom,xicom,xt,tol)
      integer n
      double precision fret,p(n),xi(n),tol
!u    uses brent,f1dim,mnbrak
      integer j
      double precision ax,bx,fa,fb,fx,xmin,xx,                          &
     &     pcom(n),xicom(n),xt(n),brent
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,pcom,xicom,xt,n)
      fret=brent(ax,xx,bx,tol,xmin,pcom,xicom,xt,n)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,pcom,xicom,xt,nmax)
      double precision ax,bx,cx,fa,fb,fc,f1dim,gold,glimit,tiny
      parameter (gold=1.618034d0, glimit=100.d0, tiny=1.d-20)
      double precision dum,fu,q,r,u,ulim
      integer nmax
      double precision pcom(nmax),xicom(nmax),xt(nmax)
      fa=f1dim(ax,pcom,xicom,xt,nmax)
      fb=f1dim(bx,pcom,xicom,xt,nmax)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+gold*(bx-ax)
      fc=f1dim(cx,pcom,xicom,xt,nmax)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
          fu=f1dim(u,pcom,xicom,xt,nmax)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+gold*(cx-bx)
          fu=f1dim(u,pcom,xicom,xt,nmax)
        else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=f1dim(u,pcom,xicom,xt,nmax)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            fb=fc
            fc=fu
            fu=f1dim(u,pcom,xicom,xt,nmax)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=f1dim(u,pcom,xicom,xt,nmax)
        else
          u=cx+gold*(cx-bx)
          fu=f1dim(u,pcom,xicom,xt,nmax)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function brent(ax,bx,cx,tol,xmin,pcom,xicom,xt,n)
      integer itmax,n
      double precision brent,ax,bx,cx,tol,xmin,f,cgold,zeps
      parameter (cgold=.3819660d0,zeps=1.0d-10)
      integer iter
      double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,                   &
     &     r,tol1,tol2,u,v,w,x,xm
      double precision pcom(n),xicom(n),xt(n), f1dim
      itmax=max(50,n)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f1dim(x,pcom,xicom,xt,n)
      fv=fx
      fw=fx
      do 11 iter=1,itmax
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+zeps
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-0.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(0.5d0*q*etemp).or.                           &
     &       p.le.q*(a-x).or.p.ge.q*(b-x))                              &
     &    goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f1dim(u,pcom,xicom,xt,n)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      write(*,'(a)') 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function f1dim(x,pcom,xicom,xt,nmax)
      double precision f1dim,func,x
      integer nmax
!u    uses func
      integer j
      double precision pcom(nmax),xicom(nmax),xt(nmax)
      do 11 j=1,nmax
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
