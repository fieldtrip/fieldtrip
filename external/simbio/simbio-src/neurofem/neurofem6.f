!  $4 19/03/2007  Anwander A.  scale difference with the dipole moment
!  $3 21/03/2003  Anwander A.  released for SIMBIO
!  $2 25.04.2001  AA, Since ipvtdi(maxnei)-Parameter in subroutine diplod
!                   is never used, it was deleted.
!  $1 21.12.2000  CW, Since srcgeo(npogeo)-Parameter in subroutine diplod
!                   is never used, it was deleted.

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem6.f :           Subroutines related to dipole sources
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

      subroutine diplod(ipoasy,indasy,mrkvec,
     &                  ipogeo,indgeo,neigeo,
     &                  xyzgeo,difvec,dipmat,diprhs,
     &                  dipsol,vecbgo,dipreg,xyzflf,
     &                  iknflf,numnei)
!
!---determine a source load configuration that generates
!---the desired component of the dipole moment vector at location xyzflf(iknflf)
!---iknflf is in the specific numbering of the relevant geometry,
!---neigeo gives the most adjacent node in the volume structure
!
      include 'neurofem.inc'
      dimension ipoasy(npogeo+1),indasy(lenasy),mrkvec(maxnei),
     &          ipogeo(npogeo),indgeo(lengeo),
     &          neigeo(npoflf)
      dimension xyzgeo(npogeo,3),difvec(3,maxnei),
     &          dipmat(ndimmt),diprhs(maxnei),
     &          dipsol(3*norder),vecbgo(npogeo),
     &          dipreg(maxnei),
     &          xyzflf(npoflf,3)
      dimension dipmom(4)
!
!---first determine the neighbouring nodes to neigeo(iknflf)
      call neibor(ipoasy,indasy,mrkvec,xyzflf,xyzgeo,
     &            neigeo(iknflf),iknflf,numnei)
!
!---then determine the coordinate difference vectors "difvec" 
!---between the neighbouring nodes and xyzflf(iknflf)
      call dipdif(mrkvec,xyzgeo,xyzflf,difvec,iknflf,numnei)
!
!---determine spatial regularizer
      call dpinmt(ipogeo,indgeo,mrkvec,
     &            dipreg,difvec,numnei)
!
!---build the system of equations for the node configuration considered
      call dipsys(dipmat,diprhs,dipsol,difvec,dipreg,
     &            numnei)
!
!---solve the system of equations
      call dppfa(dipmat,numnei,info)
      if (info.ne.0) stop 'dppfa in diplod'
      call dppsl(dipmat,numnei,diprhs)
!
!---add solution to the proper positions of the load vector
      cpumom=qscput(7,0,ierr)
      do 10 i=1,numnei
         ipos=mrkvec(i)
         vecbgo(ipos)=vecbgo(ipos)-diprhs(i)
   10 continue
      dipsum=dsum(numnei,diprhs,1)
!aa scale summe with the dipole moment
      dmax=z0
      do 15 j=1,3
         dmax=max(dmax,abs(dipsol(j)*scadip))
   15 continue
!      if (abs(dipsum-z0).gt.tol) write(*,900) iknflf,dipsum
      if (abs(dipsum/(dmax/10)-z0).gt.tol) write(*,900) iknflf,dipsum
      do 40 iorder=1,norder
         call dset(3,z0,dipmom,1)
         do 20 i=1,numnei
            do 30 idim=1,3
               dipmom(idim) = dipmom(idim)+
     &              (difvec(idim,i)*scadip)**iorder*diprhs(i)
   30       continue
   20    continue
         if (.not.linmat.and..not.linver) then
            iposol=(iorder-1)*3
            ianf=iposol+1
            iend=iposol+3
            diff=z0
            do 50 j=1,3
               idum=ianf-1+j
               diff=max(diff,abs(dipmom(j)-dipsol(idum)*scadip**iorder))
   50       continue
!aa scale difference with the dipole moment
            diff=diff/(dmax/10)
            if (diff.gt.tol) then
                  write(*,930) iorder,(dipmom(j),j=1,3),
     &                 ((dipsol(k)*scadip**iorder),k=ianf,iend)
            end if
         end if
   40 continue
!
      return
  900 format (1x,' Dipole-node ',i8,' Sum ',g12.5)
  930 format (1x,' Dipole-Order ',i3,/,
     &           ' Calculated Moments ',3(1x,g12.5),/,
     &           ' Desired    Moments ',3(1x,g12.5)   )
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      function neichk(ipoasy,npogeo,lenasy)
      dimension ipoasy(npogeo+1)
!
!---finds maximum number of neighbours for a node
      neichk=0
      do 10 i=1,npogeo
!$$$         neichk=max(neichk,ipoasy(i+1)-ipoasy(i)-1)
!---most adjacent node will also be considered
         neichk=max(neichk,ipoasy(i+1)-ipoasy(i))
   10 continue
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine neibor(ipoasy,indasy,mrkvec,xyzflf,xyzgeo,
     &           ikngeo,iknflf,numnei)
      include 'neurofem.inc'
      dimension ipoasy(npogeo+1),indasy(lenasy),mrkvec(maxnei)
      dimension xyzgeo(npogeo,3),xyzflf(npoflf,3)
!
!---find nodes adjacent to node "ikngeo" from compact asymmetric storage
!---and store them in an indexfield "mrkvec". "mrkvec" does
!---contain "ikngeo" itself, if it is not too close to iknflf
      scafak=z1/scadip
      ianf=ipoasy(ikngeo)
      iend=ipoasy(ikngeo+1)-1
      numnei=0
      do 10 j=ianf,iend
         disnod=z0
         jkngeo=indasy(j)
         do 20 idm=1,3
            dif=(xyzgeo(jkngeo,idm)-xyzflf(iknflf,idm))*scafak
            disnod=disnod+dif*dif
   20    continue
         if (sqrt(disnod).gt.tol) then
            numnei=numnei+1
            mrkvec(numnei)=indasy(j)
         else
            idum=1
         end if
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dipdif(mrkvec,xyzgeo,xyzflf,difvec,iknflf,numnei)
      include 'neurofem.inc'
      dimension mrkvec(numnei)
      dimension xyzgeo(npogeo,3),difvec(3,maxnei),
     &          xyzflf(npoflf,3)
!
      scafak=z1/scadip
      do 10 inei=1,numnei
         kno=mrkvec(inei)
         do 20 idm=1,3
            difvec(idm,inei)=(xyzgeo(kno,idm)-xyzflf(iknflf,idm))*scafak
   20    continue
   10 continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dpinmt(ipogeo,indgeo,mrkvec,
     &                  dipreg,difvec,
     &                  numnei)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo),mrkvec(maxnei)
      dimension dipreg(maxnei),difvec(3,maxnei)
!
      call dset (maxnei,z0,dipreg,1)
         do 510 j=1,numnei
            dummy=z0
            do 520 k=1,3
               dummy=dummy+difvec(k,j)*difvec(k,j)
  520       continue
            dipreg(j)=dipreg(j)+sqrt(dummy)**igladi
  510    continue
!     
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dipsys(dipmat,diprhs,dipsol,difvec,dipreg,
     &                  numnei)
      include 'neurofem.inc'
      dimension diprhs(maxnei),dipmat(ndimmt),
     &          dipsol(3*norder),difvec(3,maxnei),
     &          dipreg(maxnei)
!
!---this is generally not necessary, but just to be on the safe side:
      call dset(maxnei,z0,diprhs,1)
      call dset(ndimmt,z0,dipmat,1)
!      
      do 10 i=1,numnei
         idia=i*(i-1)/2
         do 20 j=1,i
            ipos=idia+j
!---zeroth order
            dummy=dble(3)
!---first to "norder" order
            do 30 k=1,3
               do 35 iorder=1,norder
                  dummy=dummy+(difvec(k,i)*difvec(k,j))**iorder
   35          continue
   30       continue
            dipmat(ipos)=dummy
   20    continue
!---zeroth order 
         dummy=z0
!---first to "norder" order
         do 45 iorder=1,norder
            do 40 k=1,3
               iposol=(iorder-1)*3+k
               dummy=dummy+difvec(k,i)**iorder*dipsol(iposol)
   40       continue
   45    continue
         diprhs(i)=dummy
         dipmat(idia+i)=dipmat(idia+i)+diplam*dipreg(i)
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dipana(mrkdip,ipogeo,indgeo,ibndgo,
     &                  icorno,
     &                  dipole,xyzgeo,potana,potnum,
     &                  potdif,surmat,dprodg,volmat,
     &                  rbdgeo,xyzflf)
      include 'neurofem.inc'
      dimension mrkdip(npoflf,3),ipogeo(npogeo),indgeo(lengeo),
     &          ibndgo(npogeo),icorno(npogeo)
      dimension xyzgeo(npogeo,3),potana(npogeo),potnum(npogeo),
     &          dipole(npoflf,3),potdif(npogeo),surmat(lenflf),
     &          dprodg(npogeo),volmat(lengeo),rbdgeo(npogeo),
     &          xyzflf(npoflf,3)
      dimension epsdip(3),srcloc(3),snkloc(3),srcimg(3),snkimg(3)
!
      distol=f2*disdip
      call dset (npogeo,z0,potana,1)
!
!---analytical solutions for dipole loads
!
!---analytical solution from SMYTHE: Static and Dynamic Electricity
!---Hemisphere Publishing Corporation, Third Edition, 1989, p. 259-260
         numdip=isum(3*npoflf,mrkdip,1)
         numcou=0
         call percen(0,numdip)
         do 310 ipoflf=1,npoflf
            icount=0
            do 320 j=1,3
               icount=icount+mrkdip(ipoflf,j)
  320       continue
            if (icount.eq.0) goto 310
            do 330 k=1,3
               if (mrkdip(ipoflf,k).eq.0) goto 330
!
!---initialize vectors
               call dset (3,z0,epsdip,1)
               call dset (3,z0,srcloc,1)
               call dset (3,z0,snkloc,1)
               call dset (3,z0,srcimg,1)
               call dset (3,z0,snkimg,1)
!
!--dipole moment and source strength
               dipmom=dipole(ipoflf,k)
               epsdip(k)=disdip
!$$$               source=dipmom/(z2*epsdip(k))
               source=dipmom/(epsdip(k))
               a=radius
!
!---determine location of source and sink
               a1=z0
               a2=z0
               do 340 j=1,3
                  dummy=xyzflf(ipoflf,j)-xyzcgr(j)
                  srcloc(j)=dummy+epsdip(j)*f2
                  snkloc(j)=dummy-epsdip(j)*f2
                  a1=a1+srcloc(j)*srcloc(j)
                  a2=a2+snkloc(j)*snkloc(j)
  340          continue
               a1=sqrt(a1)
               a2=sqrt(a2)
               a1=max(a1,tol)
               a2=max(a2,tol)
!
!---determine location of image source and image sink
               a1i=z0
               a2i=z0
               do 350 j=1,3
                  srcimg(j)=a*a/(a1*a1)*srcloc(j)
                  snkimg(j)=a*a/(a2*a2)*snkloc(j)
                  a1i=a1i+srcimg(j)*srcimg(j)
                  a2i=a2i+snkimg(j)*snkimg(j)
  350          continue
               a1i=sqrt(a1i)
               a2i=sqrt(a2i)
!
!---determine potentials for all nodes in the FEM mesh
               do 370 nodgeo =1,npogeo
                  r   = z0
                  r1  = z0
                  r2  = z0
                  r1i = z0
                  r2i = z0
                  do 360 j=1,3
                     distan=xyzgeo(nodgeo,j)-xyzcgr(j)
                     r  = r   + distan*distan
                     r1 = r1  + (srcloc(j)-distan)*(srcloc(j)-distan)
                     r2 = r2  + (snkloc(j)-distan)*(snkloc(j)-distan)
                     r1i= r1i + (srcimg(j)-distan)*(srcimg(j)-distan)
                     r2i= r2i + (snkimg(j)-distan)*(snkimg(j)-distan)
  360             continue
                  r   = sqrt(r)
                  r1i = sqrt(r1i)
                  r2i = sqrt(r2i)
                  if (r1.lt.distol*distol) then
                     r1inv=(z3-r1/(distol*distol))/(z2*distol)
                     r1= max(sqrt(r1),tol)
                  else
                     r1  = sqrt(r1)
                     r1inv=z1/r1
                  end if
                  if (r2.lt.distol*distol) then
                     r2inv=(z3-r2/(distol*distol))/(z2*distol)
                     r2=max(sqrt(r2),tol)
                  else
                     r2  = sqrt(r2)
                     r2inv=z1/r2
                  end if
                  if (numana.eq.1) then
!
!---smythe, 1989, insulated sphere, potentials for all nodes
                     potana(nodgeo)=potana(nodgeo)+
     &                    source/(z4*pi*permeb)
     &              * ( r1inv - r2inv + a/(a1*r1i) -a/(a2*r2i)
     &       -z1/a*log( ((a*a+a1*r1i)*(a*a+a1*r1i)-a1*a1*r*r)
     &                / ((a*a+a2*r2i)*(a*a+a2*r2i)-a2*a2*r*r) ) )
                  else if (numana.eq.2) then
!
!---fender, 1987, insulated sphere, (potentials for surface nodes only)
                     potana(nodgeo)=potana(nodgeo)+
     &               source/(z4*pi*permeb)* r1/(r*r*r)
                  else if (numana.eq.3) then
!
!--feynman, infinite space dipole formula
                     potana(nodgeo)=potana(nodgeo)+
     &               source/(z4*pi*permeb)
                  else
                     stop 'numana not defined properly'
                  end if
  370          continue

!
!---finish if all dipoles considered                     
               numcou=numcou+1
               call percen(numcou,numdip)
               if (numcou.ge.numdip) goto 2000
  330       continue
  310    continue
 2000 continue
!
!---mid potential correction
 1000 call matpro(volmat,potana,dprodg,
     &     ipogeo,indgeo,lengeo,npogeo)
      pomian=dsum(npogeo,dprodg,1)/volume
      call matpro(volmat,potnum,dprodg,
     &     ipogeo,indgeo,lengeo,npogeo)
      pominu=dsum(npogeo,dprodg,1)/volume
      write(*,910) pomian,pominu
      do 222 i=1,npogeo
         if (lmncor) then
            potana(i)=potana(i)-pomian
         end if
         potdif(i)=potana(i)-potnum(i)
  222 continue
!
!---correlation between analytical solution and FEM solution
      icount=0
      pocoan=z0
      pocofe=z0
      pocomx=z0
      euclid=z0
      eukmax=z0
      do 500 i=1,npogeo
         if(icorno(i).ne.0) then
            icount=icount+1
            pocoan=pocoan+potana(i)*potana(i)
            pocofe=pocofe+potnum(i)*potnum(i)
            pocomx=pocomx+potana(i)*potnum(i)
            euclid=euclid+potdif(i)*potdif(i)
            eukmax=eukmax+potana(i)*potana(i)
         end if
  500 continue
      correl=pocomx/(sqrt(pocoan)*sqrt(pocofe))
      write(*,900) correl
      write(*,920) sqrt(euclid)/sqrt(eukmax)
      write(*,930) icount
!
      return
  900 format('Correlation          analytical <---> FEM: ',g12.5)
  920 format('Euclidean difference analytical <---> FEM: ',g12.5)
  930 format('surface (difference) nodes ',i8)
  910 format('Mid-Potential analytical: ',g12.5,' numerical: ',g12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine inflod(mrkdip,dipole,srcgeo,srcsur,vcnflf,
     &                  iknflf,idim,icall)
      include 'neurofem.inc'
      dimension mrkdip(npoflf,3)
      dimension dipole(npoflf,3),srcgeo(npogeo),srcsur(npogeo),
     &          vcnflf(npoflf,3)
!
!---determine proper loads for the influence operator
      if (icall.gt.0) then
         if (ldipol) then
            if (lnorco) then
               do 10 i=1,3
                  mrkdip(iknflf,i)=1
                  dipole(iknflf,i)=vcnflf(iknflf,i)
   10          continue
            else
               mrkdip(iknflf,idim)=1
               dipole(iknflf,idim)=z1
            end if
         else
            if (invway.eq.1) srcgeo(iknflf)=z1
            if (invway.eq.2) srcsur(iknflf)=z1
         end if
      else if (icall.lt.0) then
!
!---zero loads from previous step
         if (ldipol) then
            if (lnorco) then
               do 20 i=1,3
                  mrkdip(iknflf,i)=0
                  dipole(iknflf,i)=z0
   20          continue
            else
               mrkdip(iknflf,idim)=0
               dipole(iknflf,idim)=z0
            end if
         else
            if (invway.eq.1) srcgeo(iknflf)=z0
            if (invway.eq.2) srcsur(iknflf)=z0
         end if
      else
         stop 'not implemented yet'
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine soldip(dipole,dipsol,iknflf)
      include 'neurofem.inc'
      dimension dipole(npoflf,3),dipsol(3*norder)
!
      call dset(3*norder,z0,dipsol,1)
!
      do 10 iorder=1,norder
         do 20 idim=1,3
            iposol=(iorder-1)*3+idim
            dipsol(iposol)=dipole(iknflf,idim)/(z2*scadip)*
     &           (disdip/(z2*scadip))**(iorder-1)*(1-(-1)**iorder)
   20    continue
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<



