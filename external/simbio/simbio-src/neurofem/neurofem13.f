!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem13.f :      	 subroutines related to the analytical 
!                            forward model
!                            WARNING: not fully tested
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

      subroutine avefey(knosur,itysur,necsur,                           &
     &                  xyzgeo,xyzflf,dipole,potsng,                    &
     &                  potint,ipoflf)
      include 'neurofem.inc'
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
      dimension xyzgeo(npogeo,ndmgeo),xyzflf(npoflf,ndmflf),            &
     &          dipole(npoflf,ndmflf),potsng(npogeo)
      dimension dipp(3),dipn(3),dippos(3)
!
      call feynma(xyzflf,xyzgeo,potsng,dipole,dipp,dipn,                &
     &           streng,ipoflf)
!
      if (logsur) then
         do 10 ielsur=1,nelsur
            call eleobf(knosur,itysur,necsur,                           &
     &           xyzgeo,ielsur,                                         &
     &           dipp,  dipn,streng,vallok)
            potint=potint+vallok
   10    continue
      end if
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine eleobf(knosur,itysur,necsur,                           &
     &                  xyzgeo,ielsur,                                  &
     &                  posp,posn,fak,vallok)
      include 'neurofem.inc'
      dimension knosur(mxkn2d,nelsur),itysur(nelsur),necsur(nelsur)
      dimension lokkno(32)
      dimension xyzgeo(npogeo,ndmgeo)
      dimension xlk(3,mxkn3d),dinpkt(4),fag(3,mxmxkn),fwe(mxmxkn),      &
     &          adjunc(3),xyzglo(3),rjacob(3,3)
      dimension rpl(3),rmi(3),posp(3),posn(3)
!
!---initialize everything
      vallok=z0
      call iset(32,0,lokkno,1)
!
!---turn to local node numbers
      do 10 i=1,necsur(ielsur)
        kno=knosur(i,ielsur)
        lokkno(i)=kno
        do 15 idim=1,ndmgeo
           xlk(idim,i) = xyzgeo(kno,idim)
           fag(idim,i)=z0
   15   continue
   10 continue
      call itpanz(itysur(ielsur),necsur(ielsur),intgrd,ninpkt,ierr)
      if (ierr.ne.0) stop 'eleobf 1'
!
      do 20 iint=1,ninpkt
!
!---determine local coordinate values and gaussian weighting factor
         call itpdat(itysur(ielsur),intgrd,iint,dinpkt,ierr)
         if (ierr.ne.0) stop 'eleobf 2'
         call fofagn(itysur(ielsur),dinpkt(1),dinpkt(2),dinpkt(3),      &
     &        lokkno(1),xlk,fwe,fag,detj,ierr)
         if (ierr.ne.0) stop 'eleobf 3'
!
!---determine global coordinates at gaussian point
         rpl2=z0
         rmi2=z0
         do 100 idim=1,ndmgeo
            xyzglo(idim)=z0
            do 90 iecke=1,necsur(ielsur)
               xyzglo(idim)=xyzglo(idim)+fwe(iecke)*xlk(idim,iecke)
   90       continue
            rpl(idim)=xyzglo(idim)-posp(idim)
            rmi(idim)=xyzglo(idim)-posn(idim)
            rpl2=rpl2+rpl(idim)*rpl(idim)
            rmi2=rmi2+rmi(idim)*rmi(idim)
  100    continue
         rpl2i=z1/sqrt(rpl2)
         rmi2i=z1/sqrt(rmi2)
!
!---determine components of the normal vector r_{n}
         call getadj(itysur(ielsur),adjunc(1),adjunc(2),adjunc(3),ierr)
         if(ierr.ne.0) stop 'eleobf 4'
!
         result=z0
!$$$         res2=z0
         do 110 k=1,ndmgeo
            result=result+adjunc(k)*(rpl(k)*rpl2i-rmi(k)*rmi2i)
!$$$            res2=res2+adjunc(k)*xyzglo(k)
  110    continue
!$$$         vallok=vallok+result*fak*f2*dinpkt(4)
         vallok=vallok+result*fak*f2*dinpkt(4)
!$$$         vallok=vallok+res2/z3
!
   20 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine feynma(xyzflf,xyzgeo,potsng,dipole,dipp,dipn,          &
     &           streng,ipoflf)
      include 'neurofem.inc'
      dimension xyzflf(npoflf,ndmflf),xyzgeo(npogeo,ndmgeo),            &
     &          potsng(npogeo),dipole(npoflf,ndmflf)
      dimension dstakt(3),dippos(3),dipakt(3),dipp(3),dipn(3)
!
!--- determines analytical potential solution for the infinite space
!--- to describe the fundamental singularity of the dipole source problem
!--- formula is valid for two (singular) monopoles
!
      distol=f2*disdip
!
      dipges=z0
      do 10 k=1,ndmflf
         dippos(k)=xyzflf(ipoflf,k)
         dipakt(k)=dipole(ipoflf,k)
         dipges=dipges+dipakt(k)*dipakt(k)
   10 continue
      dipges=sqrt(dipges)
      streng=dipges/(disdip*z4*pi*permeb)
      sc=f2*disdip/dipges
      do 15 k=1,ndmflf
         dum=sc*dipakt(k)
         dipp(k)=dippos(k)+dum
         dipn(k)=dippos(k)-dum
   15 continue
!
      do 20 i=1,npogeo
         distp=z0
         distn=z0
         do 30 k=1,ndmflf
            dsp=xyzgeo(i,k)-dipp(k)
            dsn=xyzgeo(i,k)-dipn(k)
            distp=distp+dsp*dsp
            distn=distn+dsn*dsn
   30    continue
         if (distp.lt.distol*distol) then
            distpi=(z3-distp/(distol*distol))/(z2*distol)
         else
            distpi=z1/sqrt(distp)
         end if
         if (distn.lt.distol*distol) then
            distni=(z3-distn/(distol*distol))/(z2*distol)
         else
            distni=z1/sqrt(distn)
         end if
         potsng(i)=potsng(i)+streng*(distpi-distni)
   20 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine knodef(knosur,necsur,knoobf)
      include 'neurofem.inc'
      dimension knosur(mxkn2d,nelsur),necsur(nelsur),knoobf(npogeo)
!
      call iset(npogeo,0,knoobf,1)
      do 10 iel=1,nelsur
         do 20 iecke=1,necsur(iel)
            kno=knosur(iecke,iel)
            knoobf(kno)=1
   20    continue
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine radcor(knosur,necsur,xyzgeo)
      include 'neurofem.inc'
      dimension knosur(mxkn2d,nelsur),necsur(nelsur)
      dimension xyzgeo(npogeo,ndmgeo)
!
      radmit=z0
      do 30 k=1,ndmgeo
         radmit=radmit+xyzdst(k)*f2
   30 continue
      radmit=radmit/dble(ndmgeo)
!
      do 40 iel=1,nelsur
         do 50 iecke=1,necsur(iel)
            kno=knosur(iecke,iel)
            x=xyzgeo(kno,1)
            y=xyzgeo(kno,2)
            z=xyzgeo(kno,3)
            radxyz=sqrt(x*x+y*y+z*z)
            scafac=radmit/radxyz
            xyzgeo(kno,1)=x*scafac
            xyzgeo(kno,2)=y*scafac
            xyzgeo(kno,3)=z*scafac
   50    continue
   40 continue
!
      write(*,900) radmit
      return
  900 format(' average radius of surface now ',f12.5)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine rangko(ipogeo,indgeo,knoobf,                           &
     &                  sysmat,sysrng)
      include 'neurofem.inc'
      dimension ipogeo(npogeo),indgeo(lengeo)
      dimension knoobf(npogeo)
      dimension sysmat(lengeo),sysrng(lengeo)
!     
!      if (knoobf(1).gt.0) then
!         sysrng(1)=sysmat(1)
!      end if
!      do 10 i=2,npogeo
!         if (knoobf(i).gt.0) then
!            ipoi=ipogeo(i)
!            do 20 j=ipogeo(i-1)+1,ipoi
!               sysrng(j)=sysmat(j)
!   20       continue
!            do 30 j=ipoi+1,lengeo
!               if (indgeo(j).eq.i) then
!                  sysrng(j)=sysmat(j)
!               end if
!   30       continue
!         end if
!   10 continue
      do 10 i=1,lengeo
         if (knoobf(indgeo(i)).gt.0) then
            sysrng(i)=sysmat(i)
         end if
   10 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!
