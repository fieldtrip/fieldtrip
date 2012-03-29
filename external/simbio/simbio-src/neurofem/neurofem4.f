!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     neurofem4.f :
!     					 	   -------------------
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

      subroutine poteuk(nodevl,vecref,vecsol,potmes)
      include 'neurofem.inc'
      dimension nodevl(nevlpo)
      dimension potmes(nevlpo),vecref(npogeo),vecsol(npogeo)
!
      funcnl=z0
      do 50 i=1,nevlpo
         iglo=nodevl(i)
         value=potmes(i)+vecref(i)-vecsol(iglo)
         funcnl=funcnl+value*value
         if (.not.ltrans) write(*,940) i,iglo,potmes(i)+vecref(i),      &
     &                    vecsol(iglo),value
   50 continue
      funcnl=sqrt(funcnl)
      write(*,900) funcnl
!
      return
 900  format('Euclidean norm of functional: ',e12.5)
 940  format(' Measurement point ',i6,' global node number ',i6,/,      &
     &       ' measured: ',g15.8,' calculated: ',g15.8,' diff: ',g15.8)
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine resglo(nodinf,mrkdip,                                  &
     &                  dirrhs,dipinv,solinv,                           &
     &                  srcinv,vcnflf,ikenn)
      include 'neurofem.inc'
!
      dimension nodinf(ninfpo),mrkdip(npoinv,ndminv)
      dimension dirrhs(ndirhs),dipinv(npoinv,ndminv),                   &
     &          solinv(npoinv,ndmflx),srcinv(npogeo),                   &
     &          vcnflf(npoflf,ndmflf)
!
!$$$rb      if (ldipol) then
         call iset(npoinv*ndminv, 0,mrkdip,1)
         call dset(npoinv*ndminv,z0,dipinv,1)
!$$$rb      else
         call dset(npogeo,z0,srcinv,1)
!$$$rb      end if
!
      do 40 i=1,ninfpo
         nodinv=nodinf(i)
         if (ikenn.eq.0) dummy=dirrhs(nodinv)
         if (ikenn.eq.1) dummy=solinv(nodinv,1)
         if (ldipol) then
            if (lnorco) then
               do 50 idim=1,ndminv
                  mrkdip(nodinv,idim)=1
                  dipinv(nodinv,idim) = dummy * vcnflf(nodinv,idim)
   50          continue
            else
               do 30 k=1,ndmflx
                  ipodof=(nodinv-1)*ndmflx+k
                  mrkdip(nodinv,k)=1
                  if (ikenn.eq.0) dipinv(nodinv,k)=dirrhs(ipodof)
                  if (ikenn.eq.1) dipinv(nodinv,k)=solinv(nodinv,k)
   30          continue
            end if
         else
            srcinv(nodinv)   = dummy
            dipinv(nodinv,1) = dummy
            mrkdip(nodinv,1) = 1
         end if
   40 continue
!
      return
      end
!---->---1---------2---------3---------4---------5---------6---------7--<
!$$$rb      do 40 i=1,ninfpo
!$$$rb         nodinv=nodinf(i)
!$$$rb         if (ldipol) then
!$$$rb            if (lnorco) then
!$$$rb               if (ikenn.eq.0) dummy=dirrhs(nodinv)
!$$$rb               if (ikenn.eq.1) dummy=solinv(nodinv,1)
!$$$rb               do 50 idim=1,ndminv
!$$$rb                 mrkdip(nodinv,idim)=1
!$$$rb                  dipinv(nodinv,idim) = dummy * vcnflf(nodinv,idim)
!$$$rb   50          continue
!$$$rb            else
!$$$rb               do 30 k=1,ndmflx
!$$$rb                  ipodof=(nodinv-1)*ndmflx+k
!$$$rb                  mrkdip(nodinv,k)=1
!$$$rb                  if (ikenn.eq.0) dipinv(nodinv,k)=dirrhs(ipodof)
!$$$rb                  if (ikenn.eq.1) dipinv(nodinv,k)=solinv(nodinv,k)
!$$$rb   30          continue
!$$$rb            end if
!$$$rb         else
!$$$rb            if (ikenn.eq.0) srcinv(nodinv)=dirrhs(nodinv)
!$$$rb            if (ikenn.eq.1) srcinv(nodinv)=solinv(nodinv,1)
!$$$rb            do 
!$$$rb         end if
!$$$rb   40 continue
