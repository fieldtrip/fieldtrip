!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     femutl.f :
!                              -------------------
!     begin                : Mai 2000
!
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

      subroutine bndcon(ipodia,indexj,irand,gstif,bvec,rwert,
     &                  npoin,length)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0,z1=1.d0)
      dimension ipodia(npoin),indexj(length),irand(npoin)
      dimension gstif(length),bvec(npoin),rwert(npoin)
c
c---consider boundary conditions for system A_ij x_j + b_i = 0_i
c
c--check for element 1
      if (irand(1).gt.0) then
         gstif(1)=z1
         bvec(1)=-rwert(1)
      end if
c
c--walk over the whole matrix, checking each line 
      do 10 i=2,npoin
c
c---if boundary condition shall be set for line i
         if (irand(i).gt.0) then
            do 20 j=ipodia(i-1)+1,ipodia(i)-1
c
c---modify bvec from previous lines if boundary condition was not set
               indj=indexj(j)
               if (irand(indj).eq.0) then
                  bvec(indj)=bvec(indj)+rwert(i)*gstif(j)
               end if
c
c---zero element
               gstif(j)=z0
 20         continue
c
c---set diagonal element to one and boundary value in vector bvec
            gstif(ipodia(i))=z1
            bvec(i)=-rwert(i)
         else
c
c---check for elements to be zeroed because of 
c---boundary conditions in previous lines
            do 30 j=ipodia(i-1)+1,ipodia(i)-1
               jindex=indexj(j)
               if (irand(jindex).gt.0) then
c
c--and modify vector bvec 
                  bvec(i)=bvec(i)+rwert(jindex)*gstif(j)
c
c---zero element
                  gstif(j)=z0
               end if
 30         continue
         end if
 10   continue
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine comsto(ipodia,indexj,necke,knoel,ifreej,
     &     npoin,nelem,maxkno,length,nitot,miwrk)
      implicit double precision(a-h,o-z)
      dimension ipodia(npoin),necke(nelem),ifreej(npoin),
     &          knoel(maxkno,nelem)
      dimension indexj(*)
c
c---establish compact storage
      call iset (npoin,0,ipodia,1)
      call iset (npoin,0,ifreej,1)
c
c---find out maximum number of non zero elements in each line
c---without considering diagonal elements !!!
      do 10 iel=1,nelem
         do 20 iecke=1,necke(iel)
            i=knoel(iecke,iel)
            if (i.eq.0) goto 20
            do 30 jecke=1,necke(iel)
               j=knoel(jecke,iel)
               if ( (j.lt.i).and.(j.gt.0) ) then
                  ipodia(i)=ipodia(i)+1
               end if
 30         continue
 20      continue
 10   continue
c
c---build pointer for diagonal elements under this assumption
      ipodia(1)=1
      isum=1
      do 40 ipoin=2,npoin
c
c--- "+1" for the diagonal element
         isum=isum+ipodia(ipoin)+1
c
c---pointer for diagonal elements at end of line
         ipodia(ipoin)=isum
c
c--pointer for first free position in line ipoin
         ifreej(ipoin)=ipodia(ipoin-1)+1
 40   continue
c
      length=isum
      lenold=length
c
      if ( (length+nitot).gt.miwrk ) then
         write(*,910) length+nitot
         stop 'integerfeld vergroessern'
      end if
c
c---determine indexj for the diagonal element
      do 45 ipoin=1,npoin
         indexj(ipodia(ipoin))=ipoin
 45   continue
c
c---find out indexj without diagonal elements!!!
      do 50 iel=1,nelem
         do 60 iecke=1,necke(iel)
            i=knoel(iecke,iel)
            if (i.eq.0) goto 60
            do 70 jecke=1,necke(iel)
               j=knoel(jecke,iel)
               if ( (j.lt.i).and.(j.gt.0) ) then
c
c---fill position in line i
                  indexj(ifreej(i))=j
                  ifreej(i)=ifreej(i)+1
               end if
 70         continue
 60      continue
 50   continue
      call icopy(npoin,ipodia,1,ifreej,1)
c
c---clean up lines (elimination of zeroes, double elements, sorting)
      do 80 i=2,npoin
         ianf=ifreej(i-1)+1
c
c---determine length of whole line
         lenge=ifreej(i)-ianf+1
c
c---eliminate multiple numbers 
         call elidop(indexj(ianf),lenge)
c
c---sort remaining elements
         call iasort(indexj(ianf),lenge)
         next=ipodia(i-1)+1
c
c---shift sorted line 
         do 90 j=ianf,ianf+lenge-1
            indexj(next)=indexj(j)
            next=next+1
 90      continue
         ipodia(i)=next-1
 80   continue
      length=ipodia(npoin)
c
c--clean rest
      number=lenold-length
      call iset(number,0,indexj(length+1),1)
c
c      do 333 i=2,npoin
c         write(58,900) (indexj(j),j=ipodia(i-1)+1,ipodia(i))
c 333  continue
      return
 900  format(25(i4,1x))
 910  format('Needed Integer Space: ',i8)
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine comasy(ipodia,indexj,necke,knoel,ifreej,
     &     npoin,nelem,maxkno,length,nitot,miwrk)
      implicit double precision(a-h,o-z)
      dimension ipodia(npoin+1),necke(nelem),ifreej(npoin+1),
     &          knoel(maxkno,nelem)
      dimension indexj(*)
c
c---establish compact storage
      call iset (npoin+1,0,ipodia,1)
      call iset (npoin+1,0,ifreej,1)
c
c---find out maximum number of non zero elements "ifreej(i)" in each line "i"
c---without considering diagonal elements !!!
      do 10 iel=1,nelem
         do 20 iecke=1,necke(iel)
            i=knoel(iecke,iel)
            if (i.eq.0) goto 20
            do 30 jecke=1,necke(iel)
               j=knoel(jecke,iel)
               if ( (j.gt.0).and.(j.ne.i) ) then
                  ifreej(i)=ifreej(i)+1
               end if
 30         continue
 20      continue
 10   continue
c
c---build pointer for diagonal elements under this assumption
      isum=1
      do 40 ipoin=1,npoin
c
c---pointer for diagonal elements always at beginning of line
         ipodia(ipoin)=isum
c
c--- "+1" for the next diagonal element
         isum=isum+ifreej(ipoin)+1
c
c--pointer for first free position in line ipoin
         ifreej(ipoin)=ipodia(ipoin)+1
 40   continue
c
      length=isum-1
      lenold=length
      ipodia(npoin+1)=isum
c
      if ( (length+nitot).gt.miwrk ) stop 'integerfeld vergroessern'
c
c---determine indexj for the diagonal element
      do 45 ipoin=1,npoin+1
         indexj(ipodia(ipoin))=ipoin
 45   continue
c
c---find out indexj without diagonal elements!!!
      do 50 iel=1,nelem
         do 60 iecke=1,necke(iel)
            i=knoel(iecke,iel)
            if (i.eq.0) goto 60
            do 70 jecke=1,necke(iel)
               j=knoel(jecke,iel)
               if ( (j.ne.i).and.(j.gt.0) ) then
c
c---fill position in line i
                  indexj(ifreej(i))=j
                  ifreej(i)=ifreej(i)+1
               end if
 70         continue
 60      continue
 50   continue
      call icopy(npoin+1,ipodia,1,ifreej,1)
c      do 334 i=1,npoin
c         write(59,900) (indexj(j),j=ipodia(i),ipodia(i+1)-1)
c 334  continue
c
c---clean up lines (elimination of zeroes, double elements, sorting)
      do 80 i=1,npoin
         ianf=ifreej(i)+1
c
c---determine length of line without diagonal element
         lenge=ifreej(i+1)-ianf
c
c---eliminate multiple numbers 
         call elidop(indexj(ianf),lenge)
c
c---sort remaining elements
         call iasort(indexj(ianf),lenge)
         next=ipodia(i)+1
c
c---shift whole line down to proper position
         do 90 j=ianf,ianf+lenge
            indexj(next)=indexj(j)
            next=next+1
 90      continue
         ipodia(i+1)=next-1
         indexj(next-1)=i+1
         indexj(ifreej(i+1))=0
 80   continue
      length=ipodia(npoin+1)-1
c
c--clean rest
      number=lenold-length
      call iset(number,0,indexj(length+1),1)
c
c      do 333 i=1,npoin
c         write(60,900) (indexj(j),j=ipodia(i),ipodia(i+1)-1)
c 333  continue
      return
 900  format (25(i4,1x))
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elchck(lknel,xyzglo,ityp,maxkno,ierr)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      parameter (mdim=3)
      character*3 chakri
      dimension lknel(maxkno)
      dimension xyzglo(maxkno,3)
      dimension vec12(mdim),vec23(mdim),vec14(mdim),vec15(mdim),
     &          vec45(mdim),vec56(mdim),vec67(mdim)
      dimension v12x23(mdim),v45x56(mdim),v56x67(mdim)
c
c---pruefe ob elemente richtig numeriert sind
c---momentane version: nur ebene,lineare elemente implementiert
      ierr=0
      write(chakri(1:3),'(i3)') ityp
      if (chakri(2:3).eq.'01') then
c
c---lineares Stabelement
         if(maxnum(lknel,maxkno).gt.2) goto 9001
         if ( (xyzglo(2,1)-xyzglo(1,1)).lt.z0 ) then
            call iswap (lknel(1),lknel(2))
            ierr=1
         end if
      else if (chakri(2:3).eq.'02') then
c
c---ebenes Dreieckelement
         if(maxnum(lknel,maxkno).gt.3) goto 9001
c
c---Kreuzprodukt der Vektoren von Punkt 1 nach Punkt 2 und
c                             von Punkt 2 nach Punkt 3
         diskri = ( xyzglo(2,1)-xyzglo(1,1) ) *
     &            ( xyzglo(3,2)-xyzglo(2,2) ) -
     &            ( xyzglo(3,1)-xyzglo(2,1) ) *
     &            ( xyzglo(2,2)-xyzglo(1,2) )
         if (diskri.lt.z0 ) then
            call iswap(lknel(1),lknel(3))
            ierr=1
         end if
      else if (chakri(2:3).eq.'12') then
c
c---ebenes Viereckelement
         if(maxnum(lknel,maxkno).gt.4) goto 9001
c
c---Kreuzprodukt der Vektoren von Punkt 1 nach Punkt 2 und
c                             von Punkt 2 nach Punkt 3
         diskri = ( xyzglo(2,1)-xyzglo(1,1) ) *
     &            ( xyzglo(3,2)-xyzglo(2,2) ) -
     &            ( xyzglo(3,1)-xyzglo(2,1) ) *
     &            ( xyzglo(2,2)-xyzglo(1,2) )
         if (diskri.lt.z0 ) then
c$$$ adrian 6.9.95            call iswap(lknel(1),lknel(3))
            call iswap(lknel(2),lknel(4))
            ierr=1
         end if
      else if (chakri(2:3).eq.'03') then
c
c---tetraeder
         if(maxnum(lknel,maxkno).gt.4) goto 9001
c
         do 200 i=1,ndmgeo
            vec12(i)=xyzglo(2,i)-xyzglo(1,i)
            vec23(i)=xyzglo(3,i)-xyzglo(2,i)
            vec14(i)=xyzglo(4,i)-xyzglo(1,i)
  200    continue
         call crossp(vec12,vec23,v12x23,ndmgeo)
         diskri=ddot(ndmgeo,v12x23,1,vec14,1)
         if (diskri.lt.z0) then
            call iswap(lknel(1),lknel(3))
            ierr=1
         end if
c
      else if (chakri(2:3).eq.'13') then
c
c---keil
         if(maxnum(lknel,maxkno).gt.6) goto 9001
         do 210 i=1,ndmgeo
            vec12(i)=xyzglo(2,i)-xyzglo(1,i)
            vec23(i)=xyzglo(3,i)-xyzglo(2,i)
            vec14(i)=xyzglo(4,i)-xyzglo(1,i)
            vec45(i)=xyzglo(5,i)-xyzglo(4,i)
            vec56(i)=xyzglo(6,i)-xyzglo(5,i)
  210    continue
         call crossp(vec12,vec23,v12x23,ndmgeo)
         call crossp(vec45,vec56,v45x56,ndmgeo)
         diskr1=ddot(ndmgeo,v12x23,1,vec14,1)
         diskr2=ddot(ndmgeo,v45x56,1,vec14,1)
         if (diskr1.lt.z0) then
            call iswap(lknel(1),lknel(3))
            ierr=1
         end if
         if (diskr2.lt.z0) then
            call iswap(lknel(4),lknel(6))
            ierr=1
         end if
c
      else if (chakri(2:3).eq.'23') then
c
c---quader
         if(maxnum(lknel,maxkno).gt.8) goto 9001
         do 220 i=1,ndmgeo
            vec12(i)=xyzglo(2,i)-xyzglo(1,i)
            vec23(i)=xyzglo(3,i)-xyzglo(2,i)
            vec15(i)=xyzglo(5,i)-xyzglo(1,i)
            vec56(i)=xyzglo(6,i)-xyzglo(5,i)
            vec67(i)=xyzglo(7,i)-xyzglo(6,i)
  220    continue
         call crossp(vec12,vec23,v12x23,ndmgeo)
         call crossp(vec56,vec67,v56x67,ndmgeo)
         diskr1=ddot(ndmgeo,v12x23,1,vec15,1)
         diskr2=ddot(ndmgeo,v45x56,1,vec15,1)
         if (diskr1.lt.z0) then
            call iswap(lknel(1),lknel(3))
            ierr=1
         end if
         if (diskr2.lt.z0) then
            call iswap(lknel(4),lknel(6))
            ierr=1
         end if
      else
         goto 9001
      end if
c         
      return
 9001 stop 'elchck: element nicht implementiert'
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine umspeich(wl,ww,nkno,lknel,npoin,maxkno)
      implicit double precision (a-h,o-z)
      dimension wl(maxkno),ww(npoin)
      dimension lknel(maxkno)
      do 10 i=1,nkno
         if(lknel(i).ne.0) wl(i)=ww(lknel(i))
   10 continue
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elidop(ivec,n)
c
c---eliminiert doppelt auftretende punkte in einem integer-vektor
c---und ersetzt sie durch nullen
      dimension ivec(n)
      ieff=n
      do 10 i=n,1,-1
        do 10 j=i-1,1,-1
          if ( ivec(j).eq.ivec(i) ) then
            ivec(j)=ivec(ieff)
            ivec(ieff)=0
            ieff=ieff-1
            goto 10
          end if
   10 continue
      n=ieff
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function ipoj(j,ivec,n)
      dimension ivec(n)
c
c---findet die position von j
      do 10 k=1,n
      if ( ivec(k).eq.j ) then
        ipoj=k-1
        goto 100
      end if
   10 continue
      stop 'j nicht gefunden !!!'
  100 return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine dsuch(zeil,chtext,ltext,iteil,dwert,*)
      implicit double precision (a-h,o-z)
      character zeil*80,chtext*80,fmt*7
      data fmt /'(f10.0)'/
      ifind=index(zeil(iteil:80),chtext(1:ltext))
      if (ifind.gt.0) then
         write(fmt(3:4),'(i2)') iteil-1
         read (zeil(1:iteil-1),fmt) dwert
         return 1
      else
         return
      end if
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine isuch(zeil,chtext,ltext,iteil,iwert,*)
      implicit double precision (a-h,o-z)
      character zeil*80,chtext*80,fmt*5
      data fmt /'(i10)'/
      ifind=index(zeil(iteil:80),chtext(1:ltext))
      if (ifind.gt.0) then
         write(fmt(3:4),'(i2)') iteil-1
         read (zeil(1:iteil-1),fmt) iwert
         return 1
      else
         return
      end if
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine lsuch(zeil,chtext,ltext,iteil,lwert,*)
      implicit double precision (a-h,o-z)
      character zeil*80,chtext*80
      logical lwert
      ifind=index(zeil(iteil:80),chtext(1:ltext))
      if (ifind.gt.0) then
         idum=index(zeil(1:iteil-1),'n')
           if(idum.eq.0) lwert = .true.
           if(idum.gt.0) lwert = .false.
         return 1
      else
         return
      end if
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine chsuch(zeil,chtext,ltext,iteil,chwert,len,*)
      implicit double precision (a-h,o-z)
      character zeil*80,chtext*80,chwert*80
      ifind=index(zeil(iteil:80),chtext(1:ltext))
      if (ifind.gt.0) then
        call textae(zeil(1:iteil-1),ianf,iend,len)
        if(len.gt.0) then
           chwert(1:len)=zeil(ianf:iend)
        else
           chwert(1:80)=' '
        end if
           return 1
      end if
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine elmnul(ivec,ndim)
      dimension ivec(ndim)
c
      ipos=ndim
      do 10 i=1,ndim,1
         if (ivec(i).eq.0) then
 20         ivec(i)=ivec(ipos)
            ivec(ipos)=0
            ipos=ipos-1
            if(ivec(i).eq.0) goto 20
         end if
 10   continue
      ndim=ipos
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function memint(length,nitot)
c
      memint=nitot
      nitot=nitot+length
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function memdbl(length,ndtot)
c
      memdbl=ndtot
      ndtot=ndtot+length
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function maxnum(lknel,maxkno)
      integer maxnum,lknel,maxkno
      dimension lknel(maxkno)
      do 10 i=maxkno,1,-1
         if (lknel(i).ne.0) then
            maxnum=i
            return
         end if
 10   continue
      stop 'fehler maxnum'
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function ncorn(ivec,maxkno)
      dimension ivec(maxkno)
       do 10 ikno=maxkno,1,-1
         if(ivec(ikno).ne.0) goto 20
   10 continue
   20 ncorn=ikno
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine textae(zeil,ianf,iend,lang)
c
c---bestimmt  laenge   <lang> {Output}
c             anfang   <ianf> {Output}
c               ende   <iend> {Output}  
c des characterstrings <zeil> {Input/Output}
c   
c
      character*(*) zeil
c
      do 10 jend=len(zeil),1,-1
         if(zeil(jend:jend).ne.' ') goto 20
   10 continue
      ianf=1
      iend=1
      lang=0
      return
   20 do 30 janf=1,len(zeil)
         if(zeil(janf:janf).ne.' ') goto 40
   30 continue
   40 iend=jend
      ianf=janf
      lang=iend-ianf+1
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine bndpak(irand,gstif,bvec,rwert,
     &                  npoin)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0,z1=1.d0)
      dimension irand(npoin)
      dimension gstif(npoin*(npoin+1)/2),bvec(npoin),rwert(npoin)
c
c---consider boundary conditions for system A_ij x_j = b_i 
c---in packed storage
c
c--check for element 1
      if (irand(1).gt.0) then
         gstif(1)=z1
         bvec(1)= rwert(1)
      end if
c
c--walk over the whole matrix, checking each line 
      do 10 i=2,npoin
         ianf=i*(i-1)/2+1
         iend=i*(i+1)/2-1
c
c---if boundary condition shall be set for line i
         if (irand(i).gt.0) then
            do 20 j=ianf,iend
               jindex=j-i*(i-1)/2
               if (irand(jindex).eq.0) then
                  bvec(jindex)=bvec(jindex)-rwert(i)*gstif(j)
               end if
c
c---zero element
               gstif(j)=z0
 20         continue
c
c---set diagonal element to one and boundary value in vector bvec
            gstif(iend+1)=z1
            bvec(i)= rwert(i)
         else
c
c---check for elements to be zeroed because of 
c---boundary conditions in previous lines
            do 30 j=ianf,iend
               jindex=j-i*(i-1)/2
               if (irand(jindex).gt.0) then
c
c--and modify vector bvec 
                  bvec(i)=bvec(i)-rwert(jindex)*gstif(j)
c
c---zero element
                  gstif(j)=z0
               end if
 30         continue
         end if
 10   continue
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
c-----------------------------------------------------------------------
c
c     subroutine   c o m s c a
c                  ===========
c
c Vorabaufbau der kompakten Speicherung fuer skalare Groessen
c
c EINGANGS-GROESSEN
c  iknoel(mknpel,nelem)  : Knoten-Element-Zuordnung
c  necke(nelem)          : Anzahl der Elementknoten
c  nknot                 : Anzahl der Knoten
c  nelem                 : Anzahl der Elemente
c  mknpel                : Anzahl der Knoten pro Element
c  nitot                 : bisheriges Workfeldende
c  miwrk                 : max. Workfeldgroesse
c
c AUSGANGS-GROESSEN
c  ipodia(nknot)         : Position der Diagonalelemente in der Matrix
c  indexj(*)             : Spalte zu jedem Matrixelement
c  laenge                : Laenge der Gesamtmatrix
c
c HILFSFELDER
c  ivor(*)               : Position in der Matrix, wo das vorherige 
c                          Element der gleichen Zeile steht 
c
c HINWEIS:
c  Beim Aufruf der Routine wird iw(nitot) an indexj und dw(ndtot)
c  an ivor uebergeben. Reicht der damit zur Verfuegung stehende
c  Platz nicht aus, passt auch die Matrix nicht in den Speicher.
c  ( Der Aufbau erfolgt in exakt dem Speicherbereich der spaeter
c    benoetigt wird um die Matrix aufzunehmen)
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine comsca(iknoel,ipodia,indexj,ivor ,necke,
     &                  nknot ,nelem ,mknpel,laenge,nitot,miwrk)
      implicit double precision(a-h,o-z)
      dimension iknoel(mknpel,nelem),necke(nelem),
     &          ipodia(nknot),indexj(*),ivor(*)

      do i=1,nknot
         ipodia(i)=0
      end do

c Schleife ueber alle Elemente

      ipos=0
      do iel=1,nelem
         ne=necke(iel)

c Doppelschleife (iem,jem) ueber lokale Element-Fhgr., um globale 
c besetzte Positionen (i,j) zu finden
c ipodia(i) : wird als Hilfsfeld benutzt; zeigt auf die Pos. des 
c             zuletzt gefundenen Eintrages der Zeile i
c ivor(ip)  : zeigt auf den vorherigen Eintrag der gleichen Zeile

         do iem=1,ne
            i=iknoel(iem,iel)
            do jem=1,ne
               j=iknoel(jem,iel)
               if(j.gt.i) goto 30
               
c Suche in Zeile i, ob Spalte j schon vorhanden
               
               ip=ipodia(i)
   20          continue
               if(ip.ne.0)then
                  if(indexj(ip).eq.j) goto 30
                  ip=ivor(ip)
                  goto 20
               end if

c neue Spalte fuer Zeile i gefunden

               ipos=ipos+1
               if (ipos+nitot.gt.miwrk) stop'comsto: miwrk vergroessern'
               indexj(ipos)=j
               ivor(ipos)=ipodia(i)
               ipodia(i)=ipos
 30            continue
            end do
         end do
      end do
      laenge=ipos

c indexj am Ende von ivor aufbauen und anschliessend
c nach indexj zurueckkopieren

      do i=1,nknot
         ip=ipodia(i)
         if(ip.eq.0) stop 'Leere Zeile in comsto'

         nj=0
 40      continue
         if(ip.ne.0)then
            nj=nj+1
            ivor(ipos+nj)=indexj(ip)
            ip=ivor(ip)
            goto 40
         end if

         call iasort(ivor(ipos+1),nj)
         ipos=ipos+nj
         ipodia(i)=ipos-laenge
         if(ivor(ipos).ne.i) stop 'comsto: falsches Diagonalelement'
      end do
      do ip=1,laenge
         indexj(ip)=ivor(ip+laenge)
      end do

      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine crossp(vec1,vec2,c,ndim)
      implicit double precision (a-h,o-z)
      parameter (z0=0.d0)
      parameter (mdim=3)
      dimension vec1(ndim),vec2(ndim)
      dimension a(mdim),b(mdim),c(mdim)
c
c---berechnet ein Kreuzprodukt c = a x b
      a(1)=vec1(1)
      a(2)=vec1(2)
      b(1)=vec2(1)
      b(2)=vec2(2)
      if (ndim.eq.2) then
         a(3)=z0
         b(3)=z0
      else
         a(3)=vec1(3)
         b(3)=vec2(3)
      end if
c
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=b(1)*a(3)-a(1)*b(3)
      c(3)=a(1)*b(2)-b(1)*a(2)
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine xyzdif(vec1,vec2,vec3,ndim,maxkno)
      implicit double precision (a-h,o-z)
      dimension vec1(ndim),vec2(ndim),vec3(ndim)
c
      do 10 i=1,ndim
         k=(i-1)*maxkno+1
         vec3(i)=vec2(k)-vec1(k)
   10 continue
c
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<

c***********************************************************************
c
c Sammlung von Routinen fuer den Aufbau einer kompakten Speicherung
c
c          (c) Richard Schoenen, Prof. Dr.-Ing. G. Knoll
c
c Namensgebung:  cstoXY
c
c    X: s  --> symmetrisch
c       u  --> unsymmetrisch
c    Y: s  --> skalare Groessen
c       v  --> vektorielle Groessen (mehhrere FHG pro Knoten)
c
c cstosv(iknoel,nelfhg,iknanf,ipodia,indexj,ivor  ,
c        nknot ,nelem ,mknpel,nfrhg ,laenge,mdrest,mirest,ierr)
c
c cstouv(iknoel,nelfhg,iknanf,ipodia,indexj,ivor  ,
c        nknot ,nelem ,mknpel,nfrhg ,laenge,mdrest,mirest,ierr)
c
c cstoss(iknoel,             ,ipodia,indexj,ivor  ,
c        nknot ,nelem ,mknpel,       laenge,mdrest,mirest,ierr)
c
c cstous(iknoel,             ,ipodia,indexj,ivor  ,
c        nknot ,nelem ,mknpel,       laenge,mdrest,mirest,ierr)
c
c
c                                      Aachen, Mai 94'  Richard Schoenen
c***********************************************************************
c-----------------------------------------------------------------------
c
c     subroutine   c s t o s v
c                  ===========
c
c Vorabaufbau der kompakten Speicherung
c symmetrische Matrizen
c vektorielle Groessen (mehrere FHG pro Knoten)
c
c EINGANGS-GROESSEN
c  iknoel(mknpel,nelem)  : Knoten-Element-Zuordnung
c  nelfhg(nknot)         : Anzahl der aktiven Element-FHG
c  iknanf(nknot)         : Knotenanfang in Freiheitsgradliste
c                          ( iknanf(11)+2) liefert den globalen FHG
c                            fuer den 2.lokalen FHG von Knoten 11   )
c  nknot                 : Anzahl der Knoten
c  nelem                 : Anzahl der Elemente
c  mknpel                : Anzahl der Knoten pro Element
c  nfrhg                 : Anzahl der Freiheitsgrade
c  mirest                : max. Integer-Workfeldgroesse
c  mdrest                : max. Double-Workfeldgroesse
c
c AUSGANGS-GROESSEN
c  ipodia(nfrhg)         : Position der Diagonalelemente in der Matrix
c                          (muss vorher genullt sein)
c  indexj(*)             : Spalte zu jedem Matrixelement
c  laenge                : Laenge der Gesamtmatrix
c  ierr                  : 0 --> alles ok
c                          1 --> Integer-Workfeld zu klein
c                          2 --> Double-Workfeld zu klein
c
c HILFSFELDER
c  ivor(*)               : Position in der Matrix, wo das vorherige 
c                          Element der gleichen Zeile steht 
c
c HINWEIS:
c  Beim Aufruf der Routine sollte iw(nitot) an indexj und dw(ndtot)
c  an ivor uebergeben werden. Reicht der damit zur Verfuegung stehende
c  Platz nicht aus, passt auch die Matrix nicht in den Speicher.
c  ( Der Aufbau erfolgt in exakt dem Speicherbereich der spaeter
c    benoetigt wird um die Matrix aufzunehmen)
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cstosv(iknoel,nelfhg,iknanf,ipodia,indexj,ivor ,
     &                  nknot ,nelem ,mknpel,nfrhg ,laenge,
     &                  mdrest,mirest,ierr)
      implicit double precision(a-h,o-z)
      parameter(meldim=200)
      dimension nelfhg(nelem),iknanf(nknot ),iknoel(mknpel,nelem),
     &          ipodia(nfrhg),indexj(mirest),ivor(mdrest*2)
      dimension iglof(meldim)

c Speichergrenze festlegen

      mgrenz=min(mdrest,mirest)

c ipodia nullen

      do i=1,nfrhg
         ipodia(i)=0
      end do

c Schleife ueber alle Elemente

      ipos=0
      do iel=1,nelem

c Zuordnung zwischen Element-Fhgr. und globalen Fhgr. erstellen
c nelma    : Anzahl der Freiheitsgrade des Elementes
c iglof(i) : globaler Frhgr. von Element-Frhg. i

         nelma=0
         do ilk=1,mknpel
            ikn=iknoel(ilk,iel)
            if(ikn.eq.0) goto 10
            do ilf=1,nelfhg(iel)
               nelma=nelma+1
               if(nelma.gt.meldim) stop 'cstosv: Param. meldim zu klein'
               iglof(nelma)=iknanf(ikn)+ilf
            end do
         end do
   10    continue

c Doppelschleife (iem,jem) ueber alle Element-Fhgr., um besetzte
c Positionen (i,j) zu finden
c ipodia(i) : wird als Hilfsfeld benutzt; zeigt auf die Pos. des 
c             zuletzt gefundenen Eintrages der Zeile i
c ivor(ip)  : zeigt auf den vorherigen Eintrag der gleichen Zeile

         do iem=1,nelma
            i=iglof(iem)
            do jem=1,nelma
               j=iglof(jem)
               if(j.gt.i) goto 30
               
c Suche in Zeile i, ob Spalte j schon vorhanden
               
               ip=ipodia(i)
   20          continue
               if(ip.ne.0)then
                  if(indexj(ip).eq.j) goto 30
                  ip=ivor(ip)
                  goto 20
               end if

c neue Spalte fuer Zeile i gefunden

               ipos=ipos+1
               if(ipos.gt.mgrenz) goto 99
               indexj(ipos)=j
               ivor(ipos)=ipodia(i)
               ipodia(i)=ipos

 30            continue
            end do
         end do
      end do
      laenge=ipos

c indexj am Ende von ivor aufbauen und anschliessend
c nach indexj zurueckkopieren

      do i=1,nfrhg
         ip=ipodia(i)
         if(ip.eq.0) stop 'Leere Zeile in cstosv'

         nj=0
 40      continue
         if(ip.ne.0)then
            nj=nj+1
            ivor(ipos+nj)=indexj(ip)
            ip=ivor(ip)
            goto 40
         end if

         call iasort(ivor(ipos+1),nj)
         ipos=ipos+nj
         ipodia(i)=ipos-laenge
         if(ivor(ipos).ne.i) stop 'cstosv: falsches Diagonalelement'
      end do
      do ip=1,laenge
         indexj(ip)=ivor(ip+laenge)
      end do
      ierr=0

      return

   99 continue
      if(ipos.gt.mirest) ierr=1
      if(ipos.gt.mdrest) ierr=2
      return
      end

c-----------------------------------------------------------------------
c
c     subroutine   c s t o u v
c                  ===========
c
c Vorabaufbau der kompakten Speicherung
c unsymmetrische Matrizen
c vektorielle Groessen (mehrere FHG pro Knoten)
c
c EINGANGS-GROESSEN
c  iknoel(mknpel,nelem)  : Knoten-Element-Zuordnung
c  nelfhg(nknot)         : Anzahl der aktiven Element-FHG
c  iknanf(nknot)         : Knotenanfang in Freiheitsgradliste
c                          ( iknanf(11)+2) liefert den globalen FHG
c                            fuer den 2.lokalen FHG von Knoten 11   )
c  nknot                 : Anzahl der Knoten
c  nelem                 : Anzahl der Elemente
c  mknpel                : Anzahl der Knoten pro Element
c  nfrhg                 : Anzahl der Freiheitsgrade
c  mirest                : max. Integer-Workfeldgroesse
c  mdrest                : max. Double-Workfeldgroesse
c
c AUSGANGS-GROESSEN
c  ipodia(nfrhg+1)       : Position der Diagonalelemente in der Matrix
c                          (muss vorher genullt sein)
c  indexj(*)             : Spalte zu jedem Matrixelement
c  laenge                : Laenge der Gesamtmatrix
c  ierr                  : 0 --> alles ok
c                          1 --> Integer-Workfeld zu klein
c                          2 --> Double-Workfeld zu klein
c
c HILFSFELDER
c  ivor(*)               : Position in der Matrix, wo das vorherige 
c                          Element der gleichen Zeile steht 
c
c HINWEIS:
c  Beim Aufruf der Routine sollte iw(nitot) an indexj und dw(ndtot)
c  an ivor uebergeben werden. Reicht der damit zur Verfuegung stehende
c  Platz nicht aus, passt auch die Matrix nicht in den Speicher.
c  ( Der Aufbau erfolgt in exakt dem Speicherbereich der spaeter
c    benoetigt wird um die Matrix aufzunehmen)
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cstouv(iknoel,nelfhg,iknanf,ipodia,indexj,ivor ,
     &                  nknot ,nelem ,mknpel,nfrhg ,laenge,
     &                  mdrest,mirest,ierr)
      implicit double precision(a-h,o-z)
      parameter(meldim=200)
      dimension nelfhg(nelem),iknanf(nknot ),iknoel(mknpel,nelem),
     &          ipodia(nfrhg+1),indexj(mirest),ivor(mdrest*2)
      dimension iglof(meldim)

c Speichergrenze festlegen

      mgrenz=min(mdrest,mirest)

c Schleife ueber alle Elemente

      ipos=0
      do iel=1,nelem

c Zuordnung zwischen Element-Fhgr. und globalen Fhgr. erstellen
c nelma    : Anzahl der Freiheitsgrade des Elementes
c iglof(i) : globaler Frhgr. von Element-Frhg. i

         nelma=0
         do ilk=1,mknpel
            ikn=iknoel(ilk,iel)
            if(ikn.eq.0) goto 10
            do ilf=1,nelfhg(iel)
               nelma=nelma+1
               if(nelma.gt.meldim) stop 'cstouv: Param. meldim zu klein'
               iglof(nelma)=iknanf(ikn)+ilf
            end do
         end do
   10    continue

c Doppelschleife (iem,jem) ueber alle Element-Fhgr., um besetzte
c Positionen (i,j) zu finden
c ipodia(i) : wird als Hilfsfeld benutzt; zeigt auf die Pos. des 
c             zuletzt gefundenen Eintrages der Zeile i
c ivor(ip)  : zeigt auf den vorherigen Eintrag der gleichen Zeile

         do iem=1,nelma
            i=iglof(iem)
            do jem=1,nelma
               j=iglof(jem)
               
c Suche in Zeile i, ob Spalte j schon vorhanden
               
               ip=ipodia(i)
   20          continue
               if(ip.ne.0)then
                  if(indexj(ip).eq.j) goto 30
                  ip=ivor(ip)
                  goto 20
               end if

c neue Spalte fuer Zeile i gefunden

               ipos=ipos+1
               if(ipos.gt.mgrenz) goto 99
               indexj(ipos)=j
               ivor(ipos)=ipodia(i)
               ipodia(i)=ipos

 30            continue
            end do
         end do
      end do
      laenge=ipos

c indexj am Ende von ivor aufbauen und anschliessend
c nach indexj zurueckkopieren

      do i=1,nfrhg
         ip=ipodia(i)

         nj=0
         nj=nj+1   ! fuer Diagonalelement vorgesehen
         ipodia(i)=ipos-laenge+1
 40      continue
         if(ip.ne.0)then
            ij=indexj(ip)
            if(ij.eq.i)then
               ivor(ipos+1)=ij
            else
               nj=nj+1
               ivor(ipos+nj)=ij
            end if
            ip=ivor(ip)
            goto 40
         end if

         call iasort(ivor(ipos+2),nj-1)
         ipos=ipos+nj
      end do
      ipodia(nfrhg+1)=laenge+1
      do ip=1,laenge
         indexj(ip)=ivor(ip+laenge)
      end do
      ierr=0

      return

   99 continue
      if(ipos.gt.mirest) ierr=1
      if(ipos.gt.mdrest) ierr=2
      return
      end

c-----------------------------------------------------------------------
c
c     subroutine   c s t o s s
c                  ===========
c
c Vorabaufbau der kompakten Speicherung
c symmetrische Matrizen
c skalare Groessen (1 FHG pro Knoten)
c
c EINGANGS-GROESSEN
c  iknoel(mknpel,nelem)  : Knoten-Element-Zuordnung
c  nknot                 : Anzahl der Knoten
c  nelem                 : Anzahl der Elemente
c  mknpel                : Anzahl der Knoten pro Element
c  mirest                : max. Integer-Workfeldgroesse
c  mdrest                : max. Double-Workfeldgroesse
c
c AUSGANGS-GROESSEN
c  ipodia(nfrhg)         : Position der Diagonalelemente in der Matrix
c  indexj(*)             : Spalte zu jedem Matrixelement
c  laenge                : Laenge der Gesamtmatrix
c  ierr                  : 0 --> alles ok
c                          1 --> Integer-Workfeld zu klein
c                          2 --> Double-Workfeld zu klein
c
c HILFSFELDER
c  ivor(*)               : Position in der Matrix, wo das vorherige 
c                          Element der gleichen Zeile steht 
c
c HINWEIS:
c  Beim Aufruf der Routine sollte iw(nitot) an indexj und dw(ndtot)
c  an ivor uebergeben werden. Reicht der damit zur Verfuegung stehende
c  Platz nicht aus, passt auch die Matrix nicht in den Speicher.
c  ( Der Aufbau erfolgt in exakt dem Speicherbereich der spaeter
c    benoetigt wird um die Matrix aufzunehmen)
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cstoss(iknoel,ipodia,indexj,ivor ,
     &                  nknot ,nelem ,mknpel,laenge,
     &                  mdrest,mirest,ierr)
      implicit double precision(a-h,o-z)
      dimension iknoel(mknpel,nelem),
     &          ipodia(nknot),indexj(mirest),ivor(mdrest*2)

c Speichergrenze festlegen

      mgrenz=min(mdrest,mirest)

c ipodia nullen

      do i=1,nknot
         ipodia(i)=0
      end do

c Schleife ueber alle Elemente

      ipos=0
      do iel=1,nelem

c Anzahl der Element-FHG bestimmen

         ne=0
         do ilk=1,mknpel
            ikn=iknoel(ilk,iel)
            if(ikn.eq.0) goto 10
            ne=ne+1
         end do
   10    continue

c Doppelschleife (iem,jem) ueber lokale Element-Fhgr., um globale 
c besetzte Positionen (i,j) zu finden
c ipodia(i) : wird als Hilfsfeld benutzt; zeigt auf die Pos. des 
c             zuletzt gefundenen Eintrages der Zeile i
c ivor(ip)  : zeigt auf den vorherigen Eintrag der gleichen Zeile

         do iem=1,ne
            i=iknoel(iem,iel)
            do jem=1,ne
               j=iknoel(jem,iel)
               if(j.gt.i) goto 30
               
c Suche in Zeile i, ob Spalte j schon vorhanden
               
               ip=ipodia(i)
   20          continue
               if(ip.ne.0)then
                  if(indexj(ip).eq.j) goto 30
                  ip=ivor(ip)
                  goto 20
               end if

c neue Spalte fuer Zeile i gefunden

               ipos=ipos+1
               if(ipos.gt.mgrenz) goto 99
               indexj(ipos)=j
               ivor(ipos)=ipodia(i)
               ipodia(i)=ipos

 30            continue
            end do
         end do
      end do
      laenge=ipos

c indexj am Ende von ivor aufbauen und anschliessend
c nach indexj zurueckkopieren

      do i=1,nknot
         ip=ipodia(i)
         if(ip.eq.0) stop 'Leere Zeile in cstoss'

         nj=0
 40      continue
         if(ip.ne.0)then
            nj=nj+1
            ivor(ipos+nj)=indexj(ip)
            ip=ivor(ip)
            goto 40
         end if

         call iasort(ivor(ipos+1),nj)
         ipos=ipos+nj
         ipodia(i)=ipos-laenge
         if(ivor(ipos).ne.i) stop 'cstoss: falsches Diagonalelement'
      end do
      do ip=1,laenge
         indexj(ip)=ivor(ip+laenge)
      end do
      ierr=0

      return

   99 continue
      if(ipos.gt.mirest) ierr=1
      if(ipos.gt.mdrest) ierr=2
      return
      end

c-----------------------------------------------------------------------
c
c     subroutine   c s t o u s
c                  ===========
c
c Vorabaufbau der kompakten Speicherung
c unsymmetrische Matrizen
c skalare Groessen (1 FHG pro Knoten)
c
c EINGANGS-GROESSEN
c  iknoel(mknpel,nelem)  : Knoten-Element-Zuordnung
c  nknot                 : Anzahl der Knoten
c  nelem                 : Anzahl der Elemente
c  mknpel                : Anzahl der Knoten pro Element
c  mirest                : max. Integer-Workfeldgroesse
c  mdrest                : max. Double-Workfeldgroesse
c
c AUSGANGS-GROESSEN
c  ipodia(nfrhg+1)       : Position der Diagonalelemente in der Matrix
c  indexj(*)             : Spalte zu jedem Matrixelement
c  laenge                : Laenge der Gesamtmatrix
c  ierr                  : 0 --> alles ok
c                          1 --> Integer-Workfeld zu klein
c                          2 --> Double-Workfeld zu klein
c
c HILFSFELDER
c  ivor(*)               : Position in der Matrix, wo das vorherige 
c                          Element der gleichen Zeile steht 
c
c HINWEIS:
c  Beim Aufruf der Routine sollte iw(nitot) an indexj und dw(ndtot)
c  an ivor uebergeben werden. Reicht der damit zur Verfuegung stehende
c  Platz nicht aus, passt auch die Matrix nicht in den Speicher.
c  ( Der Aufbau erfolgt in exakt dem Speicherbereich der spaeter
c    benoetigt wird um die Matrix aufzunehmen)
c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine cstous(iknoel,ipodia,indexj,ivor ,
     &                  nknot ,nelem ,mknpel,laenge,
     &                  mdrest,mirest,ierr)
      implicit double precision(a-h,o-z)
      dimension iknoel(mknpel,nelem),
     &          ipodia(nknot+1),indexj(mirest),ivor(mdrest*2)

c Speichergrenze festlegen

      mgrenz=min(mdrest,mirest)

c ipodia nullen

      do i=1,nknot
         ipodia(i)=0
      end do

c Schleife ueber alle Elemente

      ipos=0
      do iel=1,nelem

c Anzahl der Element-FHG bestimmen

         ne=0
         do ilk=1,mknpel
            ikn=iknoel(ilk,iel)
            if(ikn.eq.0) goto 10
            ne=ne+1
         end do
   10    continue

c Doppelschleife (iem,jem) ueber lokale Element-Fhgr., um globale 
c besetzte Positionen (i,j) zu finden
c ipodia(i) : wird als Hilfsfeld benutzt; zeigt auf die Pos. des 
c             zuletzt gefundenen Eintrages der Zeile i
c ivor(ip)  : zeigt auf den vorherigen Eintrag der gleichen Zeile

         do iem=1,ne
            i=iknoel(iem,iel)
            do jem=1,ne
               j=iknoel(jem,iel)
               
c Suche in Zeile i, ob Spalte j schon vorhanden
               
               ip=ipodia(i)
   20          continue
               if(ip.ne.0)then
                  if(indexj(ip).eq.j) goto 30
                  ip=ivor(ip)
                  goto 20
               end if

c neue Spalte fuer Zeile i gefunden

               ipos=ipos+1
               if(ipos.gt.mgrenz) goto 99
               indexj(ipos)=j
               ivor(ipos)=ipodia(i)
               ipodia(i)=ipos

 30            continue
            end do
         end do
      end do
      laenge=ipos

c indexj am Ende von ivor aufbauen und anschliessend
c nach indexj zurueckkopieren

      do i=1,nknot
         ip=ipodia(i)

         nj=0
         nj=nj+1   ! fuer Diagonalelement vorgesehen
         ipodia(i)=ipos-laenge+1
 40      continue
         if(ip.ne.0)then
            ij=indexj(ip)
            if(ij.eq.i)then
               ivor(ipos+1)=ij
            else
               nj=nj+1
               ivor(ipos+nj)=ij
            end if
            ip=ivor(ip)
            goto 40
         end if

         call iasort(ivor(ipos+2),nj-1)
         ipos=ipos+nj
      end do
      ipodia(nknot+1)=laenge+1
      do ip=1,laenge
         indexj(ip)=ivor(ip+laenge)
      end do
      ierr=0

      return

   99 continue
      if(ipos.gt.mirest) ierr=1
      if(ipos.gt.mdrest) ierr=2
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
c---->---1---------2---------3---------4---------5---------6---------7--<
c Shell-Sort, sortiert Integer-Feld ia in aufsteigender Reihenfolge
c Numerical Recipes, 2.Ed, S.323
c
      subroutine iasorf(ia,irf,n)
      dimension ia(n),irf(n)
      inc=1
    1 continue
      inc=3*inc+1
      if(inc.le.n) goto 1
      do i=1,n
         irf(i)=i
      end do
    2 continue
      inc=inc/3
      do i=inc+1,n
         imerk=irf(i)
         iam=ia(imerk)
         j=i
    3    continue
         itest=irf(j-inc)
         if(ia(itest).gt.iam)then
            irf(j)=itest
            j=j-inc
            if(j.gt.inc) goto 3
         end if
         irf(j)=imerk
      end do
      if(inc.gt.1) goto 2
      return
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function meminn(chaakt,chaint,istint,mxpoi,length,nitot)
      character*10 chaakt,chaint(mxpoi)
      integer      istint(mxpoi)
      save icount
      data icount /0/
c
      icount=icount+1
      if (icount.gt.mxpoi) then
         write(*,900) 
         stop 'meminn'
      end if
c
      chaint(icount)=chaakt
      istint(icount)=nitot
      meminn=nitot
      nitot=nitot+length
c
      return
  900 format(' Not enough integer pointers, please enlarge mxpoi ')
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
      function memdbn(chaakt,chadbl,istdbl,mxpoi,length,ndtot)
      character*10 chaakt,chadbl(mxpoi)
      integer      istdbl(mxpoi)
      save icount
      data icount /0/
c
      icount=icount+1
      if (icount.gt.mxpoi) then
         write(*,900)
         stop 'memdbn'
      end if
c
      chadbl(icount)=chaakt
      istdbl(icount)=ndtot
      memdbn=ndtot
      ndtot=ndtot+length
c
      return
  900 format(' Not enough double pointers, please enlarge mxpoi ')
      end
c---->---1---------2---------3---------4---------5---------6---------7--<
