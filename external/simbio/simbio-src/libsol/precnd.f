c---->---1---------2---------3---------4---------5---------6---------7--<
      subroutine precnd (gstif,ipodia,indexj,lenge,npoin,
     &                                 gstilu,ispalt,merk)
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
c
c
c     precnd:
c     =======
c
c     ( partielle l u zerlegung )
c     zerlegt die Matrix gstif unvollstndig in eine untere l
c     und in eine obere u Dreiecksmatrix,
c     die multiplikativ verknuepft gstilu ergeben !
c               gstif =~ gstilu:=l*u
c     l und u sind zusammen in gstilu abgelegt, wobei die
c     Diagonalelemente von u den wert 1 haben und in gstilu
c     nicht abgelegt werden!
c
c     hier omega=0 und nicht einstellbar!
c
c           ( mit richardscher spaltenverpointerung ! )
c     ( es sei bemerkt das die verpointerung in ispalt ans )
c     (          hauptprogramm zurueckgegeben wird !        )
c
c
c variable  typ  art  dimension   erklaerung
c---------------------------------------------------------------------
c gstif      d    i    (lenge)    matrix
c ipodia     i    i    (npoin+1)  zeiger fuer position der
c                                 diagonalelemente in gstif
c indexj     i    i    (lenge)    zeiger fuer spalte der elemente aus
c                                 gstif
c lenge      i    i               dimensionierung fuer matrixvektor
c npoin      i    i               dimensionierung fuer vektoren
c gstilu     d    o    (lenge)    matrix
c ispalt    i    o    (lenge)    zeiger fuer spalte
c
c
c merk                 (npoin)    hilfsvector
c=====================================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      logical ist
      parameter (eps=1.d-10,z0=0.d0)
      dimension gstif(lenge),gstilu(lenge)
      dimension ipodia(npoin+1),indexj(lenge)
      dimension ispalt(lenge),merk(npoin)
c
c richardsche verpointerung
c
      do 1 i=1,npoin
         merk(i)=0
 1    continue
      do 2 i=1,lenge
         if (merk(indexj(i)) .eq. 0) then
            merk(indexj(i))=i
         else
            ispalt(merk(indexj(i)))=i
            merk(indexj(i))=i
         end if
 2    continue
c
      do 3 i=1,npoin
         ispalt(merk(i))=lenge+1
 3    continue
c
c beginn von rilu.f
c
      do 10 i=1,npoin
         ilow=ipodia(i)
         ihigh=ipodia(i+1)-1
         naechst=ilow+1
         do 20 j=1,npoin
c
c diagonalelement ?
c-----------------------------
            if (i .eq. j) then
c-----------------------------
c diagonalelement !
c
c
c  s(i,j)=gstif(i,j)-summe von k=1 bis min(i,j)-1 ueber l(i,k)*u(k,j)
c
               ist= .false.
               s=z0
               do 30 k=ilow+1,ihigh
                  if (indexj(k) .ge. i) goto 23
                  if (ist) then
                     do 25 ido=1,npoin
                        if (ipodia(indexj(k))+1 .le. lhilf) goto 42
                        lhilf=ispalt(lhilf)
 25                  continue
 42                  if (ipodia(indexj(k)+1)-1 .ge. lhilf) then
                        s=s+gstilu(k)*gstilu(lhilf)
                        lhilf=ispalt(lhilf)
                     end if
                  else
                     do 40 l=ipodia(indexj(k))+1,ipodia(indexj(k)+1)-1
                        if (indexj(l) - j) 40,46,30
 46                     s=s+gstilu(k)*gstilu(l)
                        lhilf=ispalt(l)
                        ist= .true.
                        goto 30
 40                  continue
                  end if
 30            continue
 23            gstilu(ilow)=gstif(ilow)-s
c---------------
            else
c---------------
c kein diagonalelement !
c
               ijmin=min0(i,j)
               do 50 m=naechst,ihigh
                  if (indexj(m) - j) 50,52,20
c
c  s(i,j)=gstif(i,j)-summe von k=1 bis min(i,j)-1 ueber l(i,k)*u(k,j)
c
 52               naechst=m+1
                  ist= .false.
                  s=z0
                  do 130 k=ilow+1,ihigh
                     if (indexj(k) .ge. ijmin) goto 53
                     if (ist) then
                        do 70 ido=1,npoin
                           if (ipodia(indexj(k))+1 .le. lhilf) goto 99
                           lhilf=ispalt(lhilf)
 70                     continue
 99                     if (ipodia(indexj(k)+1)-1 .ge. lhilf) then
                           s=s+gstilu(k)*gstilu(lhilf)
                           lhilf=ispalt(lhilf)
                        end if
                     else
                        do 140 l=ipodia(indexj(k))+1,
     &                                      ipodia(indexj(k)+1)-1
                           if (indexj(l) - j) 140,146,130
 146                       s=s+gstilu(k)*gstilu(l)
                           lhilf=ispalt(l)
                           ist= .true.
                           goto 130
 140                    continue
                     end if
 130              continue
 53               gstilu(m)=gstif(m)-s
                  goto 20
 50            continue
c-----------------
            end if
c-----------------
 20      continue
c
c--------schummeldiagonalelement zur verhinderung von singularitaet
c
         if (dabs(gstilu(ilow)) .lt. eps)
     &                         gstilu(ilow)=dsign(eps,gstilu(ilow))
c
c--------do 60 j=i+1,npoin u(i,j)=u(i,j)/l(i,i),60 continue
c
         do 60 k=ihigh,ilow+1,-1
            if (indexj(k) .le. i) goto 10
            gstilu(k)=gstilu(k)/gstilu(ilow)
 60      continue
 10   continue
      return
      end
 
 
 
 
 
 
 
 
 
 
