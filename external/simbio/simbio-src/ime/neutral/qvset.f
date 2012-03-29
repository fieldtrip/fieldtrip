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
c@======================================================================
c Variable besetzen
c=======================================================================
      subroutine qvset(var,vwert,ierr)
      include 'imelibh'
      integer qclen,qcllen
      character*(*) var,vwert
      character*(maxlen) wert
      character*36 abc*26,abc123*36
      logical qcvgl
      save abc,abc123
      data abc123 /'abcdefghijklmnopqrstuvwxyz0123456789'/
      data abc    /'abcdefghijklmnopqrstuvwxyz'/
      
      wert=vwert
      lw=qclen(wert)
      if (lw.gt.0) wert=vwert(qcllen(vwert):lw)

      ierr=1
      lv=qclen(var)
      if (lv.le.1) return
      if (lv.gt.lenvar) then
         write(*,8001) lenvar
         lv=lenvar
      endif
      
c Test ob der Variablennamen erlaubt ist

      if (index(abc123,var(1:1)).ne.0) goto 9002
      if (index(abc,var(2:2)).eq.0) goto 9002
      do i=3,lv
         if (index(abc123,var(i:i)).eq.0) goto 9002
      enddo

c Ist die Variable schon da ?

      do i=1,maxvar
         if (lnvar(i).eq.lv) then
            if (qcvgl(nvar(i),var(1:lv))) goto 100
         endif
      enddo
      
c Freie Variable suchen

      do i=1,maxvar
         if (lnvar(i).eq.0) goto 100
      enddo

 100  continue
      if (i.gt.maxvar) goto 9001
      lw=qclen(wert)
      if (lw.gt.0) then
         nvar(i)=var
         tvar(i)=wert
         lnvar(i)=lv
         ltvar(i)=lw
      else
         lnvar(i)=0
         ltvar(i)=0
      endif
      ierr=0

      return
 8001 format(1x,'Variable ist laenger als ',i2,' Zeichen !')

 8901 format(1x,'Es koennen max. ',i4,' Variablen definiert werden !')
 8902 format(1x,'Der Variablenname "',a,'" ist nicht erlaubt !')
 9001 write(*,8901) maxvar
      return
 9002 write(*,8902) var(1:qclen(var))
      return
      end
