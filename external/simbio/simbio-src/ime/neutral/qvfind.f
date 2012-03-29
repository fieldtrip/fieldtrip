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
c Variable in einem String ersetzen
c=======================================================================
      subroutine qvfind(intext)
      include 'imelibh'
      character*(*) intext
      character*(maxlen) copy,text,kltext,chilf
      character*80 abc*36,endtxt*1,vname*(lenvar),cfuell*1
      integer qclen
      logical qcvgl
      save abc
      data abc /'abcdefghijklmnopqrstuvwxyz0123456789'/
      

      text=intext

 100  continue
      ianz=0
      ltxt=qclen(text)
      if (ltxt.le.0) return
      it=1
      kltext=text
      call qctolo(kltext(1:ltxt))
      copy=text
      lneu=0

 200  continue
      if (it.ge.ltxt) then
         if (ianz.eq.0) goto 9000
         text=copy
         goto 100
      endif

c %-Zeichen suchen

      ind=index(text(it:),'%')
      if (ind.le.0) goto 9000
      ind=ind+it

c Test auf %%

      if (text(ind:ind).eq.'%') then
         it=ind+1
         goto 200
      endif

c Formatangabe suchen

      lsub=0
      iform=0
 300  continue
      if (index('1234567890',text(ind+lsub:ind+lsub)).ne.0) then
         if (ind+lsub.eq.ltxt) goto 9000
         lsub=lsub+1
         goto 300
      endif
      lform=lsub
      if (lform.gt.0) read(text(ind:ind+lform-1),*) iform

c Variablennamen suchen

      if (text(ind+lsub:ind+lsub).eq.'{') then
         endtxt='}'
         lsub=lsub+1
      else if (text(ind+lsub:ind+lsub).eq.'(') then
         endtxt=')'
         lsub=lsub+1
      else
         endtxt=' '
      endif

      lnam=0
 400  continue
      if (index(abc,kltext(ind+lsub+lnam:ind+lsub+lnam)).ne.0) then
         lnam=lnam+1
         if (ind+lsub+lnam.le.ltxt) goto 400
      endif
      if (lnam.eq.0) then
         it=ind+lsub
         goto 200
      endif

      if (lnam.gt.lenvar) write(*,8001) lenvar
      vname=text(ind+lsub:ind+lsub+lnam-1)
      lsub=lsub+lnam
      lnam=min(lnam,lenvar)
      if (endtxt.ne.' ') then
         if (text(ind+lsub:ind+lsub).ne.endtxt) then
            it=ind+lsub
            goto 200
         endif
         lsub=lsub+1
      endif

c Variable suchen

      do ivar=1,maxvar
         if (lnvar(ivar)-1.eq.lnam) then
            if (qcvgl(nvar(ivar)(2:),vname)) goto 500
         endif
      enddo
      it=ind+lsub
      goto 200

 500  continue
      
c Wieviel Platz braucht der einzufügende Text

      chilf=copy(ind+lneu+lsub:)
      i1=ind-1+lneu

      if (lform.ne.0) then
         i2=i1+iform-1
      else
         i2=i1+ltvar(ivar)-1
      endif
      copy(i2:)=' '
      lneu=lneu+(i2-i1+1)-(lsub+1)

c Fuellzeichen ermitteln und fuellen

      if (lform.ne.0 .and. iform.gt.ltvar(ivar)) then
         cfuell='0'
         do i=1,ltvar(ivar)
            if (index('0123456789',tvar(ivar)(i:i)).le.0) cfuell='_'
         enddo
         do i=1,iform-ltvar(ivar)
            copy(i1:i1)=cfuell
            i1=i1+1
         enddo
      endif

c Variablenwert kopieren

      copy(i1:i2)=tvar(ivar)(1:ltvar(ivar))
      copy(i2+1:)=chilf
      it=ind+lsub
      ianz=ianz+1
      goto 200

c Doppelte Prozentzeichen filtern

 9000 continue
      intext=' '
      lc=qclen(copy)
      ipp=0
      do i=1,lc
         if (i+ipp.le.lc) then
            intext(i:i)=copy(i+ipp:i+ipp)
            if (copy(i+ipp:i+ipp+1).eq.'%%') ipp=ipp+1
         endif
      enddo

      return
 8001 format(1x,'Variable ist laenger als ',i2,' Zeichen !')
      end
