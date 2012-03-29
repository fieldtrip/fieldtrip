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
c Text in den Historyfile ausgeben
c
c isw=-1: Datei erstellen
c isw=0 : Text normal ausgeben
c isw=1 : Text als Kommentar ausgeben
c isw=2 : Text der mit * beginnt als Kommentar ausgeben
c=======================================================================
      subroutine qfhist(text,isw)
      include 'imelibh'
      character*(*) text
      character*(maxlen) zeil
      integer qclen,qfreeu

      l=max(1,qclen(text))
      mode=isw

      if (mode.eq.-1) then
         if (nthist.ne.0) call qfclos(nthist,0)
         nthist=qfreeu()
         call qfcrse(nthist,text,'fo',ierr)
         if (ierr.ne.0) then
            nthist=0
            iwhist=0
            call qdtext('Es wird keine History-Datei erstellt !',2)
         else
            write(zeil,8002) text(1:l)
            call qdtext(zeil,2)
         endif
         return
      endif

      if (iwhist.eq.0) return

      if (mode.eq.2) then
         mode=0
         ianf=1
 100     continue
         if (text(ianf:ianf).eq.' ') then
            if (ianf.lt.l) then
               ianf=ianf+1
               goto 100
            endif
         endif
         if (text(ianf:ianf).eq.'*') mode=1
      endif


      if (mode.eq.1) then
         write(nthist,8001) text(1:l)
      else if (mode.eq.0) then
         write(nthist,'(a)') text(1:l)
      endif

      return
 8001 format('# ',a)
 8002 format(1x,'Die History-Datei ',a,' wird erstellt.')
      end
