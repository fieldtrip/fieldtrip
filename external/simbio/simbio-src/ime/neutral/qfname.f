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
      subroutine qfname (filna, opar, text, ireran, cdfwe, queset, ierr)
c     erfragt vom benutzer einen filenamen 
c
c     ****************
c     * standard f77 *
c     ****************
c
c
c variable typ art dimension      erklaerung
c ----------------------------------------------------------------------
c filna     c   o  *80            name des geoeffneten files
c opar      c   i  *2             eroeffnungsart des files
c                                 'op' : open  'cr' : create
c text      c   i   *(*)          bei der abfrage auszugebender text
c ireran    i   i                 ausrichtung des textes
c                                     0 : linksbuendig
c                                 sonst : rechtsbuendig in spalte ireran
c cdfwe     c   i  *(*)           defaultwert der antwort
c queset    i   i                 1=ueberschreiben; 0=nicht ueberschreiben
c ierr      i   o                 fehlercode
c                                 falls ungleich null, fehler in
c                                 aufgrufenen unterprogrammen, bzw. abbruch
c
      include 'imelibh'

      character*(*) text  ,cdfwe, filna
      character*80 qichar
      character*2  opar,mode*2
      integer iunit, ireran, ieing, ierr, lfilna, queset
      logical  qflegp ,qiyeno, qfilda, ldel

      ierr=0
      mode=opar
      call qctolo(mode)

c zaehler fuer die anzahl der eingaben nullen
c (beim dritten versuch kommt ein hinweis auf abbruch)

      ieing=0
      ierr=0

c input einlesen

 10   continue
      if (ieing.eq.2) then
         ieing=0
         write (ntterm, 8000) cabbru(1:labbru)
 8000    format(1x,'abbruch der eingabe moeglich mit ', a)
      end if
      filna=qichar(text,ireran,cdfwe,lfilna,ierr)
      call qsmnam(filna,' ',filna)
      ieing=ieing+1

c check, ob eingelesener filename gleich cabbru

      if (ierr.gt.0) return

c eingabe von cr und kein default gesetzt

      if (lfilna.eq.0) then
         write (ntterm, 8010)
 8010    format(1x,'Kein Default zulaessig.')
         goto 10
      end if

c eingelesener filename ist nicht legal

      if (.not. qflegp(filna)) then
         write (ntterm, 8020)
 8020    format(1x,'Dateiname illegal.')
         goto 10
      end if

c check ob datei schon da ist

      if (mode(1:1).eq.'c' .and. queset.eq.0) then
         if (qfilda(filna)) then
            ldel=qiyeno('Alten File ueberschreiben ?',ireran,'j')
            if (.not.ldel) then
               write (ntterm, 8030)
 8030          format(1x,'Bitte anderen Dateinamen eingeben.')
               goto 10
            end if
         end if
      else if (mode(1:1).eq.'o') then
         if (.not.qfilda(filna)) then
            write (ntterm, 8040)
 8040       format(1x,'Datei existiert nicht.')
            goto 10
         end if
      end if

      return
      end
 
 
