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
      SUBROUTINE QFILAO (IUNIT, FILNA, OPAR, FMT, LREC,
     *     TEXT, IRERAN, CDFWE, QUESET, IERR)
C     Erfragt vom Benutzer einen Filenamen und oeffnet dann diesen File
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IUNIT     I   I                 FORTRAN-UNIT, auf der der File geoeffnet
C                                 werden soll
C FILNA     C   O  *80            Name des geoeffneten Files
C OPAR      C   I  *2             Eroeffnungsart des Files
C                                 'OP' : OPEN  'CR' : CREATE
C FMT       C   I  *2             Formatierung des Files
C                                 'FO' : FORMATTED  'UN' : UNFORMATTED
C LREC      I   I                 0     : Sequentiellen File oeffnen
C                                 sonst : DIRECT-ACCESS-FILE oeffnen mit
C                                         der log. Recordlaenge LREC
C TEXT      C   I   *(*)          Bei der Abfrage auszugebender Text
C IRERAN    I   I                 Ausrichtung des Textes
C                                     0 : linksbuendig
C                                 sonst : rechtsbuendig in Spalte IRERAN
C CDFWE     C   I  *(*)           Defaultwert der Antwort
C QUESET    I   I                 1=Ueberschreiben; 0=nicht Ueberschreiben
C IERR      I   O                 Fehlercode
C                                 falls ungleich Null, Fehler in
C                                 aufgrufenen Unterprogrammen, bzw. Abbruch
C
      INCLUDE 'imelibh'
C
      CHARACTER*(*) TEXT  ,CDFWE
      CHARACTER*80 FILNA ,QICHAR
      CHARACTER*2  OPAR  ,FMT, mode
      INTEGER IUNIT, LREC, IRERAN, IEING, IERR, LFILNA, QUESET
      LOGICAL  QFLEGP ,QIYENO, QFILDA

      mode=opar
      call qctolo(mode)
      call qfname(filna, opar, text, ireran, cdfwe, queset, ierr)
      if (ierr.ne.0) return

      IF (mode(1:1).EQ.'c') THEN
         IF (LREC.EQ.0) THEN
C
C     Sequentiellen File eroeffnen
C
            CALL QFCRSE(IUNIT,FILNA,FMT,IERR)
            IF (IERR .GT. 0) THEN
               RETURN
            END IF
         ELSE
C
C     DIRECT ACCESS FILE eroeffnen
C
            CALL QFCRDI(IUNIT,FILNA,FMT,LREC,IERR)
            IF (IERR .GT. 0) THEN
               RETURN
            END IF
         END IF
      ELSE
C
C     File zum Lesen oeffnen
C
         IF (LREC.EQ.0) THEN
C
C     Sequentiellen File eroeffnen
C
            CALL QFOPSE (IUNIT,FILNA,'OLD',FMT,IERR)
            IF (IERR .GT. 0) THEN
               RETURN
            END IF
         ELSE
C
C     Random Access File eroeffnen
C
            CALL QFOPDI (IUNIT,FILNA,'OLD',FMT,LREC,IERR)
            IF (IERR .GT. 0) THEN
               RETURN
            END IF
         END IF
      END IF
C
      R E T U R N
      END
 
 
